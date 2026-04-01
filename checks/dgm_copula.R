# PFS and PRO at baseline are simulated using Copula and correlated
# Then PRO observations are added to every other visit using linear model,
# adjusted for PRO at baseline and accounted for within-subject correlation

library(copula)
library(dplyr)
library(tidyr)
library(truncnorm)
library(mvtnorm)
library(purrr)

# ============================================================
# Parameters
# ============================================================
n_patients   <- 600
n_sims       <- 200
weeks        <- c(seq(3, 27, 3), seq(36, 153, 9))  # 23 post-baseline visits
rho          <- 0.5
sigma        <- 3.5
pro_mean_bl  <- 70
pro_sigma_bl <- 20
beta4        <- 0.01
scenario_names <- c("BEFORE","AT","AFTER")

# PD exponential parameters (based on medians)
lambda_control   <- log(2)/36    # control median 36 weeks
lambda_treatment <- log(2)/54    # treatment median 54 weeks

# Visit-level error covariance
pro_vcov <- diag(sigma^2, 23)
pro_vcov[pro_vcov == 0] <- rho*sigma^2

# Scenario parameters
scenario_table <- readRDS("scenario_table.rds")

# Gaussian copula correlation between baseline PRO and PFS
rho_cop <- 0.5


# ============================================================
# PRO trajectory generator
# ============================================================
pro_generator <- function(pro_bl_val, t, trt, beta1, beta3, beta4) {
  pro_bl_val + (beta1 * t) + (beta3 * t * trt) + (beta4 * pro_bl_val)
}


# ============================================================
# Helper: map continuous PD time → nearest visit week
# ============================================================
nearest_visit <- function(t) {
  wk <- weeks[weeks >= t]
  if (length(wk)==0) NA_integer_ else wk[1]
}


# ============================================================
# Helper: detect PRO event (ΔPRO <= -10)
# ============================================================
compute_PRO_event <- function(long_df) {
  long_df %>%
    group_by(patient, arm) %>%
    summarise(
      PRO_status = as.integer(any(week > 0 & PRO_change <= -10)),
      PRO_event_week = ifelse(
        PRO_status == 1,
        week[which(PRO_change <= -10 & week > 0)[1]],
        NA_integer_
      ),
      .groups = "drop"
    )
}


# ============================================================
# Data collection scenarios: full / stop / reduced
# ============================================================
collection_scenarios <- function(long_df, pro_ev_full, pd_tbl) {

  # FULL
  full <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_ev_full, by="patient") %>%
    left_join(pd_tbl, by="patient") %>%
    mutate(collection="full") %>%
    select(patient, arm, PRO_status, PRO_event_week,
           PD_status, PD_event_week, collection)


  # STOP (drop PRO after PD)
  long_stop <- long_df %>%
    left_join(pd_tbl, by=c("patient","arm")) %>%
    filter(is.na(PD_event_week) | week <= PD_event_week)

  pro_ev_stop <- compute_PRO_event(long_stop)

  pro_ev_stop <- pro_ev_stop %>%
    left_join(pd_tbl %>% select(patient, PD_event_week), by="patient") %>%
    mutate(
      PRO_event_week = ifelse(
        PRO_status == 1,
        PRO_event_week,
        PD_event_week  # censor at progression
      )
    ) %>%
    select(-PD_event_week)


  stop <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_ev_stop, by="patient") %>%
    left_join(pd_tbl,      by="patient") %>%
    mutate(collection="stop") %>%
    select(patient, arm, PRO_status, PRO_event_week,
           PD_status, PD_event_week, collection)


  # REDUCED (keep every 3rd visit after PD)
  long_reduced <- long_df %>%
    left_join(pd_tbl, by=c("patient","arm")) %>%
    group_by(patient) %>%
    mutate(
      visit_idx = row_number() - 1L,
      pd_idx = match(PD_event_week, week) - 1L,
      keep_pre  = !is.na(PD_event_week) & week <= PD_event_week,
      keep_post = !is.na(pd_idx) &
        visit_idx > pd_idx &
        ((visit_idx - pd_idx) %% 3L == 0L),
      keep = is.na(pd_idx) | keep_pre | keep_post
    ) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-visit_idx, -pd_idx, -keep_pre, -keep_post, -keep)

  pro_reduced_ev <- compute_PRO_event(long_reduced)

  # Censor at last observed PRO week if PRO event never happened
  last_weeks_reduced <- long_reduced %>%
    group_by(patient) %>%
    summarise(last_week_reduced = max(week), .groups="drop")

  pro_reduced_ev <- pro_reduced_ev %>%
    left_join(last_weeks_reduced, by="patient") %>%
    mutate(
      PRO_event_week = ifelse(
        PRO_status == 1,
        PRO_event_week,
        last_week_reduced   # censor at last reduced-visit
      )
    ) %>%
    select(-last_week_reduced)

  reduced <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_reduced_ev, by="patient") %>%
    left_join(pd_tbl,         by="patient") %>%
    mutate(collection="reduced") %>%
    select(patient, arm, PRO_status, PRO_event_week,
           PD_status, PD_event_week, collection)


  list(full=full, stop=stop, reduced=reduced)
}


# ============================================================
# SIMULATION FUNCTION
# ============================================================
simulate <- function(beta1, beta3, seed) {
  set.seed(seed)

  # ----------------------------------------------------------
  # 1) Use Gaussian copula to simulate baseline PRO + PFS time
  # ----------------------------------------------------------
  cop <- normalCopula(param = rho_cop, dim = 2)
  U <- rCopula(n_patients, cop)

  U_PRO0 <- U[,1]   # baseline PRO uniform
  U_PFS  <- U[,2]   # PFS time uniform

  # Baseline PRO from truncated normal
  PRO_baseline <- qtruncnorm(
    U_PRO0,
    a = 0, b = 100,
    mean = pro_mean_bl,
    sd = pro_sigma_bl
  )

  # Randomize arm
  arm <- sample(rep(c("control","treatment"), each = n_patients/2))

  # PFS time via exponential marginals
  PD_time <- ifelse(
    arm=="control",
    qexp(U_PFS, rate=lambda_control),
    qexp(U_PFS, rate=lambda_treatment)
  )

  PD_event_week <- sapply(PD_time, nearest_visit)
  PD_status <- ifelse(is.na(PD_event_week), 0L, 1L)

  pd_tbl <- tibble(
    patient = 1:n_patients,
    arm,
    PD_time,
    PD_event_week,
    PD_status
  )

  # ----------------------------------------------------------
  # 2) Generate longitudinal PRO from linear model
  # ----------------------------------------------------------
  pro_bl <- tibble(
    patient = 1:n_patients,
    arm     = arm,
    TRT     = ifelse(arm=="treatment",1,0),
    PRO_baseline = PRO_baseline
  )

  # Expected PRO without error
  pro_wide <- pro_bl
  for (wk in weeks) {
    pro_wide[[paste0("week", wk)]] <- pro_generator(
      pro_bl_val = pro_bl$PRO_baseline,
      t          = wk,
      trt        = pro_bl$TRT,
      beta1      = beta1,
      beta3      = beta3,
      beta4      = beta4
    )
  }

  # Apply correlated visit-level errors
  errors <- rmvnorm(n_patients, rep(0, length(weeks)), sigma=pro_vcov)
  pro_wide_err <- pro_wide[,4:ncol(pro_wide)] + errors

  # Clamp to [0,100]
  pro_wide_err[pro_wide_err < 0] <- 0
  pro_wide_err[pro_wide_err > 100] <- 100

  # Assemble final wide PRO table
  pro_wide <- cbind(
    pro_wide[,1:3],
    pro_wide_err,
    week0 = PRO_baseline
  )

  visits <- c("week0", paste0("week", weeks))
  pro_wide <- pro_wide %>% select(patient, arm, TRT, PRO_baseline, all_of(visits))

  # Long format
  pro_long <- pro_wide %>%
    pivot_longer(cols = all_of(visits),
                 names_to = "visit",
                 values_to = "PRO") %>%
    mutate(
      week = as.integer(sub("week","",visit)),
      PRO_change = PRO - PRO_baseline
    ) %>%
    select(-visit)

  # PRO event (full data)
  pro_ev_full <- compute_PRO_event(pro_long)


  pro_ev_full <- pro_ev_full %>%
    mutate(
      PRO_event_week = ifelse(
        PRO_status == 1,
        PRO_event_week,
        153L  # censor at final visit
      )
    )


  # ----------------------------------------------------------
  # 3) Build collection scenarios (full/stop/reduced)
  # ----------------------------------------------------------
  collection_scenarios(
    long_df      = pro_long,
    pro_ev_full  = pro_ev_full,
    pd_tbl       = pd_tbl
  )
}


# ============================================================
# RUN ALL SCENARIOS
# ============================================================
sim <- map(
  scenario_names,
  function(scn) {
    row <- scenario_table %>% filter(scenario == scn)
    beta1 <- row$beta1
    beta3 <- row$beta3

    map(seq_len(n_sims),
        ~ simulate(beta1, beta3, seed = 1000 + .x))
  }
)

names(sim) <- scenario_names











# Analysis and performance check for new DGM (with copula for PFS and PRO at baseline)


# Cox PH + median deterioration time
run_cox <- function(dat) {

  # Survival object
  fit <- coxph(Surv(PRO_event_week, PRO_status) ~ arm, data = dat)

  HR <- exp(coef(fit))

  # Kaplan-Meier curve to get medians
  km <- survfit(Surv(PRO_event_week, PRO_status) ~ arm, data = dat)

  median_ctrl <- summary(km)$table["arm=control","median"]
  median_trt  <- summary(km)$table["arm=treatment","median"]

  tibble(
    HR = HR,
    median_ctrl = median_ctrl,
    median_trt = median_trt
  )
}

# Risk Ratio at fixed week (27, 36, 54)
run_rr <- function(dat, cutoff_week) {

  # Define event by time cutoff
  dat2 <- dat %>%
    mutate(event_by_cutoff = ifelse(PRO_status == 1 & PRO_event_week <= cutoff_week, 1, 0))

  # Use Poisson log-link
  fit <- glm(event_by_cutoff ~ arm,
             family = poisson(link = "log"),
             data = dat2)

  tibble(RR = exp(coef(fit)["armtreatment"]))

}


# Run analysis across scenarios and simulations
analysis_copula <- map_dfr(
  scenario_names,
  function(scn) {

    map_dfr(
      seq_len(length(sim[[scn]])),
      function(i) {

        map_dfr(
          collection_types,
          function(coltype) {

            dat <- sim[[scn]][[i]][[coltype]]

            # Cox
            cox_res <- run_cox(dat)

            # RR
            rr27 <- run_rr(dat, 27)$RR
            rr36 <- run_rr(dat, 36)$RR
            rr54 <- run_rr(dat, 54)$RR

            tibble(
              scenario = scn,
              sim = i,
              collection = coltype,
              HR = cox_res$HR,
              median_ctrl = cox_res$median_ctrl,
              median_trt = cox_res$median_trt,
              RR27 = rr27,
              RR36 = rr36,
              RR54 = rr54
            )
          }
        )
      }
    )
  }
)







# Performance

# Check for CM curve for 1 scenario

dat_cop <- sim[["AFTER"]][[1]][["full"]]

surv_obj <- Surv(dat_cop$PRO_event_week, dat_cop$PRO_status)

fit_cop <- survfit(surv_obj ~ arm, data = dat_cop)

ggsurvplot(
  fit_cop,
  data = dat_cop,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "PRO progression-free survival probability",
  title = "Kaplan–Meier Curve for PRO",
  legend.title = "Arm",
  legend.labs = c("Control", "Treatment")
)


#### HR across scenarios ####
HR_cop <- analysis_copula %>%
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop"))
  )%>%
  group_by(scenario, collection) %>%
  # Geometric mean - convert HR to a log scale before calculating mean
  summarise(HR_mean = exp(mean(log(HR), na.rm = TRUE)),

            # # standard error on the log scale
            # logHR_se = sd(log(HR), na.rm = TRUE) / sqrt(sum(!is.na(HR))),
            #
            # # confidence interval for mean of log(HR), then exponentiate
            # CI_lower = exp(mean(log(HR), na.rm = TRUE) - 1.96 * logHR_se),
            # CI_upper = exp(mean(log(HR), na.rm = TRUE) + 1.96 * logHR_se),

            # 95% percentile interval for HR distribution
            p2.5  = quantile(HR, 0.025, na.rm = TRUE),
            p97.5 = quantile(HR, 0.975, na.rm = TRUE),


            # Mean median deterioration times across simulations
            median_ctrl_mean = mean(median_ctrl, na.rm = TRUE),
            median_trt_mean  = mean(median_trt,  na.rm = TRUE),

            .groups = "drop") %>%
  mutate(median_diff = median_trt_mean - median_ctrl_mean) %>%
  select(scenario, collection, HR_mean, p2.5, p97.5, median_ctrl_mean, median_trt_mean, median_diff)


# Plot for geometric mean HR with 2.5–97.5% percentile interval
ggplot(HR_cop,
       aes(x = HR_mean,
           y = interaction(scenario, collection, sep = " • "),
           color = collection)) +
  geom_segment(aes(x = p2.5, xend = p97.5,
                   yend = interaction(scenario, collection, sep = " • ")),
               linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  labs(
    title = "Geometric mean HR with 2.5–97.5% percentile interval",
    x = "Hazard Ratio (HR)",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10))


# Build summary for RR at weeks 27, 36, and 54
RR <- analysis_copula %>%
  select(scenario, collection, sim, RR27, RR36, RR54) %>%
  pivot_longer(
    cols = starts_with("RR"),
    names_to = "cutoff",
    values_to = "RR") %>%
  mutate(cutoff = recode(cutoff,
                         RR27 = "Week 27",
                         RR36 = "Week 36",
                         RR54 = "Week 54")) %>%
  group_by(scenario, collection, cutoff) %>%
  summarise(
    # Geometric mean RR - log scale
    RR_mean = exp(mean(log(RR), na.rm = TRUE)),
    # Percentile interval across 200 simulations
    p2.5  = quantile(RR, 0.025, na.rm = TRUE),
    p97.5 = quantile(RR, 0.975, na.rm = TRUE),
    .groups = "drop")

# Wide table
rr_wide <- RR %>%
  pivot_wider(
    names_from  = cutoff,
    values_from = c(RR_mean, `p2.5`, `p97.5`)) %>%
  arrange(scenario, collection)

# Plot for geometric mean RR with 2.5–97.5% percentile interval

ggplot(RR,
       aes(x = RR_mean,
           y = interaction(scenario, collection, sep = " • "),
           color = collection)) +
  geom_segment(aes(x = p2.5, xend = p97.5,
                   yend = interaction(scenario, collection, sep = " • ")),
               linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  facet_wrap(~ cutoff, scales = "fixed") +
  labs(
    title = "Geometric mean RR with 2.5–97.5% percentile interval",
    x = "Risk Ratio (RR)",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10))
