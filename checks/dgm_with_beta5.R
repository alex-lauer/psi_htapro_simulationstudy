# DGM with beta 5 (PD indicator)
# If PD happened, the slope changes and PRO deteriorates faster (no immediate jump)

library(mvtnorm)
library(truncnorm)
library(tidyr)
library(dplyr)
library(purrr)

# ============================================================
# Parameters
# ============================================================
n_patients   <- 600
n_sims       <- 200
weeks        <- c(seq(3, 27, 3), seq(36, 153, 9))
rho          <- 0.5
sigma        <- 3.5
pro_mean_bl  <- 70
pro_sigma_bl <- 20
beta4        <- 0.01


beta5_PD      <- -1     # when PD happens, slope declines


scenario_names <- c("BEFORE","AT","AFTER")
#scenario_table <- readRDS("scenario_table.rds")
#scenario_table_2 <- readRDS("calibrated_scenario_table_beta5=(-0.1).rds")
scenario_table <- readRDS("calibrated_scenario_table_beta5=(-1).rds")

# Visit-level covariance
pro_vcov <- diag(sigma^2, nrow=23)
pro_vcov[pro_vcov == 0] <- rho*sigma^2

# ============================================================
# Linear PRO generator (with β5)
# ============================================================
pro_generator <- function(pro_bl_val, t, trt, PD_indicator, time_since_PD,
                          beta1, beta3, beta4, beta5) {

  pro_bl_val +
    beta1 * t +
    beta3 * t * trt +
    beta4 * pro_bl_val +
    beta5 * PD_indicator * time_since_PD
}

# ============================================================
# Simulate PD
# ============================================================
lambda_control   <- log(2)/36
lambda_treatment <- log(2)/54

compute_PD <- function(df_patients) {
  df_patients %>%
    mutate(
      PD_time = ifelse(
        arm=="control",
        rexp(n(), lambda_control),
        rexp(n(), lambda_treatment)
      ),

      # map PD_time to next visit time
      PD_event_week_raw = sapply(PD_time, function(t) {
        w <- weeks[weeks >= t]
        if (length(w) == 0) NA_integer_ else w[1]
      }),

      # PD_status = 1 if PD observed, 0 if no PD event
      PD_status = ifelse(is.na(PD_event_week_raw), 0L, 1L),

      # if PD_status == 0 -> PD at administrative censoring = 153
      PD_event_week = ifelse(PD_status == 1,
                             PD_event_week_raw,
                             max(weeks))   # 153
    ) %>%
    select(patient, arm, PD_time, PD_event_week, PD_status)

}

# ============================================================
# PRO event detection
# ============================================================
compute_PRO_event <- function(long_df) {
  long_df %>%
    group_by(patient, arm) %>%
    summarise(
      PRO_status = as.integer(any(week > 0 & PRO_change <= -10)),
      PRO_event_week = ifelse(
        PRO_status==1,
        week[which(PRO_change <= -10 & week > 0)[1]],
        NA_integer_
      ),
      .groups="drop"
    )
}

# ============================================================
# Collection scenarios
# ============================================================
collection_scenarios <- function(long_df, pro_ev_full, pd_tbl) {

  # FULL
  full <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_ev_full, by="patient") %>%
    left_join(pd_tbl,      by="patient") %>%
    mutate(collection="full") %>%
    select(patient, arm, PRO_status, PRO_event_week,
           PD_status, PD_event_week, collection)

  # STOP
  long_stop <- long_df %>%
    left_join(pd_tbl, by=c("patient","arm")) %>%
    filter(PD_status == 0 | week <= PD_event_week)

  pro_ev_stop <- compute_PRO_event(long_stop) %>%
    left_join(pd_tbl %>% select(patient, PD_event_week), by="patient") %>%
    mutate(
      PRO_event_week = ifelse(PRO_status==1, PRO_event_week, PD_event_week)
    ) %>%
    select(-PD_event_week)

  stop <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_ev_stop, by="patient") %>%
    left_join(pd_tbl,      by="patient") %>%
    mutate(collection="stop") %>%
    select(patient, arm, PRO_status, PRO_event_week,
           PD_status, PD_event_week, collection)

  # REDUCED
  long_reduced <- long_df %>%
    left_join(pd_tbl, by=c("patient","arm")) %>%
    group_by(patient) %>%
    mutate(
      visit_idx = row_number() - 1L,
      pd_idx = ifelse(PD_status==0, NA_integer_,
                      (which(week == PD_event_week) - 1L)),
      keep = is.na(pd_idx) |
        (visit_idx <= pd_idx) |
        (visit_idx > pd_idx & ((visit_idx - pd_idx) %% 3L == 0L))
    ) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-visit_idx, -pd_idx, -keep)

  pro_ev_reduced <- compute_PRO_event(long_reduced)

  last_weeks_reduced <- long_reduced %>%
    group_by(patient) %>%
    summarise(last_week_reduced=max(week), .groups="drop")

  pro_ev_reduced <- pro_ev_reduced %>%
    left_join(last_weeks_reduced, by="patient") %>%
    mutate(
      PRO_event_week = ifelse(PRO_status==1,
                              PRO_event_week,
                              last_week_reduced)
    ) %>%
    select(-last_week_reduced)

  reduced <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_ev_reduced, by=c("patient","arm")) %>%
    left_join(pd_tbl,         by=c("patient","arm")) %>%
    transmute(
      patient, arm,
      PRO_status, PRO_event_week,
      PD_status, PD_event_week,
      collection="reduced"
    )

  list(full=full, stop=stop, reduced=reduced)
}

# ============================================================
# SIMULATE() with β5 integrated
# ============================================================
simulate <- function(beta1, beta3, beta5, seed) {

  set.seed(seed)

  # PRO at baseline
  pro_bl <- rtruncnorm(
    n = n_patients,
    mean = pro_mean_bl,
    sd = pro_sigma_bl,
    a = 0,
    b = 100) %>%
    as.data.frame()
  names(pro_bl) <- "PRO_baseline"

  pro_bl$arm     <- sample(rep(c("control","treatment"), each=300))
  pro_bl$TRT     <- as.integer(pro_bl$arm=="treatment")
  pro_bl$patient <- 1:n_patients

  # PD
  pd_tbl <- compute_PD(pro_bl)

  # PRO (wide)
  pro_wide <- pro_bl

  for (wk in weeks) {

    PD_indicator <- as.integer(pd_tbl$PD_event_week <= wk)


    # Time since PD = (wk - PD_event_week), but 0 before PD
    time_since_PD <- pmax(0, wk - pd_tbl$PD_event_week)


    pro_wide[[paste0("week", wk)]] <-
      pro_generator(
        pro_bl_val  = pro_wide$PRO_baseline,
        t           = wk,
        trt         = pro_wide$TRT,
        PD_indicator= PD_indicator,
        time_since_PD = time_since_PD,
        beta1       = beta1,
        beta3       = beta3,
        beta4       = beta4,
        beta5       = beta5
      )
  }

  # Add noise
  errors <- rmvnorm(n_patients, rep(0, length(weeks)), sigma=pro_vcov)
  pro_wide_err <- pro_wide[,5:ncol(pro_wide)] + errors

  pro_wide <- cbind(pro_wide[,1:4], pro_wide_err)
  pro_wide$week0 <- pro_wide$PRO_baseline

  visits <- c("week0", paste0("week", weeks))

  pro_long <- pro_wide %>%
    pivot_longer(all_of(visits), names_to="visit", values_to="PRO") %>%
    mutate(
      week = as.integer(sub("week","",visit)),
      PRO_change = PRO - PRO_baseline
    ) %>%
    select(-visit)

  pro_ev_full <- compute_PRO_event(pro_long) %>%
    mutate(PRO_event_week = ifelse(PRO_status==1, PRO_event_week, max(weeks)))

  collection_scenarios(pro_long, pro_ev_full, pd_tbl)
}

# Run 200 simulations for each time scenario (before, at, after)

sim_out_b5 <- map(
  scenario_names,
  function(scn) {

    row <- scenario_table %>% filter(scenario == scn)
    beta1 <- as.numeric(row$beta1)
    beta3 <- as.numeric(row$beta3)

    map(seq_len(n_sims),
        ~ simulate(
          beta1 = beta1,
          beta3 = beta3,
          beta5 = beta5_PD,
          seed  = 2000 + .x
        )
    )
  }
)

names(sim_out_b5) <- scenario_names










# Analysis

library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# Hazard Ratio from Cox PH + median PRO deterioration
run_cox_b5 <- function(dat) {

  fit <- coxph(Surv(PRO_event_week, PRO_status) ~ arm, data = dat)

  HR <- exp(coef(fit))

  km <- survfit(Surv(PRO_event_week, PRO_status) ~ arm, data = dat)

  median_ctrl <- summary(km)$table["arm=control","median"]
  median_trt  <- summary(km)$table["arm=treatment","median"]

  tibble(
    HR = HR,
    median_ctrl = median_ctrl,
    median_trt = median_trt
  )
}

# Risk Ratio at fixed timepoint
run_rr_b5 <- function(dat, cutoff_week) {

  dat2 <- dat %>%
    mutate(event_by_cutoff = ifelse(PRO_status == 1 & PRO_event_week <= cutoff_week, 1, 0))

  fit <- glm(event_by_cutoff ~ arm, family=poisson(link="log"), data=dat2)

  tibble(RR = exp(coef(fit)["armtreatment"]))
}


collection_types <- c("full", "stop", "reduced")

analysis_b5 <- map_dfr(
  scenario_names,
  function(scn) {

    sims_list <- sim_out_b5[[scn]]   # now a simple list of simulations

    map_dfr(
      seq_len(length(sims_list)),
      function(i) {

        map_dfr(
          collection_types,
          function(coltype) {

            dat <- sims_list[[i]][[coltype]]

            # HR
            cox_res <- run_cox_b5(dat)

            # RR
            rr27 <- run_rr_b5(dat, 27)$RR
            rr36 <- run_rr_b5(dat, 36)$RR
            rr54 <- run_rr_b5(dat, 54)$RR

            tibble(
              scenario   = scn,
              sim        = i,
              collection = coltype,
              HR         = cox_res$HR,
              median_ctrl = cox_res$median_ctrl,
              median_trt  = cox_res$median_trt,
              RR27        = rr27,
              RR36        = rr36,
              RR54        = rr54
            )
          }
        )
      }
    )
  }
)









# Check for CM curve for 1 scenario

dat_check_b5 <- sim_out_b5[["AFTER"]][[1]][["stop"]]

surv_obj <- Surv(dat_check_b5$PRO_event_week, dat_check_b5$PRO_status)

fit_b5 <- survfit(surv_obj ~ arm, data = dat_check_b5)

ggsurvplot(
  fit_b5,
  data = dat_check_b5,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "PRO progression-free survival probability",
  title = "Kaplan–Meier Curve for PRO",
  legend.title = "Arm",
  legend.labs = c("Control", "Treatment"))



run_cox_b5(dat_check_b5) %>%
  summary()

fit <- coxph(Surv(PRO_event_week, PRO_status) ~ arm, data = dat_check_b5)

HR <- exp(coef(fit))

km <- survfit(Surv(PRO_event_week, PRO_status) ~ arm, data = dat_check_b5)

median_ctrl <- summary(km)$table["arm=control","median"]
median_trt  <- summary(km)$table["arm=treatment","median"]

tibble(
  HR = HR,
  median_ctrl = median_ctrl,
  median_trt = median_trt
)






# Check for PD events (KM curve and HR summary)

surv_pd <- Surv(
  time  = dat_check_b5$PD_event_week,
  event = dat_check_b5$PD_status
)

fit_pd <- survfit(surv_pd ~ arm, data = dat_check_b5)

ggsurvplot(
  fit_pd,
  data = dat_check_b5,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "PFS / PD-free probability",
  title = "Kaplan–Meier Curve for PD (Scenario = AFTER, Collection = STOP)",
  legend.title = "Arm",
  legend.labs = c("Control", "Treatment")
)

cox_pd <- coxph(
  Surv(PD_event_week, PD_status) ~ arm,
  data = dat_check_b5
)

HR <- exp(coef(cox_pd))
CI <- exp(confint(cox_pd))

sf <- summary(
  survfit(Surv(PD_event_week, PD_status) ~ arm, data = dat_check_b5)
)

median_ctrl <- sf$table["arm=control",   "median"]
median_trt  <- sf$table["arm=treatment", "median"]

tibble(
  HR_PD        = HR,
#  CI_lower     = CI[1],
#  CI_upper     = CI[2],
  median_ctrl  = median_ctrl,
  median_trt   = median_trt,
#  median_diff  = median_trt - median_ctrl
)



# HR summary for PRO

HR_summary_b5 <- analysis_b5 %>%
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE","AT","AFTER")),
    collection = factor(collection, levels = c("full","reduced","stop"))
  ) %>%
  group_by(scenario, collection) %>%
  summarise(
    HR_mean = exp(mean(log(HR), na.rm = TRUE)),
    p2.5  = quantile(HR, 0.025, na.rm = TRUE),
    p97.5 = quantile(HR, 0.975, na.rm = TRUE),
    median_ctrl_mean = mean(median_ctrl, na.rm = TRUE),
    median_trt_mean  = mean(median_trt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    median_diff = median_trt_mean - median_ctrl_mean
  )


# Order scenarios for plot
HR_summary_b5 <- HR_summary_b5 %>%
  mutate(
    line_id = paste(scenario, collection, sep = " • ")
  )

# desired_order <- c(
#   "AFTER • stop • with_beta5",
#   "AFTER • stop • no_beta5",
#   "AFTER • reduced • with_beta5",
#   "AFTER • reduced • no_beta5",
#   "AFTER • full • with_beta5",
#   "AFTER • full • no_beta5",
#
#   "AT • stop • with_beta5",
#   "AT • stop • no_beta5",
#   "AT • reduced • with_beta5",
#   "AT • reduced • no_beta5",
#   "AT • full • with_beta5",
#   "AT • full • no_beta5",
#
#   "BEFORE • stop • with_beta5",
#   "BEFORE • stop • no_beta5",
#   "BEFORE • reduced • with_beta5",
#   "BEFORE • reduced • no_beta5",
#   "BEFORE • full • with_beta5",
#   "BEFORE • full • no_beta5"
# )

# HR_summary_b5 <- HR_summary_b5 %>%
#   mutate(line_id = factor(line_id, levels = desired_order))

ggplot(HR_summary_b5, aes(x = HR_mean, y = line_id, color = collection)) +

  geom_segment(aes(x = p2.5, xend = p97.5,
                   yend = line_id),
               linewidth = 1.1) +
  geom_point(size = 3) +

  geom_vline(xintercept = 1, linetype="dashed", color="grey40") +

  labs(
    title = "Geometric Mean Hazard Ratio (HR) with 95% Percentile Interval",
    subtitle = "β₅ effect: additional PRO drop after PD (-0.1 point per week)",
    x = "Hazard Ratio (HR)",
    y = "Scenario • Collection • β₅ Condition"
  ) +

  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size=9)
  )




RR_beta5 <- analysis_b5 %>%
  pivot_longer(
    cols = starts_with("RR"),
    names_to = "cutoff",
    values_to = "RR"
  ) %>%
  mutate(
    cutoff = recode(cutoff,
                    RR27 = "Week 27",
                    RR36 = "Week 36",
                    RR54 = "Week 54"),
    scenario   = factor(scenario,   levels = c("BEFORE","AT","AFTER")),
    collection = factor(collection, levels = c("full","reduced","stop"))
  ) %>%
  group_by(scenario, collection, cutoff) %>%
  summarise(
    RR_mean = exp(mean(log(RR), na.rm = TRUE)),   # geometric mean
    p2.5    = quantile(RR, 0.025, na.rm = TRUE),
    p97.5   = quantile(RR, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(line_id = paste(scenario, collection, sep = " • "))








ggplot(RR_beta5,
       aes(x = RR_mean,
           y = line_id,
           color = collection)) +

  # percentile interval
  geom_segment(aes(x = p2.5, xend = p97.5,
                   yend = line_id),
               linewidth = 1.2) +

  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype="dashed", color="grey40") +

  facet_wrap(~ cutoff, scales="fixed") +

  labs(
    title = "Risk Ratio (RR)",
    subtitle = "Geometric Mean RR with 2.5–97.5% Percentile Interval",
    x = "Risk Ratio (RR)",
    y = "Scenario • Collection",
    color = "Collection"
  ) +

  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10)
  )
