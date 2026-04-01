library(mvtnorm)
library(truncnorm)
library(tidyr)
library(dplyr)
library(purrr)

# Parameters
n_patients   <- 600         # number of patients per simulation
n_sims       <- 200         # simulations per scenario
weeks        <- c(seq(3, 27, 3), seq(36, 153, 9))  # 23 post-baseline visits
rho          <- 0.5         # correlation between visits
sigma        <- 3.5         # SD of visit-level errors
pro_mean_bl  <- 70          # baseline mean
pro_sigma_bl <- 20          # baseline SD
beta4        <- 0.01        # effect of baseline PRO on trajectory (constant across scenarios)
scenario_names <- c("BEFORE", "AT", "AFTER")

# Read scenario parameters
scenario_table <- readRDS("scenario_table.rds")

# Create a covariance matrix for errors for every visit
pro_vcov <- diag(sigma^2, nrow = 23, ncol = 23)
pro_vcov[pro_vcov == 0] <- rho*sigma^2

# Generator for expected PRO at time t
pro_generator <- function(pro_bl_val, t, trt, beta1, beta3, beta4) {
  pro_bl_val + (beta1 * t) + (beta3 * t * trt) + (beta4 * pro_bl_val)
}

# Data collection scenarios
collection_scenarios <- function(long_df, event_df, event_df_reduced) {
  merged <- long_df %>%
    distinct(patient, arm) %>%
    left_join(event_df, by = "patient")

  # Full data collection until the end of study
  full <- merged

  # STOP scenario.
  # Censor patients (set status to 0) if event happened after cutoff:
  # - week 36 for control
  # - week 54 for treatment

  stop <- merged %>%
    mutate(
      cutoff = ifelse(arm == "control", 36, 54)
    ) %>%
    mutate(status = ifelse(event_week <= cutoff, 1, 0)) %>%
    select(-cutoff)

  # REDUCED scenario
  reduced <- long_df %>%
    distinct(patient, arm) %>%
    left_join(event_df_reduced, by = "patient")


  return(list(
    full    = full,
    stop    = stop,
    reduced = reduced
  ))
}


# Simulation function

simulate <- function(beta1, beta3, seed) {
  set.seed(seed)

  # PRO at baseline
  pro_bl <- rtruncnorm(
    n = n_patients,
    mean = pro_mean_bl,
    sd = pro_sigma_bl,
    a = 0,       # lower bound
    b = 100      # upper bound
  ) %>%
    as.data.frame()

  names(pro_bl) <- "PRO_baseline"

  # Randomize patients to treatment and control 1:1
  pro_bl$arm <- sample(rep(c("control", "treatment"), each = 300))
  pro_bl$TRT <- as.integer(pro_bl$arm == "treatment")

  # Generate expected PRO values for each visit (without error)
  # Build a wide data frame
  pro_wide <- pro_bl

  for (wk in weeks) {
    pro_wide[[paste0("week", wk)]] <- pro_generator(
      pro_bl_val = pro_wide$PRO_baseline,
      t          = wk,
      trt        = pro_wide$TRT,
      beta1      = beta1,
      beta3      = beta3,
      beta4      = beta4
    )}

  # Draw errors for visits
  errors <- mvtnorm::rmvnorm(
    n = n_patients,
    mean = rep(0, nrow(pro_vcov)),
    sigma = pro_vcov) %>%
    as.data.frame()

  # Add correlated errors to generated values of PROs
  pro_wide_err <- pro_wide[, 4:ncol(pro_wide)] + errors

  # Clamp to [0, 100]
  pro_wide_err[pro_wide_err < 0] <- 0
  pro_wide_err[pro_wide_err > 100] <- 100

  # Recombine with patient metadata
  pro_wide <- cbind(pro_wide[, 1:3], pro_wide_err)

  # Add week0
  pro_wide$week0 <- pro_wide$PRO_baseline

  # Add patient column
  pro_wide$patient <- c(1:600)

  # Identify visit columns
  visits <- c("week0", paste0("week", weeks))

  # Reorder
  pro_wide <- pro_wide %>%
    dplyr::select(patient, arm, TRT, PRO_baseline, all_of(visits))

  # Build a long dataset: one row per patient-week
  pro_long <- pro_wide %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(visits),
      names_to = "visit",
      values_to = "PRO"
    ) %>%
    dplyr::mutate(
      week = as.integer(sub("week", "", visit))
    ) %>%
    dplyr::select(-c(visit))

  # Compute change from baseline (Delta PRO)
  pro_long <- pro_long %>%
    dplyr::mutate(PRO_change = PRO - PRO_baseline)

  # Event detection (first visit with deterioration: ΔPRO <= -10)
  event_df <- pro_long %>%
    group_by(patient) %>%
    summarise(
      # TRUE/FALSE if any deterioration happened
      status = any(PRO_change <= -10),

      # first week of deterioration (or NA if none)
      event_week = ifelse(
        status,
        week[which(PRO_change <= -10)[1]],
        NA_integer_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      status = as.integer(status)
    )

  # Create tables for reduced data collection scenario.
  # After the cutoff week keep only every 3rd visit:
  # For CONTROL (cutoff = 36): weeks 63, 90, 117, 144
  # For TREATMENT (cutoff = 54): weeks 81, 108, 135

  # Remove all weeks post cut-off that are not mentioned above
  long_filtered <- pro_long %>%
    mutate(
      cutoff = ifelse(arm == "control", 36, 54),
      keep_post = (week > cutoff) & ((week - cutoff) %% 27 == 0)
    ) %>%
    group_by(patient) %>%
    filter(
      week <= cutoff | keep_post
    ) %>%
    ungroup()

  # Event detection for reduced scenario
  event_df_reduced <- long_filtered %>%
    group_by(patient) %>%
    summarise(
      # TRUE/FALSE if any deterioration happened
      status = any(PRO_change <= -10),
      # first week of deterioration (or NA if none)
      event_week = ifelse(
        status,
        week[which(PRO_change <= -10)[1]],
        NA_integer_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      status = as.integer(status)   # convert TRUE/FALSE → 1/0
    )


  # Build data-collection scenario versions
  collections <- collection_scenarios(pro_long, event_df, event_df_reduced)

  return(collections)
}



# Run 200 simulations for each time scenario (before, at, after)
old_sim <- purrr::map(
  scenario_names,
  function(scn) {
    # Get betas for this scenario
    row   <- scenario_table %>% dplyr::filter(.data$scenario == scn)
    beta1 <- as.numeric(row$beta1[[1]])
    beta3 <- as.numeric(row$beta3[[1]])

    # Run n_sims simulations with different seeds
    purrr::map(
      seq_len(n_sims),
      ~ simulate(
        beta1 = beta1,
        beta3 = beta3,
        seed  = 1000 + .x
      )
    )
  }
)

names(old_sim) <- scenario_names




# Analysis

scenario_names <- names(old_sim)
collection_types <- c("full", "stop", "reduced")

# Cox PH + median deterioration time
run_cox <- function(dat) {

  # Survival object
  fit <- coxph(Surv(event_week, status) ~ arm, data = dat)

  HR <- exp(coef(fit))

  # Kaplan-Meier curve to get medians
  km <- survfit(Surv(event_week, status) ~ arm, data = dat)

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
    mutate(event_by_cutoff = ifelse(status == 1 & event_week <= cutoff_week, 1, 0))

  # Use Poisson log-link
  fit <- glm(event_by_cutoff ~ arm,
             family = poisson(link = "log"),
             data = dat2)

  tibble(RR = exp(coef(fit)["armtreatment"]))

}


# Run analysis across scenarios and simulations
old_analysis_results <- map_dfr(
  scenario_names,
  function(scn) {

    map_dfr(
      seq_len(length(old_sim[[scn]])),
      function(i) {

        map_dfr(
          collection_types,
          function(coltype) {

            dat <- old_sim[[scn]][[i]][[coltype]]

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

# KM curve check
library(survival)
library(survminer)

dat<- sim[["AFTER"]][[1]][["stop"]]

surv_obj <- Surv(dat$PRO_event_week, dat$PRO_status)

fit_2 <- survfit(surv_obj ~ arm, data = dat)


ggsurvplot(
  fit_2,
  data = dat,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "PRO progression-free survival probability",
  title = "Kaplan–Meier Curve for PRO",
  legend.title = "Treatment Arm",
  legend.labs = c("Control", "Treatment")
)



#### HR across scenarios ####
HR <- old_analysis_results %>%
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
ggplot(HR,
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



#### RR differences between scenarios ####


# Build summary for RR at weeks 27, 36, and 54
RR <- old_analysis_results %>%
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

