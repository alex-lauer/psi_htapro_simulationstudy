library(dplyr)
library(purrr)
library(tidyr)
library(survival)


# Load simulation results
sim_out <- readRDS("simulated_pro_datasets.rds")

scenario_names <- names(sim_out)
collection_types <- c("full", "stop", "reduced")

# Cox PH + median deterioration time
run_cox <- function(dat) {

  # Survival object
  fit <- coxph(Surv(PRO_event_week, PRO_status) ~ arm, data = dat)

  cox_sum <- summary(fit)

  # Hazard ratio
  HR <- cox_sum$conf.int["armtreatment", "exp(coef)"]

  # 95% confidence interval
  HR_lower <- cox_sum$conf.int["armtreatment", "lower .95"]
  HR_upper <- cox_sum$conf.int["armtreatment", "upper .95"]

  # Kaplan-Meier curve to get medians and N of events
  km <- survfit(Surv(PRO_event_week, PRO_status) ~ arm, data = dat)

  km_tbl <- as.data.frame(summary(km)$table)

  median_ctrl <- km_tbl["arm=control", "median"]
  median_trt  <- km_tbl["arm=treatment", "median"]

  median_ctrl_L <- km_tbl["arm=control", "0.95LCL"]
  median_ctrl_U <- km_tbl["arm=control", "0.95UCL"]

  median_trt_L  <- km_tbl["arm=treatment", "0.95LCL"]
  median_trt_U  <- km_tbl["arm=treatment", "0.95UCL"]

  # Number of events FROM KM fit
  n_events_ctrl <- km_tbl["arm=control", "events"]
  n_events_trt  <- km_tbl["arm=treatment", "events"]

  tibble(
    HR          = HR,
    HR_lower    = HR_lower,
    HR_upper    = HR_upper,

    n_events_ctrl = n_events_ctrl,
    n_events_trt  = n_events_trt,

    median_ctrl  = median_ctrl,
    median_ctrl_L = median_ctrl_L,
    median_ctrl_U = median_ctrl_U,

    median_trt   = median_trt,
    median_trt_L  = median_trt_L,
    median_trt_U  = median_trt_U,

    median_diff = median_trt - median_ctrl
  )
}

# Risk Ratio at fixed week (27, 36, 54)
run_rr <- function(dat, cutoff_week) {

  # Define event by time cutoff
  dat2 <- dat %>%
    mutate(event_by_cutoff = ifelse(PRO_status == 1 & PRO_event_week <= cutoff_week, 1, 0))

  # # Use Poisson log-link
  # fit <- glm(event_by_cutoff ~ arm,
  #            family = poisson(link = "log"),
  #            data = dat2)
  #
  # coefs <- summary(fit)$coefficients
  #
  # # RR and 95% CI
  # beta <- coefs["armtreatment", "Estimate"]
  # se   <- coefs["armtreatment", "Std. Error"]
  #
  # RR <- exp(beta)
  #
  # z <- qnorm(1 - 0.05 / 2)
  # RR_lower <- exp(beta - z * se)
  # RR_upper <- exp(beta + z * se)

  # Event counts and risk by arm (BY CUTOFF)
  events_tbl <- dat2 %>%
    group_by(arm) %>%
    summarise(
      n        = n(),
      n_events = sum(event_by_cutoff),
      risk     = n_events / n,
      .groups = "drop"
    )

  # N of events
  n_events_ctrl <- events_tbl %>% filter(arm == "control") %>% pull(n_events)
  n_events_trt <- events_tbl %>% filter(arm == "treatment") %>% pull(n_events)

  # N at risk
  n_ctrl  <- events_tbl %>% filter(arm == "control")   %>% pull(n)
  n_trt   <- events_tbl %>% filter(arm == "treatment") %>% pull(n)

  # Risk per arm
  risk_ctrl <- events_tbl %>% filter(arm == "control")   %>% pull(risk)
  risk_trt  <- events_tbl %>% filter(arm == "treatment") %>% pull(risk)

  # Risk CIs
  z <- qnorm(0.975)
  risk_ctrl_se <- sqrt(risk_ctrl * (1 - risk_ctrl) / n_ctrl)
  risk_trt_se  <- sqrt(risk_trt  * (1 - risk_trt)  / n_trt)

  risk_ctrl_L <- max(0, risk_ctrl - z * risk_ctrl_se)
  risk_ctrl_U <- min(1, risk_ctrl + z * risk_ctrl_se)

  risk_trt_L  <- max(0, risk_trt  - z * risk_trt_se)
  risk_trt_U  <- min(1, risk_trt  + z * risk_trt_se)

  # Risk difference
  RD <- risk_trt - risk_ctrl

  RD_se <- sqrt(
    risk_trt * (1 - risk_trt) / n_trt +
      risk_ctrl * (1 - risk_ctrl) / n_ctrl)

  RD_L <- RD - z * RD_se
  RD_U <- RD + z * RD_se

  tibble(
    cutoff_week = cutoff_week,

    # Absolute risks
    risk_ctrl   = risk_ctrl,
    risk_ctrl_L = risk_ctrl_L,
    risk_ctrl_U = risk_ctrl_U,

    risk_trt    = risk_trt,
    risk_trt_L  = risk_trt_L,
    risk_trt_U  = risk_trt_U,

    # Risk difference
    RD          = RD,
    RD_lower    = RD_L,
    RD_upper    = RD_U,

    # RR        = RR,
    # RR_lower  = RR_lower,
    # RR_upper  = RR_upper,

    n_events_ctrl = n_events_ctrl,
    n_events_trt  = n_events_trt
  )
}


# Run analysis across scenarios and simulations
analysis_results <- map_dfr(
  scenario_names,
  function(scn) {

    map_dfr(
      seq_len(length(sim_out[[scn]])),
      function(i) {

        map_dfr(
          collection_types,
          function(coltype) {

            dat <- sim_out[[scn]][[i]][[coltype]]

            # Cox
            cox_res <- run_cox(dat)

            # RR
            rr27 <- run_rr(dat, 27)
            rr36 <- run_rr(dat, 36)
            rr54 <- run_rr(dat, 54)

            tibble(
              scenario   = scn,
              sim        = i,
              collection = coltype,

              # HR + CI
              HR        = cox_res$HR,
              HR_lower  = cox_res$HR_lower,
              HR_upper  = cox_res$HR_upper,

              # KM event counts
              n_events_ctrl = cox_res$n_events_ctrl,
              n_events_trt  = cox_res$n_events_trt,

              # Medians + CI
              median_ctrl   = cox_res$median_ctrl,
              median_ctrl_L = cox_res$median_ctrl_L,
              median_ctrl_U = cox_res$median_ctrl_U,

              median_trt    = cox_res$median_trt,
              median_trt_L  = cox_res$median_trt_L,
              median_trt_U  = cox_res$median_trt_U,

              median_diff = median_trt - median_ctrl,

              # RR at week 27
              risk27_ctrl   = rr27$risk_ctrl,
              risk27_ctrl_L = rr27$risk_ctrl_L,
              risk27_ctrl_U = rr27$risk_ctrl_U,

              risk27_trt    = rr27$risk_trt,
              risk27_trt_L  = rr27$risk_trt_L,
              risk27_trt_U  = rr27$risk_trt_U,

              RD27          = rr27$RD,
              RD27_L        = rr27$RD_lower,
              RD27_U        = rr27$RD_upper,

              # RR27    = rr27$RR,
              # RR27_L  = rr27$RR_lower,
              # RR27_U  = rr27$RR_upper,

              events27_ctrl = rr27$n_events_ctrl,
              events27_trt  = rr27$n_events_trt,

              # RR at week 36
              risk36_ctrl   = rr36$risk_ctrl,
              risk36_ctrl_L = rr36$risk_ctrl_L,
              risk36_ctrl_U = rr36$risk_ctrl_U,

              risk36_trt    = rr36$risk_trt,
              risk36_trt_L  = rr36$risk_trt_L,
              risk36_trt_U  = rr36$risk_trt_U,

              RD36          = rr36$RD,
              RD36_L        = rr36$RD_lower,
              RD36_U        = rr36$RD_upper,

              # RR36    = rr36$RR,
              # RR36_L  = rr36$RR_lower,
              # RR36_U  = rr36$RR_upper,

              events36_ctrl = rr36$n_events_ctrl,
              events36_trt  = rr36$n_events_trt,

              # RR at week 54
              risk54_ctrl   = rr54$risk_ctrl,
              risk54_ctrl_L = rr54$risk_ctrl_L,
              risk54_ctrl_U = rr54$risk_ctrl_U,

              risk54_trt    = rr54$risk_trt,
              risk54_trt_L  = rr54$risk_trt_L,
              risk54_trt_U  = rr54$risk_trt_U,

              RD54          = rr54$RD,
              RD54_L        = rr54$RD_lower,
              RD54_U        = rr54$RD_upper,

              # RR54    = rr54$RR,
              # RR54_L  = rr54$RR_lower,
              # RR54_U  = rr54$RR_upper,

              events54_ctrl = rr54$n_events_ctrl,
              events54_trt  = rr54$n_events_trt
            )
          }
        )
      }
    )
  }
)

# Save analysis output
saveRDS(analysis_results, "analysis_results.rds")
