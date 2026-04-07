# Check for CM curve for 1 scenario

dat_check <- sim_out[["AFTER"]][[1]][["stop"]]

surv_obj <- Surv(dat_check$PRO_event_week, dat_check$PRO_status)

fit_2 <- survfit(surv_obj ~ arm, data = dat_check)

ggsurvplot(
  fit_2,
  data = dat_check,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "PRO progression-free survival probability",
  title = "Kaplan–Meier Curve for PRO",
  legend.title = "Arm",
  legend.labs = c("Control", "Treatment")
)

# HR and median time to deterioration for PRO
cox_fit <- coxph(
  Surv(PRO_event_week, PRO_status) ~ arm,
  data = dat_check
)

HR <- exp(coef(cox_fit))
CI <- exp(confint(cox_fit))

fit <- survfit(
  Surv(PRO_event_week, PRO_status) ~ arm,
  data = dat_check
)

median_ctrl <- summary(fit)$table["arm=control", "median"]
median_trt  <- summary(fit)$table["arm=treatment", "median"]
tibble(
  HR        = HR,
  median_ctrl  = median_ctrl,
  median_trt   = median_trt
)


# Check for KM curve for PD

# Build survival object for PD
surv_pd <- Surv(dat_check$PD_event_week, dat_check$PD_status)

# Fit Kaplan–Meier curve
fit_pd <- survfit(surv_pd ~ arm, data = dat_check)

# Plot
ggsurvplot(
  fit_pd,
  data = dat_check,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "PD-free survival probability",
  title = "Kaplan–Meier Curve for PD (Scenario = AFTER, Collection = STOP)",
  legend.title = "Arm",
  legend.labs = c("Control", "Treatment")
)


# HR and median time to deterioration for PD
cox_pd <- coxph(
  Surv(PD_event_week, PD_status) ~ arm,
  data = dat_check
)

HR <- exp(coef(cox_pd))
CI <- exp(confint(cox_pd))

sf <- summary(
  survfit(Surv(PD_event_week, PD_status) ~ arm, data = dat_check)
)

median_ctrl <- sf$table["arm=control",   "median"]
median_trt  <- sf$table["arm=treatment", "median"]

tibble(
  HR_PD        = HR,
  median_ctrl  = median_ctrl,
  median_trt   = median_trt
)
