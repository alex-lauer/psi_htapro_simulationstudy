############################################################
#  CALIBRATION SCRIPT FOR β1 AND β3  (WITH β5 SLOPE CHANGE)
#  Ensures empirical medians match protocol target medians
############################################################

library(dplyr)
library(purrr)


############################################################
# 1. Function to compute empirical medians for given betas
############################################################

get_empirical_medians <- function(beta1, beta3, beta5,
                                  target_arm = c("control","treatment"),
                                  reps = 100,
                                  n_patients = 400
){
  target_arm <- match.arg(target_arm, several.ok = TRUE)

  med_ctrl <- c()
  med_trt  <- c()

  for (i in 1:reps) {

    # simulate() must exist in the environment
    sim_res <- simulate(beta1, beta3, beta5,
                        seed = 9000 + i)[["full"]]

    # Control median
    med_ctrl[i] <- sim_res %>%
      filter(arm == "control") %>%
      summarise(med = median(PRO_event_week, na.rm=TRUE)) %>%
      pull(med)

    # Treatment median
    med_trt[i] <- sim_res %>%
      filter(arm == "treatment") %>%
      summarise(med = median(PRO_event_week, na.rm=TRUE)) %>%
      pull(med)
  }

  tibble(
    median_ctrl_emp = mean(med_ctrl, na.rm=TRUE),
    median_trt_emp  = mean(med_trt, na.rm=TRUE)
  )
}


############################################################
# 2. Loss function for optimizer
############################################################

loss_function <- function(par, beta5, targets){

  beta1 <- par[1]
  beta3 <- par[2]

  meds <- get_empirical_medians(beta1, beta3, beta5)

  (
    (meds$median_ctrl_emp - targets$ctrl)^2 +
      (meds$median_trt_emp  - targets$trt)^2
  )
}


############################################################
# 3. Optimization wrapper
############################################################

calibrate_betas <- function(
    beta1_init, beta3_init, beta5,
    target_ctrl, target_trt
){
  targets <- list(ctrl = target_ctrl, trt = target_trt)

  opt <- optim(
    par     = c(beta1_init, beta3_init),
    fn      = loss_function,
    beta5   = beta5,
    targets = targets,
    method  = "Nelder-Mead",
    control = list(maxit = 200)
  )

  tibble(
    beta1_calibrated = opt$par[1],
    beta3_calibrated = opt$par[2],
    loss             = opt$value
  )
}


############################################################
# 4. Calibrate for Each Scenario
############################################################

# Your protocol target medians:
targets_before <- list(ctrl = 27, trt = 45)
targets_at     <- list(ctrl = 36, trt = 54)
targets_after  <- list(ctrl = 45, trt = 63)

# Set β5
beta5_value <- -1

# Starting guesses for optimizer
beta1_init <- -0.30
beta3_init <-  0.10

# BEFORE scenario
calib_BEFORE <- calibrate_betas(
  beta1_init, beta3_init,
  beta5 = beta5_value,
  target_ctrl = targets_before$ctrl,
  target_trt  = targets_before$trt
)

# AT scenario
calib_AT <- calibrate_betas(
  beta1_init, beta3_init,
  beta5 = beta5_value,
  target_ctrl = targets_at$ctrl,
  target_trt  = targets_at$trt
)

# AFTER scenario
calib_AFTER <- calibrate_betas(
  beta1_init, beta3_init,
  beta5 = beta5_value,
  target_ctrl = targets_after$ctrl,
  target_trt  = targets_after$trt
)


############################################################
# 5. Combine all calibrated betas into a single table
############################################################

calibrated_betas <- tibble(
  scenario = c("BEFORE","AT","AFTER"),
  beta1    = c(calib_BEFORE$beta1_calibrated,
               calib_AT$beta1_calibrated,
               calib_AFTER$beta1_calibrated),
  beta3    = c(calib_BEFORE$beta3_calibrated,
               calib_AT$beta3_calibrated,
               calib_AFTER$beta3_calibrated),
  beta5    = beta5_value
)

print(calibrated_betas)

# Optionally save for use in simulation

#saveRDS(calibrated_betas, "calibrated_scenario_table_beta5=(-1).rds")

############################################################
# END OF SCRIPT
############################################################
