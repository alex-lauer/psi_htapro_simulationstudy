#### Parameters for PRO model
#### Simulation scenarios: PRO deterioration happens
#### before, at the same time or after PD

library(dplyr)
library(purrr)

#### Fixed parameters ####

mu_P0   <- 70          # mean baseline PRO
sd_P0   <- 20          # sd baseline PRO
beta0   <- 0           # intercept
beta2   <- 0           # no baseline treatment offset
beta4   <- 2/70        # baseline-PRO effect on slope
MCID    <- -10         # PRO decline target (MCID threshold)

#### Change from baseline: CfB(t) = β0 + β1*t + β2*Trt + β3*t*Trt + β4*P0

#### Scenario definitions: target deterioration times ####

deterioration_weeks <- list(
  BEFORE = list(ctrl = 27, trt = 45),
  AT     = list(ctrl = 36, trt = 54),
  AFTER  = list(ctrl = 45, trt = 63)
)

#### Function to compute beta1 and beta3 ####

compute_betas <- function(t_ctrl, t_trt, beta4, mu_P0) {

  # expected PRO at baseline contribution
  baseline_term <- beta4 * mu_P0

  # control arm:
  # MCID = beta1 * t_ctrl + beta4 * P0
  beta1 <- (MCID - baseline_term) / t_ctrl

  # treatment arm:
  # MCID = (beta1 + beta3) * t_trt + beta4 * P0
  beta3 <- (MCID - baseline_term) / t_trt - beta1

  return(list(beta1 = beta1, beta3 = beta3))
}


#### Build scenario table ####

scenario_table <- tibble(
  scenario = names(deterioration_weeks),
  t_ctrl   = map_dbl(deterioration_weeks, ~ .x$ctrl),
  t_trt    = map_dbl(deterioration_weeks, ~ .x$trt)
) %>%
  mutate(
    betas = map2(t_ctrl, t_trt,
                 ~ compute_betas(.x, .y, beta4 = beta4, mu_P0 = mu_P0)),
    beta1 = map_dbl(betas, "beta1"),
    beta3 = map_dbl(betas, "beta3")
  ) %>%
  select(scenario, t_ctrl, t_trt, beta1, beta3)

scenario_table

