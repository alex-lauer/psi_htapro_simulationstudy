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
  fit <- coxph(Surv(event_week, status) ~ arm, data = dat)

  # Kaplan-Meier curve to get medians
  km <- survfit(Surv(event_week, status) ~ arm, data = dat)

  median_ctrl <- summary(km)$table["arm=control","median"]
  median_trt  <- summary(km)$table["arm=treatment","median"]

  tibble(
    HR = exp(coef(fit)),
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

# Save analysis output
saveRDS(analysis_results, "analysis_results.rds")
