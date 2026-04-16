library(survival)
library(survminer)
library(patchwork)
library(dplyr)

scenario_names   <- names(sim_out)
collection_types <- c("full", "stop", "reduced")

#### PRO KM curves per scenario ####

pro_km_plots <- list()

for (scn in scenario_names) {

  sim1 <- sim_out[[scn]][[1]]

  for (coltype in collection_types) {

    dat_check <- sim1[[coltype]]

    message("Scenario = ", scn, " | Collection = ", coltype)

    # PRO Kaplan–Meier
    surv_obj_pro <- Surv(
      dat_check$PRO_event_week,
      dat_check$PRO_status
    )

    fit_pro <- survfit(
      surv_obj_pro ~ arm,
      data = dat_check
    )

    # HR
    cox_pro <- coxph(
      Surv(PRO_event_week, PRO_status) ~ arm,
      data = dat_check
    )

    HR_pro <- exp(coef(cox_pro))
    CI_pro <- exp(confint(cox_pro))

    # Medians
    sf_pro <- summary(
      survfit(
        Surv(PRO_event_week, PRO_status) ~ arm,
        data = dat_check
      )
    )

    median_ctrl_pro <- sf_pro$table["arm=control",   "median"]
    median_trt_pro  <- sf_pro$table["arm=treatment", "median"]

    # Event counts
    n_event_ctrl_pro <- sum(dat_check$PRO_status[dat_check$arm == "control"])
    n_event_trt_pro  <- sum(dat_check$PRO_status[dat_check$arm == "treatment"])

    # KM plot
    km_obj <- ggsurvplot(
      fit_pro,
      data = dat_check,
      risk.table = FALSE,
      conf.int = TRUE,
      xlab = "Weeks",
      ylab = "PRO event-free probability",
      title = paste(scn, "•", coltype),
      legend.title = "Arm",
      legend.labs = c("Control", "Treatment")
    )

    # Annotate plot
    pro_km_plots[[paste(scn, coltype, sep = "_")]] <-
      km_obj$plot +
      annotate(
        "text",
        x = Inf, y = 0.15, hjust = 1,
        label = paste0(
          "HR = ", round(HR_pro, 2),
          " (", round(CI_pro[1], 2), ", ", round(CI_pro[2], 2), ")\n",
          "Median C = ", round(median_ctrl_pro, 1),
          " | Median T = ", round(median_trt_pro, 1), "\n",
          "Events C = ", n_event_ctrl_pro,
          " | Events T = ", n_event_trt_pro
        ),
        size = 3
      )
  }
}

# Combine ALL 9 PRO KM plots
combined_PRO_KM <-
  (pro_km_plots$BEFORE_full    | pro_km_plots$BEFORE_stop    | pro_km_plots$BEFORE_reduced) /
  (pro_km_plots$AT_full        | pro_km_plots$AT_stop        | pro_km_plots$AT_reduced) /
  (pro_km_plots$AFTER_full     | pro_km_plots$AFTER_stop     | pro_km_plots$AFTER_reduced) +
  plot_annotation(
    title = "Kaplan–Meier Curves for PRO (First Replication)",
    subtitle = "Annotated with HR, medians, and number of events",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  )

print(combined_PRO_KM)






#### PD KM curve ####
# The same across scenarios

dat_check <- sim_out[["AT"]][[1]][["full"]]

surv_obj_pd <- Surv(
  dat_check$PD_event_week,
  dat_check$PD_status
)

fit_pd <- survfit(
  surv_obj_pd ~ arm,
  data = dat_check
)

# PD HR and medians
cox_pd <- coxph(
  Surv(PD_event_week, PD_status) ~ arm,
  data = dat_check
)

HR_pd <- exp(coef(cox_pd))
CI_pd <- exp(confint(cox_pd))

sf_pd <- summary(
  survfit(
    Surv(PD_event_week, PD_status) ~ arm,
    data = dat_check
  )
)

median_ctrl_pd <- sf_pd$table["arm=control",   "median"]
median_trt_pd  <- sf_pd$table["arm=treatment", "median"]


# PD event counts
n_event_ctrl_pd <- sum(dat_check$PD_status[dat_check$arm == "control"])
n_event_trt_pd  <- sum(dat_check$PD_status[dat_check$arm == "treatment"])


# PD KM plot
pd_km <- ggsurvplot(
  fit_pd,
  data = dat_check,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "PD event-free probability",
  title = paste("PD event-free probability (AT, FULL, 1st))"),
  legend.title = "Arm",
  legend.labs = c("Control", "Treatment")
)

# Annotate PD plot
print(pd_km$plot +
  annotate(
    "text",
    x = Inf, y = 0.5, hjust = 1,
    label = paste0(
      "HR = ", round(HR_pd, 2),
      " (", round(CI_pd[1], 2), ", ", round(CI_pd[2], 2), ")\n",
      "Median C = ", round(median_ctrl_pd, 1),
      " | Median T = ", round(median_trt_pd, 1), "\n",
      "Events C = ", n_event_ctrl_pd,
      " | Events T = ", n_event_trt_pd),
    size = 4))




#### PD medians accross scenarios ####
compute_PD_median <- function(dat) {

  sf_pd <- summary(
    survfit(
      Surv(PD_event_week, PD_status) ~ arm,
      data = dat
    )
  )

  tibble(
    arm = c("control", "treatment"),
    median_PD = c(
      sf_pd$table["arm=control",   "median"],
      sf_pd$table["arm=treatment", "median"]
    )
  )
}

PD_medians_long <- purrr::imap_dfr(
  sim_out,
  function(sim_list, scenario) {

    purrr::imap_dfr(
      sim_list,
      function(res, sim_id) {

        dat_check <- res$full

        compute_PD_median(dat_check) %>%
          mutate(
            scenario = scenario,
            sim      = sim_id
          )
      }
    )
  }
)

PD_medians_summary <- PD_medians_long %>%
  group_by(scenario, arm) %>%
  summarise(
    median_PD = median(median_PD, na.rm = TRUE),
    mean_PD   = mean(median_PD, na.rm = TRUE),
    .groups   = "drop"
  )
