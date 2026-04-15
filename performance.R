# Performance Measures

library(dplyr)
library(ggplot2)
library(tidyr)
library(survminer)
library(stringr)
library(scales)

#Load analysis results
analysis_results <- readRDS("analysis_results.rds")


#### HR across scenarios ####
##### HR summary table #####
HR <- analysis_results %>%
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop"))
  )%>%
  group_by(scenario, collection) %>%
  # Geometric mean - convert HR to a log scale before calculating mean
  summarise(
            # Geometric mean HR
            HR_mean = exp(mean(log(HR), na.rm = TRUE)),

            # Mean Cox PH Wald CI (log-scale averaging)
            HR_CI_lower = exp(mean(log(HR_lower), na.rm = TRUE)),
            HR_CI_upper = exp(mean(log(HR_upper), na.rm = TRUE)),

            # # 95% percentile interval for HR distribution
            # p2.5  = quantile(HR, 0.025, na.rm = TRUE),
            # p97.5 = quantile(HR, 0.975, na.rm = TRUE),

            # Event counts (from KM)
            events_ctrl = mean(n_events_ctrl, na.rm = TRUE),
            events_trt  = mean(n_events_trt,  na.rm = TRUE),

            # Mean median deterioration times across simulations
            median_ctrl = mean(median_ctrl, na.rm = TRUE),
            median_trt  = mean(median_trt,  na.rm = TRUE),
            median_diff = mean(median_diff, na.rm = TRUE),

            # Median CI
            median_ctrl_L = mean(median_ctrl_L, na.rm = TRUE),
            median_ctrl_U = mean(median_ctrl_U, na.rm = TRUE),
            median_trt_L  = mean(median_trt_L,  na.rm = TRUE),
            median_trt_U  = mean(median_trt_U,  na.rm = TRUE),

            .groups = "drop") %>%
  mutate(line_id = factor(
      paste(scenario, collection, sep = " • "),
      levels = unique(paste(scenario, collection, sep = " • ")))) %>%
  select(-median_ctrl_L, -median_ctrl_U, -median_trt_L, -median_trt_U)


##### Plot for geometric mean HR with 95%CI #####
HR_reorder <- HR %>%
  mutate(line_id = factor(line_id, levels = rev(unique(line_id))))

ggplot(HR_reorder,
       aes(x = HR_mean,
           y = line_id,
           color = collection)) +
  geom_segment(aes(x = HR_CI_lower, xend = HR_CI_upper,
                   yend = line_id),
               linewidth = 1.1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "Geometric mean HR with 95% CI",
    x = "Hazard Ratio (HR)",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(legend.position = "bottom",
    axis.text.y = element_text(size = 9))



##### Power for HR #####
# Power is defined as the probability of rejecting the null hypothesis of no treatment effect
# Rejection rule: 95%CI excludes 1.

power_HR <- analysis_results %>%
  mutate(scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop")))%>%
  mutate(reject = HR_upper < 1) %>%
  group_by(scenario, collection) %>%
  summarise(power = mean(reject),
    n_sims = n(),
    .groups = "drop")


##### Coverage for HR #####
# Coverage for a scenario = probability that the CI from that scenario contains the true HR,
# where the true HR is defined by FULL.

# Define true HR from full scenario
HR_full <-HR %>%
  filter(collection == "full") %>%
  transmute(scenario, HR_full = HR_mean)

# Does the CI for each replication cover true HR?
covered <- analysis_results %>%
  left_join(HR_full, by = "scenario")

covered <- covered %>%
  mutate(scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop")))%>%
  mutate(covered = HR_lower <= HR_full & HR_upper >= HR_full)

coverage_HR <- covered %>%
  group_by(scenario, collection) %>%
  summarise(coverage = mean(covered, na.rm = TRUE),
    n_sims   = n(),
    .groups  = "drop")

# Combine power and coverage table
power_coverage_HR <- power_HR %>%
  left_join(coverage_HR,
    by = c("scenario", "collection")) %>%
  select(-n_sims.x, -n_sims.y)

##### Combine main HR table with power and coverage #####
HR_perf <- HR %>%
  select(scenario, collection,
    HR_mean, HR_CI_lower, HR_CI_upper,
    events_ctrl, events_trt,
    median_ctrl, median_trt, median_diff) %>%
  left_join(power_coverage_HR %>%
      select(scenario, collection, power, coverage),
    by = c("scenario", "collection"))


##### Bias in HR vs full scenario ####
# Evaluate reduced/stop vs full
bias_HR <- HR_reorder %>%
  left_join(HR_full, by = "scenario") %>%
  filter(collection %in% c("reduced", "stop")) %>%
  mutate(
    # Bias of HR means on a log scale
    bias = log(HR_mean) - log(HR_full))

# Plot bias for HR means
ggplot(bias_HR, aes(
    x = bias,
    y = line_id,
    color = collection)) +
  geom_segment(aes(x = 0, xend = bias,
                   yend = line_id),
               linewidth = 1.1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title = "Bias in mean HR (log scale)",
    subtitle = "Bias = log(HR_mean(reduced/stop)) − log(HR_mean(full))",
    x = "Bias in mean HR",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10))





#### Risk differences between scenarios ####
# Build long table for RD and 95% CIs at weeks 27, 36, and 54
RD_long <- analysis_results %>%
  select(
    scenario, collection, sim,

    # Week 27
    risk27_ctrl, risk27_ctrl_L, risk27_ctrl_U,
    risk27_trt,  risk27_trt_L,  risk27_trt_U,
    RD27, RD27_L, RD27_U,
    events27_ctrl, events27_trt,

    # Week 36
    risk36_ctrl, risk36_ctrl_L, risk36_ctrl_U,
    risk36_trt,  risk36_trt_L,  risk36_trt_U,
    RD36, RD36_L, RD36_U,
    events36_ctrl, events36_trt,

    # Week 54
    risk54_ctrl, risk54_ctrl_L, risk54_ctrl_U,
    risk54_trt,  risk54_trt_L,  risk54_trt_U,
    RD54, RD54_L, RD54_U,
    events54_ctrl, events54_trt

    # RR27, RR27_L, RR27_U, events27_ctrl, events27_trt,
    # RR36, RR36_L, RR36_U, events36_ctrl, events36_trt,
    # RR54, RR54_L, RR54_U, events54_ctrl, events54_trt
  ) %>%
  pivot_longer(
    cols = -c(scenario, collection, sim),
    names_to = "name",
    values_to = "value"
  ) %>%
  mutate(
    cutoff   = as.integer(str_extract(name, "\\d+")),
    variable = str_remove(name, "\\d+")
  ) %>%
  select(-name) %>%
  pivot_wider(
    names_from  = variable,
    values_from = value
  )

##### Summary for RD #####
RD_summary <- RD_long %>%
  mutate(
    scenario   = factor(scenario, levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop"))) %>%
  group_by(scenario, collection, cutoff) %>%
  summarise(

    # Risks
    risk_ctrl_mean = mean(risk_ctrl, na.rm = TRUE),
    risk_ctrl_L    = mean(risk_ctrl_L, na.rm = TRUE),
    risk_ctrl_U    = mean(risk_ctrl_U, na.rm = TRUE),

    risk_trt_mean  = mean(risk_trt, na.rm = TRUE),
    risk_trt_L     = mean(risk_trt_L, na.rm = TRUE),
    risk_trt_U     = mean(risk_trt_U, na.rm = TRUE),

    # Risk Difference
    RD_mean = mean(RD, na.rm = TRUE),
    RD_L    = mean(RD_L, na.rm = TRUE),
    RD_U    = mean(RD_U, na.rm = TRUE),

    # Event counts
    events_ctrl = mean(events_ctrl, na.rm = TRUE),
    events_trt  = mean(events_trt,  na.rm = TRUE),

    # RR_mean = exp(mean(log(RR), na.rm = TRUE)),
    # RR_L    = exp(mean(log(RR_L), na.rm = TRUE)),
    # RR_U    = exp(mean(log(RR_U), na.rm = TRUE)),

    .groups = "drop"
  ) %>%
  mutate(
    line_id = factor(
      paste(scenario, collection, sep = " • "),
      levels = unique(paste(scenario, collection, sep = " • "))
  ))



##### Plot for absolute risks with CIs #####

risk_plot_df <- RD_summary %>%
  select(
    scenario, collection, cutoff, line_id,
    risk_ctrl_mean, risk_ctrl_L, risk_ctrl_U,
    risk_trt_mean,  risk_trt_L,  risk_trt_U
  ) %>%
  pivot_longer(
    cols = starts_with("risk_"),
    names_to = c("arm", "stat"),
    names_pattern = "risk_(ctrl|trt)_(mean|L|U)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(
    arm = recode(arm,
                 ctrl = "Control",
                 trt  = "Treatment")
  )


ggplot(
  risk_plot_df,
  aes(x = line_id, y = mean, fill = arm)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6) +
  geom_errorbar(
    aes(ymin = L, ymax = U),
    position = position_dodge(width = 0.7),
    width = 0.2,
    linewidth = 0.8) +
  geom_text(
    aes(
      label = percent(mean, accuracy = 1)
    ),
    position = position_dodge(width = 0.7),
    vjust = -2,
    size = 3.2,
    color = "black") +
  facet_grid(
    cutoff ~ .,
    labeller = labeller(cutoff = function(x) paste("Week", x))) +
  scale_fill_manual(
    values = c(
      "Control"   = "lightblue",
      "Treatment" = "darkblue"  )) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Absolute Risk of PRO Deterioration",
    subtitle = "Bars = mean risk, whiskers = 95% CI",
    x = "Scenario • Collection",
    y = "Risk",
    fill = "Arm") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"))


##### Plot for RD means with CI #####
RD_summary_reorder <- RD_summary %>%
  mutate(line_id = factor(line_id, levels = rev(unique(line_id))))


ggplot(
  RD_summary_reorder,
  aes(x = RD_mean, y = line_id, color = collection)) +
  geom_segment(
    aes(x = RD_L, xend = RD_U, yend = line_id),
    linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40") +
  facet_wrap(~ cutoff) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Mean risk difference (Treatment − Control) with 95% CI",
    x = "Risk difference",
    y = "Scenario • Collection",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9))


##### Bias for RD #####

# Define true RD = RD from the full scenario
RD_true <- RD_summary_reorder %>%
  filter(collection == "full") %>%
  select(scenario, cutoff, RD_true = RD_mean)

# Calculate RD bias for reduced and stop scenarios
RD_bias <- RD_summary_reorder %>%
  left_join(RD_true, by = c("scenario", "cutoff")) %>%
  filter(collection %in% c("reduced", "stop")) %>%
  mutate(bias = RD_mean - RD_true)

# Plot bias
ggplot(
  RD_bias,
  aes(x = bias, y = line_id, color = collection)) +
  geom_segment(
    aes(x = 0, xend = bias, yend = line_id),
    linewidth = 1.1) +
  geom_point(size = 3) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40") +
  facet_wrap(~ cutoff, scales = "fixed") +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Bias in mean Risk Difference",
    subtitle = "Bias = RD_mean(reduced/stop) − RD_mean(full)",
    x = "Bias in Risk Difference",
    y = "Scenario • Collection",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    strip.text = element_text(face = "bold"))


##### Power for RD #####
# Power = probability that the RD analysis rejects the null hypothesis RD=0
# Rejection rule: upper 95%CI is less than 0

power_RD <- RD_long %>%
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop")),
    cutoff     = factor(cutoff, levels = c(27, 36, 54)),

    reject = (RD_U < 0)) %>%
  group_by(scenario, collection, cutoff) %>%
  summarise(power = mean(reject, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(line_id = paste(scenario, collection, sep = " • "))



##### Coverage for RD #####
# Coverage = probability that the 95% CI for RD contains the true RD (FULL scenario)

covered_RD <- RD_long %>%
  left_join(RD_true, by = c("scenario", "cutoff")) %>%
  mutate(covered = (RD_L <= RD_true & RD_U >= RD_true))

coverage_RD <- covered_RD %>%
  mutate(scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop")),
    cutoff     = factor(cutoff, levels = c(27, 36, 54))) %>%
  group_by(scenario, collection, cutoff) %>%
  summarise(coverage = mean(covered, na.rm = TRUE),
    .groups  = "drop") %>%
  mutate(line_id = paste(scenario, collection, sep = " • "))

# Combined power and coverage
RD_perf <- power_RD %>%
  select(scenario, collection, cutoff, power) %>%
  left_join(coverage_RD %>%
      select(scenario, collection, cutoff, coverage),
    by = c("scenario", "collection", "cutoff")) %>%

  pivot_wider(
    names_from  = cutoff,
    values_from = c(power, coverage),
    names_glue  = "{.value}_{cutoff}"
  ) %>%
  arrange(
    factor(scenario, levels = c("BEFORE", "AT", "AFTER")),
    factor(collection, levels = c("full", "reduced", "stop"))
  )






#### Calculation for RR - not used ####

#
# ##### Plot for geometric mean RR with 95% CI #####
# ggplot(RR_summary,
#        aes(x = RR_mean, y = line_id, color = collection)) +
#   geom_segment(aes(x = RR_L, xend = RR_U, yend = line_id), linewidth = 1.2) +
#   geom_point(size = 3) +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
#   facet_wrap(~ cutoff) +
#   labs(
#     title = "Geometric mean RR with 95% CI",
#     x = "RR",
#     y = "Scenario",
#     color = "PRO data collection") +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     axis.text.y = element_text(size = 9))
#
#
# ##### Bias for RR #####
# # Define true RR = RR from the full scenario
# RR_true <- RR_summary %>%
#   filter(collection == "full") %>%
#   select(scenario, cutoff, RR_true = RR_mean)
#
# # Calculate bias
# RR_bias <- RR_summary %>%
#   left_join(RR_true, by = c("scenario", "cutoff")) %>%
#   filter(collection %in% c("reduced", "stop")) %>%
#   mutate(
#     bias = log(RR_mean) - log(RR_true)
#   )
#
# # Plot bias for RR
# ggplot(RR_bias,
#        aes(x = bias,
#            y = line_id,
#            color = collection)) +
#   geom_segment(aes(x = 0, xend = bias,
#                    yend = line_id),
#                linewidth = 1.1) +
#   geom_point(size = 3) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
#   facet_wrap(~ cutoff, scales = "fixed") +
#   labs(
#     title = "Bias in mean RR (log scale)",
#     subtitle = "Bias = log(RR_mean(reduced/stop)) − log(RR_mean(full))",
#     x = "Bias in mean RR",
#     y = "Scenario",
#     color = "PRO data collection") +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     axis.text.y = element_text(size = 10),
#     strip.text = element_text(face = "bold"))
#
# ##### Power for RR #####
# # Power = probability that the RR analysis rejects the null hypothesis RR=1
# # Rejection rule: upper 95%CI is less than 1, lower 95%CI is more than one
#
# power_RR <- RR_long %>%
#   mutate(
#     scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
#     collection = factor(collection, levels = c("full", "reduced", "stop")),
#     cutoff     = factor(cutoff, levels = c(27, 36, 54))
#   )%>%
#   mutate(
#     reject = (RR_U < 1)
#   ) %>%
#   group_by(scenario, collection, cutoff) %>%
#   summarise(
#     power = mean(reject, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     line_id = paste(scenario, collection, sep = " • ")
#   )
#
# # Plot for power
# ggplot(
#   power_RR,
#   aes(
#     x = power,
#     y = line_id,
#     color = collection
#   )
# ) +
#   geom_point(size = 3) +
#   facet_wrap(~ cutoff, scales = "fixed") +
#   geom_vline(
#     xintercept = 0.8,
#     linetype = "dashed",
#     color = "grey50"
#   ) +
#   labs(
#     title = "Power for RR",
#     subtitle = "Power = P(95% CI excludes RR = 1)",
#     x = "Power",
#     y = "Scenario • Collection",
#     color = "PRO data collection"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     axis.text.y = element_text(size = 9)
#   )
#
#
#
# ##### Coverage for RR #####
# # Coverage = probability that the 95% CI for RR contains the true RR (FULL scenario)
#
# # Does the CI for each replication cover true RR
# covered_RR <- RR_long %>%
#   left_join(RR_true, by = c("scenario", "cutoff")) %>%
#   mutate(
#     covered = (RR_L <= RR_true & RR_U >= RR_true)
#   )
#
# coverage_RR <- covered_RR %>%
#   mutate(
#     scenario   = factor(scenario, levels = c("BEFORE", "AT", "AFTER")),
#     collection = factor(collection, levels = c("full", "reduced", "stop")),
#     cutoff     = factor(
#       cutoff,
#       levels = c(27, 36, 54)
#     )
#   ) %>%
#   group_by(scenario, collection, cutoff) %>%
#   summarise(
#     coverage = mean(covered, na.rm = TRUE),
#     .groups  = "drop"
#   ) %>%
#   mutate(
#     line_id = paste(scenario, collection, sep = " • ")
#   )
#
#
# # Plot for coverage for RR
# ggplot(
#   coverage_RR,
#   aes(
#     x = coverage,
#     y = line_id,
#     color = collection
#   )
# ) +
#   geom_point(size = 3) +
#   facet_wrap(~ cutoff) +
#   geom_vline(
#     xintercept = 0.95,
#     linetype = "dashed",
#     color = "grey40"
#   ) +
#   labs(
#     title = "Coverage of RR 95% confidence intervals",
#     x = "Coverage probability",
#     y = "Scenario • Collection",
#     color = "PRO data collection"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     axis.text.y = element_text(size = 9)
#   )
