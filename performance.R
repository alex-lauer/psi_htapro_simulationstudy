# Performance Measures

library(dplyr)
library(ggplot2)
library(tidyr)
library(survminer)
library(stringr)

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
  mutate(line_id = paste(scenario, collection, sep = " • "))


##### Plot for geometric mean HR with 95%CI #####
ggplot(HR,
       aes(x = HR_mean,
           y = line_id,
           color = collection)) +
  geom_segment(aes(x = HR_CI_lower, xend = HR_CI_upper,
                   yend = line_id),
               linewidth = 1.1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  labs(
    title = "Geometric mean HR with 95% CI",
    x = "Hazard Ratio (HR)",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9))



##### Power for HR #####
# Power is defined as the probability of rejecting the null hypothesis of no treatment effect
# Rejection rule: 95%CI excludes 1.

power_HR <- analysis_results %>%
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop"))
  )%>%
  mutate(
    reject = HR_upper < 1
  ) %>%
  group_by(scenario, collection) %>%
  summarise(
    power = mean(reject),
    n_sims = n(),
    .groups = "drop"
  )


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
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop"))
  )%>%
  mutate(
    covered = HR_lower <= HR_full & HR_upper >= HR_full
  )

coverage_HR <- covered %>%
  group_by(scenario, collection) %>%
  summarise(
    coverage = mean(covered, na.rm = TRUE),
    n_sims   = n(),
    .groups  = "drop"
  )

# Combine power and coverage table

power_coverage_HR <- power_HR %>%
  left_join(
    coverage_HR,
    by = c("scenario", "collection")
  ) %>%
  select(-n_sims.x, -n_sims.y)


##### Bias in HR vs full scenario ####
# Evaluate reduced/stop vs full
bias_HR <- HR %>%
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
    legend.position = "top",
    axis.text.y = element_text(size = 10))


#### RR differences between scenarios ####
# Build long table for RR and 95% CIs at weeks 27, 36, and 54
RR_long <- analysis_results %>%
  select(
    scenario, collection, sim,
    RR27, RR27_L, RR27_U, events27_ctrl, events27_trt,
    RR36, RR36_L, RR36_U, events36_ctrl, events36_trt,
    RR54, RR54_L, RR54_U, events54_ctrl, events54_trt
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

##### Summary for RR #####
RR_summary <- RR_long %>%
  mutate(
    scenario   = factor(scenario, levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop"))) %>%
  group_by(scenario, collection, cutoff) %>%
  summarise(
    RR_mean = exp(mean(log(RR), na.rm = TRUE)),
    RR_L    = exp(mean(log(RR_L), na.rm = TRUE)),
    RR_U    = exp(mean(log(RR_U), na.rm = TRUE)),
    events_ctrl = mean(events_ctrl),
    events_trt  = mean(events_trt),
    .groups = "drop"
  ) %>%
  mutate(
    line_id = paste(scenario, collection, sep = " • ")
  )


##### Plot for geometric mean RR with 95% CI #####
ggplot(RR_summary,
       aes(x = RR_mean, y = line_id, color = collection)) +
  geom_segment(aes(x = RR_L, xend = RR_U, yend = line_id), linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  facet_wrap(~ cutoff) +
  labs(
    title = "Geometric mean RR with 95% CI",
    x = "RR",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9))


##### Bias for RR #####
# Define true RR = RR from the full scenario
RR_true <- RR_summary %>%
  filter(collection == "full") %>%
  select(scenario, cutoff, RR_true = RR_mean)

# Calculate bias
RR_bias <- RR_summary %>%
  left_join(RR_true, by = c("scenario", "cutoff")) %>%
  filter(collection %in% c("reduced", "stop")) %>%
  mutate(
    bias = log(RR_mean) - log(RR_true)
  )

# Plot bias for RR
ggplot(RR_bias,
       aes(x = bias,
           y = line_id,
           color = collection)) +
  geom_segment(aes(x = 0, xend = bias,
                   yend = line_id),
               linewidth = 1.1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ cutoff, scales = "fixed") +
  labs(
    title = "Bias in mean RR (log scale)",
    subtitle = "Bias = log(RR_mean(reduced/stop)) − log(RR_mean(full))",
    x = "Bias in mean RR",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    strip.text = element_text(face = "bold"))

##### Power for RR #####
# Power = probability that the RR analysis rejects the null hypothesis RR=1
# Rejection rule: upper 95%CI is less than 1, lower 95%CI is more than one

power_RR <- RR_long %>%
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop")),
    cutoff     = factor(cutoff, levels = c(27, 36, 54))
  )%>%
  mutate(
    reject = (RR_U < 1)
  ) %>%
  group_by(scenario, collection, cutoff) %>%
  summarise(
    power = mean(reject, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    line_id = paste(scenario, collection, sep = " • ")
  )

# Plot for power
ggplot(
  power_RR,
  aes(
    x = power,
    y = line_id,
    color = collection
  )
) +
  geom_point(size = 3) +
  facet_wrap(~ cutoff, scales = "fixed") +
  geom_vline(
    xintercept = 0.8,
    linetype = "dashed",
    color = "grey50"
  ) +
  labs(
    title = "Power for RR",
    subtitle = "Power = P(95% CI excludes RR = 1)",
    x = "Power",
    y = "Scenario • Collection",
    color = "PRO data collection"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9)
  )



##### Coverage for RR #####
# Coverage = probability that the 95% CI for RR contains the true RR (FULL scenario)

# Does the CI for each replication cover true RR
covered_RR <- RR_long %>%
  left_join(RR_true, by = c("scenario", "cutoff")) %>%
  mutate(
    covered = (RR_L <= RR_true & RR_U >= RR_true)
  )

coverage_RR <- covered_RR %>%
  mutate(
    scenario   = factor(scenario, levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop")),
    cutoff     = factor(
      cutoff,
      levels = c(27, 36, 54)
    )
  ) %>%
  group_by(scenario, collection, cutoff) %>%
  summarise(
    coverage = mean(covered, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    line_id = paste(scenario, collection, sep = " • ")
  )


# Plot for coverage for RR
ggplot(
  coverage_RR,
  aes(
    x = coverage,
    y = line_id,
    color = collection
  )
) +
  geom_point(size = 3) +
  facet_wrap(~ cutoff) +
  geom_vline(
    xintercept = 0.95,
    linetype = "dashed",
    color = "grey40"
  ) +
  labs(
    title = "Coverage of RR 95% confidence intervals",
    x = "Coverage probability",
    y = "Scenario • Collection",
    color = "PRO data collection"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9)
  )
