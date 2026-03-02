# Performance Measures

library(dplyr)
library(ggplot2)
library(tidyr)

#Load analysis results
analysis_results <- readRDS("analysis_results.rds")

#### HR across scenarios ####
analysis_results %>%
  group_by(scenario, collection) %>%
  summarise(HR_mean = mean(HR, na.rm = TRUE))

ggplot(analysis_results, aes(x = scenario, y = HR, fill = collection)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  labs(
    title = "Hazard Ratios Across Scenarios and Collection Designs",
    y = "Hazard Ratio (HR)",
    x = "Scenario"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  theme_bw()

#### Median differences across scenarios ####
analysis_results <- analysis_results %>%
  mutate(median_diff = median_trt - median_ctrl)

analysis_results %>%
  group_by(scenario, collection) %>%
  summarise(mean_median_diff = mean(median_diff, na.rm = TRUE))

#### RR differences between scenarios ####
analysis_results %>%
  group_by(scenario, collection) %>%
  summarise(
    RR27_mean = mean(RR27, na.rm = TRUE),
    RR36_mean = mean(RR36, na.rm = TRUE),
    RR54_mean = mean(RR54, na.rm = TRUE)
  )

# Long format
rr_long <- analysis_results %>%
  pivot_longer(
    cols = starts_with("RR"),
    names_to = "cutoff",
    values_to = "RR"
  ) %>%
  mutate(cutoff = factor(cutoff,
                         levels = c("RR27", "RR36", "RR54"),
                         labels = c("Week 27", "Week 36", "Week 54")))


ggplot(rr_long, aes(x = scenario, y = RR, fill = collection)) +
  geom_boxplot(position = position_dodge(width = 0.6), outlier.alpha = 0.35) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  facet_wrap(~ cutoff) +
  labs(
    title = "Risk Ratios by Scenario and Collection",
    x = "Scenario", y = "Risk Ratio (RR)") +
  theme_bw()






#### Bias relative to full design ####
bias_df <- analysis_results %>%
  group_by(scenario, sim) %>%
  mutate(
    HR_bias = HR - HR[collection == "full"],
    median_bias = median_diff - median_diff[collection == "full"],
    RR27_bias = RR27 - RR27[collection == "full"],
    RR36_bias = RR36 - RR36[collection == "full"],
    RR54_bias = RR54 - RR54[collection == "full"]
  ) %>%
  ungroup()

bias_summary <- bias_df %>%
  group_by(scenario, collection) %>%
  summarise(
    mean_HR_bias = mean(HR_bias, na.rm = TRUE),
    mean_median_bias = mean(median_bias, na.rm = TRUE),
    mean_RR27_bias = mean(RR27_bias, na.rm = TRUE),
    mean_RR36_bias = mean(RR36_bias, na.rm = TRUE),
    mean_RR54_bias = mean(RR54_bias, na.rm = TRUE),
    .groups = "drop"
  )

# Variability summary
sd_summary <- analysis_results %>%
  group_by(scenario, collection) %>%
  summarise(
    sd_HR = sd(HR, na.rm = TRUE),
    sd_median_diff = sd(median_diff, na.rm = TRUE),
    sd_RR27 = sd(RR27, na.rm = TRUE),
    sd_RR36 = sd(RR36, na.rm = TRUE),
    sd_RR54 = sd(RR54, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(bias_df, aes(x = scenario, y = HR_bias, fill = collection)) +
  geom_boxplot() +
  labs(title="HR Bias Across Scenarios", y="HR_bias", x="Scenario") +
  theme_bw()

rr_long_bias <- bias_df %>%
  pivot_longer(starts_with("RR"), names_to = "RR_type", values_to = "RR_value")

ggplot(rr_long_bias, aes(x = scenario, y = RR_value, fill = collection)) +
  geom_boxplot() +
  facet_wrap(~ RR_type, scales = "free_y") +
  theme_bw()
