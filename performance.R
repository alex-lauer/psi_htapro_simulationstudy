# Performance Measures

library(dplyr)
library(ggplot2)
library(tidyr)

#Load analysis results
analysis_results <- readRDS("analysis_results.rds")

#### HR across scenarios ####
HR <- analysis_results %>%
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

            # Average hazards per arm across simulations
            hazard_ctrl_mean = mean(hazard_ctrl, na.rm = TRUE),
            hazard_trt_mean  = mean(hazard_trt,  na.rm = TRUE),

            .groups = "drop") %>%
  mutate(median_diff = median_trt_mean - median_ctrl_mean) %>%
  select(scenario, collection, HR_mean, p2.5, p97.5, hazard_ctrl_mean, hazard_trt_mean, median_ctrl_mean, median_trt_mean, median_diff)


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



# Compute bias vs full and containment of true value by percentile intervals
# Pull full mean HR per scenario
hr_full <- HR %>%
  filter(collection == "full") %>%
  transmute(scenario, HR_full = HR_mean)

# Evaluate reduced/stop vs full
hr_eval <- HR %>%
  left_join(hr_full, by = "scenario") %>%
  filter(collection %in% c("reduced", "stop")) %>%
  mutate(
    # Bias of HR means on a log scale
    bias = log(HR_mean) - log(HR_full),
    # Containment: does reduced/stop percentile interval contain the full mean HR?
    containment = HR_full >= p2.5 & HR_full <= p97.5)

# Plot bias for HR means
ggplot(hr_eval, aes(
    x = bias,
    y = paste(scenario, collection, sep = " • "),
    color = collection)) +
  geom_segment(aes(x = 0, xend = bias,
                   yend = paste(scenario, collection, sep = " • ")),
               linewidth = 1.1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title = "Bias in mean HR (log scale)",
    subtitle = "Bias = log(HR_mean(reduced/stop)) − log(HR_mean(full)); positive = inflated HR vs full",
    x = "Bias in mean HR",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10))

# Table for containment
hr_eval %>%
  select(scenario, collection, p2.5, p97.5, HR_full, containment)

## Boxplot for HR
# ggplot(analysis_results, aes(x = scenario, y = HR, fill = collection)) +
#   geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
#   labs(
#     title = "Hazard Ratios Across Scenarios and Collection Designs",
#     y = "Hazard Ratio (HR)",
#     x = "Scenario"
#   ) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
#   theme_bw()





#### RR differences between scenarios ####


# Build summary for RR at weeks 27, 36, and 54
RR <- analysis_results %>%
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

# Bias and containment for RR

# Pull FULL's geometric mean RR per scenario
rr_full <- RR %>%
  filter(collection == "full") %>%
  select(scenario, cutoff, RR_full = RR_mean)

# Calculate bias and containment
rr_eval <- RR %>%
  left_join(rr_full, by = c("scenario", "cutoff")) %>%
  filter(collection %in% c("reduced", "stop")) %>%
  mutate(
    # Bias of RR means on a log scale
    bias = log(RR_mean) - log(RR_full),
    # Containment: does reduced/stop percentile interval contain the full mean HR?
    containment = RR_full >= p2.5 & RR_full <= p97.5)


# Plot bias for RR

ggplot(rr_eval,
       aes(x = bias,
           y = paste(scenario, collection, sep = " • "),
           color = collection)) +
  geom_segment(aes(x = 0, xend = bias,
                   yend = paste(scenario, collection, sep = " • ")),
               linewidth = 1.1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ cutoff, scales = "fixed") +
  labs(
    title = "Bias in mean RR (log scale)",
    subtitle = "Bias = log(RR_mean(reduced/stop)) − log(RR_mean(full)); positive = inflated RR vs full",
    x = "Bias in mean RR (log scale)",
    y = "Scenario",
    color = "PRO data collection") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    strip.text = element_text(face = "bold"))


# Table for containment
rr_eval %>%
  select(scenario, collection, cutoff, p2.5, p97.5, RR_full, containment)

# # Boxplot for RR
# ggplot(rr_long, aes(x = scenario, y = RR, fill = collection)) +
#   geom_boxplot(position = position_dodge(width = 0.6), outlier.alpha = 0.35) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
#   facet_wrap(~ cutoff) +
#   labs(
#     title = "Risk Ratios by Scenario and Collection",
#     x = "Scenario", y = "Risk Ratio (RR)") +
#   theme_bw()
