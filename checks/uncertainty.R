# Compare CIs empirical (2.5-97.5 percentile across 200 simulations) vs from the Cox PH model

HR_cox_CI <- analysis_results %>%
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop"))
  ) %>%
  group_by(scenario, collection) %>%
  summarise(
    # Geometric mean HR
    HR_mean = exp(mean(log(HR), na.rm = TRUE)),

    # Cox PH Wald CI aggregated across simulations
    HR_cox_L = exp(mean(log(HR_lower), na.rm = TRUE)),
    HR_cox_U = exp(mean(log(HR_upper), na.rm = TRUE)),

    .groups = "drop"
  ) %>%
  mutate(line_id = paste(scenario, collection, sep = " • "))

HR_empirical_CI <- analysis_results %>%
  mutate(
    scenario   = factor(scenario,   levels = c("BEFORE", "AT", "AFTER")),
    collection = factor(collection, levels = c("full", "reduced", "stop"))
  ) %>%
  group_by(scenario, collection) %>%
  summarise(
    HR_mean = exp(mean(log(HR), na.rm = TRUE)),
    HR_emp_L = quantile(HR, 0.025, na.rm = TRUE),
    HR_emp_U = quantile(HR, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(line_id = paste(scenario, collection, sep = " • "))


HR_compare <- HR_empirical_CI %>%
  left_join(HR_cox_CI,
            by = c("scenario", "collection", "HR_mean", "line_id"))

ggplot(HR_compare,
       aes(y = line_id)) +

  # Empirical CI (thicker, outer)
  geom_segment(
    aes(x = HR_emp_L, xend = HR_emp_U,
        yend = line_id),
    linewidth = 1.3,
    color = "steelblue"
  ) +

  # Cox PH CI (thinner, inner)
  geom_segment(
    aes(x = HR_cox_L, xend = HR_cox_U,
        yend = line_id),
    linewidth = 0.6,
    color = "darkorange"
  ) +

  # HR point
  geom_point(
    aes(x = HR_mean),
    size = 3,
    color = "black"
  ) +

  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey40") +

  labs(
    title = "HR uncertainty: empirical vs Cox PH",
    subtitle = "Blue = empirical (between simulations), Orange = Cox PH (within model)",
    x = "Hazard Ratio (HR)",
    y = "Scenario • Collection"
  ) +

  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9)
  )
