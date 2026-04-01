library(ggplot2)

# Visualization of PRO trajectory
# PRO value (y) vs time (x) for each patient
pro_over_time <- ggplot(pro_long,
  aes(x = week, y = PRO, group = patient, color = arm)) +
  geom_line(alpha = 0.15, linewidth = 0.4) +          # individual patient lines
  geom_point(alpha = 0.15, size = 0.6) +               # points at visits
  facet_wrap(~ arm, ncol = 2) +                        # separate panels per arm
  scale_x_continuous(breaks = c(0, unique(pro_long$week))) +
  labs(
    title = "Individual Patient PRO from Baseline to Each Time Point",
    x = "Week",
    y = "PRO") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(pro_over_time)


# What is the mean of CfB at week 36 and 54?
df_36 <- pro_wide %>%
  mutate(change_36 = week36 - PRO_baseline) %>%
  group_by(arm) %>%
  summarise(
    mean_change_36 = mean(change_36, na.rm = TRUE),
    sd_change_36   = sd(change_36, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

df_36

df_54 <- pro_wide %>%
  mutate(change_54 = week54 - PRO_baseline) %>%
  group_by(arm) %>%
  summarise(
    mean_change_54 = mean(change_54, na.rm = TRUE),
    sd_change_54   = sd(change_54, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

df_54
