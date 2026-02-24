library(mvtnorm)
library(truncnorm)
library(TrialSimulator)
library(tidyr)
library(dplyr)
library(ggplot2)

#set.seed(123)

# Create vector for measurement points
weeks <- c(seq(3, 27, 3), seq(36, 153, 9))

# Create a covariance matrix for errors for every visit
rho <- 0.5
sigma <- 3.5
pro_vcov <- diag(sigma^2, nrow = 23, ncol = 23)
pro_vcov[pro_vcov == 0] <- rho*sigma^2

# Draw errors for visits
errors <- rmvnorm(600, mean = rep(0, nrow(pro_vcov)), sigma = pro_vcov) |> as.data.frame()

# sd(as.numeric(errors[,1]))
# cor(errors)

# PRO at baseline
pro_mean_bl <- 70
pro_sigma_bl <- 20

pro_bl <- rtruncnorm(
  n = 600,
  mean = pro_mean_bl,
  sd = pro_sigma_bl,
  a = 0,       # lower bound
  b = 100      # upper bound
) %>%
  as.data.frame()

names(pro_bl) <- "PRO_baseline"
# hist(pro_bl$PRO_baseline)

# Randomize patients to treatment and control 1:1
pro_bl$arm <- sample(rep(c("control", "treatment"), each = 300))
pro_bl$TRT <- as.integer(pro_bl$arm == "treatment")

# Parameters for scenarios
scenario_table <- readRDS("scenario_table.rds")
beta1 <- as.numeric(scenario_table[2, "beta1"])
beta3 <- as.numeric(scenario_table[2, "beta3"])
beta4 <- 0.01


# Generator (vectorized): expected mean PRO at time t
pro_generator <- function(pro_bl_val, t, trt) {
  n <- length(pro_bl_val)
  pro_t <- pro_bl_val + (beta1 * t) + (beta3 * t * trt) + (beta4 * pro_bl_val) # + rnorm(n, mean = 0, sd = 3)
  return(pro_t)
}

# Build a wide data frame
pro_wide <- pro_bl

for (wk in weeks) {
  pro_wide[[paste0("week", wk)]] <- pro_generator(
    pro_bl_val = pro_wide$PRO_baseline,
    t          = wk,
    trt        = pro_wide$TRT)}

# Add correlated errors to generated values of PROs
pro_wide_err <- pro_wide[, 4:ncol(pro_wide)] + errors

# Clamp to [0, 100]
pro_wide_err[pro_wide_err < 0] <- 0
pro_wide_err[pro_wide_err > 100] <- 100

# Recombine with patient metadata
pro_wide <- cbind(pro_wide[, 1:3], pro_wide_err)

# Add week0
pro_wide$week0 <- pro_wide$PRO_baseline

# Add patient column
pro_wide <- pro_wide %>%
  mutate(patient = dplyr::row_number())

# Identify visit columns
visits <- c("week0", paste0("week", weeks))

# Reorder
pro_wide <- pro_wide %>%
  select(patient, arm, TRT, PRO_baseline, all_of(visits))

# Build a long dataset: one row per patient-week
pro_long <- pro_wide %>%
  pivot_longer(
    cols = all_of(visits),
    names_to = "visit",
    values_to = "PRO"
  ) %>%
  mutate(
    week = as.integer(sub("week", "", visit))  # extract numeric week
  ) %>%
  select(-c(visit))

# Compute change from baseline (Delta PRO)
pro_long <- pro_long %>%
  mutate(PRO_change = PRO - PRO_baseline)


# Visulaisation
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
