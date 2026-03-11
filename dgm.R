library(mvtnorm)
library(truncnorm)
library(tidyr)
library(dplyr)
library(purrr)

# Parameters
n_patients   <- 600         # number of patients per simulation
n_sims       <- 200         # simulations per scenario
weeks        <- c(seq(3, 27, 3), seq(36, 153, 9))  # 23 post-baseline visits
rho          <- 0.5         # correlation between visits
sigma        <- 3.5         # SD of visit-level errors
pro_mean_bl  <- 70          # baseline mean
pro_sigma_bl <- 20          # baseline SD
beta4        <- 0.01        # effect of baseline PRO on trajectory (constant across scenarios)
scenario_names <- c("BEFORE", "AT", "AFTER")

# Read scenario parameters
scenario_table <- readRDS("scenario_table.rds")

# Create a covariance matrix for errors for every visit
pro_vcov <- diag(sigma^2, nrow = 23, ncol = 23)
pro_vcov[pro_vcov == 0] <- rho*sigma^2

# Generator for expected PRO at time t
pro_generator <- function(pro_bl_val, t, trt, beta1, beta3, beta4) {
  pro_bl_val + (beta1 * t) + (beta3 * t * trt) + (beta4 * pro_bl_val)
}

# Helper to simulate PD and convert to the nearest visit week

# Exponential PFS rates from medians: lambda = log(2)/median
lambda_control   <- log(2) / 36  # control median = 36 weeks
lambda_treatment <- log(2) / 54  # treatment median = 54 weeks

compute_PD <- function(df_patients) {
  df_patients %>%
    mutate(
      PD_time = ifelse(
        arm == "control",
        rexp(n(), rate = lambda_control),
        rexp(n(), rate = lambda_treatment)
      ),
      PD_event_week = sapply(PD_time, function(t) {
        w <- weeks[weeks >= t]
        if (length(w) == 0) NA_integer_ else w[1]
      }),
      PD_status = ifelse(!is.na(PD_event_week), 1L, 0L)  # 1 = event
    ) %>%
    select(patient, arm, PD_time, PD_event_week, PD_status)
}

# Helper: PRO event detection

compute_PRO_event <- function(long_df) {
  long_df %>%
    group_by(patient, arm) %>%
    summarise(
      PRO_status = as.integer(any(week > 0 & PRO_change <= -10)),
      PRO_event_week = ifelse(
        PRO_status == 1,
        week[which(week > 0 & PRO_change <= -10)[1]],
        NA_integer_),
      .groups = "drop")
}

# Data collection scenarios

collection_scenarios <- function(long_df, pro_ev_full, pd_tbl) {

  # Full data collection until the end of study
  full <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_ev_full, by = "patient") %>%
    left_join(pd_tbl, by = "patient") %>%
    mutate(collection = "full") %>%
    select(patient, arm, PRO_status, PRO_event_week,
           PD_status, PD_event_week, collection)

  # STOP scenario
  # For patients who didn't have PD event - keep all weeks.
  # For patients who had a PD event - keep weeks up until PD event.
  long_stop <- long_df %>%
    left_join(pd_tbl, by = c("patient", "arm")) %>%
    filter(is.na(PD_event_week) | week <= PD_event_week)

  # Compute PRO events until PD event
  pro_ev_stop <- compute_PRO_event(long_stop)

  stop <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_ev_stop, by = "patient") %>%
    left_join(pd_tbl,      by = "patient") %>%
    mutate(collection = "stop") %>%
    select(patient, arm, PRO_status, PRO_event_week,
           PD_status, PD_event_week, collection)

  # REDUCED scenario

  # Collect PRO at every 3rd visit *after* PD
  # Keep all visits up to and including PD_event_week
  # After PD_event_week, keep only visits whose index distance from PD visit is a multiple of 3
  long_reduced <- long_df %>%
    left_join(pd_tbl, by = c("patient", "arm")) %>%
    group_by(patient) %>%
    mutate(
      # Visit index per patient: 0 = baseline (week0), 1 = first post-baseline, etc.
      visit_idx = dplyr::row_number() - 1L,

      # Index of the PD detection visit on the schedule; NA if no PD
      pd_idx = ifelse(
        PD_status == 0,
        NA_integer_,
        # first index where week equals PD_event_week
        (which(week == PD_event_week) - 1L)),

      # Keep rule:
      #  - If no PD, keep everything
      #  - If PD occurred, keep all rows up to PD_event_week (inclusive),
      #    and after PD keep only every 3rd visit (visit_idx - pd_idx in {3,6,9,...})
      keep = is.na(pd_idx) |
        (visit_idx <= pd_idx) |
        (visit_idx > pd_idx & ((visit_idx - pd_idx) %% 3L == 0L))
    ) %>%
    ungroup() %>%
    filter(keep)%>%
    select(-visit_idx, -pd_idx, -keep)

  # Recompute PRO event on the thinned schedule (status = 1 event, 0 no event)
  pro_reduced_ev <- compute_PRO_event(long_reduced)

  # Patient-level 'reduced' dataset (merge PRO and PD event variables)
  reduced <- long_df %>%
    distinct(patient, arm) %>%
    left_join(pro_reduced_ev, by = c("patient", "arm")) %>%
    left_join(pd_tbl, by = c("patient", "arm")) %>%
    transmute(
      patient, arm,
      PRO_status, PRO_event_week,
      PD_status,  PD_event_week,
      collection = "reduced")

  return(list(
    full    = full,
    stop    = stop,
    reduced = reduced
  ))
}


# Simulation function

simulate <- function(beta1, beta3, seed) {
  set.seed(seed)

  # PRO at baseline
  pro_bl <- rtruncnorm(
    n = n_patients,
    mean = pro_mean_bl,
    sd = pro_sigma_bl,
    a = 0,
    b = 100) %>%
    as.data.frame()
  names(pro_bl) <- "PRO_baseline"

  # Randomize patients to treatment and control 1:1
  pro_bl$arm <- sample(rep(c("control", "treatment"), each = 300))
  pro_bl$TRT <- as.integer(pro_bl$arm == "treatment")
  pro_bl$patient <- 1:n_patients

  # Simulate PD
  pd_tbl <- compute_PD(pro_bl)

  # Generate expected PRO values for each visit (without error)
  # Build a wide data frame
  pro_wide <- pro_bl
  for (wk in weeks) {
    pro_wide[[paste0("week", wk)]] <- pro_generator(
      pro_bl_val = pro_wide$PRO_baseline,
      t          = wk,
      trt        = pro_wide$TRT,
      beta1      = beta1,
      beta3      = beta3,
      beta4      = beta4
    )}

  # Draw errors for visits
  errors <- mvtnorm::rmvnorm(
    n = n_patients,
    mean = rep(0, nrow(pro_vcov)),
    sigma = pro_vcov) %>%
    as.data.frame()

  # Add correlated errors to generated values of PROs
  pro_wide_err <- pro_wide[, 5:ncol(pro_wide)] + errors

  # Clamp to [0, 100]
  pro_wide_err[pro_wide_err < 0] <- 0
  pro_wide_err[pro_wide_err > 100] <- 100

  # Recombine with patient metadata
  pro_wide <- cbind(pro_wide[, 1:4], pro_wide_err)

  # Add week0
  pro_wide$week0 <- pro_wide$PRO_baseline

  # Identify visit columns
  visits <- c("week0", paste0("week", weeks))

  # Reorder
  pro_wide <- pro_wide %>%
    dplyr::select(patient, arm, TRT, PRO_baseline, all_of(visits))

  # Build a long dataset: one row per patient-week
  pro_long <- pro_wide %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(visits),
      names_to = "visit",
      values_to = "PRO"
    ) %>%
    dplyr::mutate(
      week = as.integer(sub("week", "", visit))
    ) %>%
    dplyr::select(-c(visit))

  # Compute change from baseline (Delta PRO)
  pro_long <- pro_long %>%
    dplyr::mutate(PRO_change = PRO - PRO_baseline)

  # PRO event for FULL scenario
  pro_ev_full <- compute_PRO_event(pro_long)

  # Build data-collection scenario versions
  collections <- collection_scenarios(long_df = pro_long,
                                      pro_ev_full = pro_ev_full,
                                      pd_tbl = pd_tbl)

  return(collections)
}


# Run 200 simulations for each time scenario (before, at, after)
sim_out <- purrr::map(
  scenario_names,
  function(scn) {
    # Get betas for this scenario
    row   <- scenario_table %>% dplyr::filter(.data$scenario == scn)
    beta1 <- as.numeric(row$beta1[[1]])
    beta3 <- as.numeric(row$beta3[[1]])

    # Run n_sims simulations with different seeds
    purrr::map(
      seq_len(n_sims),
      ~ simulate(
        beta1 = beta1,
        beta3 = beta3,
        seed  = 1000 + .x
      )
    )
  }
)

names(sim_out) <- scenario_names

saveRDS(sim_out, "simulated_pro_datasets.rds")
