#!/usr/bin/env Rscript
# ============================================================
# File: process-cohort-make-distributions.R
# Description:
#   File: process-cohort-make-distributions.R
#   This script extracts core statistical distributions (gamma, mu, beta)
#   from raw cohort status and simulation trajectory data.
#   By isolating these computations from the visualization files, we ensure
#   data processing logic is cleanly decoupled and standardly formatted
#   for downstream numerical models (like incidence plotting).
# ============================================================

cat("\n\n>>> Running process-cohort-make-distributions.R ...\n\n")
start_time <- Sys.time()
suppressWarnings(suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(future.apply)
  library(here)
}))

project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

## ------------------------------------------------------------
# 1. LOAD INTERPOLATED TRAJECTORIES AND BUILD V MATRIX
# ------------------------------------------------------------
cat("\nLoading pre-interpolated monolithic simulation trajectories...\n")
interp_files <- list.files(output_dir, pattern = "cohort_sim_state_interp_.*\\.qs$", full.names = TRUE)
if (length(interp_files) == 0) stop("No interpolated trajectory files found! Run process-cohort-interpolate-solutions.R first.")
latest_interp <- interp_files[which.max(file.mtime(interp_files))]

# Extract base prefix to dynamically name output files
base_prefix <- sub("cohort_sim_state_interp_([A-Za-z0-9_-]+)\\.qs", "\\1", basename(latest_interp))

# Read the massive flat dataframe
beta_df_raw <- qs_read(latest_interp, nthreads = N_QS_THREADS)

cat("Extracting V and converting to large matrix...\n")
# The interpolation guarantees t=0.0 to 100.0 by 0.1 for every individual
time_pts <- seq(0, 100.0, by = 0.1)
n_times <- length(time_pts)
N <- nrow(beta_df_raw) / n_times

# Because beta_df_raw is stacked by ID sequentially, matrix(..., ncol=N) correctly
# spools each sequential chunk of 1001 time points into a discrete column.
V_mat <- matrix(beta_df_raw$V, nrow = n_times, ncol = N)
cohort_ids <- unique(beta_df_raw$ID)

# Clear the massive dataframe from RAM
rm(beta_df_raw)
gc()

cat(sprintf("Constructed V matrix with %d rows (time) and %d columns (individuals).\n", nrow(V_mat), ncol(V_mat)))

# ============================================================
# 2. LOAD STATUS FILE AND COMPUTE ACTIVE INDIVIDUALS
# ============================================================
cat("\nLoading status metadata...\n")
status_pattern <- "^cohort_status_P.*_xih_75_xid_85_xic_4_xir_1\\.qs$"
status_files <- list.files(output_dir, pattern = status_pattern, full.names = TRUE)
if (length(status_files) == 0) stop("No cohort_status file found matching the required thresholds!")
target_status_file <- status_files[which.max(file.mtime(status_files))]

cat("Loading status file:", basename(target_status_file), "\n")
cohort_status <- qs_read(target_status_file, nthreads = N_QS_THREADS)

cat("Computing active individuals over time...\n")
# Active individuals: have not yet died or recovered
# Therefore, an individual is active at time t if: t < tau_d AND t < tau_r
# We can find the exit time for each individual
library(tidyr)
exit_times <- pmin(
  replace_na(cohort_status$tau_d, Inf),
  replace_na(cohort_status$tau_r, Inf)
)

# Convert to active counts across all time_pts efficiently
# sapply is perfectly fine since length(time_pts) is just 1001
active_counts <- sapply(time_pts, function(t_val) {
  sum(t_val < exit_times)
})

distributions_df <- data.frame(
  time = time_pts,
  active = active_counts,
  N_total = N
)
distributions_df$empirical_survival <- distributions_df$active / distributions_df$N_total

# ============================================================
# 3. COMPUTE MU (DEATH TIME DENSITIES)
# ============================================================
cat("\nComputing mu (death time) density probabilities...\n")

dead_cohort <- cohort_status %>% filter(status == "Dead")

if (nrow(dead_cohort) < 2) {
  cat("Not enough dead patients to compute density.\n")
  distributions_df$f_d <- 0
  distributions_df$mu_P <- 0
} else {
  # Compute KDE over the time grid we established
  # Evaluate exactly over the 0 to 100 domain in 0.1 steps to perfectly match our existing time_pts
  dens <- density(dead_cohort$tau_d, from = 0, to = max(time_pts), n = length(time_pts))

  # Proportion of deaths
  p_d <- nrow(dead_cohort) / N

  distributions_df$f_d <- dens$y
  distributions_df$mu_P <- (p_d * distributions_df$f_d) / distributions_df$empirical_survival
}

# ============================================================
# 4. COMPUTE GAMMA (RECOVERY TIME DENSITIES)
# ============================================================
cat("\nComputing gamma (recovery time) density probabilities...\n")

recovered_cohort <- cohort_status %>% filter(!is.na(tau_r))

if (nrow(recovered_cohort) < 2) {
  cat("Not enough recovered patients to compute density.\n")
  distributions_df$f_r <- 0
  distributions_df$gamma_P <- 0
} else {
  dens_r <- density(recovered_cohort$tau_r, from = 0, to = max(time_pts), n = length(time_pts))
  p_r <- nrow(recovered_cohort) / N

  distributions_df$f_r <- dens_r$y
  distributions_df$gamma_P <- (p_r * distributions_df$f_r) / distributions_df$empirical_survival
}

# ============================================================
# ============================================================
# 5. COMPUTE BETA (TRANSMISSION RATES IN TIME)
# ============================================================
cat("\nIsolating active transmission periods in V matrix...\n")

# Make a copy of V_mat to apply the specific status filters
V_mat_active <- V_mat

# Ensure status thresholds map exactly to the matrix columns by matching the IDs
tau_c_v <- cohort_status$tau_c[match(cohort_ids, cohort_status$ID)]
tau_r_v <- cohort_status$tau_r[match(cohort_ids, cohort_status$ID)]
tau_d_v <- cohort_status$tau_d[match(cohort_ids, cohort_status$ID)]
tau_h_start_v <- cohort_status$tau_h_start[match(cohort_ids, cohort_status$ID)]
tau_h_end_v <- cohort_status$tau_h_end[match(cohort_ids, cohort_status$ID)]

# Apply masks iteratively to save massive amounts of RAM
V_mat_active[outer(time_pts, replace_na(tau_c_v, Inf), "<")] <- 0
V_mat_active[outer(time_pts, replace_na(tau_r_v, -Inf), ">")] <- 0
V_mat_active[outer(time_pts, replace_na(tau_d_v, Inf), ">=")] <- 0
V_mat_active[outer(time_pts, replace_na(tau_h_start_v, Inf), ">=") &
  outer(time_pts, replace_na(tau_h_end_v, -Inf), "<=")] <- 0

cat("V matrix copied and masked! 0s successfully assigned to non-transmitting phases.\n")

cat("Mapping V to beta transmission scale...\n")
alpha_v <- 16.422
k_v <- 7.49
V_mat_active <- (V_mat_active^alpha_v) / (V_mat_active^alpha_v + k_v^alpha_v)
# Where V_mat_active was 0, beta_hat correctly maps to 0.

cat("Computing row-wise beta statistics (this may take a moment)...\n")
# Using apply over rows (1)
distributions_df$beta_mean <- apply(V_mat_active, 1, mean)
distributions_df$beta_Q10 <- apply(V_mat_active, 1, quantile, probs = 0.10, names = FALSE)
distributions_df$beta_median <- apply(V_mat_active, 1, median)
distributions_df$beta_Q90 <- apply(V_mat_active, 1, quantile, probs = 0.90, names = FALSE)

# Clean up memory safely!
rm(V_mat_active)
gc()

# ============================================================
# 6. SAVE FINAL DISTRIBUTIONS DATAFRAME
# ============================================================
cat("\nSaving final distributions dataframe...\n")

out_dist <- file.path(output_dir, paste0("cohort_distributions_", base_prefix, ".qs"))
qs_save(distributions_df, out_dist, nthreads = N_QS_THREADS)
cat("Final distributions saved to", basename(out_dist), "\n")
print_end_time(start_time, "process-cohort-make-distributions.R")
