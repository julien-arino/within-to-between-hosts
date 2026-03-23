#!/usr/bin/env Rscript
# ============================================================
# File: process-cohort-make-raw-distributions.R
# Description:
#   Extracts raw statistical data (empirical survival, f_d, f_r, betas)
#   from massive simulation trajectories into lightweight intermediate files.
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("process-cohort-make-raw-distributions.R")

load_libraries(c("qs2", "dplyr"))

# Load interpolated trajectories and build V matrix
cat("\n- Loading interpolated simulation trajectories...\n")
interp_files <- list.files(output_dir, pattern = "cohort_sim_state_interp_.*\\.qs$", full.names = TRUE)
if (length(interp_files) == 0) stop("No interpolated trajectory files found! Run process-cohort-interpolate-solutions.R first.")
latest_interp <- interp_files[which.max(file.mtime(interp_files))]

base_prefix <- sub("cohort_sim_state_interp_([A-Za-z0-9_-]+)\\.qs", "\\1", basename(latest_interp))

beta_df_raw <- qs_read(latest_interp, nthreads = N_QS_THREADS)

cat("- Extracting V and converting to large matrix...\n\n")
time_pts <- seq(0, 100.0, by = 0.1)
n_times <- length(time_pts)
N <- nrow(beta_df_raw) / n_times

V_mat <- matrix(beta_df_raw$V, nrow = n_times, ncol = N)
cohort_ids <- unique(beta_df_raw$ID)

rm(beta_df_raw)
gc()

cat(sprintf("- Constructed V matrix with %d rows (time) and %d columns (individuals).\n", nrow(V_mat), ncol(V_mat)))

# Load status file and compute raw metrics
cat("\n- Executing batch raw distribution pipeline for all available status thresholds...\n")
status_pattern <- paste0("^cohort_status_", base_prefix, "_xih_75_xid_85_xic_[0-9]+_xir_[0-9]+\\.qs$")
target_status_files <- list.files(output_dir, pattern = status_pattern, full.names = TRUE)

if (length(target_status_files) == 0) stop("No cohort_status files found matching the required thresholds!")

for (file_idx in seq_along(target_status_files)) {
  target_status_file <- target_status_files[file_idx]
  cat(sprintf("\n>> Processing status file %d/%d: %s\n", file_idx, length(target_status_files), basename(target_status_file)))

  threshold_suffix <- sub(paste0("^cohort_status_", base_prefix, "(.*)\\.qs$"), "\\1", basename(target_status_file))
  cohort_status <- qs_read(target_status_file, nthreads = N_QS_THREADS)

  cat("- Computing active individuals over time...\n")
  exit_times <- pmin(
    coalesce(cohort_status$tau_d, Inf),
    coalesce(cohort_status$tau_r, Inf)
  )

  active_counts <- sapply(time_pts, function(t_val) {
    sum(t_val < exit_times)
  })

  distributions_df <- data.frame(
    time = time_pts,
    active = active_counts,
    N_total = N
  )
  distributions_df$empirical_survival <- distributions_df$active / distributions_df$N_total

  # Compute f_d
  cat("- Computing f_d (death time) density probabilities...\n")
  dead_cohort <- cohort_status %>% filter(status == "Dead")
  if (nrow(dead_cohort) < 2) {
    distributions_df$f_d <- 0
  } else {
    dens <- density(dead_cohort$tau_d, from = 0, to = max(time_pts), n = length(time_pts))
    distributions_df$f_d <- dens$y
  }

  # Compute f_r
  cat("- Computing f_r (recovery time) density probabilities...\n")
  recovered_cohort <- cohort_status %>% filter(!is.na(tau_r))
  if (nrow(recovered_cohort) < 2) {
    distributions_df$f_r <- 0
  } else {
    dens_r <- density(recovered_cohort$tau_r, from = 0, to = max(time_pts), n = length(time_pts))
    distributions_df$f_r <- dens_r$y
  }

  # Compute beta
  cat("- Isolating active transmission periods in V matrix...\n")
  V_mat_active <- V_mat

  tau_c_v <- cohort_status$tau_c[match(cohort_ids, cohort_status$ID)]
  tau_r_v <- cohort_status$tau_r[match(cohort_ids, cohort_status$ID)]
  tau_d_v <- cohort_status$tau_d[match(cohort_ids, cohort_status$ID)]
  tau_h_start_v <- cohort_status$tau_h_start[match(cohort_ids, cohort_status$ID)]
  tau_h_end_v <- cohort_status$tau_h_end[match(cohort_ids, cohort_status$ID)]

  V_mat_active[outer(time_pts, coalesce(tau_c_v, Inf), "<")] <- 0
  V_mat_active[outer(time_pts, coalesce(tau_r_v, Inf), ">")] <- 0
  V_mat_active[outer(time_pts, coalesce(tau_d_v, Inf), ">=")] <- 0
  V_mat_active[outer(time_pts, coalesce(tau_h_start_v, Inf), ">=") &
    outer(time_pts, coalesce(tau_h_end_v, Inf), "<=")] <- 0

  cat("- V matrix copied and masked! 0s successfully assigned to non-transmitting phases.\n")

  cat("- Mapping V to beta(V) = V^alpha / (V^alpha + k^alpha)...\n")
  alpha_v <- 16.422
  k_v <- 7.49
  V_mat_active <- (V_mat_active^alpha_v) / (V_mat_active^alpha_v + k_v^alpha_v)

  cat("- Computing row-wise beta statistics (this may take a moment)...\n")
  distributions_df$beta_mean <- apply(V_mat_active, 1, mean)

  calc_quant <- function(x, p) {
    active_x <- x[x > 0]
    if (length(active_x) == 0) {
      return(0)
    }
    quantile(active_x, probs = p, names = FALSE)
  }

  distributions_df$beta_mean_positive_transmission <- apply(V_mat_active, 1, function(x) {
    active_x <- x[x > 0]
    if (length(active_x) == 0) {
      return(0)
    }
    mean(active_x)
  })
  distributions_df$beta_Q10_positive_transmission <- apply(V_mat_active, 1, calc_quant, p = 0.10)
  distributions_df$beta_median_positive_transmission <- apply(V_mat_active, 1, calc_quant, p = 0.50)
  distributions_df$beta_Q90_positive_transmission <- apply(V_mat_active, 1, calc_quant, p = 0.90)

  rm(V_mat_active)
  gc()

  # Save raw distributions dataframe
  cat("- Saving raw distributions dataframe...\n")
  out_dist <- file.path(output_dir, paste0("cohort_distribution_raw_", base_prefix, threshold_suffix, ".qs"))
  qs_save(distributions_df, out_dist, nthreads = N_QS_THREADS)
  cat("Raw distributions saved to", basename(out_dist), "\n")
}

print_end_time(start_time, "process-cohort-make-raw-distributions.R")
