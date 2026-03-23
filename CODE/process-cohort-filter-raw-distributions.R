#!/usr/bin/env Rscript
# ============================================================
# File: process-cohort-filter-raw-distributions.R
# Description:
#   Loads the lightweight raw distribution components (f_d, f_r, empirical_survival)
#   and applies math smoothing algorithms (splines, rolling averages) to create
#   a structured list of data frames ($raw, $spline, $rolling_average)
#   for seamless comparative plotting.
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("process-cohort-filter-raw-distributions.R")

load_libraries(c("qs2", "dplyr", "tidyr"))

# Helper for rolling average (1-day window by default for 0.1 step = 11 pts)
smooth_rolling <- function(x, k = 31) {
  if (length(x) < k) {
    return(x)
  }
  res <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  half_k <- floor(k / 2)
  # Patch boundaries with raw data where window overflows
  res[1:half_k] <- x[1:half_k]
  res[(length(x) - half_k + 1):length(x)] <- x[(length(x) - half_k + 1):length(x)]
  return(res)
}

cat("\n- Applying smoothing filters to raw distributions...\n")
raw_pattern <- "^cohort_distribution_raw_.*\\.qs$"
raw_files <- list.files(output_dir, pattern = raw_pattern, full.names = TRUE)

if (length(raw_files) == 0) stop("No raw distribution files found!")

for (file_idx in seq_along(raw_files)) {
  raw_file <- raw_files[file_idx]
  cat(sprintf("\n>>Processing raw subset %d/%d: %s\n", file_idx, length(raw_files), basename(raw_file)))

  base_prefix <- regmatches(basename(raw_file), regexpr("P[0-9]+_DT[0-9-]+", basename(raw_file)))
  threshold_suffix <- sub(paste0("^cohort_distribution_raw_", base_prefix, "(.*)\\.qs$"), "\\1", basename(raw_file))

  raw_df <- qs_read(raw_file, nthreads = N_QS_THREADS)

  status_file <- file.path(output_dir, paste0("cohort_status_", base_prefix, threshold_suffix, ".qs"))
  if (!file.exists(status_file)) stop(paste("Matching status file not found:", status_file))

  cohort_status <- qs_read(status_file, nthreads = N_QS_THREADS)

  N <- raw_df$N_total[1]
  dead_cohort <- cohort_status %>% filter(status == "Dead")
  recovered_cohort <- cohort_status %>% filter(!is.na(tau_r))

  p_d <- nrow(dead_cohort) / N
  p_r <- nrow(recovered_cohort) / N

  # Raw
  df_raw <- data.frame(
    time = raw_df$time,
    survival = raw_df$empirical_survival,
    f_d = raw_df$f_d,
    f_r = raw_df$f_r,
    beta_mean = raw_df$beta_mean,
    beta_mean_positive_transmission = raw_df$beta_mean_positive_transmission,
    beta_median_positive_transmission = raw_df$beta_median_positive_transmission,
    beta_Q10_positive_transmission = raw_df$beta_Q10_positive_transmission,
    beta_Q90_positive_transmission = raw_df$beta_Q90_positive_transmission,
    active = raw_df$active,
    N_total = raw_df$N_total
  )

  # Compute mu_P and gamma_P on raw
  df_raw$mu_P <- ifelse(df_raw$survival > 0, (p_d * df_raw$f_d) / df_raw$survival, 0)
  df_raw$gamma_P <- ifelse(df_raw$survival > 0, (p_r * df_raw$f_r) / df_raw$survival, 0)

  # Spline
  df_spline <- data.frame(
    time = raw_df$time,
    active = raw_df$active,
    N_total = raw_df$N_total
  )
  df_spline$survival <- smooth.spline(raw_df$time, raw_df$empirical_survival)$y
  df_spline$f_d <- smooth.spline(raw_df$time, raw_df$f_d)$y
  df_spline$f_r <- smooth.spline(raw_df$time, raw_df$f_r)$y
  df_spline$beta_mean <- smooth.spline(raw_df$time, raw_df$beta_mean)$y
  df_spline$beta_mean_positive_transmission <- smooth.spline(raw_df$time, raw_df$beta_mean_positive_transmission)$y
  df_spline$beta_median_positive_transmission <- smooth.spline(raw_df$time, raw_df$beta_median_positive_transmission)$y
  df_spline$beta_Q10_positive_transmission <- smooth.spline(raw_df$time, raw_df$beta_Q10_positive_transmission)$y
  df_spline$beta_Q90_positive_transmission <- smooth.spline(raw_df$time, raw_df$beta_Q90_positive_transmission)$y

  # Compute derived hazard rates from splined components
  df_spline$mu_P <- ifelse(df_spline$survival > 0, (p_d * df_spline$f_d) / df_spline$survival, 0)
  df_spline$gamma_P <- ifelse(df_spline$survival > 0, (p_r * df_spline$f_r) / df_spline$survival, 0)

  # Rolling average
  df_rolling <- data.frame(
    time = raw_df$time,
    active = raw_df$active,
    N_total = raw_df$N_total
  )
  df_rolling$survival <- smooth_rolling(raw_df$empirical_survival)
  df_rolling$f_d <- smooth_rolling(raw_df$f_d)
  df_rolling$f_r <- smooth_rolling(raw_df$f_r)
  df_rolling$beta_mean <- smooth_rolling(raw_df$beta_mean)
  df_rolling$beta_mean_positive_transmission <- smooth_rolling(raw_df$beta_mean_positive_transmission)
  df_rolling$beta_median_positive_transmission <- smooth_rolling(raw_df$beta_median_positive_transmission)
  df_rolling$beta_Q10_positive_transmission <- smooth_rolling(raw_df$beta_Q10_positive_transmission)
  df_rolling$beta_Q90_positive_transmission <- smooth_rolling(raw_df$beta_Q90_positive_transmission)

  # Compute rolling average of the pure empirical hazard rates directly
  safe_roll <- function(x) {
    idx <- is.finite(x)
    y <- rep(0, length(x))
    if (sum(idx) > 0) y[idx] <- smooth_rolling(x[idx])
    return(y)
  }
  df_rolling$mu_P <- safe_roll(df_raw$mu_P)
  df_rolling$gamma_P <- safe_roll(df_raw$gamma_P)

  # Derived comp using rolling average
  df_comp_rolling <- df_rolling
  df_comp_rolling$mu_P <- ifelse(df_comp_rolling$survival > 0, (p_d * df_comp_rolling$f_d) / df_comp_rolling$survival, 0)
  df_comp_rolling$gamma_P <- ifelse(df_comp_rolling$survival > 0, (p_r * df_comp_rolling$f_r) / df_comp_rolling$survival, 0)

  # Package and save
  filters_list <- list(
    raw = df_raw,
    spline = df_spline,
    rolling_average = df_rolling,
    comp_using_rolling_avg = df_comp_rolling
  )

  # Save filter distributions list
  cat("- Saving filtered distributions list containing $raw, $spline, $rolling_average, and $comp_using_rolling_avg...\n")
  out_dist <- file.path(output_dir, paste0("cohort_distribution_filters_", base_prefix, threshold_suffix, ".qs"))
  qs_save(filters_list, out_dist, nthreads = N_QS_THREADS)
  cat("Filter distributions saved to", basename(out_dist), "\n")
}

print_end_time(start_time, "process-cohort-filter-raw-distributions.R")
