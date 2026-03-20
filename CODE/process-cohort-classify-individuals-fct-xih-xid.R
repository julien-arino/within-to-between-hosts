#!/usr/bin/env Rscript
# ============================================================
# File: process-cohort-classify-individuals-fct-xih-xid.R
# Description:
#   File: process-cohort-classify-individuals-fct-xih-xid.R
#   This script processes the full virtual cohort simulation results.
#   It iterates over a grid of xi_h (hospitalization) and xi_d (death)
#   thresholds, similar to process-virtual-cohort-times-vary-xih-xid.jl.
#   For each combination, it identifies individuals based on their Psi_max:
#   - Mild: Psi_max < xi_h
#   - ICU: xi_h <= Psi_max < xi_d
#   - Dead: Psi_max >= xi_d
#   For individuals classified as "Dead", their simulation trajectories are
#   truncated exactly at the time they hit their maximum damage score (t_max).
#   The cleaned and truncated data is saved for each combination of cutpoints.
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("process-cohort-classify-individuals-fct-xih-xid.R")

load_libraries(c("dplyr", "qs2", "parallel", "data.table", "here"))

# Find the latest cohort_sim_parameters_P* file
files <- list.files(output_dir, pattern = "^cohort_sim_parameters_P.*\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No cohort_sim_parameters_P... file found in ", output_dir)
}

# Sort by modification time to get the newest one if multiple exist
latest_param_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort parameters:", basename(latest_param_file), "\n")

# Infer the state file
state_file <- file.path(output_dir, sub("^cohort_sim_parameters_", "cohort_sim_state_", basename(latest_param_file)))
if (!file.exists(state_file)) {
  stop("Matching state file not found: ", state_file)
}

# Instead of a full Cartesian grid (45 pairs), explicitly define the 11 pairs
# specifically requested by downstream plotting scripts (Fig 05a, 05b, 07a, etc.)
threshold_pairs <- list(
  c(50.0, 85.0),
  c(60.0, 85.0),
  c(70.0, 75.0),
  c(70.0, 80.0),
  c(70.0, 85.0),
  c(70.0, 90.0),
  c(70.0, 95.0),
  c(75.0, 85.0), # Default value pair
  c(75.0, 90.0),
  c(75.0, 95.0),
  c(80.0, 85.0)
)

# We still need unique xi_h_vals to compute the tau_h bounds optimally just once
xi_h_vals <- unique(sapply(threshold_pairs, `[`, 1))
xi_d_vals <- unique(sapply(threshold_pairs, `[`, 2))

# ====================================================================
# PART 1: PROCESS TABLES
# ====================================================================
cat("\n[PART 1] Loading newest cohort parameters and state...\n")
cohort_params <- qs_read(latest_param_file, nthreads = N_QS_THREADS)
cohort_state <- qs_read(state_file, nthreads = N_QS_THREADS)

base_name <- sub("^cohort_sim_parameters_", "cohort_", tools::file_path_sans_ext(basename(latest_param_file)))

base_summary_df <- cohort_params %>%
  select(
    ID,
    R0_within,
    max_V,
    tau_max_V,
    max_Psi,
    tau_max_Psi
  )

# Find crossover times for hospitalisation thresholds
cat("[PART 2] Calculating hospitalisation crossovers for all xi_h thresholds...\n")

n_cores <- parallel::detectCores()
workers_to_use <- min(24, max(2, n_cores - 2))

tau_h_list <- mclapply(cohort_state, function(ind_list) {
  time <- ind_list$time
  psi <- ind_list$Psi
  res_start <- numeric(length(xi_h_vals))
  res_end <- numeric(length(xi_h_vals))

  for (i in seq_along(xi_h_vals)) {
    xi <- xi_h_vals[i]
    idx_start <- which(psi >= xi)[1]
    res_start[i] <- if (is.na(idx_start)) NA_real_ else time[idx_start]

    if (is.na(idx_start)) {
      res_end[i] <- NA_real_
    } else {
      # Look only at times after onset
      post_onset_psi <- psi[idx_start:length(psi)]
      post_onset_time <- time[idx_start:length(time)]
      idx_end <- which(post_onset_psi < xi)[1]

      res_end[i] <- if (is.na(idx_end)) NA_real_ else post_onset_time[idx_end]
    }
  }
  list(start = res_start, end = res_end)
}, mc.cores = workers_to_use)

# Separate the start and end lists
tau_h_start_list <- lapply(tau_h_list, `[[`, "start")
tau_h_end_list <- lapply(tau_h_list, `[[`, "end")

# Extremely fast O(1) allocation: unspool single vector and reshape to matrix natively
tau_h_start_mat <- matrix(unlist(tau_h_start_list, use.names = FALSE), ncol = length(xi_h_vals), byrow = TRUE)
colnames(tau_h_start_mat) <- paste0("tau_h_start_", xi_h_vals)

tau_h_end_mat <- matrix(unlist(tau_h_end_list, use.names = FALSE), ncol = length(xi_h_vals), byrow = TRUE)
colnames(tau_h_end_mat) <- paste0("tau_h_end_", xi_h_vals)

base_summary_df <- bind_cols(
  base_summary_df,
  as.data.frame(tau_h_start_mat),
  as.data.frame(tau_h_end_mat)
)

# Immediately free the state list
rm(cohort_state, tau_h_list, tau_h_start_list, tau_h_end_list, tau_h_start_mat, tau_h_end_mat)
gc()

cat("Computing and saving status tables for strictly required thresholds...\n")
for (pair in threshold_pairs) {
  xi_h <- pair[1]
  xi_d <- pair[2]

  if (xi_h >= xi_d) {
    next
  }

  tau_h_start_col <- paste0("tau_h_start_", xi_h)
  tau_h_end_col <- paste0("tau_h_end_", xi_h)

  current_summary_df <- base_summary_df %>%
    mutate(
      status = case_when(
        max_Psi < xi_h ~ "Mild",
        max_Psi < xi_d ~ "ICU",
        TRUE ~ "Dead"
      ),
      tau_d = ifelse(status == "Dead", tau_max_Psi, NA_real_),
      tau_h_start = ifelse(status %in% c("ICU", "Dead"), .data[[tau_h_start_col]], NA_real_),
      tau_h_end = ifelse(status %in% c("ICU", "Dead"), .data[[tau_h_end_col]], NA_real_)
    ) %>%
    select(-any_of(c(paste0("tau_h_start_", xi_h_vals), paste0("tau_h_end_", xi_h_vals)))) # Clear numeric trackers

  h_str <- ifelse(xi_h %% 1 == 0, as.character(as.integer(xi_h)), as.character(xi_h))
  d_str <- ifelse(xi_d %% 1 == 0, as.character(as.integer(xi_d)), as.character(xi_d))

  new_base <- sub("cohort_", "cohort_status_", base_name)
  new_base <- paste0(new_base, "_xih_", h_str, "_xid_", d_str)

  out_file <- file.path(output_dir, paste0(new_base, ".qs"))
  qs_save(current_summary_df, out_file, nthreads = N_QS_THREADS)
  cat("  -> Saved:", basename(out_file), "\n")
}


# ====================================================================
# [PART 2: PROCESS TRAJECTORIES] -> REMOVED IN OPTIMIZATION REFACTOR
# We no longer manually truncate full array sequences to write gigabytes 
# of matrices back to disk. Instead, `process-cohort-times...` computes
# tau_c and tau_r purely effectively over the lightweight arrays natively.
# ====================================================================

cat("\nAll thresholds and status files processed successfully.\n")
print_end_time(start_time, "process-cohort-classify-individuals-fct-xih-xid.R")
