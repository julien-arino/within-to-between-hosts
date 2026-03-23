#!/usr/bin/env Rscript
# ============================================================
# File: process-cohort-times-fct-xih-xid-xic-xir.R
# Description:
#   File: process-cohort-times-fct-xic-xir.R
#   This script complements zero-transmission grouping by specifically
#   identifying the exact simulation times when individuals become
#   infectious (tau_c) and when they recover/cease to be infectious (tau_r).
#   It reads the large `cohort_truncated_state` files to find V(t) crossings
#   for the given xi_c and xi_r thresholds, and appends these extracted
#   times to matching `cohort_status_` files, saving them with the
#   corresponding `_xic_X_xir_X.qs` suffixes for downstream plotting.
# ============================================================

# First things first: locate project directory and load helper functions
# and constants
project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

# Say what we are running and start the clock
start_time <- start_time_and_hello("process-cohort-times-fct-xih-xid-xic-xir.R")

# Load libraries
load_libraries(c("dplyr", "qs2", "parallel"))

# Define explicit combinations of infectiousness thresholds (xi_c, xi_r)
xi_pairs <- list(
  c(2, 1),
  c(2, 2),
  c(2, 4),
  c(4, 1),
  c(4, 2),
  c(4, 4),
  c(6, 1),
  c(6, 2),
  c(6, 4)
)

files <- list.files(output_dir, pattern = "^cohort_sim_state_P.*\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No cohort_sim_state_P... file found in ", output_dir)
}

# We process the base native simulation states, avoiding derived caches
files <- files[!grepl("_zerotrans\\.qs$", files) & !grepl("_interp", files)]

if (length(files) == 0) {
  cat("All files have already been processed.\n")
  quit(save = "no")
}

# Keep only the latest file based on the DT timestamp in the file name
dt_tags <- sapply(files, function(f) {
  m <- regmatches(basename(f), regexpr("DT[0-9]{8}-[0-9]{6}", basename(f)))
  if (length(m) > 0) m[[1]] else ""
})
if (any(dt_tags != "")) {
  files <- files[dt_tags == max(dt_tags)]
}

for (file_idx in seq_along(files)) {
  current_file <- files[file_idx]
  cat(sprintf(">>> [%d/%d] Extracting event times from: %s\n", file_idx, length(files), basename(current_file)))

  cohort_state <- qs_read(current_file, nthreads = N_QS_THREADS)

  for (pair in xi_pairs) {
    c_val <- pair[1]
    r_val <- pair[2]
    cat(sprintf("Computing tau_c/tau_r for [xi_c=%g, xi_r=%g]...\n", c_val, r_val))

    n_cores <- parallel::detectCores()
    # Using 24 cores to avoid Copy-On-Write (COW) memory overflow
    workers_to_use <- min(24, max(2, n_cores - 2))

    matches <- regmatches(basename(current_file), regexpr("P[0-9]+_DT[0-9]+-[0-9]+", basename(current_file)))
    if (length(matches) == 0) stop("Could not extract P/DT timestamp from file name.")
    timestamp_part <- matches[1]

    status_pattern <- paste0("^cohort_status_", timestamp_part, ".*\\.qs$")
    status_files <- list.files(output_dir, pattern = status_pattern, full.names = TRUE)
    status_files <- status_files[!grepl("_xic_", status_files) & !grepl("_zerotrans", status_files)]

    tau_cr_list <- mclapply(cohort_state, function(df) {
      # Find peak index
      V_max_idx <- which.max(df$V)

      # tau_c
      idx_c <- which(df$V >= c_val)[1]
      tc <- if (is.na(idx_c)) NA_real_ else df$time[idx_c]

      # tau_r
      tr <- NA_real_
      if (V_max_idx < nrow(df)) {
        post_peak_V <- df$V[(V_max_idx + 1):nrow(df)]
        r_offset <- which(post_peak_V <= r_val)[1]
        if (!is.na(r_offset)) {
          tr <- df$time[V_max_idx + r_offset]
        }
      }
      c(tau_c = tc, tau_r = tr)
    }, mc.cores = workers_to_use)

    # O(1) instantaneous memory binder (C-level vector unspooling)
    tau_cr_df <- as.data.frame(matrix(unlist(tau_cr_list, use.names = FALSE), ncol = 2, byrow = TRUE))
    colnames(tau_cr_df) <- c("tau_c", "tau_r")

    for (st_file in status_files) {
      st_df <- qs_read(st_file, nthreads = N_QS_THREADS)
      st_df$tau_c <- tau_cr_df$tau_c

      st_df$tau_r <- ifelse(
        !is.na(st_df$tau_d) & !is.na(tau_cr_df$tau_r) & tau_cr_df$tau_r >= st_df$tau_d,
        NA_real_,
        tau_cr_df$tau_r
      )

      st_base <- tools::file_path_sans_ext(basename(st_file))
      st_new_base <- paste0(st_base, "_xic_", c_val, "_xir_", r_val)
      st_out_file <- file.path(output_dir, paste0(st_new_base, ".qs"))

      qs_save(st_df, st_out_file, nthreads = N_QS_THREADS)
      cat("  -> Saved status with tau_c/tau_r:", basename(st_out_file), "\n")
    }
  }

  # Free memory between huge files
  rm(cohort_state)
  gc()
}

cat("\nAll status time extractions successfully completed!\n")
print_end_time(start_time, "process-cohort-times-fct-xih-xid-xic-xir.R")
