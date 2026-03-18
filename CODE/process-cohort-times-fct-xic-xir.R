# ==============================================================================
# File: process-cohort-times-fct-xic-xir.R
# Description:
#   This script complements zero-transmission grouping by specifically 
#   identifying the exact simulation times when individuals become 
#   infectious (tau_c) and when they recover/cease to be infectious (tau_r).
#
#   It reads the large `cohort_truncated_state` files to find V(t) crossings 
#   for the given xi_c and xi_r thresholds, and appends these extracted 
#   times to matching `cohort_status_` files, saving them with the 
#   corresponding `_xic_X_xir_X.qs` suffixes for downstream plotting.
# ==============================================================================

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))

# Define fixed infectiousness thresholds
xi_c <- seq(from = 4, to = 8, by = 1) # Start of infectious period
xi_r <- 1.0 # End of infectious period

# Find the latest truncated cohort file
if (basename(getwd()) == "CODE") {
  output_dir <- normalizePath(file.path(getwd(), "..", "OUTPUT"))
} else {
  output_dir <- normalizePath(file.path(getwd(), "OUTPUT"))
}

files <- list.files(output_dir, pattern = "^cohort_truncated_state_P.*\\.qs$", full.names = TRUE)

if (length(files) == 0) {
    stop("No cohort_truncated_state_P... file found in ", output_dir)
}

# We process the base truncated states, avoiding already padded ones
files <- files[!grepl("_zerotrans\\.qs$", files)]

if (length(files) == 0) {
    cat("All files have already been processed.\n")
    quit(save = "no")
}

for (file_idx in seq_along(files)) {
    current_file <- files[file_idx]
    cat("\n========================================\n")
    cat(sprintf("[%d/%d] Extracting event times from: %s\n", file_idx, length(files), basename(current_file)))
    cat("========================================\n")

    cohort_state <- qs_read(current_file)

    for (c_val in xi_c) {
      for (r_val in xi_r) {
        cat(sprintf("Computing tau_c/tau_r for [xi_c=%g, xi_r=%g]...\n", c_val, r_val))

        n_cores <- parallel::detectCores()
        if (n_cores >= 64) {
          workers_to_use <- max(2, round(n_cores * 2 / 3))
        } else {
          workers_to_use <- max(2, n_cores - 2)
        }
        
        matches <- regmatches(basename(current_file), regexpr("P[0-9]+_DT[0-9]+-[0-9]+", basename(current_file)))
        if (length(matches) == 0) stop("Could not extract P/DT timestamp from file name.")
        timestamp_part <- matches[1]

        status_pattern <- paste0("^cohort_status_", timestamp_part, ".*\\.qs$")
        status_files <- list.files(output_dir, pattern = status_pattern, full.names = TRUE)
        status_files <- status_files[!grepl("_xic_", status_files) & !grepl("_zerotrans", status_files)]

        cl2 <- makeCluster(workers_to_use)
        clusterExport(cl2, varlist = c("c_val", "r_val"), envir = environment())
        
        tau_cr_list <- parLapply(cl2, cohort_state, function(df) {
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
        })
        stopCluster(cl2)
        
        tau_cr_df <- as.data.frame(do.call(rbind, tau_cr_list))
        
        for (st_file in status_files) {
            st_df <- qs_read(st_file)
            st_df$tau_c <- tau_cr_df$tau_c
            
            st_df$tau_r <- ifelse(
                !is.na(st_df$tau_d) & !is.na(tau_cr_df$tau_r) & tau_cr_df$tau_r >= st_df$tau_d, 
                NA_real_, 
                tau_cr_df$tau_r
            )
            
            st_base <- tools::file_path_sans_ext(basename(st_file))
            st_new_base <- paste0(st_base, "_xic_", c_val, "_xir_", r_val)
            st_out_file <- file.path(output_dir, paste0(st_new_base, ".qs"))
            
            qs_save(st_df, st_out_file)
            cat("  -> Saved status with tau_c/tau_r:", basename(st_out_file), "\n")
        }
      }
    }
    
    # Free memory between huge files
    rm(cohort_state)
    gc()
}

cat("\nAll status time extractions successfully completed!\n")
