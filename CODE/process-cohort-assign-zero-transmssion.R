# ==============================================================================
# File: process-cohort-assign-zero-transmssion.R
# Description:
#   This script processes truncated virtual cohort simulation results.
#   It assigns V = 0 to segments of the trajectory that are considered
#   non-transmitting, based on the infectious period thresholds:
#     - xi_c (start of infectious period): V must reach this to start transmitting.
#     - xi_r (end of infectious period): V dropping below this ends transmission.
#
#   V(t) is manually set to 0.0 for all `t < tau_c` and `t > tau_r`.
#   The processed data is saved as a new `.qs` file.
# ==============================================================================

library(dplyr)
library(qs2)

# Define fixed infectiousness thresholds
xi_c <- 9.0  # Start of infectious period
xi_r <- 1.0  # End of infectious period

# Find the latest truncated cohort file
output_dir <- file.path(getwd(), "OUTPUT")
files <- list.files(output_dir, pattern = "^cohort_truncated_P.*\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No cohort_truncated_P... file found in ", output_dir)
}

# Exclude files that have already been processed to avoid infinite loops if re-run
files <- files[!grepl("_zerotrans\\.qs$", files)]

if (length(files) == 0) {
  cat("All files have already been processed.\n")
  quit(save = "no")
}

for (file_idx in seq_along(files)) {
  current_file <- files[file_idx]
  cat("\n========================================\n")
  cat(sprintf("[%d/%d] Loading truncated cohort results: %s\n", file_idx, length(files), basename(current_file)))
  cat("========================================\n")
  
  # Load the truncated dataset
  structured_output <- qs_read(current_file)


# We will modify the V trajectories within structured_output$cohort
cat(sprintf("Setting V=0 outside infectious period [xi_c=%g, xi_r=%g]...\n", xi_c, xi_r))

processed_cohort <- lapply(structured_output$cohort, function(df) {
  # Find peak index
  V_max_idx <- which.max(df$V)
  
  # Find tau_c (first time V >= xi_c)
  c_idx <- which(df$V >= xi_c)[1]
  
  # Find tau_r (first time V <= xi_r AFTER the peak)
  r_idx <- NA
  if (V_max_idx < nrow(df)) {
    post_peak_V <- df$V[(V_max_idx + 1):nrow(df)]
    r_offset <- which(post_peak_V <= xi_r)[1]
    if (!is.na(r_offset)) {
      r_idx <- V_max_idx + r_offset
    }
  }
  
  # Create a copy of the trajectory
  df_new <- df
  
  # Apply pre-infectious zeroing
  if (!is.na(c_idx) && c_idx > 1) {
    df_new$V[1:(c_idx - 1)] <- 0.0
  } else if (is.na(c_idx)) {
    # Never became infectious
    df_new$V <- 0.0
  }
  
  # Apply post-infectious zeroing
  if (!is.na(r_idx) && r_idx <= nrow(df_new)) {
    df_new$V[r_idx:nrow(df_new)] <- 0.0
  }
  
  return(df_new)
})

# Replace the old cohort data with the zeroed version
structured_output$cohort <- processed_cohort

  # Generate new filename
  base_name <- tools::file_path_sans_ext(basename(current_file))
  new_base <- paste0(base_name, "_zerotrans")
  out_file <- file.path(output_dir, paste0(new_base, ".qs"))

  # Save the zero-trans processed dataset
  qs_save(structured_output, out_file)

  cat("Saved zero-transmission cohort to:", basename(out_file), "\n")
}

cat("\nAll files successfully processed!\n")
