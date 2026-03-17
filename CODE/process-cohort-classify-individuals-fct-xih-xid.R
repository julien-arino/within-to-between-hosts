# ==============================================================================
# File: process-cohort-classify-individuals-fct-xih-xid.R
# Description:
#   This script processes the full virtual cohort simulation results.
#   It iterates over a grid of xi_h (hospitalization) and xi_d (death)
#   thresholds, similar to process-virtual-cohort-times-vary-xih-xid.jl.
#   For each combination, it identifies individuals based on their Psi_max:
#     - Mild: Psi_max < xi_h
#     - ICU: xi_h <= Psi_max < xi_d
#     - Dead: Psi_max >= xi_d
#
#   For individuals classified as "Dead", their simulation trajectories are
#   truncated exactly at the time they hit their maximum damage score (t_max).
#   The cleaned and truncated data is saved for each combination of cutpoints.
# ==============================================================================

library(dplyr)
library(qs2)

# Find the latest cohort_P* file
output_dir <- file.path(getwd(), "OUTPUT")
files <- list.files(output_dir, pattern = "^cohort_P.*\\.qs$", full.names = TRUE)
files <- files[!grepl("-times\\.qs$", files)]

if (length(files) == 0) {
  stop("No cohort_P... file found in ", output_dir)
}

# Sort by modification time to get the newest one if multiple exist
latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

# Load the non-truncated dataset
cohort_df <- qs_read(latest_file)

base_name <- tools::file_path_sans_ext(basename(latest_file))

cat("Interpolating trajectories to dt = 0.1 grid...\n")
# Interpolate each individual's time and V directly upon load
# Discard all other state variables from the cohort list
cohort_df$cohort <- lapply(seq_along(cohort_df$cohort), function(i) {
  df <- as.data.frame(cohort_df$cohort[[i]])
  
  # Interpolation grid: 0 to maximum simulated time with step 0.1
  dt <- 0.1
  # Ensure we don't go beyond the last time point
  t_interp <- seq(0, max(df$time), by = dt)
  
  # Interpolate V using linear interpolation
  # `ties = mean` handles any duplicated time points (often caused by adaptive ODE stepping)
  v_interp <- approx(df$time, df$V, xout = t_interp, ties = mean)$y
  
  data.frame(
    time = t_interp,
    V = v_interp,
    row.names = NULL
  )
})

# Define grid of thresholds to iterate over (matching the Julia script)
xi_h_vals <- c(50.0, 60.0, 70.0, 75.0, 80.0)
xi_d_vals <- c(75.0, 80.0, 85.0, 90.0, 95.0)

# Merge and basic prep of trajectories happens ONCE outside the loop to save time
# cohort_df$cohort is the list of DataFrames containing the trajectories
library(tidyr)
cat("Preparing trajectory dataframe...\n")
# WARNING: For large cohorts (e.g. N=1,000,000) binding all trajectories
# into a single dataframe can consume 60GB+ of RAM during the inner_join later.
# Ensure the machine has sufficient memory (128GB+) before running with large N.
cohort_trajectories <- bind_rows(lapply(seq_along(cohort_df$cohort), function(i) {
  cohort_df$cohort[[i]] %>%
    mutate(ID = cohort_df$parameters$ID[i])
}))

# Base summary without status
base_summary_df <- cohort_df$parameters %>%
  select(
    ID,
    Psi_max = max_Psi,
    t_max = max_Psi_t
  )

for (xi_d in xi_d_vals) {
  for (xi_h in xi_h_vals) {
    if (xi_h > xi_d) {
      next
    }

    cat("\n========================================\n")
    cat("Processing for xi_h =", xi_h, "and xi_d =", xi_d, "\n")
    cat("========================================\n")

    # Recompute status based on current thresholds
    current_summary_df <- base_summary_df %>%
      mutate(
        status = case_when(
          Psi_max < xi_h ~ "Mild",
          Psi_max < xi_d ~ "ICU",
          TRUE ~ "Dead"
        )
      )

    cohort_truncated <- cohort_trajectories %>%
      inner_join(current_summary_df, by = "ID") %>%
      filter(Psi_max < xi_d | time <= t_max)

    # Select Columns and split back into a list of data.frames
    cohort_list <- cohort_truncated %>%
      select(time, V, ID)
    
    # Split by ID into a list of dataframes, and drop the ID column from each
    cohort_list <- split(cohort_list[, c("time", "V")], cohort_list$ID)

    # Prepare parameters dataframe
    parameters_df <- current_summary_df %>%
      select(ID, status)

    # Re-assemble object to mimic Julia input
    structured_output <- list(
      cohort = cohort_list,
      parameters = parameters_df
    )

    # Format filename suffix based on thresholds
    h_str <- ifelse(xi_h %% 1 == 0, as.character(as.integer(xi_h)), as.character(xi_h))
    d_str <- ifelse(xi_d %% 1 == 0, as.character(as.integer(xi_d)), as.character(xi_d))

    new_base <- sub("cohort_", "cohort_truncated_", base_name)
    new_base <- paste0(new_base, "_xih_", h_str, "_xid_", d_str)

    out_file <- file.path(output_dir, paste0(new_base, ".qs"))

    # Save structured data
    qs_save(structured_output, out_file)

    cat("Saved:", basename(out_file), "\n")
  }
}

cat("\nAll thresholds processed successfully.\n")
