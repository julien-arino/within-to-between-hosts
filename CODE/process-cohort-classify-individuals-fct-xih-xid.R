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
library(parallel)
library(data.table)

# Find the latest cohort_P* file
if (basename(getwd()) == "CODE") {
  output_dir <- normalizePath(file.path(getwd(), "..", "OUTPUT"))
} else {
  output_dir <- normalizePath(file.path(getwd(), "OUTPUT"))
}
files <- list.files(output_dir, pattern = "^cohort_P.*\\.qs$", full.names = TRUE)
files <- files[!grepl("-times\\.qs$", files)]

if (length(files) == 0) {
  stop("No cohort_P... file found in ", output_dir)
}

# Sort by modification time to get the newest one if multiple exist
latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

# Define grid of thresholds to iterate over
xi_h_vals <- c(50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
xi_d_vals <- c(75.0, 80.0, 85.0, 90.0, 95.0)

# ====================================================================
# PART 1: PROCESS TABLES
# ====================================================================
cat("\n[PART 1] Loading newest cohort results (Tables ONLY):", basename(latest_file), "\n")
cohort_df <- qs_read(latest_file)

base_name <- tools::file_path_sans_ext(basename(latest_file))

base_summary_df <- cohort_df$parameters %>%
  select(
    ID,
    Psi_max = max_Psi,
    t_max = max_Psi_t
  )

# Set up parallel backend using PSOCK cluster to prevent memory duplication and SIGPIPEs
n_cores <- parallel::detectCores()
if (n_cores >= 64) {
  workers_to_use <- max(2, round(n_cores * 2 / 3))
} else {
  workers_to_use <- max(2, n_cores - 2)
}

cat("Extracting tau_h_start values from trajectories...\n")
cl <- makeCluster(workers_to_use)
clusterExport(cl, varlist = c("xi_h_vals"))

tau_h_list <- parLapply(cl, cohort_df$cohort, function(df) {
  psi <- df$Psi
  time <- df$time
  res <- numeric(length(xi_h_vals))
  for (i in seq_along(xi_h_vals)) {
    idx <- which(psi >= xi_h_vals[i])[1]
    res[i] <- if (is.na(idx)) NA_real_ else time[idx]
  }
  res
})
stopCluster(cl)

tau_h_mat <- do.call(rbind, tau_h_list)
colnames(tau_h_mat) <- paste0("tau_h_", xi_h_vals)

base_summary_df <- bind_cols(base_summary_df, as.data.frame(tau_h_mat))

# Immediately free the cohort dataframe
rm(cohort_df, tau_h_list, tau_h_mat)
gc()

cat("Computing and saving status tables for each threshold...\n")
for (xi_d in xi_d_vals) {
  for (xi_h in xi_h_vals) {
    if (xi_h >= xi_d) {
      next
    }

    tau_h_col <- paste0("tau_h_", xi_h)

    current_summary_df <- base_summary_df %>%
      mutate(
        status = case_when(
          Psi_max < xi_h ~ "Mild",
          Psi_max < xi_d ~ "ICU",
          TRUE ~ "Dead"
        ),
        tau_d = ifelse(status == "Dead", t_max, NA_real_),
        tau_h_start = ifelse(status %in% c("ICU", "Dead"), .data[[tau_h_col]], NA_real_)
      ) %>%
      select(-starts_with("tau_h_")) # Clear out all tau_h columns so they are not saved redundantly

    h_str <- ifelse(xi_h %% 1 == 0, as.character(as.integer(xi_h)), as.character(xi_h))
    d_str <- ifelse(xi_d %% 1 == 0, as.character(as.integer(xi_d)), as.character(xi_d))

    new_base <- sub("cohort_", "cohort_status_", base_name)
    new_base <- paste0(new_base, "_xih_", h_str, "_xid_", d_str)

    out_file <- file.path(output_dir, paste0(new_base, ".qs"))
    qs_save(current_summary_df, out_file)
    cat("  -> Saved:", basename(out_file), "\n")
  }
}


# ====================================================================
# PART 2: PROCESS TRAJECTORIES
# ====================================================================
cat("\n[PART 2] Loading newest cohort results AGAIN for trajectories ONLY...\n")
cohort_df <- qs_read(latest_file)
cohort_list_raw <- cohort_df$cohort

# Immediately free the parameters array
rm(cohort_df)
gc()

# Set up parallel backend using PSOCK cluster to prevent memory duplication and SIGPIPEs
n_cores <- parallel::detectCores()
if (n_cores >= 64) {
  workers_to_use <- max(2, round(n_cores * 2 / 3))
} else {
  workers_to_use <- max(2, n_cores - 2)
}

# Function to interpolate and truncate based purely on xi_d
interpolate_and_truncate <- function(df_element, xi_d) {
  df <- as.data.frame(df_element)

  if ("Psi" %in% names(df) && length(df$Psi) > 0) {
    psi_max <- max(df$Psi)
    if (psi_max >= xi_d) {
      t_max <- df$time[which.max(df$Psi)]
    } else {
      t_max <- max(df$time)
    }
  } else {
    t_max <- max(df$time)
  }

  t_interp <- seq(0, t_max, by = 0.1)
  v_interp <- approx(df$time, df$V, xout = t_interp, ties = mean)$y

  data.frame(
    time = t_interp,
    V = v_interp,
    row.names = NULL
  )
}

cat("\nStarting PSOCK cluster...\n")
cl <- makeCluster(workers_to_use)

for (xi_d in xi_d_vals) {
  cat("\n========================================\n")
  cat("Processing trajectories for xi_d =", xi_d, "\n")
  cat("========================================\n")

  processed_list <- parLapply(
    cl,
    cohort_list_raw,
    interpolate_and_truncate,
    xi_d = xi_d
  )

  # Mimic the previous `structured_output$cohort` format for compatibility
  # with downstream scripts (`process-cohort-assign-zero-transmssion.R`).
  # We leave `parameters` out since they are already saved separately.
  structured_output <- list(cohort = processed_list)

  d_str <- ifelse(xi_d %% 1 == 0, as.character(as.integer(xi_d)), as.character(xi_d))

  new_base <- sub("cohort_", "cohort_truncated_", base_name)
  new_base <- paste0(new_base, "_xid_", d_str)

  out_file <- file.path(output_dir, paste0(new_base, ".qs"))
  qs_save(structured_output, out_file)
  cat("Saved:", basename(out_file), "\n")
}

stopCluster(cl)
cat("\nAll thresholds and trajectories processed successfully.\n")
