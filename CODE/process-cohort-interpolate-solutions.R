#!/usr/bin/env Rscript

## process-interpolate-solutions.R
# This script reads the raw adaptive ODE output from Julia, forces it into a
# standardized 0.1 day step grid over parallel forks and saves the bound
# dataframe back to disk for plotting and distribution matrices.

suppressWarnings(suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(future.apply)
  library(here)
}))

cat("\n\n>>> Running process-interpolate-solutions.R ...\n\n")

project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

# ------------------------------------------------------------
# 1. LOAD RAW TRAJECTORIES
# ------------------------------------------------------------
sim_state_files <- list.files(output_dir, pattern = "cohort_sim_state_P.*\\.qs$", full.names = TRUE)
# Exclude any previously interpolated files if they match patterns
sim_state_files <- sim_state_files[!grepl("_interp", sim_state_files)]

if (length(sim_state_files) == 0) stop("No raw cohort_sim_state file found.")
cohort_sim_file <- sim_state_files[which.max(file.mtime(sim_state_files))]

cat("Loading newest continuous trajectories:", basename(cohort_sim_file), "\n")
cohort_list <- qs_read(cohort_sim_file, nthreads = N_QS_THREADS)

ram_size <- format(object.size(cohort_list), units = "auto")
cat(sprintf("Loaded %d trajectories into memory (Total RAM size: %s)\n", length(cohort_list), ram_size))

# ------------------------------------------------------------
# 2. RUN PARALLEL INTERPOLATION
# ------------------------------------------------------------
cat("\nComputing generalized parallel interpolation grid (0.1 steps)...\n")
n_cores <- parallel::detectCores()
# Strictly clamp worker threads to 24 for colossal high-core systems to prevent OS Out-Of-Memory (OOM) 
workers_to_use <- min(24, max(2, n_cores - 2))
plan(multicore, workers = workers_to_use)

# Adapt to max time step in the ODE
max_time <- 100.0

interpolate_trajectory <- function(lst_element) {
  if (is.null(lst_element)) return(NULL)
  df <- as.data.frame(lst_element)
  if (nrow(df) == 0) return(NULL)
  
  t_interp <- seq(0, max_time, by = 0.1)
  res <- data.frame(time = t_interp, stringsAsFactors = FALSE)
  
  if ("V" %in% names(df)) {
    res$V <- approx(df$time, df$V, xout = t_interp, rule = 2, ties = mean)$y
  } else {
    res$V <- NA
  }
  
  if ("Psi" %in% names(df)) {
    res$Psi <- approx(df$time, df$Psi, xout = t_interp, rule = 2, ties = mean)$y
  }
  
  if ("F_B" %in% names(df)) {
    res$F_B <- approx(df$time, df$F_B, xout = t_interp, rule = 2, ties = mean)$y
  }
  
  if ("F_U" %in% names(df)) {
    res$F_U <- approx(df$time, df$F_U, xout = t_interp, rule = 2, ties = mean)$y
  }
  
  return(res)
}

cohort_list_interp <- future_lapply(cohort_list, interpolate_trajectory, future.seed = TRUE)

if (!is.null(names(cohort_list))) {
  names(cohort_list_interp) <- names(cohort_list)
} else {
  names(cohort_list_interp) <- as.character(seq_along(cohort_list_interp))
}

# Free RAM 
rm(cohort_list)
gc()

cat("\nBinding interpolated trajectories into monolithic dataframe...\n")
beta_df_raw <- dplyr::bind_rows(Filter(Negate(is.null), cohort_list_interp), .id = "ID")

rm(cohort_list_interp)
gc()

# ------------------------------------------------------------
# 3. SAVE INTERPOLATED MATRIX
# ------------------------------------------------------------
out_filename <- sub("cohort_sim_state_", "cohort_sim_state_interp_", basename(cohort_sim_file))
out_interp_path <- file.path(output_dir, out_filename)

cat("\nSaving interpolated flat matrix to QS...\n")
qs_save(beta_df_raw, out_interp_path, nthreads = N_QS_THREADS)

cat("✅ Interpolation complete and saved to", out_filename, "\n")
