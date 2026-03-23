#!/usr/bin/env Rscript

# ============================================================
# File: process-cohort-prepare-Figure-C4.R
# Description: Processes massive interpolated trajectories grouped 
#              by status to produce Figure C4 aggregate timelines.
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("process-cohort-prepare-Figure-C4.R")

load_libraries(c("qs2", "dplyr", "here"))

# Automatically find the latest files
files <- list.files(output_dir, pattern = "^cohort_sim_state_interp_.*\\.qs$", full.names = TRUE)
if (length(files) == 0) stop("No cohort_sim_state_interp file found in ", output_dir)
latest_file <- files[which.max(file.mtime(files))]

# Extract standard run identifier (e.g. P1000000_DT20260320-083457)
run_id <- sub("^cohort_sim_state_interp_([P0-9_DT\\-]+)\\.qs$", "\\1", basename(latest_file))

cat("Loading newest interpolated cohort trajectories:", basename(latest_file), "\n")
cohort_df <- qs_read(latest_file, nthreads = N_QS_THREADS)

# Load target status mapping
status_pattern <- "^cohort_status_P.*_xih_75_xid_85_xic_4_xir_1\\.qs$"
status_files <- list.files(output_dir, pattern = status_pattern, full.names = TRUE)
if (length(status_files) == 0) {
  # fallback to base xih_75 if combo missing
  status_pattern <- "^cohort_status_P.*_xih_75_xid_85\\.qs$|^cohort_status_P.*_xid_85\\.qs$"
  status_files <- list.files(output_dir, pattern = status_pattern, full.names = TRUE)
}
if (length(status_files) == 0) stop("No baseline cohort_status file found")
latest_status_file <- status_files[which.max(file.mtime(status_files))]
cat("Loading newest cohort baseline statuses:", basename(latest_status_file), "\n")
cohort_status <- qs_read(latest_status_file, nthreads = N_QS_THREADS)

# Stitch on their severity status purely from the official baseline output
cols_to_grab <- intersect(names(cohort_status), c("ID", "status", "max_Psi", "tau_max_Psi", "tau_c", "tau_r", "tau_h_start", "tau_h_end"))
cohort_df <- cohort_df %>%
  inner_join(cohort_status %>% select(all_of(cols_to_grab)) %>% mutate(ID = as.character(ID)),
    by = "ID"
  )

cohort_df$status <- factor(cohort_df$status, levels = c("Mild", "ICU", "Dead"))

# Provide fallbacks if missing from interpolated sim trace
for (col in c("I", "F_U", "F_B", "Psi", "V")) {
  if (!col %in% names(cohort_df)) cohort_df[[col]] <- NA
}

# In raw data, dead individuals stop executing after xi_d (85%). Mask out post-death.
xi_d <- 85
cohort_df <- cohort_df %>%
  mutate(
    is_dead = (!is.na(max_Psi) & max_Psi >= xi_d),
    # If dead, mask variables after death time
    V = ifelse(is_dead & time > tau_max_Psi, NA, V),
    Psi = ifelse(is_dead & time > tau_max_Psi, NA, Psi),
    F_B = ifelse(is_dead & time > tau_max_Psi, NA, F_B),
    F_U = ifelse(is_dead & time > tau_max_Psi, NA, F_U),
    I = ifelse(is_dead & time > tau_max_Psi, NA, I)
  )

# Parameters
xi_c <- 4 
xi_r <- 1.0
xi_h <- 75
alpha <- 16.422
k_v <- 7.49

# Compute beta_hat
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

cohort_df$beta_hat <- compute_beta_hat(cohort_df$V)

# Build effective beta_c(a)
cohort_df <- cohort_df %>%
  mutate(
    in_window = (!is.na(tau_c) & time >= tau_c) & (is.na(tau_r) | time <= tau_r),
    in_hosp = (!is.na(tau_h_start) & time >= tau_h_start) & (is.na(tau_h_end) | time <= tau_h_end),
    beta_hat_window = case_when(
      is_dead & time > tau_max_Psi ~ NA_real_,
      in_window & !in_hosp ~ beta_hat,
      TRUE ~ 0
    )
  )

# Compute mean ± SD per status
cat("- Aggregating bounds across millions of array steps...\n")
mean_df <- cohort_df %>%
  group_by(status, time) %>%
  summarise(
    Psi_mean = mean(Psi, na.rm = TRUE),
    Psi_sd = sd(Psi, na.rm = TRUE),
    V_mean = mean(V, na.rm = TRUE),
    V_sd = sd(V, na.rm = TRUE),
    I_mean = mean(I, na.rm = TRUE),
    I_sd = sd(I, na.rm = TRUE),
    F_B_mean = mean(F_B, na.rm = TRUE),
    F_B_sd = sd(F_B, na.rm = TRUE),
    F_U_mean = mean(F_U, na.rm = TRUE),
    F_U_sd = sd(F_U, na.rm = TRUE),
    beta_mean = mean(beta_hat_window, na.rm = TRUE),
    beta_sd = sd(beta_hat_window, na.rm = TRUE),
    .groups = "drop"
  )

out_file <- file.path(output_dir, paste0("cohort_FigC4_", run_id, ".qs"))
cat("- Saving to", basename(out_file), "\n")
qs_save(mean_df, out_file, nthreads = N_QS_THREADS)

print_end_time(start_time, "process-cohort-prepare-Figure-C4.R")
