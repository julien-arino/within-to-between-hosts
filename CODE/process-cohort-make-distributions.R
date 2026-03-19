# ==============================================================================
# File: process-cohort-make-distributions.R
# Description:
#   This script extracts core statistical distributions (gamma, mu, beta) 
#   from raw cohort status and simulation trajectory data. 
# 
#   By isolating these computations from the visualization files, we ensure 
#   data processing logic is cleanly decoupled and standardly formatted 
#   for downstream numerical models (like incidence plotting).
# ==============================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(future.apply)
  library(here)
}))

cat("\n\n>>> Running process-cohort-make-distributions.R ...\n\n")

project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-all.R")) else source("functions-all.R")
}

# ------------------------------------------------------------
# 1. LOAD REQUIRED RAW COHORT FILES
# ------------------------------------------------------------
# Continuous Trajectories
sim_state_files <- list.files(output_dir, pattern = "cohort_sim_state_P.*\\.qs$", full.names = TRUE)
if (length(sim_state_files) == 0) stop("No cohort_sim_state file found.")
cohort_sim_file <- sim_state_files[which.max(file.mtime(sim_state_files))]
cat("Loading newest continuous trajectories:", basename(cohort_sim_file), "\n")
cohort_list <- qs_read(cohort_sim_file, nthreads = N_QS_THREADS)

# Extract core Base prefix string to properly name output models
base_prefix <- sub("cohort_sim_state_([A-Za-z0-9_-]+)\\.qs", "\\1", basename(cohort_sim_file))

# Status Metadata Defaults
status_files <- list.files(output_dir, pattern = "cohort_status_P.*_xih_75_xid_85(?:_xic.*)?\\.qs$", full.names = TRUE)
if (length(status_files) == 0) stop("No base cohort_status_xih_75_xid_85 file found")
base_status_file <- status_files[which.max(file.mtime(status_files))]
cat("Loading newest baseline statuses:", basename(base_status_file), "\n")
cohort_status <- qs_read(base_status_file, nthreads = N_QS_THREADS)

# Pre-extract necessary vectors and clear the massive status dataframe to save RAM!
alive_ids <- cohort_status$ID[cohort_status$status != "Dead"]

# Garbage collect cohort_status immediately
rm(cohort_status)
gc()

# ------------------------------------------------------------
# 2. COMPUTE GAMMA (RECOVERY TIME DENSITIES)
# ------------------------------------------------------------
cat("\nComputing gamma (recovery time) density probabilities via parallel workers...\n")

n_cores <- parallel::detectCores()
# Strictly clamp worker threads to 24 for colossal high-core systems to prevent OS Out-Of-Memory (OOM) 
workers_to_use <- min(24, max(2, n_cores - 2))
plan(multicore, workers = workers_to_use)

xi_h <- 75
xi_d <- 85
xi_c_values <- c("4" = 4, "6" = 6)

gamma_list <- list()

for (name in names(xi_c_values)) {
  xi_c_val <- xi_c_values[[name]]
  cat(sprintf("Processing for xi_c = %s ...\n", name))
  
  # Determine tau_r dynamically based on trajectory V falling below xi_c
  get_tau_r <- function(df_input) {
    if (is.null(df_input)) return(NA_real_)
    df <- as.data.frame(df_input)
    if (nrow(df) == 0) return(NA_real_)
    if (!"V" %in% names(df)) return(NA_real_)
    
    peak_idx <- which.max(df$V)
    if (length(peak_idx) == 0 || peak_idx == nrow(df)) return(NA_real_)
    
    post_peak_V <- df$V[(peak_idx + 1):nrow(df)]
    r_idx_relative <- which(post_peak_V <= xi_c_val)[1]
    
    if (is.na(r_idx_relative)) return(NA_real_)
    r_idx <- peak_idx + r_idx_relative
    
    return(df$time[r_idx])
  }
  
  tau_r_vals <- future_sapply(cohort_list, get_tau_r, future.seed = TRUE)
  
  if (is.null(names(tau_r_vals))) {
    names(tau_r_vals) <- as.character(seq_along(tau_r_vals))
  }
  
  # Subset to recovered individuals only (Alive) - note alive_ids was extracted at start
  alive_tau_r <- tau_r_vals[as.character(alive_ids)]
  alive_tau_r <- alive_tau_r[!is.na(alive_tau_r)]
  
  if (length(alive_tau_r) < 2) {
    cat(sprintf("Not enough recovered patients for xi_c = %s\n", name))
    next
  }
  
  # Compute KDE bounded to 35 days
  dens <- density(alive_tau_r, from = 0, to = 35)
  
  gamma_df <- tibble(
    a = dens$x,
    gamma_a = dens$y,
    xi_c = name,
    xi_label = paste0("gamma_xi_", name)
  )
  
  gamma_list[[name]] <- gamma_df
}

gamma_all <- bind_rows(gamma_list)
gamma_overall <- gamma_all %>%
  select(time = a, xi_label, gamma_a) %>%
  pivot_wider(
    id_cols     = time,
    names_from  = xi_label,
    values_from = gamma_a
  ) %>%
  arrange(time)

out_gamma <- file.path(output_dir, paste0("cohort_distributions_", base_prefix, "_gamma.qs"))
qs_save(gamma_overall, out_gamma, nthreads = N_QS_THREADS)
cat("Gamma matrix saved to", basename(out_gamma), "\n")


# ------------------------------------------------------------
# 3. COMPUTE MU (DEATH TIME DENSITIES)
# ------------------------------------------------------------
cat("\nComputing mu (death time) density probabilities...\n")

mu_list <- list()
xi_d_targets <- c("85" = "85", "95" = "95")

for (name in names(xi_d_targets)) {
  xid_val <- xi_d_targets[[name]]
  
  # Locate corresponding status file
  fname_pattern <- paste0("cohort_status_P.*_xih_75_xid_", xid_val, "\\.qs$")
  status_files_target <- list.files(output_dir, pattern = fname_pattern, full.names = TRUE)
  
  if (length(status_files_target) == 0) {
    cat(sprintf("Could not find cohort_status file for xi_d = %s\n", xid_val))
    next
  }
  
  target_file <- status_files_target[which.max(file.mtime(status_files_target))]
  target_df <- qs_read(target_file, nthreads = N_QS_THREADS)
  
  dead_cohort <- target_df %>% filter(status == "Dead")
  
  if (nrow(dead_cohort) < 2) {
    cat(sprintf("Not enough dead patients for xi_d = %s\n", xid_val))
    next
  }
  
  # Compute KDE bounded to 20 days
  dens <- density(dead_cohort$tau_d, from = 0, to = 20)
  
  mu_df <- tibble(
    a = dens$x,
    mu_a = dens$y,
    xi_d = name,
    xi_label = paste0("mu_xid_", name)
  )
  
  mu_list[[name]] <- mu_df
}

mu_all <- bind_rows(mu_list)
mu_overall <- mu_all %>%
  select(time = a, xi_label, mu_a) %>%
  pivot_wider(
    id_cols     = time,
    names_from  = xi_label,
    values_from = mu_a
  ) %>%
  arrange(time)

out_mu <- file.path(output_dir, paste0("cohort_distributions_", base_prefix, "_mu.qs"))
qs_save(mu_overall, out_mu, nthreads = N_QS_THREADS)
cat("Mu matrix saved to", basename(out_mu), "\n")

# ------------------------------------------------------------
# FREE UP RAM BEFORE PARALLELIZATION
# ------------------------------------------------------------
cat("\nClearing intermediate distribution datasets from memory before loading interpolated files...\n")
rm(list = intersect(ls(), c(
  "cohort_list",
  "gamma_list", "gamma_all", "gamma_overall", "gamma_df",
  "mu_list", "mu_all", "mu_overall", "mu_df",
  "target_df", "dead_cohort", "tau_r_vals", "alive_tau_r", "alive_ids"
)))
gc()

# ------------------------------------------------------------
# 4. COMPUTE BETA (TRANSMISSION RATES IN TIME)
# ------------------------------------------------------------
cat("\nLoading pre-interpolated monolithic simulation trajectories...\n")
interp_files <- list.files(output_dir, pattern = "cohort_sim_state_interp_.*\\.qs$", full.names = TRUE)
if (length(interp_files) == 0) stop("No interpolated trajectory files found! You must run process-interpolate-solutions.R first.")
latest_interp <- interp_files[which.max(file.mtime(interp_files))]

beta_df_raw <- qs_read(latest_interp, nthreads = N_QS_THREADS)

# Fetch Hospitalization times from standard baseline bounds AFTER interpolation
cat("Reloading cohort status to extract hospitalization times...\n")
cohort_status <- qs_read(base_status_file, nthreads = N_QS_THREADS)

tau_bind_df <- cohort_status %>%
  select(ID, tau_h_start, tau_h_end) %>%
  mutate(ID = as.character(ID))

# Free cohort_status immediately 
rm(cohort_status)
gc()

beta_df_raw <- beta_df_raw %>% left_join(tau_bind_df, by = "ID")

# Empirical compute
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}
beta_df_raw$beta_hat <- compute_beta_hat(beta_df_raw$V)

# Omit hospital periods securely
beta_df_raw <- beta_df_raw %>%
  mutate(
    tau_h_start = replace_na(tau_h_start, Inf),
    tau_h_end   = replace_na(tau_h_end, -Inf),
    beta_eff = if_else(
      time >= tau_h_start & time <= tau_h_end,
      0,
      beta_hat
    )
  )

# Transmission threshold standard filter
xi_c_transmission <- 0.001
cohort_transmitters <- beta_df_raw %>% filter(beta_eff >= xi_c_transmission)

beta_overall <- cohort_transmitters %>%
  group_by(time) %>%
  summarise(
    beta_mean   = mean(beta_eff, na.rm = TRUE),
    beta_sd     = sd(beta_eff,   na.rm = TRUE),
    beta_q10    = quantile(beta_eff, 0.10, na.rm = TRUE),
    beta_q90    = quantile(beta_eff, 0.90, na.rm = TRUE),
    beta_median = median(beta_eff, na.rm = TRUE),
    .groups     = "drop"
  )

out_beta <- file.path(output_dir, paste0("cohort_distributions_", base_prefix, "_beta.qs"))
qs_save(beta_overall, out_beta, nthreads = N_QS_THREADS)
cat("\n\tBeta matrix saved to", basename(out_beta), "\n")

cat("\n✅ All cohort distributions successfully computed and saved!\n")
