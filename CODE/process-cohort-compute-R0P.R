#!/usr/bin/env Rscript
# ============================================================
# File: process-cohort-compute-R0P.R
# Description:
#   Calculates the individual R0_i^P (R0_P2P) for the cohort by 
#   integrating the transmission rate beta(V) over the active
#   infectious periods. Elegantly replaces the legacy Julia script
#   and syncs perfectly with Figure E2.
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("process-cohort-compute-R0P.R")

load_libraries(c("qs2", "dplyr", "tidyr"))

# Between-host transmission model constants
alpha_p <- 16.422
k_v <- 7.49
S_P_0 <- 2000.0

# ------------------------------------------------------------
# 1) Locate necessary QS files
# ------------------------------------------------------------
interp_files <- list.files(output_dir, pattern = "^cohort_sim_state_interp_P.*\\.qs$", full.names = TRUE)
if (length(interp_files) == 0) stop("Cannot find interpolated cohort file.")
INTERP_FILE <- interp_files[which.max(file.mtime(interp_files))]

# Status file representing the "baseline" simulation
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xih_75_xid_85_xic_4_xir_1\\.qs$", full.names = TRUE)
if (length(status_files) == 0) {
  status_files <- list.files(output_dir, pattern = "^cohort_status_P.*\\.qs$", full.names = TRUE)
}
STATUS_FILE <- status_files[which.max(file.mtime(status_files))]

cat("Loading interpolated state:", basename(INTERP_FILE), "\n")
cohort_df <- qs_read(INTERP_FILE, nthreads = N_QS_THREADS)

cat("Loading baseline status:", basename(STATUS_FILE), "\n")
status_df <- qs_read(STATUS_FILE, nthreads = N_QS_THREADS)

# ------------------------------------------------------------
# 2) Pivot Viral Load Arrays
# ------------------------------------------------------------
cat("Pivoting viral load array...\n")
V_df <- cohort_df %>%
  select(time, ID, V) %>%
  tidyr::pivot_wider(names_from = ID, values_from = V)

time_pts <- V_df$time
cohort_ids <- as.numeric(colnames(V_df)[-1])
V_mat <- as.matrix(V_df[, -1]) # matrix structure: [Time x ID]

# Sync bounds by column order
status_df <- status_df[match(cohort_ids, status_df$ID), ]

tau_c_v <- status_df$tau_c
tau_r_v <- status_df$tau_r
tau_d_v <- status_df$tau_d
tau_h_start_v <- status_df$tau_h_start
tau_h_end_v <- status_df$tau_h_end

# ------------------------------------------------------------
# 3) Masking non-infectious / out-of-bounds periods
# ------------------------------------------------------------
cat("Applying infectious period boundaries natively...\n")
V_mat_active <- V_mat

V_mat_active[outer(time_pts, coalesce(tau_c_v, Inf), "<")] <- 0
V_mat_active[outer(time_pts, coalesce(tau_r_v, Inf), ">")] <- 0
V_mat_active[outer(time_pts, coalesce(tau_d_v, Inf), ">=")] <- 0
V_mat_active[outer(time_pts, coalesce(tau_h_start_v, Inf), ">=") &
               outer(time_pts, coalesce(tau_h_end_v, -Inf), "<=")] <- 0

# ------------------------------------------------------------
# 4) Compute transmission rates & Integrate (Trapezoidal)
# ------------------------------------------------------------
cat("Computing transmission mapping beta(V)...\n")
beta_mat <- (V_mat_active^alpha_p) / (V_mat_active^alpha_p + k_v^alpha_p)

rm(V_mat_active, V_mat, V_df, cohort_df); gc() # Aggressive RAM sweeping

cat("Integrating arrays calculating the area under beta(t)...\n")
dt <- diff(time_pts)[1] # Typically 0.1
beta_sum <- colSums(beta_mat, na.rm = TRUE)
beta_first <- beta_mat[1, ]
beta_last <- beta_mat[nrow(beta_mat), ]

integral <- dt * (beta_sum - 0.5 * beta_first - 0.5 * beta_last)
R0_P2P_values <- S_P_0 * integral

rm(beta_mat); gc() 

# ------------------------------------------------------------
# 5) Formatting output
# ------------------------------------------------------------
out_df <- data.frame(
  ID = cohort_ids,
  R0_P2P = R0_P2P_values
)

base_name <- sub("cohort_status_", "cohort_times_", basename(STATUS_FILE))
base_name <- sub("_xih.*\\.qs$", ".qs", base_name)
out_path <- file.path(output_dir, base_name)

qs_save(out_df, file = out_path, nthreads = N_QS_THREADS)
cat("Successfully computed and saved R0_P2P bindings to:\n   ", basename(out_path), "\n")

print_end_time(start_time, "process-cohort-compute-R0P.R")
