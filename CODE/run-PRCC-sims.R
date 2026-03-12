# ============================================================
# File : PRCC_parallel_categorized.R
# Parallel PRCC for within-host model with category-based plots
# ============================================================

suppressPackageStartupMessages({
  library(lhs)
  library(furrr)
  library(future)
  library(deSolve)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(latex2exp)
  library(sensitivity)
})

# ------------------------------------------------------------
# 1. Model
# ------------------------------------------------------------
source("functions-all.R") # contains extract_max_indicators()
default_params <- set_parameters()
IC <- set_IC()

OUTPUT_clo <- file.path(getwd(), "OUTPUT_clo")
FIGS_clo <- file.path(getwd(), "FIGS_clo")

dir.create(OUTPUT_clo, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGS_clo, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# 2. Parameter ranges (list format)
# ------------------------------------------------------------
param_ranges <- list(
  beta = c(0.1, 0.5),
  d_I = c(0.05, 0.15),
  d_D = c(0.001, 0.015),
  d_V = c(5, 15),
  epsilon_FI = c(1e-5, 1e-3),
  eta_FI = c(0.001, 0.05),
  k_B_F = c(0.001, 0.05),
  k_int_f = c(10, 20),
  k_lin_f = c(10, 20),
  k_U_F = c(2, 10),
  lambda_S = c(0.5, 1),
  p = c(100, 800),
  p_FI = c(0.8, 4),
  psi_F_prod = c(0.1, 0.5),
  S_max = c(0.1, 0.2),
  c_star = c(1e-5, 1e-3)
)

param_names <- names(param_ranges)
p <- length(param_ranges)

# ------------------------------------------------------------
# 3. Sampling (Latin Hypercube)
# ------------------------------------------------------------
set.seed(123)
N <- 1000000 # Increase if needed
X_unit <- randomLHS(N, p)
X <- as.data.frame(Map(function(r, u) r[1] + u * (r[2] - r[1]), param_ranges, as.data.frame(X_unit)))
colnames(X) <- param_names

# ------------------------------------------------------------
# 4. Parallelization setup
# ------------------------------------------------------------
total_cores <- parallel::detectCores()
n_cores <- max(1L, floor(0.8 * total_cores))
plan(multisession, workers = n_cores)
opts <- furrr_options(seed = TRUE, scheduling = 1)

cat("=== Parallel Execution ===\n")
cat("Detected cores :", total_cores, "\n")
cat("Cores in use   :", n_cores, "\n\n")

simulate_one <- function(i, Xdf) {
  pid <- Sys.getpid()
  if (i %% 100 == 0) {
    message(sprintf("[worker %d] simulation %d/%d", pid, i, nrow(Xdf)))
  }
  pars <- as.list(Xdf[i, , drop = FALSE])
  merged_pars <- modifyList(default_params, pars)
  extract_max_indicators(params_tmp = merged_pars, IC = IC)
}

# ------------------------------------------------------------
# 5. Run simulations
# ------------------------------------------------------------
t0 <- proc.time()[3]
res_list <- future_map(1:nrow(X), ~ simulate_one(.x, X), .progress = TRUE, .options = opts)
t1 <- proc.time()[3]

elapsed <- t1 - t0
sim_ok <- sum(!sapply(res_list, is.null))
cat("=== Simulation completed ===\n")
cat("Simulations requested :", nrow(X), "\n")
cat("Valid results         :", sim_ok, "\n")
cat(sprintf("Elapsed time          : %.2f s (%.1f sim/s)\n\n", elapsed, sim_ok / max(elapsed, 1e-9)))

res_df <- bind_rows(res_list)
stopifnot(all(c("V_max", "F_U_max", "F_B_max", "t_V_max", "t_F_U_max", "t_F_B_max") %in% names(res_df)))

# ------------------------------------------------------------
# 6. Save data for PRCC computation
# ------------------------------------------------------------
library(qs)

sim_data <- list(
  X = X,
  res_df = res_df,
  param_names = param_names
)

qsave(sim_data, file.path(OUTPUT_clo, "PRCC_sim_results.qs"))

# ------------------------------------------------------------
# 7. Clean-up
# ------------------------------------------------------------
future::plan(sequential)
cat("✅ Parallel simulations for PRCC finished. Results saved in OUTPUT_clo/PRCC_sim_results.qs.\n")
