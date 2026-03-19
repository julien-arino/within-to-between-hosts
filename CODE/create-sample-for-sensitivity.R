#!/usr/bin/env Rscript
# Set up the sample for the sensitivity analysis

cat("\n\n>>> Running create-sample-for-sensitivity.R ...\n\n")
suppressWarnings(suppressPackageStartupMessages({
    library(sensitivity)
    library(randtoolbox)
    library(here)
}))

# Set project root automatically relative to the .git tracking directory
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

# The sampling is executed inside Julia via RCall, so N is defined
# by the julia script. However, we load our utility R functions first
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

# Use the same order as in plot-sensitivity-from-julia.R
cols <- c(
  "λ_S",
  "S_max",
  "d_I",
  "d_D",
  "β_V",
  "d_V",
  "p",
  "k_U_F",
  "p_FI",
  "ψ_F_prod",
  "k_B_F",
  "k_lin_f",
  "k_int_f",
  "ε_FI",
  "η_FI",
  "c_star",
  "a_F",
  "V0"
)

pars_df <- data.frame(params = cols)

tmp <- matrix(c(
  0.5, 1, # λ_S
  0.1, 0.2, # S_max
  0.05, 0.15, # d_I
  3, 15, # d_D
  0.1, 0.5, # β_V
  5, 15, # d_V
  100, 800, # p
  2, 10, # k_U_F
  0.8, 4, # p_FI
  0.1, 0.4, # ψ_F_prod
  0.001, 0.05, # k_B_F
  10, 20, # k_lin_f
  10, 20, # k_int_f
  1e-5, 1e-3, # ε_FI
  0.001, 0.05, # η_FI
  1e-5, 1e-3, # c_star
  0.1, 0.2, # a_F
  0.5, 5 # V0
), nc = 2, byrow = TRUE)
pars_df$min <- tmp[, 1]
pars_df$max <- tmp[, 2]

# To use sensitivity::parameterSets, we need to convert the data frame to a list
pars_list <- list()
for (i in 1:dim(pars_df)[1]) {
  pars_list[[pars_df$params[i]]] <- c(pars_df$min[i], pars_df$max[i])
}
nb_pars <- length(pars_list)
pars_sobol <- parameterSets(
  par.ranges = pars_list,
  samples = N,
  method = "sobol"
)
pars_sobol <- as.data.frame(pars_sobol)
pars_sobol <- cbind(1:dim(pars_sobol)[1], pars_sobol)
colnames(pars_sobol) <- c("ID", pars_df$params)

# The initial values of the remaining static state variables:
IC <- set_IC()

# Constants
avo <- 6.02214e23
# MM_F <- 19000.0
# R_F_T <- 1000.0
# R_F_I <- 1300.0
# a_F <- (MM_F / avo) * (R_F_I + R_F_T) * (1 / 5000) * (10^9 * 1e12)

# Add these columns to pars_sobol to create individuals
individuals <- cbind(
  pars_sobol,
  avo = rep(avo, nrow(pars_sobol))
  # MM_F = rep(MM_F, nrow(pars_sobol)),
  # R_F_T = rep(R_F_T, nrow(pars_sobol)),
  # R_F_I = rep(R_F_I, nrow(pars_sobol)),
  # a_F = rep(a_F, nrow(pars_sobol))
)
# Add initial condition columns to individuals (excluding V0 since it's now in the sampled params dataframe)
individuals <- cbind(
  individuals,
  S0 = rep(IC[2], nrow(pars_sobol)),
  I0 = rep(IC[3], nrow(pars_sobol)),
  R0 = rep(IC[4], nrow(pars_sobol))
)

# Cleanup before returning to julia
rm(list = c(
  "pars_sobol",
  "pars_list",
  "pars_df",
  "IC",
  "avo"
))

cat("\n\n>>> create-sample-for-sensitivity.R successfully finished running ✅\n\n")
