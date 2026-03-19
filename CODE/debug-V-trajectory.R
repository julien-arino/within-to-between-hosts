#!/usr/bin/env Rscript
# ============================================================
# File: debug-V-trajectory.R
# Description:
# ============================================================

cat("\n\n>>> Running debug-V-trajectory.R ...\n\n")
library(qs2)
library(dplyr)
library(ggplot2)
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

latest_file <- list.files("OUTPUT", pattern = "^cohort_P.*\\.qs$", full.names = TRUE) %>% 
  .[!grepl("(_times_|_censored)", .)] %>% 
  .[which.max(file.mtime(.))]

save_data <- qs_read(latest_file, nthreads = N_QS_THREADS)
cohort <- save_data$cohort
params <- save_data$parameters

# Let's find an individual with R0 > 6.9 but max_V < 4.501
idx <- which(params$R0_within > 6.9 & sapply(cohort, function(x) x$maxima$max_V) < 4.501)[1]

cat(sprintf("Selected Individual %d for diagnostic.\n", idx))
cat(sprintf("R0_within = %.2f\n", params$R0_within[idx]))
cat(sprintf("max_V     = %.4f\n", cohort[[idx]]$maxima$max_V))

time_vec <- cohort[[idx]]$vars$time
v_vec <- cohort[[idx]]$vars$V
i_vec <- cohort[[idx]]$vars$I

# Look at the first 10 time steps
cat("\nFirst 10 timesteps:\n")
print(data.frame(time = time_vec[1:10], V = v_vec[1:10], I = i_vec[1:10]))

cat("\n\n>>> debug-V-trajectory.R successfully finished running ✅\n\n")
