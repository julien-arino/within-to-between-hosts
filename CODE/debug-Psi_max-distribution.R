#!/usr/bin/env Rscript
cat("\n\n>>> Running debug-Psi_max-distribution.R ...\n\n")
library(qs2)
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}
df <- qs_read("OUTPUT/cohort_times_P1000000_DT20260312-210843_xih_70_xid_95.qs", nthreads = N_QS_THREADS)

cat("=== Distribution of Psi_max ===\n")
print(summary(df$Psi_max))

pct_mild <- mean(df$Psi_max < 70) * 100
pct_icu <- mean(df$Psi_max >= 70 & df$Psi_max < 75) * 100
pct_dead <- mean(df$Psi_max >= 75) * 100

cat("\nPercentages for xih = 70, xid = 75:\n")
cat(sprintf("Mild (< 70%%):      %.2f%%\n", pct_mild))
cat(sprintf("ICU  (70%% to < 75%%): %.2f%%\n", pct_icu))
cat(sprintf("Dead (>= 75%%):      %.2f%%\n", pct_dead))

cat("\n\n>>> debug-Psi_max-distribution.R successfully finished running ✅\n\n")
