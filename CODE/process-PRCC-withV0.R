## process-PRCC-withV0.R
# Computes PRCC from simulation results and saves outputs specifically for the withV0 dataset.

suppressPackageStartupMessages({
    library(dplyr)
    library(sensitivity)
    library(qs2)
})

OUTPUT_dir <- file.path(getwd(), "OUTPUT")

# Find the latest sensitivity_*_withV0 file
sensitivity_files <- list.files(OUTPUT_dir, pattern = "^sensitivity_.*_withV0\\.qs$", full.names = TRUE)
if (length(sensitivity_files) == 0) {
    stop("No sensitivity 'withV0' file found in OUTPUT directory. Run run-sensitivity-analysis-sims-withV0.jl first.")
}
latest_file <- sensitivity_files[which.max(file.info(sensitivity_files)$mtime)]
cat("Loading simulation data from:", latest_file, "\n")

# Load exactly what was exported from the simulation phase
sim_data <- qs_read(latest_file)

# Extract changing variables from the individuals parameter dataframe.
# Notice "V0" is intentionally NOT excluded here because it is part of the PRCC sensitivity matrix.
exclude_cols <- c("ID", "avo", "S0", "I0", "R0", "R0_within")
param_names <- setdiff(names(sim_data$parameters), exclude_cols)
X <- sim_data$parameters[, param_names]

# Extract cohort results
res_df <- dplyr::bind_rows(sim_data$cohort)

##
## 1. PRCC computation
##
prcc_values <- list(
    V_max   = pcc(X, res_df$max_V, rank = TRUE),
    F_U_max = pcc(X, res_df$max_F_U, rank = TRUE),
    F_B_max = pcc(X, res_df$max_F_B, rank = TRUE)
)
prcc_times <- list(
    t_V_max   = pcc(X, res_df$max_V_t, rank = TRUE),
    t_F_U_max = pcc(X, res_df$max_F_U_t, rank = TRUE),
    t_F_B_max = pcc(X, res_df$max_F_B_t, rank = TRUE)
)

##
## 2. Combine results
##
PRCC_vals <- data.frame(
    names = param_names,
    V_max = prcc_values$V_max$PRCC[, 1],
    F_U_max = prcc_values$F_U_max$PRCC[, 1],
    F_B_max = prcc_values$F_B_max$PRCC[, 1]
)

PRCC_times <- data.frame(
    names = param_names,
    t_V_max = prcc_times$t_V_max$PRCC[, 1],
    t_F_U_max = prcc_times$t_F_U_max$PRCC[, 1],
    t_F_B_max = prcc_times$t_F_B_max$PRCC[, 1]
)

##
## 3. Save PRCC results for plotting
##

# Compute summed PRCC (values + times) for global summary
PRCC_vals$sum_abs <- rowSums(abs(PRCC_vals[, c("V_max", "F_U_max", "F_B_max")]))
PRCC_times$sum_abs <- rowSums(abs(PRCC_times[, c("t_V_max", "t_F_U_max", "t_F_B_max")]))

PRCC_global <- merge(
    PRCC_vals[, c("names", "sum_abs")],
    PRCC_times[, c("names", "sum_abs")],
    by = "names"
)
PRCC_global$sum_total <- PRCC_global$sum_abs.x + PRCC_global$sum_abs.y
PRCC_global <- PRCC_global[order(-PRCC_global$sum_total), ]

# Note: Exporting as _withV0 suffixes to avoid overwriting the baseline files
write.csv(PRCC_global, file.path(OUTPUT_dir, "PRCC_global_summary_withV0.csv"), row.names = FALSE)
qs_save(list(PRCC_vals = PRCC_vals, PRCC_times = PRCC_times), file.path(OUTPUT_dir, "PRCC_results_withV0.qs"))

cat("PRCC computation finished. Results saved in OUTPUT.\n")
