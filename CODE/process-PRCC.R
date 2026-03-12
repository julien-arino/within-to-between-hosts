# ============================================================
# File : process-PRCC.R
# Computes PRCC from simulation results and saves outputs
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(sensitivity)
    library(qs)
})

OUTPUT_clo <- file.path(getwd(), "OUTPUT_clo")

# Load exactly what was exported from the simulation phase
sim_data <- qread(file.path(OUTPUT_clo, "PRCC_sim_results.qs"))
X <- sim_data$X
res_df <- sim_data$res_df
param_names <- sim_data$param_names

# ------------------------------------------------------------
# 1. PRCC computation
# ------------------------------------------------------------
prcc_values <- list(
    V_max   = pcc(X, res_df$V_max, rank = TRUE),
    F_U_max = pcc(X, res_df$F_U_max, rank = TRUE),
    F_B_max = pcc(X, res_df$F_B_max, rank = TRUE)
)
prcc_times <- list(
    t_V_max   = pcc(X, res_df$t_V_max, rank = TRUE),
    t_F_U_max = pcc(X, res_df$t_F_U_max, rank = TRUE),
    t_F_B_max = pcc(X, res_df$t_F_B_max, rank = TRUE)
)

# ------------------------------------------------------------
# 2. Combine results
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 3. Save PRCC results for plotting
# ------------------------------------------------------------

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

write.csv(PRCC_global, file.path(OUTPUT_clo, "PRCC_global_summary.csv"), row.names = FALSE)
qsave(list(PRCC_vals = PRCC_vals, PRCC_times = PRCC_times), file.path(OUTPUT_clo, "PRCC_results.qs"))

cat("✅ PRCC computation finished. Results saved in OUTPUT_clo.\n")
