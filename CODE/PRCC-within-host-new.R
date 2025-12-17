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
source("functions-wb.R")  # contains extract_max_indicators()

OUTPUT_clo <- file.path(getwd(), "OUTPUT_clo")
FIGS_clo   <- file.path(getwd(), "FIGS_clo")

dir.create(OUTPUT_clo, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGS_clo,   showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# 2. Parameter ranges (list format)
# ------------------------------------------------------------
param_ranges <- list(
  beta_V       = c(0.1, 0.5),
  d_I          = c(0.05, 0.15),
  d_D          = c(0.001, 0.015),
  d_V          = c(5, 15),
  epsilon_FI   = c(1e-5, 1e-3),
  eta_FI       = c(0.001, 0.05),
  k_B_F        = c(0.001, 0.05),
  k_int_f      = c(10, 20),
  k_lin_f      = c(10, 20),
  k_U_F        = c(2, 10),
  lambda_S     = c(0.5, 1),
  p            = c(100, 800),
  p_FI         = c(0.8, 4),
  psi_F_prod   = c(0.1, 0.5),
  S_max        = c(0.1, 0.2),
  T_star       = c(1e-5, 1e-3)
  #A_F          = c(0.5e-23, 2.5e-23)
)

param_names <- names(param_ranges)
p <- length(param_ranges)

# ------------------------------------------------------------
# 3. Sampling (Latin Hypercube)
# ------------------------------------------------------------
set.seed(123)
N <- 1000000  # Increase if needed
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
  if (i %% 100 == 0)
    message(sprintf("[worker %d] simulation %d/%d", pid, i, nrow(Xdf)))
  pars <- as.list(Xdf[i, , drop = FALSE])
  extract_max_indicators_ode(pars)
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
# 6. PRCC computation
# ------------------------------------------------------------
prcc_values <- list(
  V_max   = pcc(X, res_df$V_max,   rank = TRUE),
  F_U_max = pcc(X, res_df$F_U_max, rank = TRUE),
  F_B_max = pcc(X, res_df$F_B_max, rank = TRUE)
)
prcc_times <- list(
  t_V_max   = pcc(X, res_df$t_V_max,   rank = TRUE),
  t_F_U_max = pcc(X, res_df$t_F_U_max, rank = TRUE),
  t_F_B_max = pcc(X, res_df$t_F_B_max, rank = TRUE)
)

# ------------------------------------------------------------
# 7. Combine results
# ------------------------------------------------------------
PRCC_vals <- data.frame(
  names = param_names,
  V_max   = prcc_values$V_max$PRCC[, 1],
  F_U_max = prcc_values$F_U_max$PRCC[, 1],
  F_B_max = prcc_values$F_B_max$PRCC[, 1]
)

PRCC_times <- data.frame(
  names = param_names,
  t_V_max   = prcc_times$t_V_max$PRCC[, 1],
  t_F_U_max = prcc_times$t_F_U_max$PRCC[, 1],
  t_F_B_max = prcc_times$t_F_B_max$PRCC[, 1]
)




# ------------------------------------------------------------
# 8. Plot configuration
# ------------------------------------------------------------
cat_cell   <- c("lambda_S", "S_max", "d_I", "d_D")
cat_viral  <- c("beta_V", "d_V", "p")
cat_immune <- c("k_U_F", "p_FI", "psi_F_prod", "k_B_F", "k_lin_f", "k_int_f", "epsilon_FI", "eta_FI", "T_star")

cat_bg_colors <- c(Cell = "#800080", Viral = "#A9A9A9", Immune = "#F5F5F5")
bg_alpha <- 0.15

x_labels <- c(
  "beta_V"     = TeX("$\\beta_V$"),
  "d_I"        = TeX("$d_I$"),
  "d_D"        = TeX("$d_D$"),
  "d_V"        = TeX("$d_V$"),
  "epsilon_FI" = TeX("$\\epsilon_{FI}$"),
  "eta_FI"     = TeX("$\\eta_{FI}$"),
  "k_B_F"      = TeX("$k_{BF}$"),
  "k_int_f"    = TeX("$k_{int_f}$"),
  "k_lin_f"    = TeX("$k_{lin_f}$"),
  "k_U_F"      = TeX("$k_{UF}$"),
  "lambda_S"   = TeX("$\\lambda_S$"),
  "p"          = TeX("$p$"),
  "p_FI"       = TeX("$p_{FI}$"),
  "psi_F_prod" = TeX("$\\psi_{prod}$"),
  "S_max"      = TeX("$S_{max}$"),
  "T_star"     = TeX("$T^*$")
 # "A_F"        = TeX("$A_F$")
)

indicator_colors_vals  <- c("V_max" = "#D7263D", "F_U_max" = "#1B9E77", "F_B_max" = "#0072B2")
indicator_colors_times <- c("t_V_max" = "#D7263D", "t_F_U_max" = "#1B9E77", "t_F_B_max" = "#0072B2")

# ------------------------------------------------------------
# 9. Values plot
# ------------------------------------------------------------
PRCC_vals_m <- melt(PRCC_vals, id.vars = "names")
param_ordered <- PRCC_vals_m %>% group_by(names) %>%
  summarise(max_abs = max(abs(value))) %>%
  arrange(desc(max_abs)) %>%
  pull(names)

param_bg_vals <- data.frame(
  name = param_ordered,
  xmin = seq_along(param_ordered) - 0.5,
  xmax = seq_along(param_ordered) + 0.5,
  ymin = 0, ymax = 1,
  category = ifelse(param_ordered %in% cat_cell, "Cell",
                    ifelse(param_ordered %in% cat_viral, "Viral",
                           ifelse(param_ordered %in% cat_immune, "Immune", "Immune")))
)

PRCC_vals_m$names <- factor(PRCC_vals_m$names, levels = param_ordered)

P_max <- ggplot() +
  geom_rect(data = param_bg_vals, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
            alpha = bg_alpha, color = NA, show.legend = TRUE) +
  scale_fill_manual(values = cat_bg_colors, name = "Category") +
  geom_point(data = PRCC_vals_m,
             aes(x = names, y = abs(value), color = variable, size = abs(value)),
             shape = 19, alpha = 0.7, position = position_jitter(width = 0.04, height = 0)) +
  scale_color_manual(values = indicator_colors_vals,
                     labels = c(TeX("$V_{max}$"), TeX("$F_{U_{max}}$"), TeX("$F_{B_{max}}$")),
                     name = "Indicator") +
  scale_x_discrete(labels = x_labels) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(x = "Parameter", y = "Partial Rank Correlation Coefficient") +
  scale_size_continuous(name = "|PRCC|", range = c(2.5, 8)) +
  theme(legend.position = "right", panel.grid.major.x = element_blank())

print(P_max)
ggsave(filename = file.path(FIGS_clo, "PRCC_values.png"),  
       plot = P_max, width = 25, height = 15, units = "cm", dpi = 300)

# ------------------------------------------------------------
# 10. Time-to-max plot
# ------------------------------------------------------------
PRCC_times_m <- melt(PRCC_times, id.vars = "names")
param_ordered_t <- PRCC_times_m %>% group_by(names) %>%
  summarise(max_abs = max(abs(value))) %>%
  arrange(desc(max_abs)) %>%
  pull(names)

param_bg_times <- data.frame(
  name = param_ordered_t,
  xmin = seq_along(param_ordered_t) - 0.5,
  xmax = seq_along(param_ordered_t) + 0.5,
  ymin = 0, ymax = 1,
  category = ifelse(param_ordered_t %in% cat_cell, "Cell",
                    ifelse(param_ordered_t %in% cat_viral, "Viral",
                           ifelse(param_ordered_t %in% cat_immune, "Immune", "Immune")))
)

PRCC_times_m$names <- factor(PRCC_times_m$names, levels = param_ordered_t)

P_tmax <- ggplot() +
  geom_rect(data = param_bg_times, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
            alpha = bg_alpha, color = NA, show.legend = TRUE) +
  scale_fill_manual(values = cat_bg_colors, name = "Category") +
  geom_point(data = PRCC_times_m,
             aes(x = names, y = abs(value), color = variable, size = abs(value)),
             shape = 19, alpha = 0.7, position = position_jitter(width = 0.04, height = 0)) +
  scale_color_manual(values = indicator_colors_times,
                     labels = c(TeX("$t_{V_{max}}$"), TeX("$t_{F_{U_{max}}}$"), TeX("$t_{F_{B_{max}}}$")),
                     name = "Indicator") +
  scale_x_discrete(labels = x_labels) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(x = "Parameter", y = "PRCC of time to max") +
  scale_size_continuous(name = "|PRCC|", range = c(2.5, 8)) +
  theme(legend.position = "right", panel.grid.major.x = element_blank())

print(P_tmax)
ggsave(filename = file.path(FIGS_clo, "PRCC_values_times.png"),     
       plot = P_tmax, width = 25, height = 15, units = "cm", dpi = 300)



# ------------------------------------------------------------
# 10b. Global summed PRCC plot (values + times)
# ------------------------------------------------------------

# Compute summed PRCC (values + times)
PRCC_vals$sum_abs  <- rowSums(abs(PRCC_vals[, c("V_max", "F_U_max", "F_B_max")]))
PRCC_times$sum_abs <- rowSums(abs(PRCC_times[, c("t_V_max", "t_F_U_max", "t_F_B_max")]))

PRCC_global <- merge(
  PRCC_vals[, c("names", "sum_abs")],
  PRCC_times[, c("names", "sum_abs")],
  by = "names"
)
PRCC_global$sum_total <- PRCC_global$sum_abs.x + PRCC_global$sum_abs.y
PRCC_global <- PRCC_global[order(-PRCC_global$sum_total), ]

# Plot simple global influence figure (LaTeX labels)
P_global <- ggplot(PRCC_global, aes(x = reorder(names, -sum_total), y = sum_total)) +
  geom_line(group = 1, color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +
  scale_x_discrete(labels = x_labels) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Parameter",
    y = TeX("Global $\\Sigma|PRCC|$ (values + times)"),
    title = "Global Influence of Parameters"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(P_global)

ggsave(
  filename = file.path(FIGS_clo, "PRCC_global.png"),
  plot = P_global,
  width = 25, height = 15, units = "cm", dpi = 300
)

write.csv(PRCC_global, file.path(OUTPUT_clo, "PRCC_global_summary.csv"), row.names = FALSE)
cat("\n✅ Global PRCC plot saved:", file.path(FIGS_clo, "PRCC_global.pdf"), "\n")


# ------------------------------------------------------------
# 11. Clean-up
# ------------------------------------------------------------
future::plan(sequential)
cat("✅ Parallel PRCC finished. Plots saved.\n")
