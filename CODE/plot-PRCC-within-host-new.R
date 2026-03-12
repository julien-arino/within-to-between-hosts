# ============================================================
# File : plot-PRCC-within-host-new.R
# category-based plots for PRCC values
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(latex2exp)
    library(qs)
})

OUTPUT_clo <- file.path(getwd(), "OUTPUT_clo")
FIGS_clo <- file.path(getwd(), "FIGS_clo")

dir.create(FIGS_clo, showWarnings = FALSE, recursive = TRUE)

# Load data
prcc_data <- qread(file.path(OUTPUT_clo, "PRCC_results.qs"))
PRCC_vals <- prcc_data$PRCC_vals
PRCC_times <- prcc_data$PRCC_times
PRCC_global <- read.csv(file.path(OUTPUT_clo, "PRCC_global_summary.csv"))

# ------------------------------------------------------------
# 8. Plot configuration
# ------------------------------------------------------------
cat_cell <- c("lambda_S", "S_max", "d_I", "d_D")
cat_viral <- c("beta", "d_V", "p")
cat_immune <- c("k_U_F", "p_FI", "psi_F_prod", "k_B_F", "k_lin_f", "k_int_f", "epsilon_FI", "eta_FI", "c_star")

cat_bg_colors <- c(Cell = "#800080", Viral = "#A9A9A9", Immune = "#F5F5F5")
bg_alpha <- 0.15

x_labels <- c(
    "beta"       = TeX("$\\beta$"),
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
    "c_star"     = TeX("$c^*$")
)

indicator_colors_vals <- c("V_max" = "#D7263D", "F_U_max" = "#1B9E77", "F_B_max" = "#0072B2")
indicator_colors_times <- c("t_V_max" = "#D7263D", "t_F_U_max" = "#1B9E77", "t_F_B_max" = "#0072B2")

# ------------------------------------------------------------
# 9. Values plot
# ------------------------------------------------------------
PRCC_vals_m <- melt(PRCC_vals, id.vars = "names")
param_ordered <- PRCC_vals_m %>%
    group_by(names) %>%
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
            ifelse(param_ordered %in% cat_immune, "Immune", "Immune")
        )
    )
)

PRCC_vals_m$names <- factor(PRCC_vals_m$names, levels = param_ordered)

P_max <- ggplot() +
    geom_rect(
        data = param_bg_vals, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
        alpha = bg_alpha, color = NA, show.legend = TRUE
    ) +
    scale_fill_manual(values = cat_bg_colors, name = "Category") +
    geom_point(
        data = PRCC_vals_m,
        aes(x = names, y = abs(value), color = variable, size = abs(value)),
        shape = 19, alpha = 0.7, position = position_jitter(width = 0.04, height = 0)
    ) +
    scale_color_manual(
        values = indicator_colors_vals,
        labels = c(TeX("$V_{max}$"), TeX("$F_{U_{max}}$"), TeX("$F_{B_{max}}$")),
        name = "Indicator"
    ) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    labs(x = "Parameter", y = "Partial Rank Correlation Coefficient") +
    scale_size_continuous(name = "|PRCC|", range = c(2.5, 8)) +
    theme(legend.position = "right", panel.grid.major.x = element_blank())

print(P_max)
ggsave(
    filename = file.path(FIGS_clo, "PRCC_values.png"),
    plot = P_max, width = 25, height = 15, units = "cm", dpi = 300
)

# ------------------------------------------------------------
# 10. Time-to-max plot
# ------------------------------------------------------------
PRCC_times_m <- melt(PRCC_times, id.vars = "names")
param_ordered_t <- PRCC_times_m %>%
    group_by(names) %>%
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
            ifelse(param_ordered_t %in% cat_immune, "Immune", "Immune")
        )
    )
)

PRCC_times_m$names <- factor(PRCC_times_m$names, levels = param_ordered_t)

P_tmax <- ggplot() +
    geom_rect(
        data = param_bg_times, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
        alpha = bg_alpha, color = NA, show.legend = TRUE
    ) +
    scale_fill_manual(values = cat_bg_colors, name = "Category") +
    geom_point(
        data = PRCC_times_m,
        aes(x = names, y = abs(value), color = variable, size = abs(value)),
        shape = 19, alpha = 0.7, position = position_jitter(width = 0.04, height = 0)
    ) +
    scale_color_manual(
        values = indicator_colors_times,
        labels = c(TeX("$t_{V_{max}}$"), TeX("$t_{F_{U_{max}}}$"), TeX("$t_{F_{B_{max}}}$")),
        name = "Indicator"
    ) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    labs(x = "Parameter", y = "PRCC of time to max") +
    scale_size_continuous(name = "|PRCC|", range = c(2.5, 8)) +
    theme(legend.position = "right", panel.grid.major.x = element_blank())

print(P_tmax)
ggsave(
    filename = file.path(FIGS_clo, "PRCC_values_times.png"),
    plot = P_tmax, width = 25, height = 15, units = "cm", dpi = 300
)

# ------------------------------------------------------------
# 10b. Global summed PRCC plot (values + times)
# ------------------------------------------------------------
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

cat("\n✅ Global PRCC plot saved:", file.path(FIGS_clo, "PRCC_global.png"), "\n")
