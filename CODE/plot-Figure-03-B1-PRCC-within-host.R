#!/usr/bin/env Rscript
# ============================================================
# File : plot-PRCC-within-host-new.R
# category-based plots for PRCC values
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(latex2exp)
    library(qs2)
    library(here)
})

# Set project root automatically relative to the .git tracking directory
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

OUTPUT <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  source(file.path(project_dir, "CODE", "functions-and-definitions.R"))
}
FIGS <- file.path(project_dir, "FIGS")

dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

# Load data
prcc_result_files <- list.files(OUTPUT, pattern = "^sensitivity_P.*-PRCC-results\\.qs$", full.names = TRUE)
if (length(prcc_result_files) == 0) {
    stop("Cannot find PRCC results qs file in ", OUTPUT)
}
LATEST_PRCC_FILE <- prcc_result_files[which.max(file.mtime(prcc_result_files))]
cat("Loading PRCC data from:", basename(LATEST_PRCC_FILE), "\n")

prcc_data <- qs_read(LATEST_PRCC_FILE, nthreads = N_QS_THREADS)
PRCC_vals <- prcc_data$PRCC_vals
PRCC_times <- prcc_data$PRCC_times

PRCC_summary_file <- sub("-PRCC-results\\.qs$", "-PRCC-summary.csv", LATEST_PRCC_FILE)
PRCC_global <- read.csv(PRCC_summary_file)

# ------------------------------------------------------------
# 8. Plot configuration
# ------------------------------------------------------------
cat_cell <- c("λ_S", "S_max", "d_I", "d_D")
cat_viral <- c("β_V", "d_V", "p", "V0")
cat_immune <- c("k_U_F", "p_FI", "ψ_F_prod", "k_B_F", "k_lin_f", "k_int_f", "ε_FI", "η_FI", "c_star", "a_F")

cat_bg_colors <- c(Cell = "#800080", Viral = "#A9A9A9", Immune = "#F5F5F5")
bg_alpha <- 0.15

x_labels <- c(
    "β_V"        = "beta[V]",
    "d_I"        = "d[I]",
    "d_D"        = "d[D]",
    "d_V"        = "d[V]",
    "ε_FI"       = "epsilon[FI]",
    "η_FI"       = "eta[FI]",
    "k_B_F"      = "k[BF]",
    "k_int_f"    = "k[int[f]]",
    "k_lin_f"    = "k[lin[f]]",
    "k_U_F"      = "k[UF]",
    "λ_S"        = "lambda[S]",
    "p"          = "p",
    "p_FI"       = "p[FI]",
    "ψ_F_prod"   = "psi[prod]",
    "S_max"      = "S[max]",
    "c_star"     = "c^'*'",
    "a_F"        = "a[F]",
    "V0"         = "V[0]"
)

indicator_colors_vals <- c("V_max" = "#D7263D", "F_U_max" = "#1B9E77", "F_B_max" = "#0072B2")
indicator_colors_times <- c("t_V_max" = "#D7263D", "t_F_U_max" = "#1B9E77", "t_F_B_max" = "#0072B2")

# ------------------------------------------------------------
# 9. Values plot
# ------------------------------------------------------------
PRCC_vals_m <- melt(PRCC_vals, id.vars = "names") %>%
    filter(variable %in% names(indicator_colors_vals))
param_ordered <- PRCC_vals_m %>%
    group_by(names) %>%
    summarise(max_abs = max(abs(value), na.rm = TRUE)) %>%
    arrange(desc(max_abs)) %>%
    pull(names)

param_bg_vals <- data.frame(
    name = param_ordered,
    xmin = seq_along(param_ordered) - 0.5,
    xmax = seq_along(param_ordered) + 0.5,
    ymin = -1, ymax = 1,
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
        aes(x = names, y = value, color = variable, size = abs(value)),
        shape = 19, alpha = 0.7, position = position_jitter(width = 0.04, height = 0)
    ) +
    scale_color_manual(
        values = indicator_colors_vals,
        labels = c(parse(text = "V[max]"), parse(text = "F[U[max]]"), parse(text = "F[B[max]]")),
        name = "Indicator", na.translate = FALSE
    ) +
    scale_x_discrete(labels = function(x) parse(text = x_labels[x])) +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_minimal(base_size = 14) +
    labs(x = "Parameter", y = "Partial Rank Correlation Coefficient") +
    scale_size_continuous(name = "|PRCC|", range = c(2.5, 8), guide = "none") +
    theme(legend.position = "right", panel.grid.major.x = element_blank())

print(P_max)
ggsave(
    filename = file.path(FIGS, "Figure-B1a-PRCC-values.png"),
    plot = P_max, width = 25, height = 15, units = "cm", dpi = 300
)
ggsave(
    filename = file.path(FIGS, "Figure-B1a-PRCC-values.pdf"),
    plot = P_max, width = 25, height = 15, units = "cm", dpi = 300, device = cairo_pdf
)

# ------------------------------------------------------------
# 10. Time-to-max plot
# ------------------------------------------------------------
PRCC_times_m <- melt(PRCC_times, id.vars = "names") %>%
    filter(variable %in% names(indicator_colors_times))
param_ordered_t <- PRCC_times_m %>%
    group_by(names) %>%
    summarise(max_abs = max(abs(value), na.rm = TRUE)) %>%
    arrange(desc(max_abs)) %>%
    pull(names)

param_bg_times <- data.frame(
    name = param_ordered_t,
    xmin = seq_along(param_ordered_t) - 0.5,
    xmax = seq_along(param_ordered_t) + 0.5,
    ymin = -1, ymax = 1,
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
        aes(x = names, y = value, color = variable, size = abs(value)),
        shape = 19, alpha = 0.7, position = position_jitter(width = 0.04, height = 0)
    ) +
    scale_color_manual(
        values = indicator_colors_times,
        labels = c(parse(text = "t[V[max]]"), parse(text = "t[F[U[max]]]"), parse(text = "t[F[B[max]]]")),
        name = "Indicator", na.translate = FALSE
    ) +
    scale_x_discrete(labels = function(x) parse(text = x_labels[x])) +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_minimal(base_size = 14) +
    labs(x = "Parameter", y = "PRCC of time to max") +
    scale_size_continuous(name = "|PRCC|", range = c(2.5, 8), guide = "none") +
    theme(legend.position = "right", panel.grid.major.x = element_blank())

print(P_tmax)
ggsave(
    filename = file.path(FIGS, "Figure-B1b-PRCC-values-times.png"),
    plot = P_tmax, width = 25, height = 15, units = "cm", dpi = 300
)
ggsave(
    filename = file.path(FIGS, "Figure-B1b-PRCC-values-times.pdf"),
    plot = P_tmax, width = 25, height = 15, units = "cm", dpi = 300, device = cairo_pdf
)

# ------------------------------------------------------------
# 10b. Global summed PRCC plot (values + times)
# ------------------------------------------------------------
param_ordered_g <- PRCC_global$names

# Background categories (same style as others)
param_bg_global <- data.frame(
    name = param_ordered_g,
    xmin = seq_along(param_ordered_g) - 0.5,
    xmax = seq_along(param_ordered_g) + 0.5,
    ymin = 0,
    ymax = max(PRCC_global$sum_total, na.rm = TRUE) * 1.05,
    category = ifelse(param_ordered_g %in% cat_cell, "Cell",
        ifelse(param_ordered_g %in% cat_viral, "Viral",
            ifelse(param_ordered_g %in% cat_immune, "Immune", "Immune")
        )
    )
)

PRCC_global$names <- factor(PRCC_global$names, levels = param_ordered_g)

P_global <- ggplot() +
    geom_rect(
        data = param_bg_global, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
        alpha = bg_alpha, color = NA, show.legend = TRUE
    ) +
    scale_fill_manual(values = cat_bg_colors, name = "Category") +
    geom_point(
        data = PRCC_global,
        aes(x = names, y = sum_total, size = sum_total),
        color = "blue"
    ) +
    scale_size_continuous(name = TeX("$\\Sigma|PRCC|$"), range = c(2.5, 8), guide = "none") +
    scale_x_discrete(labels = function(x) parse(text = x_labels[x])) +
    theme_minimal(base_size = 14) +
    labs(
        x = "Parameter",
        y = parse(text = '"Global "*Sigma*"|PRCC| (values + times)"')
    ) +
    theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold")
    )

print(P_global)

ggsave(
    filename = file.path(FIGS, "Figure-03c-PRCC-global.png"),
    plot = P_global,
    width = 25, height = 15, units = "cm", dpi = 300
)
ggsave(
    filename = file.path(FIGS, "Figure-03c-PRCC-global.pdf"),
    plot = P_global,
    width = 25, height = 15, units = "cm", dpi = 300, device = cairo_pdf
)

cat("\nFigures saved in", FIGS, ":\n")
cat("  - Figure-B1a-PRCC-values (.png, .pdf)\n")
cat("  - Figure-B1b-PRCC-values-times (.png, .pdf)\n")
cat("  - Figure-03c-PRCC-global (.png, .pdf)\n\n")
