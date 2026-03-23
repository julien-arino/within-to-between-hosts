#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-C4-time-evolution-quantities.R
# Description:
#   File: plot-Figure-C4-time-evolution-quantities.R
#   6-panel plot showing time evolution of V, beta, I, Psi, FB, FU
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-C4-time-evolution-quantities.R")

load_libraries(c("qs2", "dplyr", "ggplot2", "patchwork", "here", "future.apply"))

# Load INTERPOLATED sim_state
# Automatically find the latest pre-computed aggregate file
files <- list.files(output_dir, pattern = "^cohort_FigC4_P.*\\.qs$", full.names = TRUE)
if (length(files) == 0) stop("No pre-computed cohort_FigC4_P...qs file found. Please run process-cohort-prepare-Figure-C4.R first.")
latest_file <- files[which.max(file.mtime(files))]

cat("Loading aggregated Figure C4 matrix:", basename(latest_file), "\n")
mean_df <- qs_read(latest_file, nthreads = N_QS_THREADS)

status_cols <- c(
  Mild = "dodgerblue4",
  ICU  = "orange",
  Dead = "red"
)

# Helper function
make_plot <- function(y_mean, y_sd, ylab, xlim = NULL, ylim = NULL, legend = FALSE) {
  ggplot(mean_df, aes(
    x = time,
    y = !!sym(y_mean),
    color = status,
    fill = status
  )) +
    geom_ribbon(
      aes(
        ymin = pmax(!!sym(y_mean) - !!sym(y_sd), 0),
        ymax = !!sym(y_mean) + !!sym(y_sd)
      ),
      alpha = 0.22, color = NA
    ) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = status_cols) +
    scale_fill_manual(values = status_cols) +
    labs(x = "Time (days)", y = ylab) +
    {
      if (!is.null(ylim)) coord_cartesian(ylim = ylim) else NULL
    } +
    {
      if (!is.null(xlim)) coord_cartesian(xlim = xlim)
    } +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = if (legend) "top" else "none",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Create the 6 panels
p_V <- make_plot(
  "V_mean", "V_sd",
  expression("Viral load" ~ log[10] * "(copies/ml)")
)

p_beta <- make_plot("beta_mean", "beta_sd",
  expression("Transmission rate"),
  xlim = c(0, 30)
)

p_I <- make_plot(
  "I_mean", "I_sd",
  expression("Infected cells" ~ 10^9 * " cells/ml")
)

p_Psi <- make_plot("Psi_mean", "Psi_sd",
  expression("Lung tissue damage" ~ Psi(t) ~ "%"),
  ylim = c(0, 100),
  legend = TRUE
)

p_FB <- make_plot(
  "F_B_mean", "F_B_sd",
  expression("Bound IFN" ~ 10^{
    -5
  } ~ "pg/ml")
)

p_FU <- make_plot(
  "F_U_mean", "F_U_sd",
  expression("Unbound IFN (pg/ml)")
)

# Combine 2 rows x 3 columns
combined_all <-
  (p_Psi | p_V | p_beta) /
    (p_I | p_FB | p_FU)

combined_all <- combined_all +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

suppressWarnings(print(combined_all))

# Save single-page figure

out_pdf <- file.path(figs_dir, "Figure-C4-time-evolution-quantities.pdf")
out_png <- file.path(figs_dir, "Figure-C4-time-evolution-quantities.png")

suppressWarnings({
  ggsave(
    out_pdf,
    plot = combined_all,
    width = 26, height = 16, units = "cm", dpi = 300
  )

  ggsave(
    out_png,
    plot = combined_all,
    width = 26, height = 16, units = "cm", dpi = 300
  )
})

cat("6-panel evolution plot saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
print_end_time(start_time, "plot-Figure-C4-time-evolution-quantities.R")
