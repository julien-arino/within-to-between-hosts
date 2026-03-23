#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-E2-R0P-vs-R0within.R
# Description:
#   Plots the difference between R0_P2P (R_0i^P) and R0_within 
#   against time to maximum tissue damage (tau_max_Psi).
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-E2-R0P-vs-R0within.R")

load_libraries(c("ggplot2", "dplyr", "qs2", "tidyr", "patchwork"))

# 1a) Find latest processed status file (base file)
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*\\.qs$", full.names = TRUE)
status_files <- status_files[!grepl("_xic_|_xir_", status_files)]
if (length(status_files) == 0) {
  stop("Cannot find cohort_status qs file in ", output_dir)
}
STATUS_FILE <- status_files[which.max(file.mtime(status_files))]

# 1b) Find latest cohort_times file
times_files <- list.files(output_dir, pattern = "^cohort_times_P.*\\.qs$", full.names = TRUE)
if (length(times_files) == 0) {
  stop("Cannot find cohort_times qs file in ", output_dir)
}
TIMES_FILE <- times_files[which.max(file.mtime(times_files))]

cat("Loading status:", basename(STATUS_FILE), "\n")
status_df <- qs_read(STATUS_FILE, nthreads = N_QS_THREADS)

cat("Loading times:", basename(TIMES_FILE), "\n")
times_df <- qs_read(TIMES_FILE, nthreads = N_QS_THREADS)

# Join by ID
summary_df <- status_df %>%
  inner_join(times_df %>% select(ID, R0_P2P), by = "ID")

# Formatting columns for the plot
summary_df <- summary_df %>%
  mutate(
    R0 = R0_within,
    t_Psi_max = tau_max_Psi,
    R0_diff = (R0_P2P - R0_within) / R0_within
  ) %>%
  filter(!is.na(R0) & !is.na(R0_P2P)) %>%
  filter(t_Psi_max < 100)

summary_df$status <- factor(summary_df$status, levels = c("Mild", "ICU", "Dead"))

status_colors <- c(
  Mild = "dodgerblue4",
  ICU  = "orange",
  Dead = "red"
)

# ---------------------------
# Plotting
# ---------------------------
p2 <- ggplot(summary_df, aes(x = t_Psi_max, y = R0_diff)) +

  # Draw a zero-line across the difference axis
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +

  # Patients with R0 >= 1 (colored by status)
  geom_point(
    data = subset(summary_df, R0 >= 1),
    aes(color = status),
    alpha = 0.18,
    size = 0.6
  ) +

  # Patients with R0 < 1 (forced explicitly to DARK GREEN)
  geom_point(
    data = subset(summary_df, R0 < 1),
    color = "darkgreen",
    alpha = 0.4,
    size = 0.8
  ) +

  scale_color_manual(
    values = status_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(-1000, -100, -10, 0, 10, 100, 1000, 10000)
  ) +
  coord_cartesian(
    xlim = c(0, 100)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    x = expression("Time to maximum tissue damage" ~ tau[i]^{Psi[max]} ~ "(days)"),
    y = expression("Normalized Difference" ~ (R["0,i"]^P - R["0,within"]) / R["0,within"])
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 2, alpha = 1)
    )
  )

print(p2)

# ---------------------------
# Save figure
# ---------------------------
out_png <- file.path(figs_dir, "Figure-E2-R0P-vs-R0within.png")
out_pdf <- file.path(figs_dir, "Figure-E2-R0P-vs-R0within.pdf")

ggsave(filename = out_png, plot = p2, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(filename = out_pdf, plot = p2, width = 25, height = 15, units = "cm", dpi = 300, device = cairo_pdf)
cat("\nFigure E2 saved to:\n  -", out_png, "\n  -", out_pdf, "\n")

cat("\n--- Statistics ---\n")
probs <- seq(0, 1, by = 0.1)
perc_names <- paste0(seq(0, 100, by = 10), "%")

stats_df <- data.frame(
  Metric = c("R0_within", "R0_P2P"),
  Mean = c(
    mean(summary_df$R0_within, na.rm = TRUE),
    mean(summary_df$R0_P2P, na.rm = TRUE)
  )
)

perc_within <- quantile(summary_df$R0_within, probs = probs, na.rm = TRUE)
perc_P2P    <- quantile(summary_df$R0_P2P, probs = probs, na.rm = TRUE)

for (i in seq_along(probs)) {
  stats_df[[perc_names[i]]] <- c(perc_within[i], perc_P2P[i])
}

print(stats_df, row.names = FALSE)
cat("------------------\n")

print_end_time(start_time, "plot-Figure-E2-R0P-vs-R0within.R")
