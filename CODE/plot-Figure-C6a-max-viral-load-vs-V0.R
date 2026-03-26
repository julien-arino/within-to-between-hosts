#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-C6a-max-viral-load-vs-V0.R
# Description:
#   Plot maximum viral load per individual against initial viral load (V0), colored by status
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-C6a-max-viral-load-vs-V0.R")

load_libraries(c("ggplot2", "dplyr", "qs2", "here"))

# ------------------------------------------------------------
# 1. Automatically find the latest cohort files
# ------------------------------------------------------------
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xid_[0-9]+\\.qs$", full.names = TRUE)
if (length(status_files) == 0) {
  stop("No base cohort_status file found in ", output_dir)
}
latest_status_file <- status_files[which.max(file.mtime(status_files))]
cat("Loading newest cohort status results:", basename(latest_status_file), "\n")
cohort_df <- qs_read(latest_status_file, nthreads = N_QS_THREADS)

param_files <- list.files(output_dir, pattern = "^cohort_sim_parameters_P.*\\.qs$", full.names = TRUE)
if (length(param_files) == 0) {
  stop("No base cohort parameters file found in ", output_dir)
}
latest_param_file <- param_files[which.max(file.mtime(param_files))]
cat("Loading newest cohort parameters results:", basename(latest_param_file), "\n")
params_df <- qs_read(latest_param_file, nthreads = N_QS_THREADS)

cat("Files loaded. Joining on ID...\n")

# ------------------------------------------------------------
# 2. Extract and join data
# ------------------------------------------------------------
viral_max_df <- cohort_df %>%
  select(ID, status, V_max = max_V)

v0_df <- params_df %>%
  select(ID, V0)

plot_df <- inner_join(viral_max_df, v0_df, by = "ID")

cat("Joined dataframe created for", nrow(plot_df), "individuals\n")

# ------------------------------------------------------------
# 3. Define custom colors and factor levels
# ------------------------------------------------------------
status_colors <- c(Mild = "dodgerblue4", ICU = "orange", Dead = "red")

# Update factor levels for status
plot_df$status <- factor(plot_df$status, levels = c("Mild", "ICU", "Dead"))

# ------------------------------------------------------------
# 4. Create the plot
# ------------------------------------------------------------
p_scatter <- ggplot(plot_df, aes(x = V0, y = V_max, color = status)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_color_manual(values = status_colors) +
  labs(
    x = expression("log"[10]*"(Initial Viral load (V"[0]*") [copies/ml])"),
    y = expression("log"[10]*"(Max Viral load [copies/ml])"),
    color = "Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top"
  )

# ------------------------------------------------------------
# 5. Print the plot to screen (or save directly)
# ------------------------------------------------------------
print(p_scatter)

# ------------------------------------------------------------
# 6. Save plot to PDF and PNG
# ------------------------------------------------------------
out_pdf <- file.path(figs_dir, "Figure-C6a-max-viral-load-vs-V0.pdf")
out_png <- file.path(figs_dir, "Figure-C6a-max-viral-load-vs-V0.png")

ggsave(
  filename = out_png, 
  plot = p_scatter, width = 9, height = 7, units = "in", dpi = 300
)
ggsave(
  filename = out_pdf, 
  plot = p_scatter, width = 9, height = 7, units = "in", dpi = 300
)

cat("Figures saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
print_end_time(start_time, "plot-Figure-C6a-max-viral-load-vs-V0.R")
