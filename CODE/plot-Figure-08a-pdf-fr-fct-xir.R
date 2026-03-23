#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-08a-pdf-fr-fct-xic.R
# Description:
# ============================================================

# First things first: locate project directory and load helper functions
# and constants
project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

# Say what we are running and start the clock
start_time <- start_time_and_hello("plot-Figure-08a-pdf-fr-fct-xic.R")

# Load libraries
load_libraries(c("ggplot2", "dplyr", "qs2", "tidyr"))

# Load f_r from distribution files
xir_vals <- c(1, 4)
data_list <- list()

for (xir in xir_vals) {
  dist_pattern <- paste0("^cohort_distribution_filters_P.*_xih_75_xid_85_xic_4_xir_", xir, "\\.qs$")
  dist_file <- list.files(output_dir, pattern = dist_pattern, full.names = TRUE)
  if (length(dist_file) == 0) stop(paste("No distribution file found for xir =", xir))

  latest_file <- dist_file[which.max(file.mtime(dist_file))]
  dist_list <- qs_read(latest_file, nthreads = N_QS_THREADS)

  df_xir <- data.frame(
    a = dist_list$rolling_average$time,
    f_r_a = dist_list$rolling_average$f_r,
    xi_r = factor(as.character(xir), levels = as.character(xir_vals))
  )
  data_list[[as.character(xir)]] <- df_xir
}
gamma_all <- dplyr::bind_rows(data_list)

# Plot
p <- ggplot(
  gamma_all,
  aes(
    x = a,
    y = f_r_a,
    color = xi_r,
    fill = xi_r
  )
) +
  geom_line(linewidth = 0.8) +
  geom_area(alpha = 0.2, position = "identity") +
  labs(
    x = "Age of infection (days)",
    y = expression(f[r](a)),
    color = expression(xi^r),
    fill = expression(xi^r)
  ) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_color_manual(
    values = c("1" = "firebrick3", "4" = "mistyrose3"),
    labels = c("1" = expression(1), "4" = expression(4)),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c("1" = "firebrick3", "4" = "mistyrose3"),
    labels = c("1" = expression(1), "4" = expression(4)),
    drop = FALSE
  ) +
  theme_minimal(base_size = 14)

print(p)

# Save results
out_pdf <- file.path(figs_dir, "Figure-08a-pdf-fr-fct-xir.pdf")
out_png <- file.path(figs_dir, "Figure-08a-pdf-fr-fct-xir.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm", dpi = 300)

cat("\nFigures saved to", out_pdf, "and", out_png, "\n")

cat("\nSummary Statistics of Recovery Time:\n")
suppressWarnings({
  stats_df <- gamma_all %>%
    group_by(xi_r) %>%
    arrange(a) %>%
    mutate(da = c(a[2] - a[1], diff(a))) %>%
    summarise(
      mean_val = sum(a * f_r_a * da, na.rm = TRUE),
      median_val = a[which(cumsum(f_r_a * da) >= 0.5)[1]],
      .groups = 'drop'
    )
})

for (i in seq_len(nrow(stats_df))) {
  cat(sprintf("  - xi^r = %s: mean = %.2f days, median = %.2f days\n", 
              as.character(stats_df$xi_r[i]), 
              stats_df$mean_val[i], 
              stats_df$median_val[i]))
}

print_end_time(start_time, "plot-Figure-08a-pdf-fr-fct-xic.R")
