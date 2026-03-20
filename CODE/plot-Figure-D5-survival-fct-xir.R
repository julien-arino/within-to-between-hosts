#!/usr/bin/env Rscript

# ============================================================
# File: plot-Figure-D5-survival-fct-xic.R
# Description: Plot the survival function for two values of xi_r (1 and 4), fixing xi_c = 4
# ============================================================

# First things first: locate project directory and load helper functions
# and constants
project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

# Say what we are running and start the clock
start_time <- start_time_and_hello("plot-Figure-D5-survival-fct-xic.R")

# Load libraries
load_libraries(c("ggplot2", "dplyr", "qs2", "tidyr"))

# Load empirical_survival from distribution files
xir_vals <- c(1, 4)
data_list <- list()

for (xir in xir_vals) {
  dist_pattern <- paste0("^cohort_distribution_filters_P.*_xih_75_xid_85_xic_4_xir_", xir, "\\.qs$")
  dist_file <- list.files(output_dir, pattern = dist_pattern, full.names = TRUE)
  if (length(dist_file) == 0) stop(paste("No distribution file found for xir =", xir))

  latest_file <- dist_file[which.max(file.mtime(dist_file))]
  dist_list <- qs_read(latest_file, nthreads = N_QS_THREADS)

  df_xir <- data.frame(
    a = dist_list$spline$time,
    survival = dist_list$spline$survival,
    xi_r = factor(as.character(xir), levels = as.character(xir_vals))
  )
  data_list[[as.character(xir)]] <- df_xir
}
survival_all <- dplyr::bind_rows(data_list) %>%
  dplyr::filter(is.finite(survival))

# Plot
p <- ggplot(
  survival_all,
  aes(
    x = a,
    y = survival,
    color = xi_r,
    fill = xi_r
  )
) +
  geom_line(linewidth = 0.8) +
  geom_area(alpha = 0.2, position = "identity") +
  labs(
    x = "Age of infection (days)",
    y = "Probability of survival",
    color = expression(xi^r),
    fill = expression(xi^r)
  ) +
  coord_cartesian(xlim = c(0, 80)) +
  scale_color_manual(
    values = c("1" = "#e7298a", "2" = "#1b9e77", "4" = "#d95f02"),
    labels = c("1" = expression(1), "2" = expression(2), "4" = expression(4)),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c("1" = "#e7298a", "2" = "#1b9e77", "4" = "#d95f02"),
    labels = c("1" = expression(1), "2" = expression(2), "4" = expression(4)),
    drop = FALSE
  ) +
  theme_minimal(base_size = 14)

print(p)

# Save results
out_pdf <- file.path(figs_dir, "Figure-D5-survival-fct-xir.pdf")
out_png <- file.path(figs_dir, "Figure-D5-survival-fct-xir.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm", dpi = 300)

cat("Figures saved to", out_pdf, "and", out_png, "\n")

print_end_time(start_time, "plot-Figure-D5-survival-fct-xic.R")
