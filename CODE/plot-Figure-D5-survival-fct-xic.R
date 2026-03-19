#!/usr/bin/env Rscript

# ============================================================
# File: plot-Figure-D5-survival-fct-xic.R
# Description: Plot the survival function for two values of xi_c (4 and 6)
# ============================================================

# First things first: locate project directory and load helper functions
# and constants
project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

# Say what we are running and start the clock
start_time <- start_time_and_hello("plot-Figure-D5-survival-fct-xic.R")

# Load libraries
load_libraries(c("ggplot2", "dplyr", "qs2", "tidyr"))

# Compute empirical_survival locally for transmitting individuals explicitly
xic_vals <- c(4, 6)
data_list <- list()

for (xic in xic_vals) {
  status_pattern <- paste0("^cohort_status_P.*_xih_75_xid_85_xic_", xic, "_xir_1\\.qs$")
  status_file <- list.files(output_dir, pattern = status_pattern, full.names = TRUE)
  if (length(status_file) == 0) stop(paste("No status file found for xic =", xic))

  latest_file <- status_file[which.max(file.mtime(status_file))]
  cohort_status <- qs_read(latest_file, nthreads = N_QS_THREADS)

  # Cohort constraint: Transmitters ONLY
  cohort_status <- cohort_status %>% filter(!is.na(tau_c))

  N <- nrow(cohort_status)
  time_pts <- seq(0, 100.0, by = 0.1)

  exit_times <- pmin(
    tidyr::replace_na(cohort_status$tau_d, Inf),
    tidyr::replace_na(cohort_status$tau_r, Inf)
  )
  active_counts <- sapply(time_pts, function(t_val) {
    sum(t_val < exit_times)
  })
  empirical_survival <- active_counts / N

  df_xic <- data.frame(
    a = time_pts,
    survival = empirical_survival,
    xi_c = factor(as.character(xic), levels = c("4", "6"))
  )
  data_list[[as.character(xic)]] <- df_xic
}
survival_all <- dplyr::bind_rows(data_list) %>%
  dplyr::filter(is.finite(survival))

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(
  survival_all,
  aes(
    x = a,
    y = survival,
    color = xi_c,
    fill = xi_c
  )
) +
  geom_line(linewidth = 0.8) +
  geom_area(alpha = 0.2, position = "identity") +
  labs(
    x = "Age of infection (days)",
    y = "Probability of survival",
    color = expression(xi^c),
    fill = expression(xi^c)
  ) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_color_manual(
    values = c("4" = "#e7298a", "6" = "#1b9e77"),
    labels = c("4" = expression(4), "6" = expression(6)),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c("4" = "#e7298a", "6" = "#1b9e77"),
    labels = c("4" = expression(4), "6" = expression(6)),
    drop = FALSE
  ) +
  theme_minimal(base_size = 14)

print(p)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------
out_pdf <- file.path(figs_dir, "Figure-D5-survival-fct-xic.pdf")
out_png <- file.path(figs_dir, "Figure-D5-survival-fct-xic.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm", dpi = 300)

cat("Figures saved to", out_pdf, "and", out_png, "\n")

print_end_time(start_time, "plot-Figure-D5-survival-fct-xic.R")
