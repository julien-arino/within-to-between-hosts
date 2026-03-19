#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-08a-pdf-fr-fct-xic.R
# Description:
# ============================================================

cat("\n\n>>> Running plot-Figure-08a-pdf-fr-fct-xic.R ...\n\n")
start_time <- Sys.time()
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

# ------------------------------------------------------------
# 1. Setup paths and load the full continuous state array
# ------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

# ------------------------------------------------------------
# 2. Compute f_r locally for transmitting individuals explicitly
# ------------------------------------------------------------
xic_vals <- c(4, 6)
data_list <- list()

for (xic in xic_vals) {
  status_pattern <- paste0("^cohort_status_P.*_xih_75_xid_85_xic_", xic, "_xir_1\\.qs$")
  status_file <- list.files(output_dir, pattern = status_pattern, full.names = TRUE)
  if (length(status_file) == 0) stop(paste("No status file found for xic =", xic))
  
  # grab the newest
  latest_file <- status_file[which.max(file.mtime(status_file))]
  cohort_status <- qs_read(latest_file, nthreads = N_QS_THREADS)
  
  # Cohort constraint: Transmitters ONLY (individuals who entered the transmitting stage)
  cohort_status <- cohort_status %>% filter(!is.na(tau_c))
  
  N <- nrow(cohort_status)
  time_pts <- seq(0, 100.0, by = 0.1)
  
  recovered_cohort <- cohort_status %>% filter(!is.na(tau_r))
  dens_r <- density(recovered_cohort$tau_r, from = 0, to = max(time_pts), n = length(time_pts))
  
  df_xic <- data.frame(
    a = time_pts,
    f_r_a = dens_r$y,
    xi_c = factor(as.character(xic), levels = c("4", "6"))
  )
  data_list[[as.character(xic)]] <- df_xic
}
gamma_all <- dplyr::bind_rows(data_list)

# ------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------
p <- ggplot(
  gamma_all,
  aes(
    x = a,
    y = f_r_a,
    color = xi_c,
    fill = xi_c,
    linetype = xi_c
  )
) +
  geom_line(linewidth = 0.8) +
  geom_area(alpha = 0.2, position = "identity") +
  labs(
    x = "Age of infection (days)",
    y = "Probability density function",
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
  scale_linetype_manual(
    values = c("4" = "solid", "6" = "dashed"),
    labels = c("4" = expression(4), "6" = expression(6)),
    drop = FALSE
  ) +
  theme_minimal(base_size = 14)

print(p)

# ------------------------------------------------------------
# 3. Save results
# ------------------------------------------------------------
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-08a-pdf-fr-fct-xic.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-08a-pdf-fr-fct-xic.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm", dpi = 300)

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")

print_end_time(start_time, "plot-Figure-08a-pdf-fr-fct-xic.R")
