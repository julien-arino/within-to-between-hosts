#!/usr/bin/env Rscript
# ============================================================
# File: debug-Figure-08b-gammaP-fct-xir.R
# Description: Plot gamma_P for raw, spline, and rolling average
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("debug-Figure-08b-gammaP-fct-xir.R")

load_libraries(c("ggplot2", "dplyr", "qs2", "tidyr"))

# Fixed test configuration
xic <- 4
xir <- 2

cat("Loading isolated filtering test data for xic=4, xir=4...\n")
dist_pattern <- "^cohort_distribution_filters_P.*_xih_75_xid_85_xic_4_xir_4\\.qs$"
dist_files <- list.files(output_dir, pattern = dist_pattern, full.names = TRUE)

if (length(dist_files) == 0) stop("No distribution file found for xic=4, xir=4!")

latest_file <- dist_files[which.max(file.mtime(dist_files))]
filters_list <- qs_read(latest_file, nthreads = N_QS_THREADS)

cat("Combining filter subsets for comparative rendering...\n")

df_raw <- data.frame(
  a = filters_list$raw$time,
  gamma_a = filters_list$raw$gamma_P,
  filter_type = "raw"
)
df_spline <- data.frame(
  a = filters_list$spline$time,
  gamma_a = filters_list$spline$gamma_P,
  filter_type = "spline"
)
df_rolling <- data.frame(
  a = filters_list$rolling_average$time,
  gamma_a = filters_list$rolling_average$gamma_P,
  filter_type = "rolling_average"
)
df_comp_rolling <- data.frame(
  a = filters_list$comp_using_rolling_avg$time,
  gamma_a = filters_list$comp_using_rolling_avg$gamma_P,
  filter_type = "comp_using_rolling_avg"
)

gamma_all <- dplyr::bind_rows(df_raw, df_rolling, df_comp_rolling, df_spline) %>%
  dplyr::filter(is.finite(gamma_a))
gamma_all$filter_type <- factor(gamma_all$filter_type, levels = c("raw", "rolling_average", "comp_using_rolling_avg", "spline"))

cat("- Rendering Comparative Filter Plot...\n")

p <- ggplot(
  gamma_all,
  aes(
    x = a,
    y = gamma_a,
    color = filter_type,
    fill = filter_type
  )
) +
  geom_line(linewidth = 0.8) +
  geom_area(alpha = 0.2, position = "identity") +
  labs(
    x = "Age of infection (days)",
    y = expression("Hazard rate" ~ gamma[P]),
    color = "Filter Type",
    fill = "Filter Type"
  ) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_color_manual(
    values = c("raw" = "#e41a1c", "rolling_average" = "#377eb8", "comp_using_rolling_avg" = "#984ea3", "spline" = "#4daf4a"),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c("raw" = "#e41a1c", "rolling_average" = "#377eb8", "comp_using_rolling_avg" = "#984ea3", "spline" = "#4daf4a"),
    drop = FALSE
  ) +
  theme_minimal(base_size = 14) +
  facet_wrap(~filter_type, ncol = 1)

out_pdf <- file.path(figs_dir, "debug-Figure-08b-gammaP-fct-xir.pdf")
out_png <- file.path(figs_dir, "debug-Figure-08b-gammaP-fct-xir.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm", dpi = 300)

cat("\nFigures saved to", out_pdf, "and", out_png, "\n")
print_end_time(start_time, "debug-Figure-08b-gammaP-fct-xir.R")
