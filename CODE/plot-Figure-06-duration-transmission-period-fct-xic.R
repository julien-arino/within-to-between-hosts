#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-06-duration-transmission-period-fct-xic.R
# Description:
#   Generates a grouped boxplot showing the duration of the transmission period
#   (tau_r - tau_c) for varying transmissionness onset thresholds (xi^c),
#   color-coded by recovery threshold (xi^r), while holding hospitalisation
#   and death thresholds fixed.
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-06-duration-transmission-period-fct-xic.R")

load_libraries(c("ggplot2", "dplyr", "qs2"))

# USER SETTINGS
xi_h_target <- 75
xi_d_target <- 85

# Load all available files for dynamic tracking
all_files <- list.files(output_dir, pattern = "^cohort_status_P.*\\.qs$", full.names = TRUE)
if (length(all_files) == 0) {
  stop("No cohort_status files found in ", output_dir)
}

results <- list()
percent_df <- data.frame()

# Compile all valid combinations mathematically
valid_files <- data.frame(file = all_files, stringsAsFactors = FALSE) %>%
  filter(grepl("_xih_", file), grepl("_xid_", file), grepl("_xic_", file), grepl("_xir_", file)) %>%
  mutate(
    xih = as.numeric(gsub(".*_xih_([0-9]+).*", "\\1", basename(file))),
    xid = as.numeric(gsub(".*_xid_([0-9]+).*", "\\1", basename(file))),
    xic = as.numeric(gsub(".*_xic_([0-9.]+).*", "\\1", basename(file))),
    xir = as.numeric(gsub(".*_xir_([0-9.]+).*", "\\1", basename(file)))
  )

valid_subset <- valid_files %>% filter(xih == xi_h_target, xid == xi_d_target)
if (nrow(valid_subset) == 0) stop(paste("No files matched the target xih =", xi_h_target, "and xid =", xi_d_target))

# Sort logically to track progress cleanly
valid_subset <- valid_subset %>% arrange(xic, xir)

for (i in seq_len(nrow(valid_subset))) {
  f <- valid_subset$file[i]
  xi_c <- valid_subset$xic[i]
  xi_r <- valid_subset$xir[i]

  cat("Loading parsed metrics for xi_c =", xi_c, "| xi_r =", xi_r, "\n")
  cohort_df <- qs_read(f, nthreads = N_QS_THREADS)

  total_patients <- nrow(cohort_df)

  sim_end_time <- ceiling(max(c(cohort_df$tau_r, cohort_df$tau_d, cohort_df$tau_max_V), na.rm = TRUE) / 10) * 10

  tau_df <- cohort_df %>%
    mutate(
      transmission_duration = case_when(
        is.na(tau_c) ~ NA_real_,
        !is.na(tau_d) ~ tau_d - tau_c,
        is.na(tau_r) ~ sim_end_time - tau_c,
        TRUE ~ tau_r - tau_c
      ),
      xi_c = xi_c,
      xi_r = xi_r
    ) %>%
    filter(!is.na(transmission_duration) & transmission_duration > 0)

  results[[f]] <- tau_df

  percent <- 100 * nrow(tau_df) / total_patients

  percent_df <- rbind(
    percent_df,
    data.frame(
      xi_c = xi_c,
      xi_r = xi_r,
      percent = percent
    )
  )
}

cat("\n==== Transmitter Percentages ===\n")
cat(" xi_c | xi_r |  % Transmitters\n")
cat("-------------------------------\n")
for (i in seq_len(nrow(percent_df))) {
  cat(sprintf(" %4s | %4s |  %6.2f%%\n", percent_df$xi_c[i], percent_df$xi_r[i], percent_df$percent[i]))
}
cat("===============================\n\n")

duration_df <- bind_rows(results)

# Deduplicate percentages so we only get one text label per xi_c group
percent_df <- percent_df %>% distinct(xi_c, percent)

# Compute the highest visible point across all dodging geometries 
# (either the top whisker or the red mean point) to tightly crop the empty y-axis space
visible_extents <- duration_df %>%
  group_by(xi_c, xi_r) %>%
  summarise(
    mean_val = mean(transmission_duration, na.rm = TRUE),
    q3 = quantile(transmission_duration, 0.75, na.rm = TRUE),
    iqr = IQR(transmission_duration, na.rm = TRUE),
    top_whisker = suppressWarnings(max(transmission_duration[transmission_duration <= q3 + 1.5 * iqr], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(max_vis = pmax(mean_val, top_whisker, na.rm = TRUE))

max_y_plot <- max(visible_extents$max_vis, na.rm = TRUE)

# Labels for percentages 
percent_df$label <- paste0(round(percent_df$percent, 1), "%")
# Set the bottom of the info ~15% above the highest extent of the boxes/means
percent_df$ypos <- max_y_plot * 1.15

# Define custom color palette depending on what xi_r values exist
unique_r <- as.character(sort(unique(duration_df$xi_r)))
color_palette <- setNames(
  c("lightblue", "steelblue", "royalblue4", "navy")[seq_along(unique_r)],
  unique_r
)

# Plot
p <- ggplot(
  duration_df,
  aes(
    x = factor(xi_c),
    y = transmission_duration,
    fill = factor(xi_r)
  )
) +
  geom_boxplot(
    alpha = 0.8,
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "red",
    size = 2,
    position = position_dodge(width = 0.8),
    show.legend = FALSE
  ) +
  geom_text(
    data = percent_df,
    aes(
      x = factor(xi_c),
      y = ypos,
      label = label
    ),
    size = 3.5,
    fontface = "bold",
    vjust = -0.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = expression(xi^r),
    values = color_palette
  ) +
  labs(
    x = expression(xi^c ~ "(log"[10] ~ "viral load threshold)"),
    y = "Transmission period duration (days)"
  ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right") +
    annotate(
      "text",
      x = length(unique(duration_df$xi_c)) / 2 + 0.5,
      y = max_y_plot * 1.25,
      label = "Percentage of transmitters",
      size = 4,
      fontface = "italic"
    ) +
    # Gracefully crop the camera viewport (to cut out invisible y=100 outliers) 
    # without destroying the underlying math calculations that `ylim()` would trigger
    coord_cartesian(ylim = c(0, max_y_plot * 1.35))

print(p)

# Save
out_pdf <- file.path(figs_dir, "Figure-06-duration-transmission-period-fct-xic.pdf")
out_png <- file.path(figs_dir, "Figure-06-duration-transmission-period-fct-xic.png")

ggsave(out_pdf, plot = p, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 25, height = 15, units = "cm", dpi = 300)

cat("\nSaved to:\n  -", out_pdf, "\n  -", out_png, "\n")

print_end_time(start_time, "plot-Figure-06-duration-transmission-period-fct-xic.R")
