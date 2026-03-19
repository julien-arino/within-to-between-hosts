## plot-Figure-05a-length-of-hospital-stay.R
# Generates a boxplot showing the length of stay in hospital for varying 
# hospitalisation thresholds (xi^h) while holding the death threshold (xi^d) fixed.
#
# IMPORTANT: This script explicitly conditions on survival through the 
# hospitalisation period. Individuals who ultimately die (tau_d is not NA) 
# are excluded from the duration calculations.

suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

cat("\n\n>>> Running plot-Figure-05a-length-of-hospital-stay.R ...\n\n")

# USER SETTINGS
xi_d_target <- 85
xi_h_values <- c(50, 60, 70, 75, 80)

# Set project root automatically relative to the .git tracking directory
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-all.R")) else source("functions-all.R")
}

# Function to load threshold data and compute hospitalization
results <- list()
percent_df <- data.frame()

for (xi_h in xi_h_values) {
  
  cat("Computing hospitalisation for xi_h =", xi_h, "and fixed xi_d =", xi_d_target, "\n")
  
  # Find latest file for this specific xih and xid combination
  pattern_str <- sprintf("^cohort_status_P.*_xih_%d_xid_%d\\.qs$", xi_h, xi_d_target)
  files <- list.files(output_dir, pattern = pattern_str, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No cohort_status file found matching pattern: ", pattern_str)
  }
  
  latest_file <- files[which.max(file.mtime(files))]
  cohort_df <- qs_read(latest_file, nthreads = N_QS_THREADS)
  
  total_patients <- nrow(cohort_df)
  
  # Keep hospitalised only AND condition on survival (exclude those who die)
  hosp_df <- cohort_df %>%
    filter(!is.na(tau_h_start) & is.na(tau_d)) %>%
    mutate(
      duration = tau_h_end - tau_h_start,
      xi_h = xi_h
    ) %>%
    filter(is.finite(duration))
  
  results[[as.character(xi_h)]] <- hosp_df
  
  percent <- 100 * nrow(hosp_df) / total_patients
  
  percent_df <- rbind(
    percent_df,
    data.frame(
      xi_h = xi_h,
      percent = percent
    )
  )
}

hospital_df <- bind_rows(results)

# Print Summary Table for Legend Generation
cat("\n==== Hospitalization Percentages (xi_d = ", xi_d_target, ") ===\n", sep="")
cat(" xi_h    % Hospitalized\n")
cat("-------------------------\n")
for (i in seq_len(nrow(percent_df))) {
  cat(sprintf("  %2d       %6.2f%%\n", percent_df$xi_h[i], percent_df$percent[i]))
}
cat("=========================\n\n")

# Position for percentage labels
percent_df$ypos <- max(hospital_df$duration, na.rm=TRUE) * 1.05
percent_df$label <- paste0(round(percent_df$percent, 1), "%")

# Plot
p <- ggplot(
  hospital_df,
  aes(x = factor(xi_h),
      y = duration)
) +
  
  geom_boxplot(
    fill = "steelblue",
    alpha = 0.7
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    color = "deeppink",
    size = 3
  ) +
  
  geom_text(
    data = percent_df,
    aes(
      x = factor(xi_h),
      y = ypos,
      label = label
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4
  ) +
  
  labs(
    x = expression(xi^h~"(Hospitalisation threshold %)"),
    y = "Length of stay in hospital (days)"
  ) +
  
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = 3, # center of 5 categories
    y = max(hospital_df$duration, na.rm=TRUE) * 1.15,
    label = "Surviving hospitalised patients (%)",
    size = 4,
    fontface = "italic"
  )

print(p)

# Save
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-05a-length-of-hospital-stay.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-05a-length-of-hospital-stay.png")

ggsave(
  out_pdf,
  plot = p,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 300
)

ggsave(
  out_png,
  plot = p,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 300
)
cat("\nSaved to:\n  -", out_pdf, "\n  -", out_png, "\n")
