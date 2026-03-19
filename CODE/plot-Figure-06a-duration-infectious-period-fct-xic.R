#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-06a-duration-infectious-period-fct-xic.R
# Description:
#   Generates a boxplot showing the duration of the infectious period
#   (tau_r - tau_c) for varying infectiousness onset thresholds (xi^c),
#   while holding hospitalisation and death thresholds fixed.
# ============================================================

cat("\n\n>>> Running plot-Figure-06a-duration-infectious-period-fct-xic.R ...\n\n")
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

# USER SETTINGS
xi_h_target <- 75
xi_d_target <- 85
xi_c_values <- 4:8

# Set project root automatically relative to the .git tracking directory
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

# Load all available files for dynamic tracking
all_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xir_1\\.qs$", full.names = TRUE)
if (length(all_files) == 0) {
  stop("No cohort_status files found with xir_1 in ", output_dir)
}

results <- list()
percent_df <- data.frame()

# First compile all valid combinations found
valid_files <- data.frame(file = all_files, stringsAsFactors = FALSE) %>%
  mutate(
    xih = as.numeric(gsub(".*_xih_([0-9]+).*", "\\1", basename(file))),
    xid = as.numeric(gsub(".*_xid_([0-9]+).*", "\\1", basename(file))),
    xic = as.numeric(gsub(".*_xic_([0-9.]+).*", "\\1", basename(file)))
  )

# If the targeted base isn't in the dataset, fallback to the clearest baseline
if (!(xi_h_target %in% unique(valid_files$xih))) {
  xi_h_target <- valid_files$xih[1]
  cat("Target xi_h=75 not found. Falling back to plotting xi_h =", xi_h_target, "\n")
}
if (!(xi_d_target %in% unique(valid_files$xid))) {
  xi_d_target <- valid_files$xid[1]
  cat("Target xi_d=85 not found. Falling back to plotting xi_d =", xi_d_target, "\n")
}

valid_subset <- valid_files %>% filter(xih == xi_h_target, xid == xi_d_target)
if(nrow(valid_subset) == 0) stop("No files matched the combined fallback targets!")

# Process only those within the subset
for (i in seq_len(nrow(valid_subset))) {
  f <- valid_subset$file[i]
  xi_c <- valid_subset$xic[i]
  
  cat("Loading parsed metrics for xi_c =", xi_c, "\n")
  cohort_df <- qs_read(f, nthreads = N_QS_THREADS)
  
  total_patients <- nrow(cohort_df)
  
  # Establish the simulation end time (usually 30 or 60 days) to bound NA tau_r values
  sim_end_time <- ceiling(max(c(cohort_df$tau_r, cohort_df$tau_d, cohort_df$tau_max_V), na.rm = TRUE) / 10) * 10
  
  # Calculate valid infectious durations exactly as mathematically defined in the manuscript
  tau_df <- cohort_df %>%
    mutate(
      infectious_duration = case_when(
        # If they never became infectious, they have no duration
        is.na(tau_c) ~ NA_real_,
        
        # For those who die, infectious period is from tau_c to tau_d
        !is.na(tau_d) ~ tau_d - tau_c,
        
        # For those who don't die, it is from tau_c to tau_r (or maximum time if tau_r is NA)
        is.na(tau_r) ~ sim_end_time - tau_c,
        
        # Else, standard clearance
        TRUE ~ tau_r - tau_c
      ),
      xi_c = xi_c
    ) %>%
    filter(!is.na(infectious_duration) & infectious_duration > 0)
  
  results[[as.character(xi_c)]] <- tau_df
  
  percent <- 100 * nrow(tau_df) / total_patients
  
  percent_df <- rbind(
    percent_df,
    data.frame(
      xi_c = xi_c,
      percent = percent
    )
  )
}

cat("\n==== Transmitter Percentages ===\n")
cat(" xi_c  |  % Transmitters\n")
cat("-------------------------------\n")
for(i in seq_len(nrow(percent_df))) {
  cat(sprintf(" %4s  |  %6.2f%%\n", percent_df$xi_c[i], percent_df$percent[i]))
}
cat("===============================\n\n")

cat("==== Cohort Viral Load Stats ===\n")
cat("Summary of Maximum Viral Load (max_V) across the cohort:\n")
print(summary(cohort_df$max_V))
cat(sprintf("Fraction of individuals with max_V >= 5.0: %.2f%%\n", 100 * mean(cohort_df$max_V >= 5.0, na.rm=TRUE)))
cat(sprintf("Fraction of individuals with max_V >= 4.5: %.2f%%\n", 100 * mean(cohort_df$max_V >= 4.5, na.rm=TRUE)))
cat("=================================\n\n")

duration_df <- bind_rows(results)

# Labels for percentages
percent_df$label <- paste0(round(percent_df$percent, 1), "%")
percent_df$ypos <- max(duration_df$infectious_duration, na.rm=TRUE) * 1.08

# Plot
p <- ggplot(
  duration_df,
  aes(
    x = factor(xi_c),
    y = infectious_duration
  )
) +
  geom_boxplot(
    fill = "steelblue",
    alpha = 0.7
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "red",
    size = 3
  ) +
  geom_text(
    data = percent_df,
    aes(
      x = factor(xi_c),
      y = ypos,
      label = label
    ),
    size = 4,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  labs(
    x = expression(xi^c~"(log"[10]~"viral load threshold)"),
    y = "Infectious period duration (days)"
  ) +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = length(unique(duration_df$xi_c)) / 2 + 0.5,
    y = max(duration_df$infectious_duration, na.rm=TRUE) * 1.18,
    label = "Percentage of transmitters",
    size = 4,
    fontface = "italic"
  )

print(p)

# Save
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-06a-duration-infectious-period-fct-xic.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-06a-duration-infectious-period-fct-xic.png")

ggsave(out_pdf, plot = p, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 25, height = 15, units = "cm", dpi = 300)

cat("\nSaved to:\n  -", out_pdf, "\n  -", out_png, "\n")

cat("\n\n>>> plot-Figure-06a-duration-infectious-period-fct-xic.R successfully finished running ✅\n\n")
