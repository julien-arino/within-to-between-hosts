## plot-Figure-06a-duration-infectious-period-fct-xic.R
# Generates a boxplot showing the duration of the infectious period 
# (tau_r - tau_c) for varying infectiousness onset thresholds (xi^c),
# while holding hospitalisation and death thresholds fixed.

library(qs2)
library(dplyr)
library(ggplot2)

# USER SETTINGS
xi_h_target <- 75
xi_d_target <- 85
xi_c_values <- 1:10
output_dir <- file.path(getwd(), "OUTPUT")

# Load and compute duration
results <- list()
percent_df <- data.frame()

for (xi_c in xi_c_values) {
  
  cat("Computing infectious duration for xi_c =", xi_c, "\n")
  
  # Format string for file matching (handling integers from Julia output)
  pattern_str <- sprintf("^cohort_times_.*_xih_%d_xid_%d_xic_%s\\.qs$", 
                         xi_h_target, xi_d_target, xi_c)
  files <- list.files(output_dir, pattern = pattern_str, full.names = TRUE)
  
  if (length(files) == 0) {
    # It might have been saved as .0 by Julia
    pattern_str <- sprintf("^cohort_times_.*_xih_%d_xid_%d_xic_%s\\.0\\.qs$", 
                           xi_h_target, xi_d_target, xi_c)
    files <- list.files(output_dir, pattern = pattern_str, full.names = TRUE)
  }
  
  if (length(files) == 0) {
    warning("No cohort_times file found for xi_c = ", xi_c)
    next
  }
  
  latest_file <- files[which.max(file.mtime(files))]
  cohort_df <- qs_read(latest_file)
  
  total_patients <- nrow(cohort_df)
  
  # Establish the simulation end time (usually 30 or 60 days) to bound NA tau_r values
  sim_end_time <- ceiling(max(c(cohort_df$tau_r, cohort_df$tau_d, cohort_df$max_V_t), na.rm = TRUE) / 10) * 10
  
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
    x = expression(xi^c~"(Log10 viral load threshold)"),
    y = "Infectious duration (days)"
  ) +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = length(xi_c_values) / 2 + 0.5,
    y = max(duration_df$infectious_duration, na.rm=TRUE) * 1.18,
    label = "Percentage transmitters",
    size = 4,
    fontface = "italic"
  )

print(p)

# Save
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/Figure-06a-duration-infectious-period-fct-xic.pdf"
out_png <- "FIGS/Figure-06a-duration-infectious-period-fct-xic.png"

ggsave(out_pdf, plot = p, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 25, height = 15, units = "cm", dpi = 300)

cat("\nSaved to", out_pdf, "and", out_png, "\n")
