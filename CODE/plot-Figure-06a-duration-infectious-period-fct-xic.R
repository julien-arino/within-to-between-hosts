library(qs2)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
output_dir <- file.path(getwd(), "OUTPUT")
files <- list.files(output_dir, pattern = "cohort_censored_.*\\.qs$|cohort-censored_.*\\.qs$|cohort_results_truncated\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No truncated/censored cohort file found in ", output_dir)
}

latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

# Load cohort
cohort_df <- qs_read(latest_file)

# ------------------------------------------------------------
# Compute beta_hat
# ------------------------------------------------------------
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

cohort_df$beta_hat <- compute_beta_hat(cohort_df$V)

# ------------------------------------------------------------
# Sort data
# ------------------------------------------------------------
cohort_df <- cohort_df %>%
  arrange(individual_id, time)

# ------------------------------------------------------------
# xi_c values
# ------------------------------------------------------------
xi_seq <- c(1e-06, 1e-05, 1e-04, 0.001, 0.01)

results <- list()
percent_results <- data.frame()

total_patients <- length(unique(cohort_df$individual_id))

# ------------------------------------------------------------
# Loop over xi_c
# ------------------------------------------------------------
for (xi_val in xi_seq) {
  
  cat("Computing infectious duration for xi_c =", xi_val, "\n")
  
  tau_df <- cohort_df %>%
    group_by(individual_id) %>%
    summarise(
      
      tau_c = {
        idx <- which(beta_hat >= xi_val)
        if (length(idx) > 0) min(time[idx]) else Inf
      },
      
      tau_r = {
        idx <- which(beta_hat >= xi_val)
        if (length(idx) > 0) max(time[idx]) else Inf
      },
      
      .groups = "drop"
    )
  
  tau_df <- tau_df %>%
    mutate(
      infectious_duration = tau_r - tau_c,
      xi_c = xi_val
    )
  
  # patients who actually transmit
  tau_df <- tau_df %>%
    filter(is.finite(infectious_duration))
  
  results[[as.character(xi_val)]] <- tau_df
  
  percent <- 100 * nrow(tau_df) / total_patients
  
  percent_results <- rbind(
    percent_results,
    data.frame(
      xi_c = xi_val,
      percent = percent
    )
  )
}

duration_df <- bind_rows(results)

# ------------------------------------------------------------
# Labels for percentages
# ------------------------------------------------------------
percent_results$label <- paste0(round(percent_results$percent,1), "%")
percent_results$ypos <- max(duration_df$infectious_duration, na.rm=TRUE) * 1.08

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
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
    data = percent_results,
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
    x = expression(xi^c~"(Log10 values)"),
    y = "Infectious duration (days)"
  ) +
  
  theme_minimal(base_size = 14)

# indication of percentage
p <- p +
  annotate(
    "text",
    x = 3,
    y = max(duration_df$infectious_duration, na.rm=TRUE) * 1.18,
    label = "Percentage transmitters",
    size = 4,
    fontface = "italic"
  )

print(p)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/plot_Figure_06a_duration_infectious_period_fct_xic.pdf"
out_png <- "FIGS/plot_Figure_06a_duration_infectious_period_fct_xic.png"

ggsave(
  filename = out_pdf,
  plot = p,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = out_png,
  plot = p,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 300
)
cat("\nSaved to", out_pdf, "and", out_png, "\n")
