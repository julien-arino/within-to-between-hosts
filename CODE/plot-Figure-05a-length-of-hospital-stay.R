library(qs2)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Load cohort
# ------------------------------------------------------------
cohort_df <- qs_read("OUTPUT/cohort_results_truncated.qs")

# ------------------------------------------------------------
# Hospitalisation thresholds
# ------------------------------------------------------------
xi_h_values <- c(50, 60, 70, 80, 90)

results <- list()
percent_df <- data.frame()

total_patients <- length(unique(cohort_df$individual_id))

for (xi_h in xi_h_values) {
  
  cat("Computing hospitalisation for xi_h =", xi_h, "\n")
  
  hosp_df <- cohort_df %>%
    group_by(individual_id) %>%
    summarise(
      
      tau_h_start = {
        idx <- which(Psi >= xi_h)
        if (length(idx) > 0) min(time[idx]) else NA
      },
      
      tau_h_end = {
        idx <- which(Psi >= xi_h)
        if (length(idx) > 0) max(time[idx]) else NA
      },
      
      status = first(status),
      .groups = "drop"
    )
  
  # Keep hospitalised only
  hosp_df <- hosp_df %>%
    filter(!is.na(tau_h_start))
  
  hosp_df <- hosp_df %>%
    mutate(
      duration = tau_h_end - tau_h_start,
      xi_h = xi_h
    )
  
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

# ------------------------------------------------------------
# Position for percentage labels
# ------------------------------------------------------------
percent_df$ypos <- max(hospital_df$duration, na.rm=TRUE) * 1.1
percent_df$label <- paste0(round(percent_df$percent,1), "%")

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
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
  
  theme_minimal(base_size = 14)+
  annotate(
    "text",
    x = 3,
    y = max(hospital_df$duration, na.rm=TRUE) * 1.25,
    label = "Hospitalised patients (%)",
    size = 4,
    fontface = "italic"
  )

print(p)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/plot_Figure_05a_hospitalisation_duration_xih.pdf"
out_png <- "FIGS/plot_Figure_05a_hospitalisation_duration_xih.png"

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
cat("\nSaved to", out_pdf, "and", out_png, "\n")
