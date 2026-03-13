# ============================================================
# File: plot-Figure-B1-max-viral-load-fct-status.R
# Purpose: Plot maximum viral load per patient (colored by status)
# ============================================================

library(ggplot2)
library(dplyr)
library(qs2)

# ------------------------------------------------------------
# 1. Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
output_dir <- file.path(getwd(), "OUTPUT")
files <- list.files(output_dir, pattern = "cohort_censored_.*\\.qs$|cohort-censored_.*\\.qs$|cohort_results_truncated\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No truncated/censored cohort file found in ", output_dir)
}

latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

cohort_df <- qs_read(latest_file)

cat("✅ Cohort data loaded with", nrow(cohort_df), "rows\n")

# ------------------------------------------------------------
# 2️⃣ Compute maximum Viral Load per patient
# ------------------------------------------------------------
viral_max_df <- cohort_df %>%
  group_by(individual_id, status) %>%
  summarise(
    V_max = max(V, na.rm = TRUE),
    .groups = "drop"
  )

cat("✅ Computed maximum viral load for", nrow(viral_max_df), "patients\n")

# ------------------------------------------------------------
# 3️⃣ Define custom colors
# ------------------------------------------------------------
status_colors <- c(Mild = "dodgerblue4", ICU = "orange", Dead = "red")

# Update factor levels for status
viral_max_df$status <- factor(viral_max_df$status, levels = c("Mild", "ICU", "Dead"))

# ------------------------------------------------------------
# 4️⃣ Create the plot
# ------------------------------------------------------------
p_Vmax <- ggplot(viral_max_df, aes(x = status, y = V_max, color = status)) +
  geom_jitter(width = 0.25, alpha = 0.3, size = 1) +                # scatter points
  geom_boxplot(outlier.shape = NA, fill=NA, alpha = 0.4, width = 0.5, color="black") +      # overlay boxplot
  scale_color_manual(values = status_colors) +
  labs(
    x = "Severity group",
    y = expression("log"[10]*"(Max Viral load [copies/ml])"),
    color = "Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    # plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ------------------------------------------------------------
# 5️⃣ Print the plot to screen
# ------------------------------------------------------------
print(p_Vmax)

# ------------------------------------------------------------
# 6️⃣ Save plot to PDF and PNG
# ------------------------------------------------------------

dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/plot_Figure_B1_max_viral_load_fct_status.pdf"
out_png <- "FIGS/plot_Figure_B1_max_viral_load_fct_status.png"

ggsave(
  filename = out_png, 
  plot = p_Vmax, width = 9, height = 5, units = "in", dpi = 300
)
ggsave(
  filename = out_pdf, 
  plot = p_Vmax, width = 9, height = 5, units = "in", dpi = 300
)

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")
