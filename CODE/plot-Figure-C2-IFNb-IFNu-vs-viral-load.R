# ============================================================
# File: plot-Figure-C2-IFNb-IFNu-vs-viral-load.R
# Purpose: Plot (V, F_B) and (V, F_U) side by side (no log scale)
# ============================================================

library(ggplot2)
library(dplyr)
library(qs2)
library(patchwork)

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
cat("✅ Cohort loaded with", nrow(cohort_df), "rows\n")

# ------------------------------------------------------------
# 2️⃣ Optional sampling (for clarity)
# ------------------------------------------------------------
set.seed(123)
cohort_sample <- cohort_df %>%
  group_by(status) %>%
  sample_n(size = min(3000, n()/3)) %>%
  ungroup()

# ------------------------------------------------------------
# 3️⃣ Define colors
# ------------------------------------------------------------
status_colors <- c(Mild = "dodgerblue4", ICU = "orange", Dead = "red")

# Capitalize factor levels to match colors and expected labels
cohort_sample$status <- factor(cohort_sample$status, levels = c("Mild", "ICU", "Dead"))

# ------------------------------------------------------------
# 4️⃣ Plot V vs F_B
# ------------------------------------------------------------
p_FB <- ggplot(cohort_sample, aes(x = V, y = F_B, color = status)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = status_colors) +
  labs(
    x = "Viral load",
    y = "Bound IFN",
    color = "Status"
  ) +
  theme_minimal(base_size = 14) 

# ------------------------------------------------------------
# 5️⃣ Plot V vs F_U
# ------------------------------------------------------------
p_FU <- ggplot(cohort_sample, aes(x = V, y = F_U, color = status)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = status_colors) +
  labs(
    x = "Viral load",
    y = "Unbound IFN",
    color = "Status"
  ) +
  theme_minimal(base_size = 14)

# ------------------------------------------------------------
# 6️⃣ Combine (2 columns, 1 row)
# ------------------------------------------------------------
combined <- (p_FB | p_FU) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top"
  )
print(combined)

# ------------------------------------------------------------
# 7️⃣ Save to PDF & PNG
# ------------------------------------------------------------
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/Figure-C2-IFNb-IFNu-vs-viral-load.pdf"
out_png <- "FIGS/Figure-C2-IFNb-IFNu-vs-viral-load.png"

ggsave(
  filename = out_pdf, 
  plot = combined, width = 9, height = 4, units = "in", dpi = 300
)

ggsave(
  filename = out_png, 
  plot = combined, width = 9, height = 4, units = "in", dpi = 300
)

cat("✅ Plots saved to:", out_pdf, "and", out_png, "\n")
