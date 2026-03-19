#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-C1-max-viral-load-fct-status.R
# Purpose: Plot maximum viral load per individual (colored by status)
# ============================================================

cat("\n\n>>> Running plot-Figure-C1-max-viral-load-fct-status.R ...\n\n")
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(here)))

# ------------------------------------------------------------
# 1. Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
# Set project root automatically relative to the .git tracking directory
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}
files <- list.files(output_dir, pattern = "^cohort_status_P.*_xid_[0-9]+\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No base cohort_status file found in ", output_dir)
}

latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort status results:", basename(latest_file), "\n")

cohort_df <- qs_read(latest_file, nthreads = N_QS_THREADS)

cat("Cohort data loaded with", nrow(cohort_df), "individuals\n")

# ------------------------------------------------------------
# 2️⃣ Compute maximum Viral Load per individual
# ------------------------------------------------------------
cat("Extracting V_max per individual from status dataframe...\n")
viral_max_df <- cohort_df %>%
  select(ID, status, V_max = max_V)

cat("Computed maximum viral load for", nrow(viral_max_df), "individuals\n")

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

dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-C1-max-viral-load-fct-status.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-C1-max-viral-load-fct-status.png")

ggsave(
  filename = out_png, 
  plot = p_Vmax, width = 9, height = 5, units = "in", dpi = 300
)
ggsave(
  filename = out_pdf, 
  plot = p_Vmax, width = 9, height = 5, units = "in", dpi = 300
)

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")

cat("\n\n>>> plot-Figure-C1-max-viral-load-fct-status.R successfully finished running ✅\n\n")
