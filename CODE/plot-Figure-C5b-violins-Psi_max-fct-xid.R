# ============================================================
# File: plot-Figure-C5b-violins-Psi_max-fct-xid.R
# Purpose: Violin plots of Psi_max fct xid
# ============================================================

suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))

# ------------------------------------------------------------
# 1. Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-C5b-violins-Psi_max-fct-xid.R ...\n\n")

# Set project root automatically relative to the .git tracking directory
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

output_dir <- file.path(project_dir, "OUTPUT")
files <- list.files(output_dir, pattern = "^cohort_truncated_state_.*\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No truncated/censored cohort file found in ", output_dir)
}

latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

# --- Load data ---
cohort_df <- qs_read(latest_file)

# --- One Psi_max per patient ---
patient_summary <- cohort_df %>%
  group_by(individual_id, status) %>%
  summarise(Psi_max = max(Psi, na.rm=TRUE), .groups = "drop")

# --- Thresholds we want to compare ---
xi_values <- c(85, 90, 95)

# --- Build cumulative groups ---
psi_groups <- bind_rows(
  lapply(xi_values, function(th) {
    patient_summary %>%
      filter(Psi_max >= th) %>%
      mutate(xi_group = factor(th, levels = xi_values))
  })
)

# Check counts
print(psi_groups %>% count(xi_group))

# --- Colors (one time only) ---
col_map <- c("85"="red1", "90"="#ef8a62", "95"="darkred")

# ------------------------------------------------------------
# Violin version
# ------------------------------------------------------------
p_violin <- ggplot(psi_groups,
                   aes(x = xi_group, y = Psi_max, fill = xi_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.alpha = 0.2, fill = "white") +
  scale_fill_manual(values = col_map, labels = xi_values) +
  scale_x_discrete(labels = xi_values) +
  labs(
    title = "",
    x = expression(xi^d~"(%)"),
    y = expression(Psi[max]~"(%)")
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_violin)

# --- Save both PNG and PDF versions ---
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-C5b-violins-Psi_max-fct-xid.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-C5b-violins-Psi_max-fct-xid.png")

ggsave(
  filename = out_pdf,
  plot = p_violin,
  width = 25, height = 15, units = "cm", dpi = 300
)

ggsave(
  filename = out_png,
  plot = p_violin,
  width = 25, height = 15, units = "cm", dpi = 300
)

cat("✅ Plot saved to", out_pdf, "and", out_png, "\n")
