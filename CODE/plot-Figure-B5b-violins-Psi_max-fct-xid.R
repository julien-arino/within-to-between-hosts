# ============================================================
# File: plot-Figure-B5b-violins-Psi_max-fct-xid.R
# Purpose: Violin plots of Psi_max fct xid
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
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/plot_Figure_B5b_violins_Psi_max_fct_xid.pdf"
out_png <- "FIGS/plot_Figure_B5b_violins_Psi_max_fct_xid.png"

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
