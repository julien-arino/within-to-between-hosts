library(ggplot2)
library(dplyr)
library(qs2)

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

cohort_df <- qs_read(latest_file)

# ------------------------------------------------------------
# Extract maximum Psi and time of maximum per patient
# ------------------------------------------------------------
death_summary <- cohort_df %>%
  group_by(individual_id, status) %>%
  summarise(
    psi_max = max(Psi, na.rm = TRUE),
    tau_d   = time[which.max(Psi)],
    .groups = "drop"
  )

# ------------------------------------------------------------
# Build dataset for thresholds 85 and 95
# ------------------------------------------------------------
tau_violin <- bind_rows(
  
  death_summary %>%
    filter(status == "Dead" | psi_max >= 85) %>%
    mutate(threshold = "85"),
  
  death_summary %>%
    filter(status == "Dead" | psi_max >= 95) %>%
    mutate(threshold = "95")
  
) %>%
  mutate(threshold = factor(threshold, levels = c("85","95")))

# ------------------------------------------------------------
# Plot horizontal violins
# ------------------------------------------------------------
p_violin <- ggplot(
  tau_violin,
  aes(x = tau_d, y = threshold, fill = threshold)
) +
  
  geom_violin(alpha = 0.6, color = NA) +
  
  geom_boxplot(
    width = 0.15,
    fill = "white",
    alpha = 0.7,
    outlier.shape = NA
  ) +
  
  coord_cartesian(xlim = c(0,20)) +
  
  labs(
    x = expression("Time to death"~tau[i]^d~"(days)"),
    y = expression("Damage threshold"~xi^d),
    fill = NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "85" = "mistyrose3",
      "95" = "firebrick3"
    ),
    labels = c(
      expression(xi^d == 85),
      expression(xi^d == 95)
    )
  ) +
  
  theme_minimal(base_size = 14)

print(p_violin)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/plot_Figure_07a_death_time_fct_xid.pdf"
out_png <- "FIGS/plot_Figure_07a_death_time_fct_xid.png"

ggsave(
  out_pdf,
  p_violin,
  width = 20,
  height = 8,
  units = "cm"
)

ggsave(
  out_png,
  p_violin,
  width = 20,
  height = 8,
  units = "cm"
)
cat("\nSaved to", out_pdf, "and", out_png, "\n")
