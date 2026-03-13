library(dplyr)
library(ggplot2)
library(qs2)
library(tidyr)

# ------------------------------------------------------------
# Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
output_dir <- file.path(getwd(), "OUTPUT")
files <- list.files(output_dir, pattern = "cohort_censored_.*\\.qs$|cohort-censored_.*\\.qs$|cohort_results_truncated\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No truncated/censored cohort file found in ", output_dir)
}

# Sort by modification time to get the newest one if multiple exist
latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

# Load data
cohort_df <- qs_read(latest_file)

# Compute Psi_max per patient
patient_summary <- cohort_df %>%
  group_by(individual_id, status) %>%
  summarise(Psi_max = max(Psi, na.rm=TRUE), .groups = "drop")

# Thresholds ξᵈ to evaluate
xi_values <- c(75, 80, 85, 90, 95)

# Function to classify outcomes for each threshold
classify_outcomes <- function(df, xi) {
  df %>%
    mutate(
      outcome = case_when(
        Psi_max < 75 ~ "Mild",          # toujours < 75%
        Psi_max >= xi ~ "Dead",         # au-dessus de xi_d → mort
        TRUE ~ "ICU"                    # entre 75 et xi_d → ICU
      ),
      xi_d = xi
    )
}

# Apply to all thresholds
classified <- bind_rows(lapply(xi_values, classify_outcomes, df = patient_summary))

# Compute proportions
prop_df <- classified %>%
  count(xi_d, outcome) %>%
  group_by(xi_d) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

# Outcome colors
col_outcomes <- c(
  Mild = "dodgerblue4",
  ICU  = "orange",
  Dead = "red"
)

# Capitalize factor levels to match colors
prop_df$outcome <- factor(prop_df$outcome, levels = c("Mild", "ICU", "Dead"))

# Plot: stacked bar chart
p_outcomes <- ggplot(prop_df, aes(x = factor(xi_d), y = percent, fill = outcome)) +
  geom_bar(stat = "identity", width = 0.75, color = "white") +
  scale_fill_manual(values = col_outcomes) +
  labs(
    x = expression(xi^d~"(%)"),
    y = "Percentage of individuals (%)",
    fill = "Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_outcomes)

# --- Save both PNG and PDF versions ---
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/plot_Figure_05b_pct_outcomes_fct_xid.pdf"
out_png <- "FIGS/plot_Figure_05b_pct_outcomes_fct_xid.png"

ggsave(
  filename = out_pdf,
  plot = p_outcomes,
  width = 25, height = 15, units = "cm", dpi = 300
)

ggsave(
  filename = out_png,
  plot = p_outcomes,
  width = 25, height = 15, units = "cm", dpi = 300
)

cat("\nSaved to", out_pdf, "and", out_png, "\n")
