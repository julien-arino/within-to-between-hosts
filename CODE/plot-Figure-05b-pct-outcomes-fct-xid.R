library(dplyr)
library(ggplot2)
library(qs2)
library(tidyr)

# USER SETTINGS
xi_h_target <- 70
xi_values <- c(75, 80, 85, 90, 95)
output_dir <- file.path(getwd(), "OUTPUT")

# Function to load threshold data and count outcomes
get_outcome_counts <- function(xi_d_val) {
  # Find latest file for this specific xih and xid combination
  pattern_str <- sprintf("^cohort_times_.*_xih_%d_xid_%d\\.qs$", xi_h_target, xi_d_val)
  files <- list.files(output_dir, pattern = pattern_str, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No cohort_times file found matching pattern: ", pattern_str)
  }
  
  latest_file <- files[which.max(file.mtime(files))]
  cat("Loading:", basename(latest_file), "\n")
  
  cohort_df <- qs_read(latest_file)
  
  # Based on strict thresholds applied to the actual maximum tissue damage
  class_df <- cohort_df %>%
    mutate(
      outcome = case_when(
        Psi_max >= xi_d_val ~ "Dead",
        Psi_max >= xi_h_target ~ "ICU",
        TRUE ~ "Mild"
      ),
      xi_d = xi_d_val
    )
  
  counts <- class_df %>% count(xi_d, outcome)
  return(counts)
}

# Apply to all xi_d thresholds
cat("Aggregating outcomes across xi_d values, holding xi_h =", xi_h_target, "...\n")
counts_all <- bind_rows(lapply(xi_values, get_outcome_counts))

# Compute proportions
prop_df <- counts_all %>%
  group_by(xi_d) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

# Outcome colors
col_outcomes <- c(
  Mild = "dodgerblue4",
  ICU  = "orange",
  Dead = "red"
)

# Capitalize factor levels to match colors and standardise ordering
# To put Mild at the bottom, ICU in the middle, Dead at the top, we need factor levels c("Dead", "ICU", "Mild")
prop_df$outcome <- factor(prop_df$outcome, levels = c("Dead", "ICU", "Mild"))

# Print summary table
cat("\n==== Percentages per Outcome (xi_h = ", xi_h_target, ") ===\n")
summary_table <- prop_df %>%
  select(xi_d, outcome, percent) %>%
  pivot_wider(names_from = outcome, values_from = percent) %>%
  arrange(xi_d) %>%
  # We can reorder columns for display: xi_d, Mild, ICU, Dead
  select(xi_d, Mild, ICU, Dead) %>%
  mutate(across(c(Mild, ICU, Dead), ~sprintf("%.2f%%", .x)))

print(as.data.frame(summary_table), row.names = FALSE)
cat("===============================================\n\n")

# Plot: stacked bar chart
p_outcomes <- ggplot(prop_df, aes(x = factor(xi_d), y = percent, fill = outcome)) +
  geom_bar(stat = "identity", width = 0.75, color = "white") +
  scale_fill_manual(
    values = col_outcomes,
    breaks = c("Mild", "ICU", "Dead")
  ) +
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

# Save both PNG and PDF versions
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/Figure-05b-pct-outcomes-fct-xid.pdf"
out_png <- "FIGS/Figure-05b-pct-outcomes-fct-xid.png"

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

cat("\nSaved Figure 05b to", out_pdf, "and", out_png, "\n")
