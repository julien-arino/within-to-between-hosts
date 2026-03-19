#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

cat("\n\n>>> Running plot-Figure-05b-pct-outcomes-fct-xid.R ...\n\n")

# USER SETTINGS
xi_h_target <- 70
xi_values <- c(75, 80, 85, 90, 95)

# Set project root automatically relative to the .git tracking directory
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

# Load all available base files (excluding the _xic_..._xir_... extensions)
all_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xid_[0-9]+\\.qs$", full.names = TRUE)
if (length(all_files) == 0) {
  stop("No base cohort_status files found directly ending in _xid_XX.qs in ", output_dir)
}

results <- list()

for (f in all_files) {
  # Extract threshold digits from filename
  match <- regexpr("_xih_([0-9]+)_xid_([0-9]+)", basename(f))
  if (match == -1) next
  
  # Only add if it's the latest generation (avoid double counting old runs)
  # The filenames differentiate themselves by timestamp, so we group by the prefixes
  # Note: A simpler approach is to read it, since we're just plotting aggregates of whatever latest data exists!
  
  cohort_df <- qs_read(f, nthreads = N_QS_THREADS)
  
  # If the file didn't actually contain xi_h and xi_d inside, we extract it from the path string
  xih_val <- as.numeric(gsub(".*_xih_([0-9]+).*", "\\1", basename(f)))
  xid_val <- as.numeric(gsub(".*_xid_([0-9]+).*", "\\1", basename(f)))
  
  # Skip processing if we only want one specific xi_h line (e.g. tracking purely xi_d curves)
  # BUT gracefully fallback to whatever exists if 70 isn't available
  
  class_df <- cohort_df %>%
    mutate(
      outcome = case_when(
        max_Psi >= xid_val ~ "Dead",
        max_Psi >= xih_val ~ "ICU",
        TRUE ~ "Mild"
      ),
      xi_d = xid_val,
      xi_h = xih_val
    )
  
  counts <- class_df %>% count(xi_h, xi_d, outcome)
  results[[basename(f)]] <- counts
}

counts_all <- bind_rows(results)

# If the targeted xi_h isn't in the dataset, fallback to the clearest baseline
if (!(xi_h_target %in% unique(counts_all$xi_h))) {
  xi_h_target <- counts_all$xi_h[1]
  cat("\nTarget xi_h=70 not found. Falling back to plotting xi_h =", xi_h_target, "\n")
}

counts_all <- counts_all %>% filter(xi_h == xi_h_target)

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
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-05b-pct-outcomes-fct-xid.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-05b-pct-outcomes-fct-xid.png")

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

cat("\nSaved Figure 05b to:\n  -", out_pdf, "\n  -", out_png, "\n")
