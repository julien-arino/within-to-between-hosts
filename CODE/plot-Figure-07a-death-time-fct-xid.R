## plot-Figure-07a-death-time-fct-xid.R
# Generates a horizontal violin plot showing the time to death (tau_d)
# for two different death thresholds (xi^d = 85 and 95).

library(ggplot2)
library(dplyr)
library(qs2)

# USER SETTINGS
xi_h_target <- 70 # Fixed xi_h just to find a valid file (tau_d is independent of xi_h)
xi_d_values <- c(85, 95)
output_dir <- file.path(getwd(), "OUTPUT")

# Load and process threshold data
results <- list()

for (xi_d in xi_d_values) {
  cat("Loading data for xi_d =", xi_d, "...\n")
  
  pattern_str <- sprintf("^cohort_times_.*_xih_%d_xid_%d\\.qs$", xi_h_target, xi_d)
  files <- list.files(output_dir, pattern = pattern_str, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No cohort_times file found matching pattern: ", pattern_str)
  }
  
  latest_file <- files[which.max(file.mtime(files))]
  cohort_df <- qs_read(latest_file)
  
  # Filter for individuals who died (tau_d is not NA)
  death_df <- cohort_df %>%
    filter(!is.na(tau_d)) %>%
    select(tau_d) %>%
    mutate(threshold = as.character(xi_d))
    
  results[[as.character(xi_d)]] <- death_df
}

tau_violin <- bind_rows(results) %>%
  mutate(threshold = factor(threshold, levels = as.character(xi_d_values)))

# Plot horizontal violins
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
  
  coord_cartesian(xlim = c(0, 20)) +
  
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

# Save
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/Figure-07a-death-time-fct-xid.pdf"
out_png <- "FIGS/Figure-07a-death-time-fct-xid.png"

ggsave(out_pdf, plot = p_violin, width = 20, height = 8, units = "cm", dpi = 300)
ggsave(out_png, plot = p_violin, width = 20, height = 8, units = "cm", dpi = 300)

cat("\nSaved to", out_pdf, "and", out_png, "\n")
