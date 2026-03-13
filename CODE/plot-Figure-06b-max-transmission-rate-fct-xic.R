library(qs2)
library(dplyr)
library(ggplot2)

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

# Load cohort
cohort_df <- qs_read(latest_file)

# ------------------------------------------------------------
# 2 Compute beta_hat
# ------------------------------------------------------------
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

cohort_df$beta_hat <- compute_beta_hat(cohort_df$V)

# ------------------------------------------------------------
# 3 Compute beta_max per patient
# ------------------------------------------------------------
beta_max_df <- cohort_df %>%
  group_by(individual_id) %>%
  summarise(
    beta_max = max(beta_hat, na.rm = TRUE),
    status = first(status),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4 Values of xi_c
# ------------------------------------------------------------
xi_seq <- c(1e-06, 1e-05, 1e-04, 0.001, 0.01)

# ------------------------------------------------------------
# 5 Collect transmitters and percentages
# ------------------------------------------------------------
results <- list()
percent_results <- data.frame()

total_patients <- nrow(beta_max_df)

for (xi_val in xi_seq) {
  
  cat("Processing transmitters for xi_c =", xi_val, "\n")
  
  transmitters <- beta_max_df %>%
    filter(beta_max >= xi_val)
  
  transmitters$xi_c <- xi_val
  
  results[[as.character(xi_val)]] <- transmitters
  
  percent <- 100 * nrow(transmitters) / total_patients
  
  percent_results <- rbind(
    percent_results,
    data.frame(
      xi_c = xi_val,
      percent = percent
    )
  )
}

transmitters_df <- bind_rows(results)

# ------------------------------------------------------------
# 6 Labels for percentages
# ------------------------------------------------------------
percent_results$label <- paste0(round(percent_results$percent,1), "%")
percent_results$ypos <- max(transmitters_df$beta_max, na.rm=TRUE) * 1.05

# ------------------------------------------------------------
# 7 Plot
# ------------------------------------------------------------
p <- ggplot(
  transmitters_df,
  aes(
    x = factor(xi_c),
    y = beta_max
  )
) +
  
  geom_boxplot(
    fill = "gray62",
    alpha = 0.7
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    color = "deeppink1",
    size = 3
  ) +
  
  geom_text(
    data = percent_results,
    aes(
      x = factor(xi_c),
      y = ypos,
      label = label
    ),
    size = 4,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  labs(
    x = expression("Threshold"~xi^c),
    y = expression("Maximum transmission rate"~beta[max])
  ) +
  
  theme_minimal(base_size = 14)+
  annotate(
    "text",
    x = 3,
    y = max(transmitters_df$beta_max, na.rm=TRUE) * 1.12,
    label = "Percentage of transmitters",
    size = 4,
    fontface = "italic"
  )

print(p)

# ------------------------------------------------------------
# 8 Save figure
# ------------------------------------------------------------
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/plot_Figure_06b_max_transmission_rate_fct_xic.pdf"
out_png <- "FIGS/plot_Figure_06b_max_transmission_rate_fct_xic.png"

ggsave(
  out_pdf,
  plot = p,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 300
)

ggsave(
  out_png,
  plot = p,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 300
)

cat("\nSaved to", out_pdf, "and", out_png, "\n")
