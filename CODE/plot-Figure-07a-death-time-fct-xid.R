## plot-Figure-07a-death-time-fct-xid.R
# Generates a horizontal violin plot showing the time to death (tau_d)
# for two different death thresholds (xi^d = 85 and 95).

suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(here)))

cat("\n\n>>> Running plot-Figure-07a-death-time-fct-xid.R ...\n\n")

# USER SETTINGS
xi_h_target <- 70 # Fixed xi_h just to find a valid file (tau_d is independent of xi_h)
# USER SETTINGS
xi_d_values <- c(85, 95)

# Set project root automatically relative to the .git tracking directory
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")

# Load all available files to map available parameters dynamically
all_files <- list.files(output_dir, pattern = "^cohort_status_P.*\\.qs$", full.names = TRUE)
if (length(all_files) == 0) {
  stop("No cohort_status files found with any parameters in ", output_dir)
}

# Compile all combinations
valid_files <- data.frame(file = all_files, stringsAsFactors = FALSE) %>%
  mutate(
    xih = as.numeric(gsub(".*_xih_([0-9]+).*", "\\1", basename(file))),
    xid = as.numeric(gsub(".*_xid_([0-9]+).*", "\\1", basename(file)))
  )

# If the targeted xi_h isn't in the dataset, fallback to the clearest baseline
if (!(xi_h_target %in% unique(valid_files$xih))) {
  xi_h_target <- valid_files$xih[1]
  cat("Target xi_h=70 not found. Falling back to extracting death times using baseline xi_h =", xi_h_target, "\n")
}

results <- list()

for (xi_d in xi_d_values) {
  cat("Loading data for xi_d =", xi_d, "...\n")
  
  target_file <- valid_files %>% filter(xih == xi_h_target, xid == xi_d)
  
  if (nrow(target_file) == 0) {
    warning("No status data available for targeted xi_d = ", xi_d)
    next
  }
  
  # If there happen to be multiple matching this pattern exactly, grab the newest
  latest_file <- target_file$file[which.max(file.mtime(target_file$file))]
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
      "85" = "firebrick3",
      "95" = "mistyrose3"
    ),
    labels = c(
      expression(xi^d == 85),
      expression(xi^d == 95)
    )
  ) +
  
  theme_minimal(base_size = 14)

print(p_violin)

# Save
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-07a-death-time-fct-xid.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-07a-death-time-fct-xid.png")

ggsave(out_pdf, plot = p_violin, width = 20, height = 8, units = "cm", dpi = 300)
ggsave(out_png, plot = p_violin, width = 20, height = 8, units = "cm", dpi = 300)

cat("\nSaved to:\n  -", out_pdf, "\n  -", out_png, "\n")
