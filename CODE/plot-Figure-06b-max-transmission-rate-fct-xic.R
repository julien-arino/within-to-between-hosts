suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

# ------------------------------------------------------------
# Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
# Set project root automatically relative to the .git tracking directory
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

# USER SETTINGS
xi_h_target <- 75
xi_d_target <- 85
xi_r_target <- 1

output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-all.R")) else source("functions-all.R")
}

pattern_str <- sprintf("^cohort_status_P.*_xih_%s_xid_%s_xic_.*_xir_%s\\.qs$", 
                       xi_h_target, xi_d_target, xi_r_target)

files <- list.files(output_dir, pattern = pattern_str, full.names = TRUE)

if (length(files) == 0) {
  stop("No matching cohort_status file found in ", output_dir)
}

cat("Found", length(files), "status files. Processing...\n")

compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

results <- list()
percent_results <- data.frame()

for (f in files) {
  # Extract xic from filename
  xic_val <- as.numeric(gsub(".*_xic_([0-9.]+)_.*", "\\1", basename(f)))
  
  cat("Processing transmitters for xi_c =", xic_val, "\n")
  
  cohort_df <- qs_read(f, nthreads = N_QS_THREADS)
  total_patients <- nrow(cohort_df)
  
  # A transmitter is someone who became infectious (tau_c is not NA)
  transmitters <- cohort_df %>%
    filter(!is.na(tau_c)) %>%
    mutate(
      beta_max = compute_beta_hat(max_V),
      xi_c = xic_val
    ) %>%
    select(ID, beta_max, status, xi_c)
    
  results[[as.character(xic_val)]] <- transmitters
  
  percent <- 100 * nrow(transmitters) / total_patients
  
  percent_results <- rbind(
    percent_results,
    data.frame(
      xi_c = xic_val,
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
    x = expression(xi^c~"(log"[10]~"viral load threshold)"),
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
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-06b-max-transmission-rate-fct-xic.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-06b-max-transmission-rate-fct-xic.png")

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
