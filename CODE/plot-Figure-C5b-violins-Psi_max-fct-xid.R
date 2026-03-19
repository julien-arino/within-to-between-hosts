#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-C5b-violins-Psi_max-fct-xid.R
# Purpose: Violin plots of Psi_max fct xid
# ============================================================

suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))

# ------------------------------------------------------------
# 1. Setup paths and load dynamically matched cohort_status files
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-C5b-violins-Psi_max-fct-xid.R ...\n\n")

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

# --- Setup target xi_h and load all status files ---
xi_h_target <- 75
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*\\.qs$", full.names = TRUE)
# Exclude files with xic and xir changes
status_files <- status_files[!grepl("_xic_|_xir_", status_files)]

if (length(status_files) == 0) {
  stop("No base cohort_status files found in ", output_dir)
}

valid_files <- data.frame(file = status_files, stringsAsFactors = FALSE) %>%
  mutate(
    xih = as.numeric(gsub(".*_xih_([0-9]+).*", "\\1", basename(file))),
    xid = as.numeric(gsub(".*_xid_([0-9]+).*", "\\1", basename(file)))
  )

# --- Thresholds we want to compare ---
xi_values <- c(85, 90, 95)

# --- Build cumulative groups ---
psi_groups <- bind_rows(
  lapply(xi_values, function(th) {
    
    cat("Processing xi_d =", th, "\n")
    # 1. Load the specific status file for this xi_d
    target_file_df <- valid_files %>% filter(xih == xi_h_target, xid == th)
    if (nrow(target_file_df) == 0) {
      warning("No status data available for targeted xi_h=75 and xi_d=", th)
      return(NULL)
    }
    
    latest_status_file <- target_file_df$file[which.max(file.mtime(target_file_df$file))]
    status_df <- qs_read(latest_status_file, nthreads = N_QS_THREADS)
    
    # 2. Extract Psi_max directly from status
    patient_summary <- status_df %>%
      mutate(Psi_max = max_Psi)
    
    # 3. Filter and annotate groups
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
