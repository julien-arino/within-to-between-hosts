# ============================================================
# File: plot-Figure-C5a-max-tissue-damage-fct-time-to-max-damage.R
# Purpose: Mean Psi-max vs Time to maximum tissue damage
# ============================================================

suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))

SHOW_TITLES <- FALSE

# ------------------------------------------------------------
# 1. Setup paths and load dynamically matched cohort_status files
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-C5a-max-tissue-damage-fct-time-to-max-damage.R ...\n\n")

# Set project root automatically relative to the .git tracking directory
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")

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

# --- Thresholds to plot ---
xi_values <- c(85, 90, 95)

# --- Compute mean ± SD psi_max vs tau_d for each xi ---
psi_tau_list <- lapply(xi_values, function(xi){
  
  cat("Processing xi_d =", xi, "\n")
  # 1. Load the specific status file for this xi_d
  target_file_df <- valid_files %>% filter(xih == xi_h_target, xid == xi)
  if (nrow(target_file_df) == 0) {
    warning("No status data available for targeted xi_h=75 and xi_d=", xi)
    return(NULL)
  }
  
  latest_status_file <- target_file_df$file[which.max(file.mtime(target_file_df$file))]
  status_df <- qs_read(latest_status_file)
  
  # 2. Extract psi_max and tau_d directly from status
  dead_summary_xi <- status_df %>%
    mutate(
      psi_max = max_Psi,
      tau_d   = tau_max_Psi
    ) %>%
    filter(psi_max >= xi)
  
  if(nrow(dead_summary_xi)==0) return(NULL)
  
  # 3. Bin tau_d to integers to compute mean and sd for the envelope
  psi_tau_xi <- dead_summary_xi %>%
    mutate(tau_d_bin = round(tau_d)) %>%
    group_by(tau_d_bin) %>%
    summarise(
      psi_mean = mean(psi_max, na.rm=TRUE),
      psi_sd   = sd(psi_max, na.rm=TRUE),
      n        = n(),
      .groups="drop"
    ) %>%
    # Use the bin as the x-axis value (tau_d), and only keep points with at least 2 patients to have a valid sd
    rename(tau_d = tau_d_bin) %>%
    filter(n > 1) %>%
    mutate(xi_d = factor(xi))
  
  psi_tau_xi
})

psi_tau_all <- bind_rows(psi_tau_list)

# Colors
cols_line <- c("85"="red1", "90"="#ef8a62", "95"="darkred")
cols_fill <- c("85"="red1", "90"="#ef8a62", "95"="darkred")

# --- PLOT ---
p_multi <- ggplot(psi_tau_all,
                  aes(x=tau_d, y=psi_mean,
                      color=xi_d, fill=xi_d)) +
  
  geom_ribbon(aes(ymin = pmax(psi_mean - psi_sd, 0),
                  ymax = pmin(psi_mean + psi_sd, 100)),
              alpha=0.20, color=NA) +
  
  geom_line(linewidth=1.2) +
  
  # thresholds
  geom_hline(data=data.frame(xi_d=factor(xi_values), y=xi_values),
             aes(yintercept=y, color=xi_d),
             linetype="dashed", linewidth=0.9) +
  
  scale_color_manual(values=cols_line) +
  scale_fill_manual(values=cols_fill) +
  
  labs(
    title = if(SHOW_TITLES) "Mean Psi-max vs Time of Death for Each Threshold" else NULL,
    x = expression("Time to maximum tissue damage"~tau[Psi[max]]^d),
    y = expression("Maximum tissue damage"~Psi[max]),
    color = expression(xi^d~"(%)"),
    fill  = expression(xi^d~"(%)")
  ) +
  coord_cartesian(xlim=c(0,20), ylim=c(80,100)) +
  theme_minimal(base_size=14) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    legend.position="right"
  )

print(p_multi)

# --- Save both PNG and PDF versions ---
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-C5a-max-tissue-damage-fct-time-to-max-damage.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-C5a-max-tissue-damage-fct-time-to-max-damage.png")

ggsave(
  filename = out_pdf,
  plot = p_multi,
  width = 25, height = 15, units = "cm", dpi = 300
)

ggsave(
  filename = out_png,
  plot = p_multi,
  width = 25, height = 15, units = "cm", dpi = 300
)

cat("✅ Plot saved to", out_pdf, "and", out_png, "\n")
