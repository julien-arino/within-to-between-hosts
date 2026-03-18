# ============================================================
# File: plot-Figure-C5a-max-tissue-damage-fct-time-to-max-damage.R
# Purpose: Mean Psi-max vs Time to maximum tissue damage
# ============================================================

suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))

SHOW_TITLES <- FALSE

# ------------------------------------------------------------
# 1. Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-C5a-max-tissue-damage-fct-time-to-max-damage.R ...\n\n")

# Set project root automatically relative to the .git tracking directory
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

output_dir <- file.path(project_dir, "OUTPUT")
files <- list.files(output_dir, pattern = "^cohort_truncated_state_.*\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No truncated/censored cohort file found in ", output_dir)
}

latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

# --- Load your cohort data ---
cohort_df <- qs_read(latest_file)

# --- Compute Psi_max and tau_d per patient ---
death_summary <- cohort_df %>%
  group_by(individual_id, status) %>%
  summarise(
    psi_max = max(Psi, na.rm = TRUE),
    tau_d   = time[which.max(Psi)],
    .groups = "drop"
  )

# --- Thresholds to plot ---
xi_values <- c(85, 90, 95)

# --- Compute mean ± SD psi_max vs tau_d for each xi ---
psi_tau_list <- lapply(xi_values, function(xi){
  
  # We KEEP ONLY patients that exceed the threshold
  dead_summary_xi <- death_summary %>%
    filter(psi_max >= xi)
  
  if(nrow(dead_summary_xi)==0) return(NULL)
  
  psi_tau_xi <- dead_summary_xi %>%
    group_by(tau_d) %>%
    summarise(
      psi_mean = mean(psi_max),
      psi_sd   = sd(psi_max),
      .groups="drop"
    ) %>%
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
