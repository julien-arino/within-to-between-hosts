#!/usr/bin/env Rscript
# ============================================================
# File: plot-Psi_max-vs-Psi_max_time-Figure-04.R
# Goal:
#   1) Load truncated cohort results to get Psi_max & status per patient
#   2) Load process-virtual-cohort-results.csv to get R0_within and tau_Psi_max
#   3) Join and make figure (Psi_max vs t_Psi_max, split by R0_within<1 vs >=1)
# Outputs:
#   - FIGS/plot_Psi_max_vs_Psi_max_time_Figure_04.png and .pdf
# ============================================================

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(here)))

cat("\n\n>>> Running plot-Figure-04-Psi_max-vs-Psi_max_time.R ...\n\n")

# Set project root automatically relative to the .git tracking directory
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

dir.create(file.path(project_dir, "OUTPUT"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)

output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

# 1a) Find latest processed status file (xi_h=75, xi_d=85)
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xih_75_xid_85\\.qs$", full.names = TRUE)
if (length(status_files) == 0) stop("Cannot find status qs file in ", output_dir)
STATUS_FILE <- status_files[which.max(file.mtime(status_files))]

# ---------------------------
# 1) Load data
# ---------------------------
cat("Loading processed status:", basename(STATUS_FILE), "\n")
status_df <- qs_read(STATUS_FILE, nthreads = N_QS_THREADS)

# We don't need to merge; status_df already has all parameter metrics natively embedded!
summary_df <- status_df

# ---------------------------
# 2) Formatting for plotting (rename columns to match plot)
# ---------------------------
summary_df <- summary_df %>%
  mutate(
    R0 = R0_within,
    t_Psi_max = tau_max_Psi,
    Psi_max = max_Psi
  ) %>%
  filter(!is.na(R0))

summary_df$status <- factor(summary_df$status, levels = c("Mild", "ICU", "Dead"))

status_colors <- c(
  Mild = "dodgerblue4",
  ICU  = "orange",
  Dead = "red"
)

# ---------------------------
# Count and filter out values hitting the simulation time boundary (x >= 100)
# ---------------------------
dropped_count <- sum(summary_df$t_Psi_max >= 100, na.rm = TRUE)
summary_df <- summary_df %>% filter(t_Psi_max < 100)

# ---------------------------
# 5) Create plot
# ---------------------------
p2 <- ggplot(summary_df, aes(x = t_Psi_max, y = Psi_max)) +
  
  # Patients with R0 >= 1 (colored by status)
  geom_point(
    data = subset(summary_df, R0 >= 1),
    aes(color = status),
    alpha = 0.18,
    size = 0.6
  ) +
  
  # Patients with R0 < 1 (DARK GREEN)
  geom_point(
    data = subset(summary_df, R0 < 1),
    color = "darkgreen",
    alpha = 0.4,
    size = 0.8
  ) +
  
  coord_cartesian(xlim=c(0,100), ylim=c(0,100)) +
  
  scale_color_manual(
    values = status_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    x = expression("Time to maximum tissue damage"~tau[Psi[max]]~"(days)"),
    y = expression("Maximum tissue damage"~Psi[max]),
    parse = TRUE
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 2, alpha = 1)
    )
  )

p2 <- p2 +
  annotate("text",
           x = 1, y = 15,
           label = "R[0] < 1",
           size = 4,
           color = "darkgreen",
           parse = TRUE,
           fontface = "bold") +
  annotate("text",
           x = 45, y = 95,
           label = "R[0] >= 1",
           size = 4,
           color = "black",
           parse = TRUE,
           fontface = "bold")
print(p2)

# ---------------------------
# 6) Save figure
# ---------------------------
out_png <- file.path(project_dir, "FIGS", "Figure-04-Psi_max-vs-Psi_max_time.png")
out_pdf <- file.path(project_dir, "FIGS", "Figure-04-Psi_max-vs-Psi_max_time.pdf")

ggsave(filename = out_png, plot = p2, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(filename = out_pdf, plot = p2, width = 25, height = 15, units = "cm", dpi = 300)

cat("\nFigure 4 saved to:\n  -", out_png, "\n  -", out_pdf, "\n")
cat("Note:", dropped_count, "individuals who hit the simulation horizon (t >= 100 days) were omitted from this plot.\n")
