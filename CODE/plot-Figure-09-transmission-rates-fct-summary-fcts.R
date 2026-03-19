#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(future.apply)
}))

# ------------------------------------------------------------
# 1. Automatically find the latest files
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-09-transmission-rates-fct-summary-fcts.R ...\n\n")

project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

# --- Find and load pre-computed distributions ---
dist_files <- list.files(output_dir, pattern = "^cohort_distributions_P.*_beta\\.qs$", full.names = TRUE)
if (length(dist_files) == 0) stop("No beta distributions file found in OUTPUT")
latest_dist <- dist_files[which.max(file.mtime(dist_files))]
cat("Loading newest beta distributions:", basename(latest_dist), "\n")
beta_overall <- qs_read(latest_dist, nthreads = N_QS_THREADS)

# ------------------------------------------------------------
# 7. Plot
# ------------------------------------------------------------
p <- ggplot(beta_overall, aes(x=time)) +
  
  geom_ribbon(
    aes(
      ymin = pmax(beta_mean-beta_sd,0),
      ymax = beta_mean+beta_sd,
      fill="mean ± SD"
    ),
    alpha=0.15,
    color=NA
  )+
  
  geom_line(
    aes(y=beta_mean,color="mean"),
    linewidth=1.3
  )+
  
  geom_line(
    aes(y=beta_q10,color="Q10"),
    linetype="22",
    linewidth=0.7
  )+
  
  geom_line(
    aes(y=beta_median,color="median"),
    linewidth=1
  )+
  
  geom_line(
    aes(y=beta_q90,color="Q90"),
    linetype="22",
    linewidth=0.7
  )+
  
  scale_color_manual(
    breaks = c("mean","Q10","median","Q90"),
    values = c(
      "mean"   = "dodgerblue4",
      "Q10"    = "black",
      "median" = "magenta4",
      "Q90"    = "orange"
    )
  )+
  
  scale_fill_manual(
    values=c("mean ± SD"="dodgerblue4")
  )+
  
  coord_cartesian(xlim=c(0,80))+
  
  labs(
    x="Time (days)",
    y="Transmission rate",
    color="",
    fill=""
  )+
  
  theme_minimal(base_size=14)

print(p)

# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-09-transmission-rates-fct-summary-fcts.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-09-transmission-rates-fct-summary-fcts.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm")
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm")

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")
