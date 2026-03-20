#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-09-transmission-rates-fct-summary-fcts.R
# Description:
# ============================================================

# First things first: locate project directory and load helper functions
# and constants
project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

# Say what we are running and start the clock
start_time <- start_time_and_hello("plot-Figure-09-transmission-rates-fct-summary-fcts.R")

# Load libraries
load_libraries(c("ggplot2", "dplyr", "qs2", "tidyr"))

# ------------------------------------------------------------
# 1. Automatically find the latest files
# ------------------------------------------------------------

# --- Find and load pre-computed distributions ---
pattern_str <- "^cohort_distribution_filters_P.*_xic_[0-9]+_.*\\.qs$"
dist_files <- list.files(
  output_dir,
  pattern = pattern_str,
  full.names = TRUE
)

if (length(dist_files) == 0) {
  stop("No beta distributions file found in OUTPUT")
}

latest_dist <- dist_files[which.max(file.mtime(dist_files))]
cat("Loading newest distributions:", basename(latest_dist), "\n")

beta_list <- qs_read(latest_dist, nthreads = N_QS_THREADS)

if (!is.null(beta_list$rolling_average$beta_mean)) {
  beta_overall <- beta_list$rolling_average
} else {
  beta_overall <- beta_list$raw
}

# ------------------------------------------------------------
# 7. Plot
# ------------------------------------------------------------
p <- ggplot(beta_overall, aes(x=time)) +
  
  geom_ribbon(
    aes(
      ymin = pmax(beta_Q10, 0),
      ymax = beta_Q90,
      fill="Q10 to Q90"
    ),
    alpha=0.15,
    color=NA
  )+
  
  geom_line(
    aes(y=beta_mean,color="mean"),
    linewidth=1.3
  )+
  
  geom_line(
    aes(y=beta_Q10,color="Q10"),
    linetype="22",
    linewidth=0.7
  )+
  
  geom_line(
    aes(y=beta_median,color="median"),
    linewidth=1
  )+
  
  geom_line(
    aes(y=beta_Q90,color="Q90"),
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
    values=c("Q10 to Q90"="dodgerblue4")
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
out_pdf <- file.path(figs_dir, "Figure-09-transmission-rates-fct-summary-fcts.pdf")
out_png <- file.path(figs_dir, "Figure-09-transmission-rates-fct-summary-fcts.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm")
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm")

cat("\nFigures saved:\n  - ", out_pdf, "\n  - ", out_png, "\n")

print_end_time(start_time, "plot-Figure-09-transmission-rates-fct-summary-fcts.R")
