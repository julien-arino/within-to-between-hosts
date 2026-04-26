#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-C2b-unboundIFN-vs-viral-load.R
# Description:
#   Plot max(F_U) against max(V) for the cohort
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-C2b-unboundIFN-vs-viral-load.R")

load_libraries(c("ggplot2", "dplyr", "qs2", "here", "ggrastr"))

# 1a. Load simulation parameters (contains the max metrics)
param_files <- list.files(output_dir, pattern = "^cohort_sim_parameters_P.*\\.qs$", full.names = TRUE)
if (length(param_files) == 0) stop("No cohort_sim_parameters file found in ", output_dir)
latest_param_file <- param_files[which.max(file.mtime(param_files))]
cat("Loading newest cohort parameters:", basename(latest_param_file), "\n")
cohort_params <- qs_read(latest_param_file, nthreads = N_QS_THREADS)

# 1b. Load corresponding patient status summary
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xih_75_xid_85\\.qs$|^cohort_status_P.*_xid_85\\.qs$", full.names = TRUE)
if (length(status_files) == 0) stop("No baseline cohort_status file found in ", output_dir)
latest_status_file <- status_files[which.max(file.mtime(status_files))]
cat("Loading newest cohort baseline statuses:", basename(latest_status_file), "\n")
cohort_status <- qs_read(latest_status_file, nthreads = N_QS_THREADS)

# 2. Merge parameters and statuses
cohort_merged <- cohort_status %>%
  select(ID, status) %>%
  inner_join(cohort_params, by = "ID")

cat("Merged dataset contains", nrow(cohort_merged), "individuals\n")

# 3. Define colors
status_colors <- c(Mild = "dodgerblue4", ICU = "orange", Dead = "red")
cohort_merged$status <- factor(cohort_merged$status, levels = c("Mild", "ICU", "Dead"))

# 4. Plot max(V) vs max(F_U)
p_FU <- ggplot(cohort_merged, aes(x = max_V, y = max_F_U, color = status)) +
  geom_point_rast(alpha = 0.4, size = 1) +
  scale_color_manual(values = status_colors) +
  labs(
    x = "Maximum viral load",
    y = "Maximum unbound IFN",
    color = "Status"
  ) +
  theme_minimal(base_size = 20) + 
  theme(legend.position = "top")

print(p_FU)

# 5. Save to PDF & PNG
out_pdf <- file.path(figs_dir, "Figure-C2b-unboundIFN-vs-viral-load.pdf")
out_png <- file.path(figs_dir, "Figure-C2b-unboundIFN-vs-viral-load.png")

ggsave(filename = out_pdf, plot = p_FU, width = 6, height = 5, units = "in", dpi = 300)
ggsave(filename = out_png, plot = p_FU, width = 6, height = 5, units = "in", dpi = 300)

cat("Plots saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
print_end_time(start_time, "plot-Figure-C2b-unboundIFN-vs-viral-load.R")
