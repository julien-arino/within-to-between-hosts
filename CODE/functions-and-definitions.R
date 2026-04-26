#!/usr/bin/env Rscript
## functions-and-definitions.R
#
# This file contains shared configurations and initial conditions used in the R scope.
#

# Load libraries suppressing the output messages
# Beware: this will suppress ALL output messages from the loaded packages, so if you have
# any warnings or messages from the loaded packages, you will not see them.
# Run the install-required-packages.R script in order to see warnings, if need be.
load_libraries <- function(packages) {
  for (pkg in packages) {
    suppressWarnings(suppressPackageStartupMessages(library(pkg, character.only = TRUE)))
  }
}

# Start timer and print message
start_time_and_hello <- function(script_name) {
  cat(sprintf("\n\n>>> Starting %s ...\n\n", script_name))
  return(Sys.time())
}

# Timing formatting
print_end_time <- function(start_time, script_name) {
  end_time <- Sys.time()
  time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))
  mins <- floor(time_diff / 60)
  secs <- round(time_diff %% 60, 1)
  cat(sprintf("\n\n>>> Finished %s in %d minutes and %.1f seconds ✅\n\n", script_name, mins, secs))
}

set_IC <- function() {
  csv_path <- "data-initial-conditions.csv"
  if (exists("project_dir")) {
    csv_path <- file.path(project_dir, "CODE", "data-initial-conditions.csv")
  }

  if (!file.exists(csv_path)) {
    csv_path <- "data-initial-conditions.csv"
  }

  ic_data <- read.csv(csv_path, stringsAsFactors = FALSE)
  IC <- setNames(ic_data$value, ic_data$variable)
  return(IC)
}


# Define project preferred base size (often 14 or 15 in this project)
base_font_size <- 14

load_libraries("ggplot2")

# Set standard minimal theme with baseline customization
# commonly used throughout overleaf-within-to-between-hosts-new/CODE/CODE_clo/
theme_set(
  theme_minimal(base_size = base_font_size) +
    theme(
      # Bold and centered titles
      plot.title = element_text(hjust = 0.5, face = "bold"),

      # Clean axis titles
      legend.position = "right"
    )
)

# ==============================================================================
# Standard Project Colors
# ==============================================================================

# Commonly used discrete colors for consistency across plots
project_colors <- c(
  blue   = "dodgerblue4",
  orange = "orange",
  red    = "red3",
  grey   = "gray60",
  black  = "black"
)

# Standard mapping for Clinical Status
status_cols <- c(
  "Mild" = "dodgerblue3",
  "ICU" = "orange",
  "Dead" = "red3"
)


###
### CONSTANTS
###
# Output path
output_dir <- "/mnt/NAS-small-OUTPUT/within-to-between-hosts"
# Figure path
figs_dir <- file.path(project_dir, "FIGS")
# Number of cores for qs2 (parallel)
n_cores_logical <- parallel::detectCores(logical = FALSE)
if (is.na(n_cores_logical)) n_cores_logical <- parallel::detectCores()
N_QS_THREADS <- min(16, n_cores_logical)
