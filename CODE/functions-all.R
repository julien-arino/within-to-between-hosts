## functions-all.R
#
# This file contains shared configurations and initial conditions used in the R scope.
#

suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

n_cores_logical <- parallel::detectCores(logical = FALSE)
if (is.na(n_cores_logical)) n_cores_logical <- parallel::detectCores()
N_QS_THREADS <- min(16, n_cores_logical)

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
    "ICU"  = "orange",
    "Dead" = "red3"
)
