# ==============================================================================
# File: plot-theme.R
# Description:
#   Sets a standardized ggplot2 theme for the project to ensure consistent
#   aesthetics across all plots. This script should be sourced at the
#   beginning of any plotting script. Matches the style used in CODE_clo.
# ==============================================================================

library(ggplot2)

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

cat("✅ Loaded standard ggplot2 theme (plot-theme.R)\n")
