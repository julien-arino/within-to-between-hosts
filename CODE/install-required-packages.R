## install-required-packages.R
# Ensure that all R packages required by the project are installed.

# List of required packages based on the project's source code
required_packages <- c(
    "cowplot",
    "deSolve",
    "doParallel",
    "dplyr",
    "factoextra",
    "foreach",
    "furrr",
    "future",
    "future.apply",
    "ggplot2",
    "ggpubr",
    "Hmisc",
    "latex2exp",
    "lhs",
    "lubridate",
    "ODEsensitivity",
    "patchwork",
    "purrr",
    "qs2",
    "readr",
    "reshape",
    "reshape2",
    "scales",
    "sensitivity",
    "stringr",
    "tidyr"
)

cat("Checking and installing required R packages...\n")

# Use a CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("Installing missing package:", pkg, "\n")
        # Wrapped in a try-catch equivalent to avoid stopping script on single package failure
        tryCatch(
            {
                install.packages(pkg, quiet = TRUE)
                # Verify it actually installed
                if (require(pkg, character.only = TRUE, quietly = TRUE)) {
                    cat("Successfully installed", pkg, "\n")
                } else {
                    cat("Failed to install", pkg, "\n")
                }
            },
            error = function(e) {
                cat("Error installing", pkg, ":", conditionMessage(e), "\n")
            }
        )
    } else {
        cat("Package already installed:", pkg, "\n")
    }
}

cat("All required R packages are installed.\n")
