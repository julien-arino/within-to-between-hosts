# Create the virtual cohort: set parameters for all individuals in the cohort,
# save the list of parameters in a (potentially) shared directory for all
# nodes to be able to process the same list.

# Needed for sim date
library(lubridate)

# Read in file with all functions used
source("functions_all.R")

# Set directory where the results are located, typically not synced to GitHub
# but available on local NAS for sharing between nodes.
NAS_small_OUTPUT = "/home/jarino/NAS-small-OUTPUT/within-to-between-hosts"

# Set general parameters and (common) initial conditions
params = set_parameters()
IC = set_IC()

# Number of patients in the virtual cohort
N = 100000
# Generate virtual cohort. Do it in one go, even if we split sims, so that any
# scheme (LHS, etc.) applies to the entire cohort, not each batch of patients.
patients_df = generate_params_patients(n = N, params = params)

# Record date-time at start to have common file name
date_time_start = 
  format(now(tzone = "UTC"), "D%Y%m%d-T%H%M%SZ")

# Start preparing the save variable
SAVE = list()
SAVE$parameters = patients_df
# Add IC and results to save variable
SAVE$IC = IC

# Save the results
saveRDS(SAVE,
        file = sprintf("%s/params_P%07d_%s.Rds",
                       NAS_small_OUTPUT,
                       N, date_time_start))
