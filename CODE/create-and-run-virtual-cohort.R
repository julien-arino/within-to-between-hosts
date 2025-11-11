library(deSolve)
library(lubridate)
library(future.apply)

source("functions-all.R")

OUTPUT_NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
OUTPUT_LOCAL = "OUTPUT/"

# Run parallel?
PARALLEL = TRUE

# Weight of sims in RAM may lead to explosion of RAM usage (or swapping). Set
# a maximum number of individuals to be simulated at once.
max_patients_per_batch = 5000
# Number of patients in the virtual cohort
N = 1000000
# Number of sims needed to reach N
nb_batches = ceiling(N / max_patients_per_batch)

# Set general parameters and (common) initial conditions
params = set_parameters()
IC = set_IC()

# Generate virtual cohort. Do it in one go, even if we split sims, so that any
# scheme (LHS, etc.) applies to the entire cohort, not each batch of patients.
tmp_df = generate_params_patients(n = N, params = params)
# Needed for the parallel call
tmp_idx = 1:N
# Now prepare batches
patients_df = list()
patients_idx = list()
for (b in 1:nb_batches) {
  start_current_batch = (b-1)*max_patients_per_batch+1
  end_current_batch = b*max_patients_per_batch
  patients_df[[b]] = tmp_df [start_current_batch:end_current_batch,]
  patients_idx[[b]] = tmp_idx[start_current_batch:end_current_batch]
}

writeLines("Starting computations")
if (PARALLEL) {
  # Prepare parallel processing environment
  if (availableCores() >= 64) {
    # We're on the cluster.. Be reasonable to avoid overheating
    nb_cores_to_use = 32
  } else {
    # We're on something smaller, reserve 2 cores to avoid sluggish
    nb_cores_to_use = availableCores()-2
  }
  if (Sys.getenv("RSTUDIO") == "1") {
    # If running from within Rstudio
    writeLines("Running from Rstudio, using plan(multisession)")
    plan(multisession, workers = nb_cores_to_use)
  } else {
    # If running from the command line
    writeLines("Running from the command line, using plan(multicore)")
    plan(multicore, workers = nb_cores_to_use)
  }
} else {
  plan(sequential)
}

# Record date-time at start to have common file name
date_time_start = 
  format(now(tzone = "UTC"), "%Y%m%d-%H%M%S")
# Run for all the batches
for (b in 1:nb_batches) {
  writeLines(paste0("Starting batch ", b, " out of ", nb_batches))
  tictoc::tic()
  # Take the patients for this batch
  patients = patients_df[[b]]
  # Indices for this batch
  patients_idx_this_batch = patients_idx[[b]]
  # Run computation
  COHORT = 
    future_lapply(
      X = patients_idx_this_batch, 
      FUN = function(x) 
        run_one_patient(idx = x,
                        patients = patients,
                        IC = IC),
      future.seed = TRUE
    )

  tictoc::toc()
  
  # Start preparing the save variable
  SAVE = list()
  SAVE$parameters = patients
  # Add IC and results to save variable
  SAVE$IC = IC
  SAVE$cohort = COHORT
  
  writeLines("Saving results")
  saveRDS(SAVE, 
          file = sprintf("%s/sim_P%07d_DT%s_part-%03d-of-%03d.Rds",
                         OUTPUT_LOCAL,
                         N, 
                         date_time_start,
                         b, nb_batches))
  # Clean up to avoid a run with previous run in RAM
  rm(COHORT)
  rm(SAVE)
}

