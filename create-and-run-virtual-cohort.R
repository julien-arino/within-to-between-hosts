library(deSolve)
library(parallel)
library(lubridate)

source("functions_all.R")

OUTPUT_NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
OUTPUT_LOCAL = "/home/jarino/OUTPUT_local/within-to-between-hosts"

# Run parallel?
PARALLEL = TRUE

# Weight of sims in RAM may lead to explosion of RAM usage (or swapping). Set
# a maximum number of individuals to be simulated at once.
max_patients_per_batch = 20000
# Number of patients in the virtual cohort
N = 100000
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
# Set up the cluster if needed
if (PARALLEL) {
  # RUN IN PARALLEL
  # Detect number of cores
  no_cores <- detectCores()
  if (no_cores > 124) {
    # Detect rich person's problem. 
    # (Could also recompile R setting the number of sockets higher than the 
    # default... not done here.)
    no_cores = 124
  }
  # Initiate cluster with parts shared across all batches
  cl <- makeCluster(no_cores)
  # Export needed variables
  clusterEvalQ(cl,{
    library(deSolve)
  })
  clusterExport(cl,
                c("run_one_patient",
                  "add_IC_to_params",
                  "rhs_within_host_deSolve",
                  "IC"),
                envir = .GlobalEnv)
} else {
  writeLines("We're running sequentially!")
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
  if (PARALLEL) {
    # RUN IN PARALLEL, set parameters specific to this batch.
    clusterExport(cl,
                  c("patients"),
                  envir = .GlobalEnv)
    # Run computation
    COHORT = 
      parLapply(cl = cl, 
                X = patients_idx_this_batch,
                fun = function(x) 
                  run_one_patient(idx = x,
                                  patients = patients,
                                  IC = IC))
  } else {
    # RUN SEQUENTIALLY
    writeLines("Going old school (debugging, probably)")
    COHORT = lapply(X = patients_idx_this_batch,
                    FUN = function(x) 
                      run_one_patient(idx = x,
                                      patients = patients,
                                      IC = IC))
  }
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

# Stop cluster
stopCluster(cl)
