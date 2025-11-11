# Run the virtual cohort

library(deSolve)
library(parallel)

source("functions_all.R")

# Directory where the parameters are located, typically not synced to GitHub
# but on local NAS for sharing between nodes.
NAS_small_OUTPUT = "/home/jarino/NAS-small-OUTPUT/within-to-between-hosts"

# Directory where the local results are located, typically not synced to GitHub.
LOCAL_OUTPUT = "/home/jarino/LOCAL-OUTPUT/within-to-between-hosts"

# Run parallel?
PARALLEL = FALSE

# Find a parameter file to process on the NAS
params_files = data.frame(
  file_name = list.files(path = NAS_small_OUTPUT,
                          pattern = glob2rx("params*.Rds"))
)
tmp = strsplit(params_files$file_name, "-")
params_files$cohort_size_str = 
  unlist(lapply(tmp, 
                function(x) x[2]))
params_files$cohort_size = 
  as.numeric(gsub("P", "", params_files$cohort_size_str))
params_files$date = 
  unlist(lapply(tmp, 
                function(x) x[3]))
params_files$time = 
  unlist(lapply(tmp, 
                function(x) x[4]))
params_files$time = gsub(".Rds", "", params_files$time)

# Load the parameters from the NAS
# For now, for simplicity, do the most recent
params_tmp = readRDS(
  sprintf("%s/%s",
          NAS_small_OUTPUT,
          params_files$file_name[dim(params_files)[1]])
)
# Set general parameters and (common) initial conditions
params_df_tmp = params_tmp$parameters
IC = params_tmp$IC

# Weight of sims in RAM may lead to explosion of RAM usage (or swapping). Set
# a maximum number of individuals to be simulated at once.
max_patients_per_batch = 500
# Number of patients in the virtual cohort
N = params_files$cohort_size[dim(params_files)[1]]
# Number of sims needed to reach N
nb_batches = ceiling(N / max_patients_per_batch)

# Needed for the parallel call
tmp_idx = 1:N
# Now prepare batches
patients_df = list()
patients_idx = list()
for (b in 1:nb_batches) {
  start_current_batch = (b-1)*max_patients_per_batch+1
  end_current_batch = b*max_patients_per_batch
  patients_df[[b]] = params_df_tmp[start_current_batch:end_current_batch,]
  patients_idx[[b]] = tmp_idx[start_current_batch:end_current_batch]
}

writeLines("Starting computations")
# Set up the cluster if needed
if (PARALLEL) {
  # RUN IN PARALLEL
  # Detect number of cores
  no_cores <- detectCores()-1
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


# Are we resuming an ongoing (interrupted) computation
# Pattern of the result files, if any
pattern = sprintf("sim-%s-%s-%s-part-",
                  params_files$cohort_size_str[dim(params_files)[1]], 
                  params_files$date[dim(params_files)[1]], 
                  params_files$time[dim(params_files)[1]])
# Get list of stuff we have already done (saved locally)
done_files = data.frame(
  files = list.files(LOCAL_OUTPUT, pattern)
)
if (dim(done_files)[1]>0) {
  tmp = strsplit(done_files$files, "-")
  done_files$part = as.numeric(unlist(lapply(tmp, function(x) x[6])))
  if (max(done_files$part)<nb_batches) {
    start = max(done_files$part)+1
    RUN_SIMS = TRUE
  } else {
    RUN_SIMS = FALSE
  }
} else {
  start = 1
  RUN_SIMS = TRUE
}
if (RUN_SIMS) {
  # Run for all the batches
  for (b in start:nb_batches) {
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
            file = sprintf("%s/sim-%s-%s-%s-part-%03d-of-%03d.Rds",
                           LOCAL_OUTPUT,
                           params_files$cohort_size_str[dim(params_files)[1]], 
                           params_files$date[dim(params_files)[1]], 
                           params_files$time[dim(params_files)[1]],
                           b, nb_batches))
    # Clean up to avoid a run with previous run in RAM
    rm(COHORT)
    rm(SAVE)
  }
  
  # Stop cluster
  if (PARALLEL) {
    stopCluster(cl)
  }
}