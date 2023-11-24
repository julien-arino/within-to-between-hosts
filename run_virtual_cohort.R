library(deSolve)
library(dde)
library(parallel)
library(lubridate)

source("functions_all.R")

OUTPUT_NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
OUTPUT_LOCAL = "/home/jarino/OUTPUT_local/within-to-between-hosts"

params = set_parameters()
IC = set_IC()

N = 50000
params.df = generate_params_patients(n = N, params = params)
patient_idx = 1:N

writeLines("Starting computations")
tictoc::tic()
if (TRUE) {
  # RUN IN PARALLEL
  # Detect number of cores
  no_cores <- detectCores()
  if (no_cores > 124) {
    no_cores = 120
  }
  # Initiate cluster
  cl <- makeCluster(no_cores)
  # Export needed variables
  clusterEvalQ(cl,{
    library(deSolve)
  })
  clusterExport(cl,
                c("run_one_patient",
                  "add_IC_to_params",
                  "rhs_within_host_deSolve",
                  "params.df",
                  "IC"),
                envir = .GlobalEnv)
  # Run computation
  COHORT = 
    parLapply(cl = cl, 
              X = patient_idx,
              fun = function(x) run_one_patient(idx = x,
                                                params.df = params.df,
                                                IC = IC))
  # Stop cluster
  stopCluster(cl)
} else {
  # RUN SEQUENTIALLY
  COHORT = lapply(X = patient_idx,
                  FUN = function(x) run_one_patient(idx = x,
                                                    params.df = params.df,
                                                    IC = IC))
}
tictoc::toc()

SAVE = list()
SAVE$parameters = params.df
SAVE$IC = IC
SAVE$cohort = COHORT

writeLines("Saving results")
saveRDS(SAVE, 
        file = sprintf("%s/sim_P%07d_DT%s.Rds",
                       OUTPUT_NAS,
                       N, 
                       format(now(tzone = "UTC"), "%Y_%m_%d-%H_%M_%S")))
