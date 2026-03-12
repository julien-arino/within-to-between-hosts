library(deSolve)
library(lubridate)
library(future.apply)
source("functions-all.R")

OUTPUT_NAS <- "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
OUTPUT_LOCAL <- "../OUTPUT/"

# Run parallel?
PARALLEL <- TRUE

# Weight of sims in RAM may lead to explosion of RAM usage (or swapping). Set
# a maximum number of individuals to be simulated at once.
max_individuals_per_batch <- 5
# Number of individuals in the virtual cohort
N <- 10 # Testing number
# Number of sims needed to reach N
nb_batches <- ceiling(N / max_individuals_per_batch)

# Set general parameters and (common) initial conditions
params <- set_parameters()
IC <- set_IC()

# Generate virtual cohort by sourcing the central sample generation script.
# The script relies on 'N' being defined in the environment.
source("create-sample-for-sensitivity.R")
individuals_df <- individuals
# Convert DataFrame to a list of lists for future_lapply
individuals_list <- split(individuals_df, seq(nrow(individuals_df)))
individuals_list <- lapply(individuals_list, as.list)
writeLines("Starting computations")
if (PARALLEL) {
  # Prepare parallel processing environment
  if (availableCores() >= 64) {
    # We're on the cluster.. Be reasonable to avoid overheating
    nb_cores_to_use <- 32
  } else {
    # We're on something smaller, reserve 2 cores to avoid sluggish
    nb_cores_to_use <- availableCores() - 2
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
date_time_start <-
  format(now(tzone = "UTC"), "%Y%m%d-%H%M%S")

tictoc::tic()

# Run computation via future_lapply
COHORT <-
  future_lapply(
    X = individuals_list,
    FUN = function(p) {
      run_one_individual(
        params_tmp = p,
        IC = IC
      )
    },
    future.seed = TRUE
  )

tictoc::toc()

# Compute R0 for the batch individuals
individuals_df <- compute_R0_cohort(individuals_df)

# Start preparing the save variable
SAVE <- list()
SAVE$parameters <- individuals_df
# Add IC and results to save variable
SAVE$IC <- IC
SAVE$cohort <- COHORT

writeLines("Saving results")
saveRDS(SAVE,
  file = sprintf(
    "%s/sim_P%07d_DT%s.Rds",
    OUTPUT_LOCAL,
    N,
    date_time_start
  )
)
writeLines("Done!")
