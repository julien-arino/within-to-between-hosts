library(parallel)

source("functions_all.R")

NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
# Load the file. Let's do the latest one for now..
cohort = readRDS(sprintf("%s/sim_10000invididuals_2023_07_17-21_10_34.Rds",
                      NAS))
COHORT = cohort$cohort

# Make a list with each entry the state variables. Assumes all 
# individuals have the same time output and state variables
# Shift time to start at tau_I, because up to then, we are in the initial
# time delay part.
STATE_VARS = list()
for (n in colnames(COHORT[[1]])) {
  if (n == "time") {
    STATE_VARS[[n]] = COHORT[[1]][,"time"]
    idx_tauI = which(STATE_VARS$time==cohort$parameters$tau_I[1])[1]
    STATE_VARS[[n]] = 
      STATE_VARS[[n]][idx_tauI:length(STATE_VARS[[n]])] -
      cohort$parameters$tau_I[1]
  } else { #if ((n != "A_I") && (n != "A_R")) {
    STATE_VARS[[n]] = mat.or.vec(nr = dim(COHORT[[1]])[1], 
                                 nc = length(COHORT))
    for (i in 1:length(COHORT)) {
      STATE_VARS[[n]][,i] = COHORT[[i]][,n]
    }
    STATE_VARS[[n]] = STATE_VARS[[n]][idx_tauI:dim(STATE_VARS[[n]])[1],]
  }
}

# Create summary view for state variables
names_variables = setdiff(colnames(COHORT[[1]]), 
                          c("time", "A_I", "A_R"))
names_variables = setdiff(colnames(COHORT[[1]]), 
                          c("time"))
SUMMARIES = list()
for (n in names_variables) {
  SUMMARIES[[n]] = t(apply(STATE_VARS[[n]], 1, 
                           function(x) 
                             quantile(x, 
                                      probs = c(0.025,0.25,0.5,0.75,0.975))))
  SUMMARIES[[n]] = as.data.frame(SUMMARIES[[n]])
  SUMMARIES[[n]]$mean = apply(STATE_VARS[[n]], 1, mean)
}
