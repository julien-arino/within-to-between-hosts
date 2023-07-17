library(parallel)

source("functions_all.R")

# Make a list with each entry the state variables. Assumes all 
# individuals have the same time output and state variables
STATE_VARS = list()
for (n in colnames(COHORT[[1]])) {
  if (n == "time") {
    STATE_VARS[[n]] = COHORT[[1]][,"time"]
  } else if ((n != "A_I") && (n != "A_R")) {
    STATE_VARS[[n]] = mat.or.vec(nr = dim(COHORT[[1]])[1], 
                                 nc = length(COHORT))
    for (i in 1:length(COHORT)) {
      STATE_VARS[[n]][,i] = COHORT[[i]][,n]
    }
  }
}

# Create summary view for state variables
names_variables = setdiff(colnames(COHORT[[1]]), 
                          c("time", "A_I", "A_R"))
SUMMARIES = list()
for (n in names_variables) {
  SUMMARIES[[n]] = t(apply(STATE_VARS[[n]], 1, quantile))
  SUMMARIES[[n]] = as.data.frame(SUMMARIES[[n]])
  SUMMARIES[[n]]$mean = apply(STATE_VARS[[n]], 1, mean)
}
