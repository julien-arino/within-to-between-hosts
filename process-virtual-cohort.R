source("functions_all.R")

NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
# Load the file. Let's do the latest one for now..
cohort = readRDS(sprintf("%s/sim_50000invididuals_2023_07_18-07_44_19.Rds",
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

# Make a matrix with S_max values as columns to be able to compute Phi
S_max_tmp = mat.or.vec(nr = dim(STATE_VARS$S)[1],
                       nc = dim(STATE_VARS$S)[2])
for (c in 1:dim(STATE_VARS$S)[2]) {
  S_max_tmp[,c] = rep(cohort$parameters$S_max[c],
                      dim(STATE_VARS$S)[1])
}
STATE_VARS$Phi = (STATE_VARS$S+STATE_VARS$R-S_max_tmp) / S_max_tmp * 100

# Create summary view for state variables
names_variables = setdiff(c(colnames(COHORT[[1]]),"Phi"), 
                          c("time", "A_I", "A_R"))
SUMMARIES = list()
for (n in names_variables) {
  SUMMARIES[[n]] = t(apply(STATE_VARS[[n]], 1, 
                           function(x) 
                             quantile(x, 
                                      probs = c(0.025,0.25,0.5,0.75,0.975))))
  SUMMARIES[[n]] = as.data.frame(SUMMARIES[[n]])
  SUMMARIES[[n]]$mean = apply(STATE_VARS[[n]], 1, mean)
}
# Which of the patients are in the different zones?
SUMMARIES$lung_loss = data.frame(max = apply(STATE_VARS$Phi, 2, min))
SUMMARIES$lung_loss$outcome = 4-findInterval(SUMMARIES$lung_loss$max,
                                             vec = c(-200,-85,-75,0))


# Prepare right censored data: for individuals who die, write down time 
# of death this happens (>=85% lung cell loss)
STATE_VARS$idx_censored = rep(NA, dim(STATE_VARS$V)[2])
STATE_VARS$time_censored = rep(NA, dim(STATE_VARS$V)[2])
for (i in 1:dim(STATE_VARS$V)[2]) {
  if (SUMMARIES$lung_loss$outcome[i] == 3) {
    STATE_VARS$idx_censored[i] = which(STATE_VARS$Phi[,i]<=-85)[1]
    STATE_VARS$time_censored[i] = STATE_VARS$time[STATE_VARS$idx_censored[i]]
  }
}
# Write down death times
SUMMARIES$time_of_death = 
  sort(STATE_VARS$time_censored[!is.na(STATE_VARS$time_censored)])
# Do the same for hospitalisations
STATE_VARS$idx_hosp = rep(NA, dim(STATE_VARS$V)[2])
STATE_VARS$time_hosp = rep(NA, dim(STATE_VARS$V)[2])
for (i in 1:dim(STATE_VARS$V)[2]) {
  if (SUMMARIES$lung_loss$outcome[i] == 2) {
    STATE_VARS$idx_hosp[i] = which(STATE_VARS$Phi[,i]<=-75)[1]
    STATE_VARS$time_hosp[i] = STATE_VARS$time[STATE_VARS$idx_hosp[i]]
  }
}
# Write down hospitalisation times
SUMMARIES$time_of_hospitalisation = 
  sort(STATE_VARS$time_hosp[!is.na(STATE_VARS$time_hosp)])
# Now recoveries. Suppose Phi gets larger than some threshold theta
threshold = -15
idx_uncensored = which(is.na(STATE_VARS$idx_censored))
min_phi = apply(STATE_VARS$Phi, 2, min)
idx_min_phi = c()
idx_cross_threshold = c()
for (i in 1:dim(STATE_VARS$Phi)[2]) {
  idx_min_phi = c(idx_min_phi, which(STATE_VARS$Phi[,i] == min_phi[i]))
  idx_tmp = idx_min_phi[i]:dim(STATE_VARS$Phi)[1]
  idx_cross_threshold = c(idx_cross_threshold,
                          which(STATE_VARS$Phi[idx_tmp,i] >= threshold)[1])
}
SUMMARIES$time_of_recovery = 
  sort(STATE_VARS$time[idx_cross_threshold[idx_uncensored]])
#SUMMARIES$time_of_recovery = 
#  SUMMARIES$time_of_recovery[which(SUMMARIES$time_of_recovery > 0)]
