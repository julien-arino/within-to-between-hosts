library(parallel)
library(deSolve)
library(reshape)
library(ggplot2)

source("functions_sim.R")

params = set_parameters()
IC = set_IC_2()

N = 10000
params.df = generate_params_patients(n = N, params = params)
patient_idx = 1:N

tictoc::tic()
if (TRUE) {
  # RUN IN PARALLEL
  # Detect number of cores, use all but 1
  no_cores <- detectCores()
  if (no_cores > 124) {
    no_cores = 124
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

saveRDS(SAVE, file = "sim.Rds")

# t = COHORT[[1]][,"time"]
# V = mat.or.vec(nr = length(t), nc = N)
# for (i in 1:N) {
#   V[,i] = COHORT[[i]][,"V"]
# }
# 
# df <- data.frame(x = t, V)
# # Long format
# df <- melt(df, id.vars = "x")
# 
# ggplot(df, aes(x = x, y = value, color = variable)) +
#   geom_line() + 
#   theme(legend.position = "none")
# ggsave(filename = "FIGS/viral-load.png",
#        width = 20, height = 15, units = "cm")
