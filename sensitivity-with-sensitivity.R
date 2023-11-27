library(sensitivity)

source("functions_all.R")

rhs_within_host_deSolve = function(t, x, p) {
  # Set stuff that doesn't change
  avo=6.02214e23
  MM_F = 19000
  R_F_T = 1000
  R_F_I = 1300
  A_F = as.numeric((MM_F/avo) *
                     (R_F_I+R_F_T) *
                     (1/5000)*(10^9*1e12))
  V0 = 1
  S0 = 0.16
  I0 = 0
  R0 = 0 
  # Rest of the right hand side
  with(as.list(c(x, p)),{
    # Variables are (in order) V,S,I,R,D,F_U,F_B
    # So to get the lagvalue values...
    if (t<tau_I) {
      V_t = V0
      S_t = S0
      I_t = 0
      R_t = 0
    } else {
      V_t = lagvalue(t-tau_I,1)
      S_t = lagvalue(t-tau_I,2)
      I_t = lagvalue(t-tau_I,3)
      R_t = lagvalue(t-tau_I,4)
    }
    dV = p*I-d_V*V
    dS = lambda_S*(1-(S+I+D+R)/S_max)*S-beta*S*V
    dI = beta*S_t*V_t*(1-F_B/(epsilon_FI+F_B))*A_I-d_I*I
    dR = lambda_S*(1-(S+I+D+R)/S_max)*R +
      beta*S_t*V_t*(F_B/(epsilon_FI+F_B))*A_R
    dD = d_I*I-d_D*D
    dF_U = psi_F_prod+p_FI*I/(I+eta_FI) - 
      k_lin_f*F_U-k_B_F*((T_star+I)*A_F-F_B)*F_U+k_U_F*F_B
    dF_B = -k_int_f*F_B+k_B_F*((T_star+I)*A_F-F_B)*F_U-k_U_F*F_B
    dA_I = delta*A_I*(I_t-I)
    dA_R = delta*A_R*(R_t-R)
    return(list(c(dV, dS, dI, dR, dD, dF_U, dF_B,dA_I,dA_R)))
  })
}

# Given a patient index idx, a parameters data frame patients and 
# initial conditions IC, run the simulation of the within host model for
# this patient
# Return only the indicators of that patient (V_max, etc.)
run_one_patient_indicators = function(idx = 1, 
                                      patients, 
                                      IC) {
  writeLines(paste0("patient index = ", idx))
  params_tmp = patients[which(patients$ID == idx),]
  params_tmp = add_IC_to_params(params_tmp, IC)
  times <- c(seq(0, ceiling(params_tmp$tau_I), by = 0.01), 
             seq(ceiling(params_tmp$tau_I), 200, by = 0.1))
  tmp <- dede(y = IC, 
               times = times, 
               func = rhs_within_host_deSolve, 
               parms = params_tmp)
  OUT = list()
  OUT$V_max = max(tmp[, "V"])
  O
  return(OUT)
}

# Create response for all computed indicators, returning a data frame for each
# different point in parameter space provided as argument.
value_indicators = function(params_change, 
                            params_fixed, 
                            t_f = 200,
                            ncpus = 60,
                            parallel = TRUE) {
  writeLines("Creating list with parameters")
  col_names = colnames(params_change)
  PARS = list()
  for (i in 1:dim(params_change)[1]) {
    PARS[[i]] = list()
    for (j in 1:dim(params_change)[2]) {
      PARS[[i]][[col_names[j]]] = params_change[i,j]
    }
  }
  writeLines("Setting up cluster")
  if (parallel) {
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
  writeLines("Starting simulations")
  writeLines("Assembling results")
}


pars = data.frame(params = 
                    c("beta",
                      "d_D",
                      "d_I",
                      "d_V",
                      "delta",
                      "epsilon_FI",
                      "eta_FI",
                      "k_B_F",
                      "k_int_f",
                      "k_lin_f",
                      "k_U_F",
                      "lambda_S",
                      "p",
                      "p_FI",
                      "psi_F_prod",
                      "S_max",
                      "tau_I",
                      "T_star"))
tmp = matrix(c( 
  0.1, 0.5, # beta
  3, 15, # d_D
  0.05, 0.15, # d_I
  5, 15, # d_V
  0.05, 0.15, # delta
  1e-5, 1e-3, # epsilon_FI
  0.001, 0.05, # eta_FI
  0.001, 0.05, # k_B_F
  10, 20, # k_int_f
  10, 20, # k_lin_f
  2, 10, # k_U_F
  0.5, 1, # lambda_S
  200, 500, # p
  1, 4, # p_FI
  0.1, 0.5, # psi_F_prod
  0.1, 0.2, # S_max
  0.1, 0.3, # tau_I
  1e-5, 1e-3), # T_star)
  nc = 2, byrow = TRUE)
pars$min = tmp[,1]
pars$max = tmp[,2]

# Set all parameters (including ones we don't change)
params_all = set_parameters()
# The initial values of the state variables:
IC = set_IC()

# a 100-sample with X1 ~ U(0.5, 1.5)
# X2 ~ U(1.5, 4.5)
# X3 ~ U(4.5, 13.5)
# library(boot)
# n <- 1000
# X <- data.frame(X1 = runif(n, 0.5, 2.5),
#                 X2 = runif(n, 1.5, 4.5),
#                 X3 = runif(n, 4.5, 13.5))
# # linear model : Y = X1^2 + X2 + X3
# y <- with(X, X1^2 + X2 + X3)
# # sensitivity analysis
# x <- pcc(X, y, nboot = 100)
# print(x)
# plot(x)
# library(ggplot2)
# ggplot(x)
# ggplot(x, ylim = c(-1.5,1.5))
# x <- pcc(X, y, semi = TRUE, nboot = 100)
# print(x)
# plot(x)
