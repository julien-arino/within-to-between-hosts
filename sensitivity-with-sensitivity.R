library(deSolve)
library(parallel)
library(lubridate)
library(sensitivity)
library(ggplot2)
library(randtoolbox)
library(latex2exp)


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
  OUT$F_U_max = max(tmp[, "F_U"])
  OUT$F_B_max = max(tmp[, "F_B"])
  OUT$A_I_max = max(tmp[, "A_I"])
  OUT$A_R_max = max(tmp[, "A_R"])
  OUT$T_max_V = times[which.max(tmp[, "V"])]
  OUT$T_max_F_U = times[which.max(tmp[, "F_U"])]
  OUT$T_max_F_B = times[which.max(tmp[, "F_B"])]
  OUT$T_max_A_I = times[which.max(tmp[, "A_I"])]
  OUT$T_max_A_R = times[which.max(tmp[, "A_R"])]
  return(OUT)
}

# Create response for all computed indicators, returning a data frame for each
# different point in parameter space provided as argument.
value_indicators = function(params_change, 
                            params_fixed, 
                            t_f = 200,
                            ncpus = 60,
                            parallel = TRUE) {
  writeLines("Finalising parameters data frame")
  col_names = colnames(params_change)[2:dim(params_change)[2]]
  col_names_fixed = names(params_fixed)
  col_names_fixed = setdiff(col_names_fixed, col_names)
  for (i in 1:length(col_names_fixed)) {
    params_change[[col_names_fixed[i]]] = params_fixed[[col_names_fixed[i]]]
  }
  # Indices for this batch
  patients_idx = 1:dim(params_change)[1]
  if (parallel) {
    # RUN IN PARALLEL
    writeLines("Setting up cluster")
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
                  c("run_one_patient_indicators",
                    "add_IC_to_params",
                    "rhs_within_host_deSolve",
                    "IC"),
                  envir = .GlobalEnv)
    # Run computation in parallel
    writeLines("Starting simulations in parallel")
    COHORT = 
      parLapply(cl = cl, 
                X = patients_idx,
                fun = function(x) 
                  run_one_patient_indicators(idx = x,
                                             patients = params_change,
                                             IC = IC))
    # Stop cluster
    stopCluster(cl)
  } else {
    # Run computation sequentially
    writeLines("Going old school, running sequentially (debugging, probably)")
    COHORT = lapply(X = patients_idx,
                    FUN = function(x) 
                      run_one_patient_indicators(idx = x,
                                                 patients = params_change,
                                                 IC = IC))
  }
  # We now need to transform the cohort (list) into a data frame
  writeLines("Assembling results")
  COHORT = as.data.frame(as.matrix(do.call(rbind, COHORT)))
}

pars.df = data.frame(params = 
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
  100, 800, # p
  0.8, 4, # p_FI
  0.1, 0.5, # psi_F_prod
  0.1, 0.2, # S_max
  0.1, 0.3, # tau_I
  1e-5, 1e-3), # T_star)
  nc = 2, byrow = TRUE)
pars.df$min = tmp[,1]
pars.df$max = tmp[,2]

# Number of patients in the virtual cohort (sample size for the sensitivity)
N = 100000

# To use sensitivity::parameterSets, we need to convert the data frame to a list
pars.list = list()
for (i in 1:dim(pars.df)[1]) {
  pars.list[[pars.df$params[i]]] = c(pars.df$min[i], pars.df$max[i])
}
nb_pars = length(pars.list)
# pars.grid = parameterSets(par.ranges = pars.list, 
#                           samples = nb_pars, 
#                           method = "grid")
pars.sobol = parameterSets(par.ranges = pars.list, 
                           samples = N, 
                           method = "sobol")
pars.sobol = as.data.frame(pars.sobol)
pars.sobol = cbind(1:dim(pars.sobol)[1], pars.sobol)
colnames(pars.sobol) = c("ID", pars.df$params)

# Set all parameters (including ones we don't change)
params_all = set_parameters()
# The initial values of the state variables:
IC = set_IC()
# Generate virtual cohort. Do it in one go, even if we split sims, so that any
# scheme (LHS, etc.) applies to the entire cohort, not each batch of patients.
# tmp_df = generate_params_patients(n = N, params = params_all)
# Compute indicators for each patient
vals = value_indicators (params_change = pars.sobol, params_fixed = params_all)
# Compute the correlation between the indicators and the parameters
writeLines("Computing partial rank correlations - V_max")
tmp = as.numeric(vals$V_max)
x_V_max= pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
             rank = TRUE, semi = FALSE)
ggplot(x_V_max)
writeLines("Computing partial rank correlations - F_U_max")
tmp = as.numeric(vals$F_U_max)
x_F_U_max = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                rank = TRUE, semi = FALSE)
ggplot(x_F_U_max)
writeLines("Computing partial rank correlations - F_B_max")
tmp = as.numeric(vals$F_B_max)
x_F_B_max = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                rank = TRUE, semi = FALSE)
ggplot(x_F_B_max)
writeLines("Computing partial rank correlations - T_max_V")
tmp = as.numeric(vals$T_max_V)
x_T_max_V = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
               rank = TRUE, semi = FALSE)
ggplot(x_T_max_V)
writeLines("Computing partial rank correlations - T_max_F_U")
tmp = as.numeric(vals$T_max_F_U)
x_T_max_F_U = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                  rank = TRUE, semi = FALSE)
ggplot(x_T_max_F_U)
writeLines("Computing partial rank correlations - T_max_F_B")
tmp = as.numeric(vals$T_max_F_B)
x_T_max_F_B = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                  rank = TRUE, semi = FALSE)
ggplot(x_T_max_F_B)
# Combine values for aggregate plots. All values, then by type
PRCC = data.frame(
  names = names(pars.list),
  V_max = x_V_max$PRCC,
  F_U_max = x_F_U_max$PRCC,
  F_B_max = x_F_B_max$PRCC,
  T_max_V = x_T_max_V$PRCC,
  T_max_F_U = x_T_max_F_U$PRCC,
  T_max_F_B = x_T_max_F_B$PRCC
)
colnames(PRCC) = c("names", "V_max", "F_U_max", "F_B_max", "T_max_V", "T_max_F_U", "T_max_F_B")
PRCC_values = PRCC[,c("names", "V_max", "F_U_max", "F_B_max")]
PRCC_times = PRCC[,c("names", "T_max_V", "T_max_F_U", "T_max_F_B")]
# Order the parameters by decreasing absolute value of PRCC
order_PRCC_values = order(apply(abs(PRCC_values[2:4]), 1, max), 
                          decreasing = TRUE)
order_PRCC_times = order(apply(abs(PRCC_times[2:4]), 1, max), 
                         decreasing = TRUE)
PRCC_values = PRCC_values[order_PRCC_values,]
PRCC_times = PRCC_times[order_PRCC_times,]
# Reshape for plotting 
PRCC_values = reshape2::melt(PRCC_values, id.vars = "names")
PRCC_times = reshape2::melt(PRCC_times, id.vars = "names")
# Needed for x axis labels
x_labels = c(TeX("$\\beta$"), 
             TeX("$d_D$"), 
             TeX("$d_I$"), 
             TeX("$d_V$"), 
             TeX("$\\delta$"), 
             TeX("$\\epsilon_{FI}$"), 
             TeX("$\\eta_{FI}$"), 
             TeX("$k_{B_F}$"), 
             TeX("$k_{int_f}$"), 
             TeX("$k_{lin_f}$"), 
             TeX("$k_{U_f}$"), 
             TeX("$\\lambda_{S}$"),
             TeX("$p$"),
             TeX("$p_{FI}$"),
             TeX("$\\psi_{F_{prod}}$"),
             TeX("$S_{max}$"),
             TeX("$\\tau_{I}$"),
             TeX("$T_*$"))
# Do some scaling for the size of points
scaled_size_values = PRCC_values$value
max_value_values = max(abs(PRCC_values$value))
half_saturation_values = 0.01
for (i in 1:dim(PRCC_values)[1]) {
  if (PRCC_values$value[i] < 0) {
    scaled_size_values[i] = 
      PRCC_values$value[i]/(-half_saturation_values-max_value_values)
  } else {
    scaled_size_values[i] = 
      PRCC_values$value[i]/(half_saturation_values+max_value_values)
  }
}
scaled_size_times = PRCC_times$value
max_value_times = max(abs(PRCC_times$value))
half_saturation_times = 0.01
for (i in 1:dim(PRCC_times)[1]) {
  if (PRCC_times$value[i] < 0) {
    scaled_size_times[i] = 
      PRCC_times$value[i]/(-half_saturation_times-max_value_times)
  } else {
    scaled_size_times[i] = 
      PRCC_times$value[i]/(half_saturation_times+max_value_times)
  }
}
# Plot the values
ggplot(data = PRCC_values, aes(names, value, col = variable)) + 
  geom_point(size = scaled_size_values*12) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  ylim(-1,1) +
  scale_color_discrete(labels = c(TeX("$V_{max}$"), 
                                  TeX("$F_{U_{max}}$"), 
                                  TeX("$F_{B_{max}}$"))) +
  labs(color="Indicator") +
  scale_x_discrete(labels = x_labels[order_PRCC_values],
                   limits = names(pars.list)[order_PRCC_values]) +
  theme(legend.position = c(0.9,0.9))
ggsave(filename = "FIGS/sensitivities-PRCC-values.png",
       width = 20, height = 15, units = "cm")
# Plot the times
ggplot(data = PRCC_times, aes(names, value, col = variable)) + 
  geom_point(size = scaled_size_times*12) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  ylim(-1,1) +
  scale_color_discrete(labels = c(TeX("$t_{V_{max}}$"), 
                                  TeX("$t_{F_{U_{max}}}$"), 
                                  TeX("$t_{F_{B_{max}}}$"))) +
  labs(color="Indicator") +
  scale_x_discrete(labels = x_labels[order_PRCC_times],
                   limits = names(pars.list)[order_PRCC_times]) +
  theme(legend.position = c(0.9,0.9))
ggsave(filename = "FIGS/sensitivities-PRCC-times.png",
       width = 20, height = 15, units = "cm")

