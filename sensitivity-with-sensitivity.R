library(deSolve)
library(parallel)
library(lubridate)
library(sensitivity)
library(randtoolbox)

source("functions-all.R")

OUTPUT_NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
OUTPUT_LOCAL = "/home/jarino/OUTPUT_local/within-to-between-hosts"

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
N = 500

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

# Set all parameters (including ones we don't change, if any)
params_all = set_parameters()
# The initial values of the state variables:
IC = set_IC()

# Compute indicators for each patient in the virtual cohort.
# Can be split runs to avoid overloading RAM if too many patients.
# Weight of sims in RAM may lead to explosion of RAM usage (or swapping). 
# Set a maximum number of individuals to be simulated at once.
max_patients_per_batch = 300
# Number of batches needed to reach cohort size
nb_batches = ceiling(N / max_patients_per_batch)
# The current time, so files have a common name
curr_TD = format(now(tzone = "UTC"), "%Y_%m_%d-%H_%M_%S")
# Run each batch, save it, clean up..
for (b in 1:nb_batches) {
  writeLines(paste0("Starting batch ", b, " out of ", nb_batches))
  tictoc::tic()
  # The current batch
  idx_start_batch = (b - 1) * max_patients_per_batch + 1
  idx_end_batch = min(b * max_patients_per_batch, N)
  tmp.sobol = pars.sobol[idx_start_batch:idx_end_batch,]
  # Run current batch
  vals = value_indicators(params_change = tmp.sobol, params_fixed = params_all)
  tictoc::toc()
  # Save current batch
  writeLines("Saving results")
  saveRDS(vals, 
          file = sprintf("%s/sensi_P%07d_DT%s_part-%03d-of-%03d.Rds",
                         OUTPUT_LOCAL,
                         N, 
                         curr_TD,
                         b, nb_batches))
}


writeLines("Computing partial rank correlations")
x_V_max = compute_PRCC(vals$V_max, pars.sobol[,2:dim(pars.sobol)[2]])
x_F_U_max = compute_PRCC(vals$F_U_max, pars.sobol[,2:dim(pars.sobol)[2]])
x_F_B_max = compute_PRCC(vals$F_B_max, pars.sobol[,2:dim(pars.sobol)[2]])
x_T_max_V = compute_PRCC(vals$T_max_V, pars.sobol[,2:dim(pars.sobol)[2]])
x_T_max_F_U = compute_PRCC(vals$T_max_F_U, pars.sobol[,2:dim(pars.sobol)[2]])
x_T_max_F_B = compute_PRCC(vals$T_max_F_B, pars.sobol[,2:dim(pars.sobol)[2]])
# Combine values for aggregate plots. All values, then by type
PRCC = list()
PRCC$all = data.frame(
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
             TeX("$T^*$"))
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
  theme_minimal() + 
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
  theme_minimal() +
  theme(legend.position = c(0.9,0.9))
ggsave(filename = "FIGS/sensitivities-PRCC-times.png",
       width = 20, height = 15, units = "cm")

colnames(PRCC$all) = c("names", "V_max", "F_U_max", "F_B_max", "T_max_V", "T_max_F_U", "T_max_F_B")
PRCC$values = PRCC$all[,c("names", "V_max", "F_U_max", "F_B_max")]
PRCC$times = PRCC$all[,c("names", "T_max_V", "T_max_F_U", "T_max_F_B")]
# Save the results
