library(ODEsensitivity)

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
# The timepoints of interest:
times <- c(0.01, seq(1, 50, by = 1))
# Morris screening:
set.seed(7292)
# Warning: The following code might take very long!
LVres_morris <- DDEmorris(mod = rhs_within_host_deSolve,
                          pars = pars$params,
                          state_init = IC,
                          times = times,
                          binf = pars$min,
                          bsup = pars$max,
                          r = 500,
                          design = list(type = "oat",
                                        levels = 10, grid.jump = 1),
                          scale = TRUE,
                          ode_method = "lsoda",
                          parallel_eval = TRUE,
                          parallel_eval_ncores = 6)

set.seed(59281)
# Sobol' sensitivity analysis (here only with n = 500, but n = 1000 is
# recommended):
# Warning: The following code might take very long!
LVres_sobol <- ODEsobol(mod = LVmod,
                        pars = LVpars,
                        state_init = LVinit,
                        times = LVtimes,
                        n = 500,
                        rfuncs = "runif",
                        rargs = paste0("min = ", LVbinf,
                                       ", max = ", LVbsup),
                        sobol_method = "Martinez",
                        ode_method = "lsoda",
                        parallel_eval = TRUE,
                        parallel_eval_ncores = 6)


##### FitzHugh-Nagumo equations (Ramsay et al., 2007) #####
FHNmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
    dCurrent <- - 1 / s *(Voltage - a + b * Current)
    return(list(c(dVoltage, dCurrent)))
  })
}
# Warning: The following code might take very long!
FHNres_morris <- ODEmorris(mod = FHNmod,
                           pars = c("a", "b", "s"),
                           state_init = c(Voltage = -1, Current = 1),
                           times = seq(0.1, 50, by = 5),
                           binf = c(0.18, 0.18, 2.8),
                           bsup = c(0.22, 0.22, 3.2),
                           r = 500,
                           design = list(type = "oat",
                                         levels = 50, grid.jump = 1),
                           scale = TRUE,
                           ode_method = "adams",
                           parallel_eval = TRUE,
                           parallel_eval_ncores = 6)
  
# Warning: The following code might take very long!

FHNres_sobol <- ODEsobol(mod = FHNmod,
                         pars = c("a", "b", "s"),
                         state_init = c(Voltage = -1, Current = 1),
                         times = seq(0.1, 50, by = 5),
                         n = 500,
                         rfuncs = "runif",
                         rargs = c(rep("min = 0.18, max = 0.22", 2),
                                   "min = 2.8, max = 3.2"),
                         sobol_method = "Martinez",
                         ode_method = "adams",
                         parallel_eval = TRUE,
                         parallel_eval_ncores = 6)

# Just for demonstration purposes: The use of different distributions for the 
# parameters (here, the distributions and their arguments are chosen 
# completely arbitrarily):
# Warning: The following code might take very long!

demo_dists <- ODEsobol(mod = FHNmod,
                       pars = c("a", "b", "s"),
                       state_init = c(Voltage = -1, Current = 1),
                       times = seq(0.1, 50, by = 5),
                       n = 500,
                       rfuncs = c("runif", "rnorm", "rexp"),
                       rargs = c("min = 0.18, max = 0.22",
                                 "mean = 0.2, sd = 0.2 / 3",
                                 "rate = 1 / 3"),
                       sobol_method = "Martinez",
                       ode_method = "adams",
                       parallel_eval = TRUE,
                       parallel_eval_ncores = 6)
