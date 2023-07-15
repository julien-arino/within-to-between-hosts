library(deSolve)
library(dde)

rhs_within_host_deSolve = function(t, x, p) {
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

# rhs_within_host_dde = function(t, x, p) {
#   with(as.list(c(x, p)),{
#     # Variables are (in order) V,S,I,R,D,F_U,F_B
#     # So to get the lagvalue values...
#     V_t = ylag(tau_I, 1L)
#     S_t = ylag(tau_I, 2L)
#     I_t = ylag(tau_I, 3L)
#     R_t = ylag(tau_I, 4L)
#     # Now the RHS itself
#     dV = p*I-d_V*V
#     dS = lambda_S*(1-(S+I+D+R)/S_max)*S-beta*S*V
#     dI = beta*S_t*V_t*(1-F_B/(epsilon_FI+F_B))*A_I-d_I*I
#     dR = lambda_S*(1-(S+I+D+R)/S_max)*R +
#       beta*S_t*V_t*(F_B/(epsilon_FI+F_B))*A_R
#     dD = d_I*I-d_D*D
#     dF_U = psi_F_prod+p_FI*I/(I+eta_FI) - 
#       k_lin_f*F_U-k_B_F*((T_star+I)*A_F-F_B)*F_U+k_U_F*F_B
#     dF_B = -k_int_f*F_B+k_B_F*((T_star+I)*A_F-F_B)*F_U-k_U_F*F_B
#     dA_I = delta*A_I*(I_t-I)
#     dA_R = delta*A_R*(R_t-R)
#     return(list(c(dV, dS, dI, dR, dD, dF_U, dF_B,dA_I,dA_R)))
#   })
# }

set_IC_1 = function() {
  # | Variable | Definition    | Value  | Unit |
  # |----------|---------------|--------|------|
  # | V   | Viral load         | 4.5    | log_10(copies/ml) |
  # | S   | Susceptible cell   | 0.16   | 10^9 cells/ml |
  # | I   | Infected cell      | 0      | 10^9 cells/ml |
  # | R   | Resistant cell     | 0      | 10^9 cells/ml |
  # | D   | Dead cell          | 0      | 10^9 cells/ml |
  # | F_U | Unbound interferon | 0.015  | pg/ml |
  # | F_B | Bound interferon   | 1.1e-8 | pg/ml |
  IC = c(V = 4.5, S = 0.16, I = 0, R = 0, D = 0, 
         F_U = 0.015, F_B = 1.1e-8, A_I = 1, A_R = 1)
  return(IC)
}

set_IC_2 = function() {
  IC = c(V = 1, S = 0.16, I = 0, R = 0, D = 0, 
         F_U = 0.015, F_B = 1.1e-8, A_I = 1, A_R = 1)
  return(IC)
}

set_parameters = function() {
  # | Parameter  | Definition | Value | Unit |
  # | lambda_S   |  Proliferation of susceptible cells   | 0.74  |  $day^{-1}$  |
  # | S_max      |  Target cell concentration   | 0.16  |  $10^9$cells/ml  |
  # | d_I        |  Death rate of infected cells   | 0.1  |  $day^{-1}$  |
  # | d_D        |  Degradation rate of dead cells   | 8  |  $day^{-1}$  |
  # | tau_I      |  Eclipse time   | 0.17  |  $day$ |
  # | beta       |  Viral infection rate   | 0.3 (SD:0.1994)  |  $day^{-1}cop/ml$  |
  # | d_V        | Viral decay rate    |  8.4 (SD:0.67)  |  $day^{-1}$ |
  # | p          |  Viral production rate  |   394 (SD:158.65) |  $day^{-1}(cop/10^9cells)$
  # | K_U_F      |  IFN unbinding rate   |  6.072  | $day^{-1}$  |
  # | p_FI       |  IFN production by infected cells   |  2.8235 (SD:1.8741) |  $day^{-1}(pg/ml)$  |
  # | psi_F_prod | IFN production by macrophages and monocytes  |  0.25 |  $day^{-1} (pg/ml)$  |
  # | k_B_F      | IFN binding rate  |  0.0107 |  $day^{-1} (ml/pg)$  |
  # | k_lin_f    |  IFN renal clearance rate  | 16.635 (SD:2.49)  |  $day^{-1}$  |
  # | k_int_f    | IFN internalization rate  |  16.968 (SD:2.54) |  $day^{-1}$  |
  # | epsilon_FI |  Half maximal response   | 2E-4 |  $10^9 cell/ml$  |
  # | eta_FI     | Half-maximal response  | 0.022  |  $10^9$cells/ml  |
  # | T_star     | Initial CD8+ T cells |  1.104E-4 |  $10^9$cells/ml |
  params = c(
    lambda_S = 0.74,
    S_max = 0.16,
    d_I = 0.1,
    d_D = 8,
    tau_I = 0.17,
    beta = 0.3,
    beta_stddev = 0.1994,
    d_V = 8.4,
    d_V_stddev = 0.67,
    p = 394,
    p_stddev = 158.65,
    k_U_F = 6.072,
    p_FI = 2.8235,
    p_FI_stddev = 1.8741,
    psi_F_prod = 0.25,
    k_B_F = 0.0107,
    k_lin_f = 16.635,
    k_lin_f_stddev = 2.49,
    k_int_f = 16.968,
    k_int_f_stddev = 2.54,
    epsilon_FI = 2e-4,
    eta_FI = 0.022328,
    T_star = 1.104e-4,
    delta = 0.1,
    avo=6.02214e23,
    MM_F = 19000,
    R_F_T = 1000,
    R_F_I = 1300
  )
  params = c(params,
             A_F = as.numeric(
               (params["MM_F"]/params["avo"]) *
                 (params["R_F_I"]+params["R_F_T"]) *
                 (1/5000)*(10^9*1e12)
             ))
  return(params)
}

add_IC_to_params = function(params, IC) {
  params = c(params, 
             V0 = as.numeric(IC["V"]), 
             S0 = as.numeric(IC["S"]),
             I0 = as.numeric(IC["I"]), 
             R0 = as.numeric(IC["R"]))
return(params)
}

# Given parameters, find those with a standard deviation
# given and generate a table with sampled values of these parameters, 
# regular values of the others, as well as initial conditions. 
# This way, we have all that's needed for the cohort.
generate_params_patients = function(params, n = 1000) {
  # The whole list of parameters, including values of std dev for some
  names_params = names(params)
  # Which are the parameters that contain std dev information
  idx_stddev = grep("stddev", names_params)
  # Which parameters is that std dev info for
  params_with_stddev = gsub("_stddev", 
                            "", 
                            names_params[idx_stddev])
  idx_params_with_stddev = which(names_params %in% params_with_stddev)
  OUT = 
    data.frame(
      mat.or.vec(nr = n,
                 nc = length(names_params)-length(idx_params_stddev))
    )
  colnames(OUT) = names_params[setdiff(1:length(names_params),
                                       idx_stddev)]
  for (curr_col in colnames(OUT)) {
    if (!(curr_col %in% params_with_stddev)) {
      # This is a "regular" parameter, we just replicate it n times
      OUT[[curr_col]] = rep(params[curr_col], n)
    } else {
      # This a parameter we must sample. As indicated, we sample from a 
      # normal distribution with mean the value given and std dev 3 times
      # the given std dev
      OUT[[curr_col]] = rnorm(n = n,
                              mean = params[curr_col],
                              sd = 3*params[sprintf("%s_stddev",
                                                    curr_col)])
      # This can give us negative values. If so, just take the mean 
      # (for now)
      OUT[[curr_col]][which(OUT[[curr_col]]<0)] = params[curr_col]
    }
  }
  return(OUT)
}

# Given a patient index idx, a parameters data frame params.df and 
# initial conditions IC, run the simulation of the within host model for
# this patient
run_one_patient = function(idx = 1, 
                           params.df, 
                           IC) {
  writeLines(paste0("patient index = ", idx))
  params_tmp = params.df[idx,]
  params_tmp = add_IC_to_params(params_tmp, IC)
  times <- seq(0, 200, by = 0.1)
  yout <- dede(y = IC, 
               times = times, 
               func = rhs_within_host_deSolve, 
               parms = params_tmp)
  return(yout)
}