library(deSolve)

rhs_within_host = function(t, x, p) {
  with(as.list(c(x, p)),{
    # Variables are (in order) V,S,I,R,D,F_U,F_B
    # So to get the lagvalue values...
    if (t<tau_I) {
      V_t = V0
      S_t = S0
    } else {
      V_t = lagvalue(t-tau_I,1)
      S_t = lagvalue(t-tau_I,2)
    }
    dV = p*I-d_V*V
    dS = lambda_S*(1-(S+I+D+R)/S_max)*S-beta*S*V
    dI = beta*S_t*V_t*(1-F_B/(epsilon_FI+F_B))-d_I*I
    dR = lambda_S*(1-(S+I+D+R)/S_max)*R +
      beta*S_t*V_t*(F_B/(epsilon_FI+F_B))
    dD = d_I*I-d_D*D
    dF_U = psi_F_prod+p_FI*I/(I+eta_FI) - 
      k_lin_f*F_U-k_B_F*((T_star+I)*A_F-F_B)*F_U+k_U_F*F_B
    dF_B = -k_int_f*F_B+k_B_F*((T_star+I)*A_F-F_B)*F_U-k_U_F*F_B
    return(list(c(dV, dS, dI, dR, dD, dF_U, dF_B)))
  })
}

set_IC = function() {
  # | Variable | Definition    | Value  | Unit |
  # |----------|---------------|--------|------|
  # | S   | Susceptible cell   | 0.16   | 10^9 cells/ml |
  # | I   | Infected cell      | 0      | 10^9 cells/ml |
  # | R   | Resistant cell     | 0      | 10^9 cells/ml |
  # | D   | Dead cell          | 0      | 10^9 cells/ml |
  # | V   | Viral load         | 4.5    | log_10(copies/ml) |
  # | F_U | Unbound interferon | 0.015  | pg/ml |
  # | F_B | Bound interferon   | 1.1e-8 | pg/ml |
  IC = c(S = 0.16, I = 0, R = 0, D = 0, V = 4.5, 
         F_U = 0.015, F_B = 1.1e-8)
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
    d_V = 8.4,
    p = 394,
    k_U_F = 6.072,
    p_FI = 2.8235,
    psi_F_prod = 0.25,
    k_B_F = 0.0107,
    k_lin_f = 16.635,
    k_int_f = 16.968,
    epsilon_FI = 2e-4,
    eta_FI = 0.022,
    T_star = 1.104e-4,
    tau_I = 0.1,
    A_F = 0.01
  )
  return(params)
}

add_IC_to_params = function(params, IC) {
  params = c(params, 
             V0 = as.numeric(IC["V"]), 
             S0 = as.numeric(IC["S"]))
  return(params)
}
