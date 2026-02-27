# ============================================================
# File : model_ode_functions.R
# Within-host ODE model (no infection delay)
# ============================================================

library(deSolve)

# ---- Default parameter set ----
default_params <- list(
  beta_V     = 0.3,
  d_I        = 0.1,
  d_D        = 8,
  d_V        = 8.4,
  epsilon_FI = 2e-4,
  eta_FI     = 0.022328,
  k_B_F      = 0.0107,
  k_int_f    = 16.968,
  k_lin_f    = 16.635,
  k_U_F      = 6.072,
  lambda_S   = 0.74,
  p          = 394,
  p_FI       = 2.8235,
  psi_F_prod = 0.25,
  S_max      = 0.16,
  T_star     = 1.104e-4
)


# ---- Initial conditions + time grid ----
yini  <- c(V = 1, S = 0.16, I = 0, R = 0, D = 0,
           F_U = 0.015, F_B = 1.1e-8)
temps <- seq(0, 80, by = 0.1)

# ---- ODE system ----
within_host_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    avo  <- 6.02214e23
    MM_F <- 19000
    R_F_T <- 1000
    R_F_I <- 1300
    A_F <- as.numeric((MM_F / avo) * (R_F_I + R_F_T) *
                        (1 / 5000) * (10^9 * 1e12))
    
    N <- S + I + R + D
    
    dV  <- p * I - d_V * V
    dS  <- lambda_S * (1 - N / S_max) * S - beta_V * S * V
    dI  <- beta_V * S * V * (1 - F_B / (epsilon_FI + F_B)) - d_I * I
    dR  <- lambda_S * (1 - N / S_max) * R +
      beta_V * S * V * (F_B / (epsilon_FI + F_B))
    dD  <- d_I * I - d_D * D
    dF_U <- psi_F_prod + (p_FI * I) / (I + eta_FI) -
      k_lin_f * F_U -
      k_B_F * ((T_star + I) * A_F - F_B) * F_U +
      k_U_F * F_B
    dF_B <- -k_int_f * F_B +
      k_B_F * ((T_star + I) * A_F - F_B) * F_U -
      k_U_F * F_B
    
    list(c(dV, dS, dI, dR, dD, dF_U, dF_B))
  })
}

# ---- Run ODE simulation ----
simulate_model_ode <- function(params) {
  ode(y = yini, times = temps, func = within_host_ode,
      parms = as.list(params), atol = 1e-8, rtol = 1e-8)
}

# ---- Extract maxima and time-to-max indicators ----
extract_max_indicators_ode <- function(params) {
  out <- tryCatch(simulate_model_ode(params), error = function(e) NULL)
  if (is.null(out)) {
    return(data.frame(
      V_max = NA_real_, F_U_max = NA_real_, F_B_max = NA_real_,
      t_V_max = NA_real_, t_F_U_max = NA_real_, t_F_B_max = NA_real_
    ))
  }
  
  out <- as.data.frame(out)
  
  V_max   <- max(out$V,   na.rm = TRUE)
  F_U_max <- max(out$F_U, na.rm = TRUE)
  F_B_max <- max(out$F_B, na.rm = TRUE)
  
  t_V_max   <- out$time[which.max(out$V)]
  t_F_U_max <- out$time[which.max(out$F_U)]
  t_F_B_max <- out$time[which.max(out$F_B)]
  
  data.frame(
    V_max   = as.numeric(V_max),
    F_U_max = as.numeric(F_U_max),
    F_B_max = as.numeric(F_B_max),
    t_V_max   = as.numeric(t_V_max),
    t_F_U_max = as.numeric(t_F_U_max),
    t_F_B_max = as.numeric(t_F_B_max)
  )
}

# ---- Quick test ----
if (sys.nframe() == 0) {
  cat("Running quick test for the within-host ODE model...\n")
  result <- extract_max_indicators_ode(default_params)
  print(result)
}


# ---- Compute disease severity (Ψ_i) and classification ----
compute_severity_ode <- function(out, S_max) {
  # Ensure correct format
  out <- as.data.frame(out)
  
  # 1️⃣ Compute % of tissue damage
  Psi <- 100 * (S_max - (out$S + out$R)) / S_max
  
  # 2️⃣ Find maximum damage
  Psi_max <- max(Psi, na.rm = TRUE)
  
  # 3️⃣ Classify patient outcome
  if (Psi_max >= 85) {
    status <- "dead"
  } else if (Psi_max >= 75) {
    status <- "ICU"
  } else {
    status <- "rest"
  }
  
  # 4️⃣ Record times when thresholds are crossed (optional)
  t_ICU   <- ifelse(any(Psi >= 75), out$time[which(Psi >= 75)[1]], NA)
  t_death <- ifelse(any(Psi >= 85), out$time[which(Psi >= 85)[1]], NA)
  
  # 5️⃣ Return all useful information
  return(list(
    Psi = Psi,           # trajectory over time
    Psi_max = Psi_max,   # max damage %
    status = status,     # ICU / dead / rest
    t_ICU = t_ICU,       # time to ICU threshold
    t_death = t_death    # time to death threshold
  ))
}

if (sys.nframe() == 0) {
  cat("Running quick test for severity classification...\n")
  out <- as.data.frame(simulate_model_ode(default_params))
  sev <- compute_severity_ode(out, default_params$S_max)
  print(paste("Max tissue damage:", round(sev$Psi_max, 2), "%"))
  print(paste("Status:", sev$status))
}
