# ============================================================
# File: plot-Figure-D4-incidence-fct-time-same-R0P-different-summaries.R
# Purpose: BETWEEN-HOST SIMULATOR – Compare 4 transmission profiles
# beta_mean, beta_median, beta_Q10, beta_Q90
# All rescaled to the same R0^P = 2.5
# ============================================================

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))

# ------------------------------------------------------------
# USER PARAMETERS
# ------------------------------------------------------------
S0   <- 2000
U0   <- 1
d_P  <- 0
b_P  <- 0
Tmax <- 110

R0_target <- 2.5

# ------------------------------------------------------------
# Load data dynamically
# ------------------------------------------------------------
beta_df       <- qs_read("OUTPUT/beta_overall_transmitters.qs")
gamma_overall <- qs_read("OUTPUT/gamma_overall.qs")
mu_overall    <- qs_read("OUTPUT/mu_overall.qs")

gamma_col <- "gamma_xi_1e5"
mu_col    <- "mu_xid_85"

# ------------------------------------------------------------
# Align datasets
# ------------------------------------------------------------
n0 <- min(nrow(beta_df), nrow(gamma_overall), nrow(mu_overall))

beta_df       <- beta_df[1:n0,]
gamma_overall <- gamma_overall[1:n0,]
mu_overall    <- mu_overall[1:n0,]

a_vals <- beta_df$time

beta_mean   <- beta_df$beta_mean
beta_q10    <- beta_df$beta_q10
beta_q90    <- beta_df$beta_q90
beta_median <- beta_df$beta_median

gamma_a <- gamma_overall[[gamma_col]]
mu_a    <- mu_overall[[mu_col]]

gamma_a[is.na(gamma_a)] <- 0
mu_a[is.na(mu_a)] <- 0

dt <- mean(diff(a_vals))

t_vals <- seq(0, Tmax, by = dt)
nT <- length(t_vals)
nA <- length(a_vals)

# ------------------------------------------------------------
# Kernel (age-dependent hazards)
# ------------------------------------------------------------
haz_var <- d_P + gamma_a + mu_a
cumhaz_var <- cumsum(haz_var * dt)
kernel_var <- exp(-cumhaz_var)

# ------------------------------------------------------------
# Simulator
# ------------------------------------------------------------
simulate_case <- function(beta_age){
  
  S <- numeric(nT)
  U <- numeric(nT)
  
  S[1] <- S0
  U[1] <- U0
  
  for(t in 2:nT){
    
    a_len <- min(nA, t-1)
    
    U_past <- U[(t-a_len):(t-1)]
    
    integral <- sum(beta_age[1:a_len] *
                      rev(U_past) *
                      kernel_var[1:a_len]) * dt
    
    U[t] <- S[t-1] * integral
    
    S[t] <- S[t-1] + dt*(b_P - d_P*S[t-1] - U[t])
    
    S[t] <- max(S[t],0)
  }
  
  tibble(time=t_vals, U=U)
}

# ------------------------------------------------------------
# R0 computation
# ------------------------------------------------------------
compute_R0 <- function(beta_age){
  S0 * sum(beta_age * kernel_var) * dt
}

# ------------------------------------------------------------
# Base R0
# ------------------------------------------------------------
R0_mean   <- compute_R0(beta_mean)
R0_q10    <- compute_R0(beta_q10)
R0_q90    <- compute_R0(beta_q90)
R0_median <- compute_R0(beta_median)

# ------------------------------------------------------------
# Scaling factors
# ------------------------------------------------------------
k_mean   <- R0_target / R0_mean
k_q10    <- R0_target / R0_q10
k_q90    <- R0_target / R0_q90
k_median <- R0_target / R0_median

beta_mean_scaled   <- k_mean   * beta_mean
beta_q10_scaled    <- k_q10    * beta_q10
beta_q90_scaled    <- k_q90    * beta_q90
beta_median_scaled <- k_median * beta_median

# ------------------------------------------------------------
# Simulations
# ------------------------------------------------------------
res_mean   <- simulate_case(beta_mean_scaled)   %>% mutate(profile="mean")
res_q10    <- simulate_case(beta_q10_scaled)    %>% mutate(profile="Q10")
res_q90    <- simulate_case(beta_q90_scaled)    %>% mutate(profile="Q90")
res_median <- simulate_case(beta_median_scaled) %>% mutate(profile="median")

res_all <- bind_rows(res_mean,res_q10,res_median,res_q90) %>%
  mutate(
    profile = factor(
      profile,
      levels = c("mean","Q10","median","Q90")
    )
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(res_all, aes(time,U,color=profile,linetype=profile))+
  
  geom_line(linewidth=1)+
  
  scale_color_manual(
    values=c(
      "mean"="dodgerblue4",
      "Q10"="black",
      "median"="magenta4",
      "Q90"="orange"
    ),
    breaks=c("mean","Q10","median","Q90"),
    labels=c(
      expression(beta[mean]),
      expression(beta[Q10]),
      expression(beta[median]),
      expression(beta[Q90])
    )
  )+
  
  scale_linetype_manual(
    values=c(
      "mean"="solid",
      "Q10"="21",    # Actually maps to dashed lines in ggplot2 due to a custom override in guides below
      "median"="solid",
      "Q90"="21"
    )
  )+
  
  guides(
    linetype="none",
    color=guide_legend(override.aes=list(
      linetype=c("solid","21","solid","21"),
      linewidth=1
    ))
  )+
  
  labs(
    x="Time (days)",
    y=expression("Incidence"~U[P](t)),
    color=bquote(R[0]^P==.(R0_target))
  )+
  
  coord_cartesian(xlim=c(0,50))+
  
  theme_minimal(base_size=14)

print(p)

# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/Figure-D4-incidence-fct-time-same-R0P-different-summaries.pdf"
out_png <- "FIGS/Figure-D4-incidence-fct-time-same-R0P-different-summaries.png"

ggsave(
  filename = out_pdf,
  plot = p,
  width = 20,
  height = 8,
  units = "cm"
)

ggsave(
  filename = out_png,
  plot = p,
  width = 20,
  height = 8,
  units = "cm"
)

cat("✅ Plot saved to", out_pdf, "and", out_png, "\n")
