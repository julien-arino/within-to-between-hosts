# ============================================================
# File: plot-Figure-D1-incidence-fct-time-various-scenarios-mean.R
# Purpose: 8 STRUCTURAL SCENARIOS (mean beta)
# 2 R0 VALUES PER PANEL
# SAME R0 WITHIN EACH PANEL (calibrated)
# ============================================================

library(dplyr)
library(ggplot2)
library(qs2)
library(tidyr)

# ------------------------------------------------------------
# USER CONTROLS
# ------------------------------------------------------------
R0_targets <- c(2.5, 5)
S0    <- 2000
U0    <- 1
d_P   <- 0
b_P   <- 0
Tmax  <- 200
Tplot <- 85

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
beta_df       <- qs_read("OUTPUT/beta_overall_transmitters.qs")
gamma_overall <- qs_read("OUTPUT/gamma_overall.qs")
mu_overall    <- qs_read("OUTPUT/mu_overall.qs")

xi_r_target <- "1e5"
xi_d_target <- "85"

gamma_col <- paste0("gamma_xi_", xi_r_target)
mu_col    <- paste0("mu_xid_",   xi_d_target)

# ------------------------------------------------------------
# ALIGN DATA
# ------------------------------------------------------------
n0 <- min(nrow(beta_df), nrow(gamma_overall), nrow(mu_overall))
beta_df       <- beta_df[1:n0, ]
gamma_overall <- gamma_overall[1:n0, ]
mu_overall    <- mu_overall[1:n0, ]

a_vals  <- beta_df$time
beta_a  <- beta_df$beta_mean
gamma_a <- gamma_overall[[gamma_col]]
mu_a    <- mu_overall[[mu_col]]

ok <- complete.cases(a_vals,beta_a,gamma_a,mu_a)

a_vals  <- a_vals[ok]
beta_a  <- beta_a[ok]
gamma_a <- gamma_a[ok]
mu_a    <- mu_a[ok]

ord <- order(a_vals)
a_vals  <- a_vals[ord]
beta_a  <- beta_a[ord]
gamma_a <- gamma_a[ord]
mu_a    <- mu_a[ord]

dt <- mean(diff(a_vals))
nA <- length(a_vals)

# ------------------------------------------------------------
# TIME GRID
# ------------------------------------------------------------
t_vals <- seq(0,Tmax,by=dt)
nT     <- length(t_vals)

# ------------------------------------------------------------
# CONSTANTS
# ------------------------------------------------------------
gamma_const <- mean(gamma_a)
mu_const    <- mean(mu_a)

beta_const_scalar <- mean(beta_a)
beta_const <- rep(beta_const_scalar,nA)

# ------------------------------------------------------------
# KERNELS
# ------------------------------------------------------------
kernel_gc_mc <- exp(-(d_P + gamma_const + mu_const)*a_vals)
kernel_ga_mc <- exp(-cumsum((d_P+gamma_a+mu_const)*dt))
kernel_gc_ma <- exp(-cumsum((d_P+gamma_const+mu_a)*dt))
kernel_ga_ma <- exp(-cumsum((d_P+gamma_a+mu_a)*dt))

# ------------------------------------------------------------
# R0 FUNCTION
# ------------------------------------------------------------
compute_R0 <- function(beta,kernel){
  S0*sum(beta*kernel)*dt
}

# ------------------------------------------------------------
# SIMULATOR
# ------------------------------------------------------------
simulate_case <- function(beta_raw,kernel,panel_name){
  
  res_list <- list()
  
  for(R0_target in R0_targets){
    
    R0_base <- compute_R0(beta_raw,kernel)
    k <- R0_target/R0_base
    beta_scaled <- k*beta_raw
    
    S <- numeric(nT)
    U <- numeric(nT)
    
    S[1] <- S0
    U[1] <- U0
    
    for(t in 2:nT){
      
      a_len <- min(nA,t-1)
      U_past <- U[(t-a_len):(t-1)]
      
      integral <- sum(beta_scaled[1:a_len]*
                        rev(U_past)*
                        kernel[1:a_len])*dt
      
      U[t] <- S[t-1]*integral
      S[t] <- S[t-1]+dt*(b_P-d_P*S[t-1]-U[t])
      S[t] <- max(S[t],0)
    }
    
    res_list[[as.character(R0_target)]] <-
      tibble(time=t_vals,
             U=U,
             panel=panel_name,
             R0=factor(paste0("R0=",R0_target)))
  }
  
  bind_rows(res_list)
}

# ------------------------------------------------------------
# RUN 8 PANELS
# ------------------------------------------------------------
res_all <- bind_rows(
  
  simulate_case(beta_const,kernel_gc_mc,"beta_c_gamma_c_mu_c"),
  simulate_case(beta_const,kernel_ga_mc,"beta_c_gamma_a_mu_c"),
  simulate_case(beta_const,kernel_gc_ma,"beta_c_gamma_c_mu_a"),
  simulate_case(beta_const,kernel_ga_ma,"beta_c_gamma_a_mu_a"),
  
  simulate_case(beta_a,kernel_gc_mc,"beta_a_gamma_c_mu_c"),
  simulate_case(beta_a,kernel_ga_mc,"beta_a_gamma_a_mu_c"),
  simulate_case(beta_a,kernel_gc_ma,"beta_a_gamma_c_mu_a"),
  simulate_case(beta_a,kernel_ga_ma,"beta_a_gamma_a_mu_a")
  
) %>% filter(time<=Tplot)

# ------------------------------------------------------------
# FACET LABELS
# ------------------------------------------------------------
panel_levels <- unique(res_all$panel)

panel_labels <- c(
  "beta[P]^c*','*gamma[P]^c*','*mu[P]^c",
  "beta[P]^c*','*gamma[P](a)*','*mu[P]^c",
  "beta[P]^c*','*gamma[P]^c*','*mu[P](a)",
  "beta[P]^c*','*gamma[P](a)*','*mu[P](a)",
  "beta[P](a)*','*gamma[P]^c*','*mu[P]^c",
  "beta[P](a)*','*gamma[P](a)*','*mu[P]^c",
  "beta[P](a)*','*gamma[P]^c*','*mu[P](a)",
  "beta[P](a)*','*gamma[P](a)*','*mu[P](a)"
)

res_all$panel <- factor(res_all$panel,
                        levels=panel_levels,
                        labels=panel_labels)

# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------
p <- ggplot(res_all,aes(time,U,color=R0))+
  geom_line(linewidth=1.2)+
  facet_wrap(~panel,nrow=2,labeller=label_parsed)+
  scale_color_manual(values=c("R0=2.5"="navy",
                              "R0=5"="darkred"),
                     labels=c(expression(R[0]==2.5),
                              expression(R[0]==5)))+
  labs(x="Time (days)",
       y=expression("Incidence"~U[P](t)),
       color=NULL)+
  theme_minimal(base_size=14)+
  theme(legend.position="right",
        strip.text=element_text(face="bold"))

print(p)

dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
ggsave("FIGS/Figure-D1-incidence-fct-time-various-scenarios-mean.pdf",
       p, width=22, height=12, units="cm")
ggsave("FIGS/Figure-D1-incidence-fct-time-various-scenarios-mean.png",
       p, width=22, height=12, units="cm")
cat("✅ Plot saved to FIGS/Figure-D1-incidence-fct-time-various-scenarios-mean.pdfn")
