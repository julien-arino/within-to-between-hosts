#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-D1-incidence-fct-time-various-scenarios-mean.R
# Description:
#   File: plot-Figure-D1-incidence-fct-time-various-scenarios-mean.R
#   8 STRUCTURAL SCENARIOS (mean beta)
#   2 R0 VALUES PER PANEL
#   SAME R0 WITHIN EACH PANEL (calibrated)
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-D1-incidence-fct-time-various-scenarios-mean.R")

load_libraries(c("dplyr", "ggplot2", "qs2", "tidyr", "here"))

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


get_latest_dist <- function(pattern) {
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found for pattern: ", pattern)
  files[which.max(file.mtime(files))]
}

xi_r_target <- "4"
xi_d_target <- "85"

pattern_str <- paste0("^cohort_distribution_filters_P.*_xid_", xi_d_target, "_.*_xir_", xi_r_target, "\\.qs$")
latest_dist <- get_latest_dist(pattern_str)
dist_list <- qs_read(latest_dist, nthreads = N_QS_THREADS)

if (!is.null(dist_list$rolling_average$beta_mean)) {
  dist_df <- dist_list$rolling_average
} else {
  dist_df <- dist_list$raw
}

# ------------------------------------------------------------
# ALIGN DATA
# ------------------------------------------------------------
combined_df <- dist_df %>%
  select(time, beta_mean, gamma_P, mu_P) %>%
  mutate(
    gamma_P = tidyr::replace_na(gamma_P, 0),
    mu_P = tidyr::replace_na(mu_P, 0)
  ) %>%
  arrange(time)

gamma_col <- "gamma_P"
mu_col    <- "mu_P"

a_vals  <- combined_df$time
beta_a  <- combined_df$beta_mean
gamma_a <- combined_df[[gamma_col]]
mu_a    <- combined_df[[mu_col]]

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

out_pdf <- file.path(figs_dir, "Figure-D1-incidence-fct-time-various-scenarios-mean.pdf")
out_png <- file.path(figs_dir, "Figure-D1-incidence-fct-time-various-scenarios-mean.png")

ggsave(out_pdf, p, width=22, height=12, units="cm")
ggsave(out_png, p, width=22, height=12, units="cm")

cat("Plot saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
print_end_time(start_time, "plot-Figure-D1-incidence-fct-time-various-scenarios-mean.R")
