#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-D3-incidence-fct-time-various-scenarios-median-AUC.R
# Description:
#   File: plot-Figure-D3-incidence-fct-time-various-scenarios-median-AUC.R
#   8 STRUCTURAL SCENARIOS (median beta) WITH AUC
#   2 R0 VALUES PER PANEL
#   SAME R0 WITHIN EACH PANEL (calibrated)
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-D3-incidence-fct-time-various-scenarios-median-AUC.R")

load_libraries(c("dplyr", "ggplot2", "qs2", "tidyr", "here"))

R0_targets <- c(2.5, 5)
S0    <- 2000
U0    <- 1
d_P   <- 0
b_P   <- 0
Tmax  <- 300
Tplot <- 300



get_latest_dist <- function(pattern) {
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found for pattern: ", pattern)
  files[which.max(file.mtime(files))]
}

beta_df       <- qs_read(get_latest_dist("^cohort_distributions_P.*_beta\\.qs$"), nthreads = N_QS_THREADS)
gamma_overall <- qs_read(get_latest_dist("^cohort_distributions_P.*_gamma\\.qs$"), nthreads = N_QS_THREADS)
mu_overall    <- qs_read(get_latest_dist("^cohort_distributions_P.*_mu\\.qs$"), nthreads = N_QS_THREADS)

xi_r_target <- "4"
xi_d_target <- "85"

gamma_col <- paste0("gamma_xi_", xi_r_target)
mu_col    <- paste0("mu_xid_",   xi_d_target)

combined_df <- beta_df %>%
  select(time, beta_median) %>%
  left_join(gamma_overall %>% select(time, all_of(gamma_col)), by = "time") %>%
  left_join(mu_overall %>% select(time, all_of(mu_col)), by = "time") %>%
  mutate(
    across(c(all_of(gamma_col), all_of(mu_col)), ~ replace_na(.x, 0))
  ) %>%
  arrange(time)

a_vals  <- combined_df$time
beta_a  <- combined_df$beta_median
gamma_a <- combined_df[[gamma_col]]
mu_a    <- combined_df[[mu_col]]

dt <- median(diff(a_vals))
nA <- length(a_vals)

t_vals <- seq(0,Tmax,by=dt)
nT     <- length(t_vals)

gamma_const <- median(gamma_a)
mu_const    <- median(mu_a)

beta_const_scalar <- median(beta_a)
beta_const <- rep(beta_const_scalar,nA)

kernel_gc_mc <- exp(-(d_P + gamma_const + mu_const)*a_vals)
kernel_ga_mc <- exp(-cumsum((d_P+gamma_a+mu_const)*dt))
kernel_gc_ma <- exp(-cumsum((d_P+gamma_const+mu_a)*dt))
kernel_ga_ma <- exp(-cumsum((d_P+gamma_a+mu_a)*dt))

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

summary_stats <- res_all %>%
  group_by(panel, R0) %>%
  summarise(
    AUC  = sum(U)*dt,
    .groups="drop"
  ) %>%
  group_by(panel) %>%
  mutate(
    label = paste0("AUC = ", round(AUC,0)),
    y_pos = max(res_all$U) * 0.85 -
      (as.numeric(factor(R0)) - 1) * max(res_all$U) * 0.08
  )

p <- ggplot(res_all,aes(time,U,color=R0))+
  geom_line(linewidth=1.2)+
  facet_wrap(~panel,nrow=2,labeller=label_parsed)+
  geom_text(data=summary_stats,
            aes(x=Tplot*0.6,
                y=y_pos,
                label=label,
                color=R0),
            inherit.aes=FALSE,
            size=2,
            show.legend=FALSE)+
  scale_color_manual(values=c("R0=2.5"="navy",
                              "R0=5"="darkred"),
                     labels=c(expression(R[0]^P==2.5),
                              expression(R[0]^P==5)))+
  labs(x="Time (days)",
       y=expression("Incidence"~U[P](t)),
       color=NULL)+
  theme_minimal(base_size=14)+
  theme(legend.position="right",
        strip.text=element_text(face="bold"))

print(p)

out_pdf <- file.path(figs_dir, "Figure-D3-incidence-fct-time-various-scenarios-median-AUC.pdf")
out_png <- file.path(figs_dir, "Figure-D3-incidence-fct-time-various-scenarios-median-AUC.png")

ggsave(out_pdf, p, width=22, height=12, units="cm")
ggsave(out_png, p, width=22, height=12, units="cm")

cat("Plot saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
print_end_time(start_time, "plot-Figure-D3-incidence-fct-time-various-scenarios-median-AUC.R")
