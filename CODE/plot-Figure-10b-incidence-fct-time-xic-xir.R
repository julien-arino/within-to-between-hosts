#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-10b-incidence-fct-time-xic-xir.R
# Description:
#   BETWEEN-HOST SIMULATOR – 1 FIGURE
#   Compare R0=2.5 using beta_mean across varying (xi_c, xi_r) pairs
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-10b-incidence-fct-time-xic-xir.R")

load_libraries(c("ggplot2", "dplyr", "qs2"))

# USER CONTROLS
R0_t <- 2.5
S0 <- 2000
U0 <- 1
d_P <- 1 / (60 * 365) # 0.00
b_P <- 2000 * d_P # 0.00
Tmax <- 110

xi_pairs <- list(
  c(4, 1),
  c(4, 4)
)

# Simulator (returns U only)
simulate_case_U <- function(beta_age, kernel_age, dt, nT, nA) {
  S <- numeric(nT)
  U <- numeric(nT)
  S[1] <- S0
  U[1] <- U0
  
  for (t in 2:nT) {
    a_len <- min(nA, t - 1)
    U_past <- U[(t - a_len):(t - 1)]
    
    integral <- sum(beta_age[1:a_len] * rev(U_past) * kernel_age[1:a_len]) * dt
    U[t] <- S[t - 1] * integral
    
    S[t] <- S[t - 1] + dt * (b_P - d_P * S[t - 1] - U[t])
    S[t] <- max(S[t], 0)
  }
  
  U
}

res_list <- list()

for (pair in xi_pairs) {
  xic <- pair[1]
  xir <- pair[2]
  
  # Base pattern string targeting the combinations
  pattern_str <- paste0("^cohort_distribution_filters_P.*_xic_", xic, "_xir_", xir, "\\.qs$")
  
  # find latest dist
  files <- list.files(output_dir, pattern = pattern_str, full.names = TRUE)
  if (length(files) == 0) {
    cat(sprintf("No files found for xic=%g xir=%g. Skipping...\n", xic, xir))
    next
  }
  latest_dist <- files[which.max(file.mtime(files))]
  cat("Loading:", basename(latest_dist), "\n")
  dist_list <- qs_read(latest_dist, nthreads = N_QS_THREADS)
  
  dist_df <- dist_list$rolling_average
  if (is.null(dist_df$beta_mean)) dist_df <- dist_list$raw
  
  combined_df <- dist_df %>%
    select(time, beta_mean, gamma_P, mu_P) %>%
    filter(time <= 50) %>%
    mutate(
      beta_mean = coalesce(beta_mean, 0),
      gamma_P   = coalesce(gamma_P, 0),
      mu_P      = coalesce(mu_P, 0)
    )
  
  a_vals <- combined_df$time
  beta_base <- combined_df$beta_mean
  gamma_a <- combined_df$gamma_P
  mu_a <- combined_df$mu_P
  
  gamma_a[is.na(gamma_a)] <- 0
  mu_a[is.na(mu_a)] <- 0
  
  ord <- order(a_vals)
  a_vals <- a_vals[ord]
  beta_base <- beta_base[ord]
  gamma_a <- gamma_a[ord]
  mu_a <- mu_a[ord]
  
  dt <- mean(diff(a_vals))
  if (!is.finite(dt) || dt <= 0) stop("dt is not valid. Check a_vals grid.")
  
  t_vals <- seq(0, Tmax, by = dt)
  nT <- length(t_vals)
  nA <- length(a_vals)
  
  haz_var <- d_P + gamma_a + mu_a
  cumhaz_var <- cumsum(haz_var * dt)
  kernel_var <- exp(-cumhaz_var)
  
  R0_base_var <- S0 * sum(beta_base * kernel_var, na.rm = TRUE) * dt
  
  # scale beta(a) to hit R0_t exactly
  k_var <- R0_t / R0_base_var
  beta_var <- k_var * beta_base
  
  # simulate
  U_var <- simulate_case_U(beta_var, kernel_var, dt, nT, nA)
  
  AUC_val <- sum(U_var) * dt
  cat(sprintf("[R0=2.5] AUC for xic=%g, xir=%g: %f\n", xic, xir, AUC_val))
  
  pair_label <- paste0("xic", xic, "_xir", xir)
  
  out_var <- tibble(
    time = t_vals,
    U = U_var,
    xic = xic,
    xir = xir,
    pair_label = pair_label
  )
  
  res_list[[length(res_list) + 1]] <- out_var
}

if (length(res_list) == 0) {
  stop("No resulting simulations could be run. Check pattern generation or files.")
}

res_all <- bind_rows(res_list) %>% arrange(time)

# Plot
y_max <- max(res_all$U, na.rm = TRUE)

p <- ggplot(res_all, aes(x = time, y = U, color = pair_label)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(
    values = c(
      "xic4_xir1" = "darkblue",
      "xic4_xir4" = "darkred"
    ),
    labels = c(
      "xic4_xir1" = expression(xi^c == 4 ~ "," ~ xi^r == 1),
      "xic4_xir4" = expression(xi^c == 4 ~ "," ~ xi^r == 4)
    )
  ) +
  labs(
    x = "Time (days)", 
    y = expression("Incidence"~U[P](t)), 
    color = "Thresholds"
  ) +
  coord_cartesian(xlim = c(0, Tmax), ylim = c(0, y_max)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

print(p)

out_pdf <- file.path(figs_dir, "Figure-10b-incidence-fct-time-xic-xir.pdf")
out_png <- file.path(figs_dir, "Figure-10b-incidence-fct-time-xic-xir.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm")
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm")

cat("Plot saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
print_end_time(start_time, "plot-Figure-10b-incidence-fct-time-xic-xir.R")
