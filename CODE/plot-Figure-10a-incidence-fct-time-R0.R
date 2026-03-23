#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-10a-incidence-fct-time-R0.R
# Description:
#   BETWEEN-HOST SIMULATOR – 1 FIGURE
#   Compare R0=2.5 vs R0=5 using ONLY beta_mean
#   - Solid: age-dependent (beta(a), gamma(a), mu(a))
#   - Dashed: constants (beta_c, gamma_c, mu_c)
# ============================================================

project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

start_time <- start_time_and_hello("plot-Figure-10a-incidence-fct-time-R0.R")

load_libraries(c("ggplot2", "dplyr", "qs2"))

# Parameters and simulation control
R0_targets <- c(2.5, 5.0)
S0 <- 2000
U0 <- 1
d_P <- 1 / (60 * 365) # 0.00
b_P <- 2000 * d_P # 0.00
Tmax <- 110

# Load data
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

# Kernels: constant vs age-dependent
gamma_const <- mean(gamma_a, na.rm = TRUE)
mu_const <- mean(mu_a, na.rm = TRUE)

kernel_const <- exp(-(d_P + gamma_const + mu_const) * a_vals)

haz_var <- d_P + gamma_a + mu_a
cumhaz_var <- cumsum(haz_var * dt)
kernel_var <- exp(-cumhaz_var)

# Simulator (returns U only)
simulate_case_U <- function(beta_age, kernel_age, label = "") {
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

  tibble(time = t_vals, U = U)
}

compute_R0P <- function(beta_age, kernel_age) {
  S0 * sum(beta_age * kernel_age, na.rm = TRUE) * dt
}

# Build curves for each R0 target
denom_const <- S0 * sum(kernel_const, na.rm = TRUE) * dt
R0_base_var <- compute_R0P(beta_base, kernel_var)

res_list <- list()

for (R0_target in R0_targets) {
  # scale beta(a) to hit R0_target exactly
  k_var <- R0_target / R0_base_var
  beta_var <- k_var * beta_base

  # pick constant beta_c to hit same R0_target exactly under const hazards
  beta_c <- R0_target / denom_const
  beta_const <- rep(beta_c, nA)

  # simulate
  out_var <- simulate_case_U(beta_var, kernel_var, paste0("var_R0_", R0_target)) %>%
    mutate(R0 = R0_target, structure = "var_age")

  out_const <- simulate_case_U(beta_const, kernel_const, paste0("const_R0_", R0_target)) %>%
    mutate(R0 = R0_target, structure = "const")

  res_list[[paste0("var_", R0_target)]] <- out_var
  res_list[[paste0("const_", R0_target)]] <- out_const
}

res_all <- bind_rows(res_list) %>%
  mutate(
    R0 = factor(R0, levels = R0_targets, labels = c("R0=2.5", "R0=5")),
    structure = factor(structure,
      levels = c("const", "var_age"),
      labels = c("constants", "age-dependent")
    )
  )

# Plot
y_max <- max(res_all$U, na.rm = TRUE)

p <- ggplot(res_all, aes(x = time, y = U, color = R0, linetype = structure)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(
    values = c(
      "R0=2.5" = "darkblue",
      "R0=5" = "darkred"
    ),
    labels = c(
      expression(R[0]^P == 2.5),
      expression(R[0]^P == 5)
    )
  ) +
  scale_linetype_manual(values = c("constants" = "21", "age-dependent" = "solid")) +
  labs(x = "Time (days)", y = expression("Incidence" ~ U[P](t)), color = NULL, linetype = NULL) +
  coord_cartesian(xlim = c(0, Tmax), ylim = c(0, y_max)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

print(p)

out_pdf <- file.path(figs_dir, "Figure-10a-incidence-fct-time-R0.pdf")
out_png <- file.path(figs_dir, "Figure-10a-incidence-fct-time-R0.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm")
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm")

cat("Plot saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
print_end_time(start_time, "plot-Figure-10a-incidence-fct-time-R0.R")
