# ============================================================
# File: plot-Figure-C1-incidence-fct-time.R
# Purpose: BETWEEN-HOST SIMULATOR – 2 MAIN CASES
#   Case 1: βc, γc, μc  (all constant)
#   Case 2: β(a), γ(a), μ(a)  (all age-dependent)
# ============================================================

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

# Set project root automatically relative to the .git tracking directory
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

# ------------------------------------------------------------
# 1. Load beta, gamma(a), mu(a) dynamically from OUTPUT
# ------------------------------------------------------------
output_dir <- file.path(project_dir, "OUTPUT")

get_latest_dist <- function(pattern) {
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found for pattern: ", pattern)
  files[which.max(file.mtime(files))]
}

beta_df       <- qs_read(get_latest_dist("^cohort_distributions_P.*_beta\\.qs$"))
gamma_overall <- qs_read(get_latest_dist("^cohort_distributions_P.*_gamma\\.qs$"))
mu_overall    <- qs_read(get_latest_dist("^cohort_distributions_P.*_mu\\.qs$"))

# ------------------------------------------------------------
# 2. Choose gamma_xi_* and mu_xid_* columns
# ------------------------------------------------------------
xi_r_target <- "4"
xi_d_target <- "85"

gamma_col <- paste0("gamma_xi_", xi_r_target)
mu_col    <- paste0("mu_xid_",   xi_d_target)

# Check if columns exist
if (!gamma_col %in% colnames(gamma_overall)) {
  stop("Column ", gamma_col, " not found in gamma_overall.")
}
if (!mu_col %in% colnames(mu_overall)) {
  stop("Column ", mu_col, " not found in mu_overall.")
}

# ------------------------------------------------------------
# 3. Align lengths and substitute missing density bounds
# ------------------------------------------------------------

combined_df <- beta_df %>%
  select(time, beta_mean, beta_q90) %>%
  left_join(gamma_overall %>% select(time, all_of(gamma_col)), by = "time") %>%
  left_join(mu_overall %>% select(time, all_of(mu_col)), by = "time") %>%
  filter(time <= 50) %>%
  mutate(
    across(c(all_of(gamma_col), all_of(mu_col)), ~ replace_na(.x, 0))
  ) %>%
  arrange(time)

a_vals         <- combined_df$time
beta_mean_vals <- combined_df$beta_mean
beta_q90_vals  <- combined_df$beta_q90
gamma_a        <- combined_df[[gamma_col]]
mu_a           <- combined_df[[mu_col]]

# Time step
dt <- mean(diff(a_vals))

# ------------------------------------------------------------
# 4. β constant vs β(a)
# ------------------------------------------------------------

beta_const_mean_scalar <- mean(beta_mean_vals, na.rm = TRUE)
beta_const_q90_scalar  <- mean(beta_q90_vals,  na.rm = TRUE)

beta_const_mean <- rep(beta_const_mean_scalar, length(a_vals))
beta_const_q90  <- rep(beta_const_q90_scalar,  length(a_vals))

# ------------------------------------------------------------
# 5. Demography parameters and kernels for (γ, μ)
# ------------------------------------------------------------

# --- Demography (modifiable) ---
S0  <- 500           # initial susceptible population
d_P <- 0.00          # natural death rate (per day)
b_P <- d_P * S0      # births so that DFE has S ≈ S0 (here 0)

# --- Constant hazards for γ, μ ---
gamma_const <- mean(gamma_a, na.rm = TRUE)
mu_const    <- mean(mu_a,    na.rm = TRUE)

# Case 1: γc, μc : constant hazard
kernel_gc_mc <- exp(-(d_P + gamma_const + mu_const) * a_vals)

# Case 2: γ(a), μ(a)
haz_ga_ma    <- d_P + gamma_a + mu_a
cumhaz_ga_ma <- cumsum(haz_ga_ma * dt)
kernel_ga_ma <- exp(-cumhaz_ga_ma)

# ------------------------------------------------------------
# 6. Simulation function  (with demography structure)
# ------------------------------------------------------------
simulate_case <- function(beta_vals, kernel) {
  n   <- length(a_vals)
  S   <- numeric(n)
  U   <- numeric(n)
  FOI <- numeric(n)
  
  # Initial conditions
  S[1]   <- S0
  U[1]   <- 1
  FOI[1] <- 0
  
  for (t in 2:n) {
    a_len  <- min(t, n)
    U_hist <- c(rep(0, n - a_len), U[1:a_len])
    
    # Force of infection kernel: ∑ β(a) U(t-a) K(a) dt
    integral <- sum(
      beta_vals[1:a_len] *
        rev(U_hist[(n - a_len + 1):n]) *
        kernel[1:a_len]
    ) * dt
    
    FOI[t] <- integral
    
    # Incidence U_P(t) = S_P(t) * ∫ β(a) U_P(t-a) K(a) da
    U[t] <- S[t - 1] * integral
    
    # Susceptibles: dS_P/dt = b_P - d_P S_P - U_P(t)
    S[t] <- S[t - 1] + dt * (b_P - d_P * S[t - 1] - U[t])
    S[t] <- max(S[t], 0)
  }
  
  tibble(time = a_vals, U = U)
}

# ------------------------------------------------------------
# 7. Run simulations ONLY for:
#    - Case 1: βc, γc, μc
#    - Case 2: β(a), γ(a), μ(a)
# ------------------------------------------------------------

# Case 1: βc, γc, μc
res_const_mean <- simulate_case(beta_const_mean, kernel_gc_mc)
res_const_q90  <- simulate_case(beta_const_q90,  kernel_gc_mc)
res_const_mean$beta_scenario <- "mean"
res_const_q90$beta_scenario  <- "q90"
res_const <- bind_rows(res_const_mean, res_const_q90) %>%
  mutate(panel = "beta_c_gamma_c_mu_c")

# Case 2: β(a), γ(a), μ(a)
res_var_mean <- simulate_case(beta_mean_vals, kernel_ga_ma)
res_var_q90  <- simulate_case(beta_q90_vals,  kernel_ga_ma)
res_var_mean$beta_scenario <- "mean"
res_var_q90$beta_scenario  <- "q90"
res_var <- bind_rows(res_var_mean, res_var_q90) %>%
  mutate(panel = "beta_a_gamma_a_mu_a")

# Combine only these two
res_all <- bind_rows(res_const, res_var) %>%
  filter(time <= 10)

# ------------------------------------------------------------
# 8. Factors and labels for 2 facets
# ------------------------------------------------------------

res_all$beta_scenario <- factor(res_all$beta_scenario,
                                levels = c("mean", "q90"))

panel_levels <- c(
  "beta_c_gamma_c_mu_c",
  "beta_a_gamma_a_mu_a"
)

panel_labels <- c(
  "beta[c]*', '*gamma[c]*', '*mu[c]",   # Case 1: all constant
  "beta(a)*', '*gamma(a)*', '*mu(a)"    # Case 2: all variable
)

res_all$panel <- factor(res_all$panel, levels = panel_levels)

res_all$panel_lab <- factor(
  res_all$panel,
  levels = panel_levels,
  labels = panel_labels
)

y_max <- max(res_all$U, na.rm = TRUE)

# ------------------------------------------------------------
# 9. Plot: 2 panels (1 row x 2 columns)
# ------------------------------------------------------------
p <- ggplot(res_all,
            aes(x = time, y = U, color = beta_scenario)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(
    ~ panel_lab,
    nrow = 1,
    labeller = label_parsed
  ) +
  scale_color_manual(
    values = c("mean" = "navy", "q90" = "darkred"),
    labels = c(
      expression(beta[mean]),
      expression(beta[Q90])
    )
  ) +
  labs(
    x = "Time (days)",
    y = "Incidence U(t)",
    color = NULL
  ) +
  coord_cartesian(ylim = c(0, y_max)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

print(p)

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

fig_dir <- file.path(project_dir, "FIGS")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(fig_dir, "Figure-10-incidence-fct-time.pdf")
out_png <- file.path(fig_dir, "Figure-10-incidence-fct-time.png")

ggsave(
  filename = out_pdf,
  plot = p,
  width = 10,
  height = 5,
  units = "in"
)

ggsave(
  filename = out_png,
  plot = p,
  width = 10,
  height = 5,
  units = "in"
)

cat("✅ Extracted between-host simulation plot saved to", out_pdf, "and", out_png, "\n")
