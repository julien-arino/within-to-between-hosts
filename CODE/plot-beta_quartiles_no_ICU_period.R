# ============================================================
# File: compute_beta_clean_ICU_isolation.R
# Plot with mean ± SD ribbon + q10 and q90 lines + LEGEND
# ICU do not contribute when Psi > 75
# ============================================================
library(qs2)
library(ggplot2)
library(dplyr)

SHOW_TITLES <- FALSE

# 1) Load the cohort data
cohort_df <- qs_read("OUTPUT/cohort_results_truncated.qs")

# 2) Hill function
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

# 3) Compute β̂(t)
cohort_df$beta_hat <- compute_beta_hat(cohort_df$V)

# 3b) Apply ICU isolation rule:
# - Mild, Dead: beta_hat unchanged
# - ICU: contribution = 0 when Psi > 75, beta_hat otherwise
cohort_df <- cohort_df %>%
  mutate(
    beta_hat_eff = if_else(
      status == "ICU" & Psi > 75,
      0,          # ICU isolated when Psi > 75
      beta_hat    # otherwise normal contribution
    )
  )

# 4) Aggregate by time (using beta_hat_eff)
beta_overall <- cohort_df %>%
  group_by(time) %>%
  summarise(
    beta_mean = mean(beta_hat_eff, na.rm = TRUE),
    beta_sd   = sd(beta_hat_eff,   na.rm = TRUE),
    beta_q10  = quantile(beta_hat_eff, probs = 0.10, na.rm = TRUE),
    beta_q90  = quantile(beta_hat_eff, probs = 0.90, na.rm = TRUE),
    .groups   = "drop"
  )

# Save for future plots
qs_save(beta_overall, "OUTPUT/beta_overall_ICU_isolation.qs")

# 5) Plot
p_beta_hat <- ggplot(beta_overall, aes(x = time)) +
  # ribbon: mean ± SD
  geom_ribbon(aes(ymin = pmax(beta_mean - beta_sd, 0),
                  ymax = beta_mean + beta_sd,
                  fill = "mean ± SD"),
              alpha = 0.15,
              color = NA) +
  # mean line
  geom_line(aes(y = beta_mean, color = "mean"), linewidth = 1.3) +
  # q10 and q90
  geom_line(aes(y = beta_q10, color = "q10"), linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = beta_q90, color = "q90"), linetype = "dashed", linewidth = 0.6) +
  labs(
    title = if (SHOW_TITLES) "Transmissibility (ICU isolated when Psi > 75)" else NULL,
    x = "Time (days)",
    y = "Density",
    color = "",
    fill  = ""
  ) +
  scale_color_manual(
    values = c(
      "mean" = "dodgerblue4",
      "q10"  = "black",
      "q90"  = "orange"
    ),
    labels = c(
      "mean" = "mean",
      "q10"  = "Q10",
      "q90"  = "Q90"
    )
  ) +
  scale_fill_manual(
    values = c("mean ± SD" = "dodgerblue4"),
    labels = c("mean ± SD" = "mean ± SD")
  ) +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 0.55)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# 6) Save
ggsave("FIGS/beta_hat_overall_ICU_iso.png",
       plot = p_beta_hat, width = 15, height = 10, units = "cm", dpi = 300)

ggsave(
  filename = "/home/cdjuikem/overleaf-within-to-between-hosts-new/FIGS/beta_hat_quartiles_ICU_iso.png",
  plot = p_beta_hat,
  width = 25, height = 15, units = "cm", dpi = 300
)

ggsave(
  filename = "/home/cdjuikem/overleaf-within-to-between-hosts-new/FIGS/beta_hat_quartiles_ICU_iso.pdf",
  plot = p_beta_hat,
  width = 25, height = 15, units = "cm", dpi = 300
)

print(p_beta_hat)
