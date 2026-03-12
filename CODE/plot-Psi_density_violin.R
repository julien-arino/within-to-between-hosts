library(ggplot2)
library(dplyr)
library(qs)

# --- Load data ---
cohort_df <- qread("OUTPUT_clo/cohort_results_truncated.qs")

# --- One Psi_max per individual ---
individual_summary <- cohort_df %>%
  group_by(individual_id, status) %>%
  summarise(Psi_max = max(Psi), .groups = "drop")

# --- Thresholds we want to compare ---
xi_values <- c(85, 90, 95)

# --- Build cumulative groups ---
psi_groups <- bind_rows(
  lapply(xi_values, function(th) {
    individual_summary %>%
      filter(Psi_max >= th) %>%
      mutate(xi_group = factor(th, levels = xi_values))   # << simplicité
  })
)

# Check counts
print(psi_groups %>% count(xi_group))

# --- Colors (one time only) ---
col_map <- c("85"="red1", "90"="#ef8a62", "95"="darkred")

# ------------------------------------------------------------
# Density curves ------------------------------------------------------------
p_density <- ggplot(psi_groups,
                    aes(x = Psi_max, color = xi_group, fill = xi_group)) +
  geom_density(alpha = 0.25, linewidth = 1.2) +
  scale_color_manual(values = col_map, labels = xi_values) +
  scale_fill_manual(values = col_map, labels = xi_values) +
  labs(
    title = "",
    x = expression(Psi[max]~"(%)"),
    y = "Density",
    color = expression(xi^d~"(%)"),
    fill  = expression(xi^d~"(%)")
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_density)

ggsave("/home/cdjuikem/overleaf-within-to-between-hosts-new/FIGS/Psi_density_thresholds.pdf", p_density,
       width = 20, height = 12, units = "cm", dpi = 300)
# ------------------------------------------------------------
# Violin version ------------------------------------------------------------
p_violin <- ggplot(psi_groups,
                   aes(x = xi_group, y = Psi_max, fill = xi_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.alpha = 0.2, fill = "white") +
  scale_fill_manual(values = col_map, labels = xi_values) +
  scale_x_discrete(labels = xi_values) +
  labs(
    title = "",
    x = expression(xi^d~"(%)"),
    y = expression(Psi[max]~"(%)")
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_violin)

ggsave("/home/cdjuikem/overleaf-within-to-between-hosts-new/FIGS/Psi_violin_thresholds.pdf", p_violin,
       width = 20, height = 12, units = "cm", dpi = 300)
