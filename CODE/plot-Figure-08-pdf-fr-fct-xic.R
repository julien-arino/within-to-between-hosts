library(ggplot2)
library(dplyr)
library(qs2)
library(tidyr)

# ------------------------------------------------------------
# 1. Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
output_dir <- file.path(getwd(), "OUTPUT")
files <- list.files(output_dir, pattern = "cohort_censored_.*\\.qs$|cohort-censored_.*\\.qs$|cohort_results_truncated\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No truncated/censored cohort file found in ", output_dir)
}

latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

cohort_df <- qs_read(latest_file)

# ------------------------------------------------------------
# 2. Compute beta_hat
# ------------------------------------------------------------
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

cohort_df$beta_hat <- compute_beta_hat(cohort_df$V)

# ------------------------------------------------------------
# 3. Sort + derivative
# ------------------------------------------------------------
cohort_df <- cohort_df %>%
  arrange(individual_id, time) %>%
  group_by(individual_id) %>%
  mutate(d_beta = beta_hat - lag(beta_hat)) %>%
  ungroup()

# ------------------------------------------------------------
# 4. Thresholds
# ------------------------------------------------------------
xi_c_values <- c(1e-5, 1e-3)

gamma_list <- list()

for (xi_c in xi_c_values) {
  
  cat("\nComputing for xi_c =", xi_c, "\n")
  
  tau_df <- cohort_df %>%
    group_by(individual_id) %>%
    summarise(
      
      tau_c = {
        idx <- which(beta_hat >= xi_c)
        if (length(idx) > 0) min(time[idx]) else NA
      },
      
      tau_r = {
        idx_c <- which(beta_hat >= xi_c)
        
        if (length(idx_c) == 0) {
          NA
        } else {
          
          t_c <- min(time[idx_c])
          
          idx_r <- which(
            time > t_c &
              beta_hat < xi_c &
              d_beta < 0
          )
          
          if (length(idx_r) > 0) min(time[idx_r]) else NA
        }
      },
      
      status = first(status),
      .groups = "drop"
    )
  
  tau_df <- tau_df %>%
    filter(status != "Dead") %>% # match new capitalization
    filter(!is.na(tau_c) & !is.na(tau_r))
  
  cat("Recovered patients:", nrow(tau_df), "\n")
  
  if (nrow(tau_df) > 0) {
    summary_stats <- tau_df %>%
      summarise(
        mean_tau_r   = mean(tau_r),
        median_tau_r = median(tau_r),
        n = n()
      ) %>%
      mutate(xi_c = xi_c)
    
    print(summary_stats)
    
    dens <- density(tau_df$tau_r, from = 0, to = 30)
    
    gamma_df <- data.frame(
      a = dens$x,
      gamma_a = dens$y,
      xi_c = xi_c
    )
    
    gamma_list[[as.character(xi_c)]] <- gamma_df
  }
}

gamma_all <- bind_rows(gamma_list)

# ------------------------------------------------------------
# 5. Factor ordering
# ------------------------------------------------------------
gamma_all$xi_c <- factor(gamma_all$xi_c, levels = c(1e-5, 1e-3))

# ------------------------------------------------------------
# 5bis. Put γ in wide format, one column per ξʳ
# ------------------------------------------------------------
gamma_overall <- gamma_all %>%
  mutate(
    xi_label = case_when(
      xi_c == "1e-05" | xi_c == "1e-5" ~ "gamma_xi_1e5",
      xi_c == "0.001" | xi_c == "1e-3" ~ "gamma_xi_1e3",
      TRUE ~ paste0("gamma_", xi_c)
    )
  ) %>%
  select(time = a, xi_label, gamma_a) %>%
  pivot_wider(
    id_cols   = time,
    names_from  = xi_label,
    values_from = gamma_a
  ) %>%
  arrange(time)

# ------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------
p <- ggplot(
  gamma_all,
  aes(x = a,
      y = gamma_a,
      color = xi_c,
      fill = xi_c)
) +
  
  geom_line(linewidth = 0.7) +
  geom_area(alpha = 0.2, position = "identity") +
  
  labs(
    x = "Time of infection (days)",
    y = "Density",
    color = expression(xi^c),
    fill  = expression(xi^c)
  ) +
  
  coord_cartesian(xlim = c(0,30)) +
  
  scale_color_manual(
    values = c("#e7298a", "#1b9e77"),
    labels = c(expression(10^-5), expression(10^-3))
  ) +
  
  scale_fill_manual(
    values = c("#e7298a", "#1b9e77"),
    labels = c(expression(10^-5), expression(10^-3))
  ) +
  
  theme_minimal(base_size = 14)

print(p)

# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/Figure-08-pdf-fr-fct-xic.pdf"
out_png <- "FIGS/Figure-08-pdf-fr-fct-xic.png"

ggsave(
  out_pdf,
  plot = p,
  width = 20,
  height = 8,
  units = "cm",
  dpi = 300
)

ggsave(
  out_png,
  plot = p,
  width = 20,
  height = 8,
  units = "cm",
  dpi = 300
)

qs_save(gamma_overall, "OUTPUT/gamma_overall.qs")

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")
cat("✅ Wide γ saved to OUTPUT/gamma_overall.qs\n")
