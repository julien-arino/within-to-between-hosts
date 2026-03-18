suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

# ------------------------------------------------------------
# 1. Load cohort data
# ------------------------------------------------------------
cat("Loading cohort data...\n")
suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
files <- list.files(output_dir, pattern = "cohort_censored_.*\\.qs$|cohort-censored_.*\\.qs$|cohort_results_truncated\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No truncated/censored cohort file found in ", output_dir)
}

latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")
cohort_df <- qs_read(latest_file)

# ------------------------------------------------------------
# 2. Function to compute death time τᵢᵈ for a given Ψ threshold
# ------------------------------------------------------------
compute_tau_d <- function(time, Psi, threshold = 85) {
  dPsi <- diff(Psi) / diff(time)
  dPsi <- c(dPsi, tail(dPsi, 1))  # align length
  idx <- which(Psi >= threshold & dPsi > 0)
  if (length(idx) == 0) return(NA_real_)
  return(time[idx[1]])
}

# ------------------------------------------------------------
# 3. Loop over several Ψ thresholds
# ------------------------------------------------------------
xi_d_values <- c(85, 95)
mu_list <- list()
summary_df <- data.frame()

cat("Computing τᵢᵈ for Ψ thresholds:", xi_d_values, "%\n")

for (xi in xi_d_values) {
  cat("\n--- Threshold Ψ ≥", xi, "% ---\n")
  
  death_times <- cohort_df %>%
    filter(status == "Dead") %>%   # Match capitalized 'Dead' from new cohort format
    group_by(individual_id) %>%    # Match individual_id
    summarise(tau_d = compute_tau_d(time, Psi, threshold = xi), .groups = "drop") %>%
    filter(!is.na(tau_d))
  
  n_dead <- nrow(death_times)
  if (n_dead > 0) {
    mean_d <- mean(death_times$tau_d)
    med_d  <- median(death_times$tau_d)
    cat("✅ Found", n_dead, "deaths.\n")
    cat("   Mean:", round(mean_d, 2), "days | Median:", round(med_d, 2), "days\n")
    
    dens <- density(death_times$tau_d, from = 0, to = 20)
    mu_df <- data.frame(a = dens$x, mu_a = dens$y, xi_d = paste0(xi))
    mu_list[[as.character(xi)]] <- mu_df
    
    summary_df <- rbind(summary_df,
                        data.frame(xi_d = xi,
                                   mean_tau = mean_d,
                                   median_tau = med_d,
                                   n = n_dead))
  } else {
    cat("⚠️  No deaths detected for threshold", xi, "%\n")
  }
}

# ------------------------------------------------------------
# 4. Combine densities and set factor order (high → low)
# ------------------------------------------------------------

mu_all <- bind_rows(mu_list)
mu_all$xi_d <- factor(mu_all$xi_d, levels = c("95", "85"))

# ------------------------------------------------------------
# 4bis. Put μ in wide format, one column per Ψ threshold
#       → one row per time, like beta_overall / gamma_overall
# ------------------------------------------------------------
mu_overall <- mu_all %>%
  mutate(
    # column names like mu_xid_90, mu_xid_85, ...
    xi_label = paste0("mu_xid_", xi_d)
  ) %>%
  select(time = a, xi_label, mu_a) %>%
  tidyr::pivot_wider(
    id_cols   = time,
    names_from  = xi_label,
    values_from = mu_a
  ) %>%
  arrange(time)


# ------------------------------------------------------------
# 5. Plot all μ_P(a) densities together
# ------------------------------------------------------------
p_mu <- ggplot(mu_all, aes(x = a, y = mu_a, color = xi_d, fill = xi_d)) +
  geom_area(alpha = 0.15, position = "identity") +
  geom_line(linewidth = 0.8) +
  labs(
    x = expression("Time since infection (days)"),
    y = "Density",
    color = expression(xi^d ~ "(%)"),
    fill  = expression(xi^d ~ "(%)")
  ) +
  scale_color_manual(
    values = c("95" = "black", "85" = "#D73027"),
    labels = c(
      expression(xi^d == 95),
      expression(xi^d == 85)
    )) +
  scale_fill_manual(
    values = c("95" = "black", "85" = "#D73027"),
    labels = c(
      expression(xi^d == 95),
      expression(xi^d == 85)
    ))+
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 0.5)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "italic"),
    legend.text  = element_text(size = 11)
  )

print(p_mu)

# ------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------

dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-07b-pdf-fd-fct-xid.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-07b-pdf-fd-fct-xid.png")

ggsave(out_pdf, plot = p_mu, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(out_png, plot = p_mu, width = 25, height = 15, units = "cm", dpi = 300)

# Densities μ(a) in wide format (like beta_overall / gamma_overall)
qs_save(mu_overall, file.path(project_dir, "OUTPUT", "mu_overall.qs"))

# Summary of death times (mean etc.)
saveRDS(summary_df, file.path(project_dir, "OUTPUT", "death_summary_thresholds.rds"))

# Optional: keep long format too
saveRDS(mu_all, file.path(project_dir, "OUTPUT", "mu_density_thresholds.rds"))

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")
cat("✅ Wide μ saved to OUTPUT/mu_overall.qs\n")
cat("✅ Long μ saved to OUTPUT/mu_density_thresholds.rds\n")
cat("✅ Summary (means/medians) saved to OUTPUT/death_summary_thresholds.rds\n")
