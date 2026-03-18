suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

# ------------------------------------------------------------
# 1. Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
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

xi_h <- 75
xi_c <- 0.001

# ------------------------------------------------------------
# 2. Compute beta
# ------------------------------------------------------------
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

cohort_df$beta_hat <- compute_beta_hat(cohort_df$V)

# ------------------------------------------------------------
# 3. Hospitalisation times
# ------------------------------------------------------------
tau_df <- cohort_df %>%
  group_by(individual_id) %>% # Updated from patient_id
  summarise(
    tau_h_start = {
      idx <- which(Psi >= xi_h)
      if(length(idx)>0) min(time[idx]) else Inf
    },
    tau_h_end = {
      idx <- which(Psi >= xi_h)
      if(length(idx)>0) max(time[idx]) else -Inf
    },
    .groups="drop"
  )

cohort_df <- left_join(cohort_df, tau_df, by="individual_id")

# ------------------------------------------------------------
# 4. Remove hospitalisation period
# ------------------------------------------------------------
cohort_df <- cohort_df %>%
  mutate(
    beta_eff = if_else(
      time >= tau_h_start & time <= tau_h_end,
      0,
      beta_hat
    )
  )

# ------------------------------------------------------------
# 5. Keep transmitters at each time
# ------------------------------------------------------------
cohort_transmitters <- cohort_df %>%
  filter(beta_eff >= xi_c)

# ------------------------------------------------------------
# 6. Compute statistics
# ------------------------------------------------------------
beta_overall <- cohort_transmitters %>%
  group_by(time) %>%
  summarise(
    beta_mean   = mean(beta_eff, na.rm = TRUE),
    beta_sd     = sd(beta_eff,   na.rm = TRUE),
    beta_q10    = quantile(beta_eff, 0.10, na.rm = TRUE),
    beta_q90    = quantile(beta_eff, 0.90, na.rm = TRUE),
    beta_median = median(beta_eff, na.rm = TRUE),
    .groups     = "drop"
  )

qs_save(
  beta_overall,
  file.path(project_dir, "OUTPUT", "beta_overall_transmitters.qs")
)

# ------------------------------------------------------------
# 7. Plot
# ------------------------------------------------------------
p <- ggplot(beta_overall, aes(x=time)) +
  
  geom_ribbon(
    aes(
      ymin = pmax(beta_mean-beta_sd,0),
      ymax = beta_mean+beta_sd,
      fill="mean ± SD"
    ),
    alpha=0.15,
    color=NA
  )+
  
  geom_line(
    aes(y=beta_mean,color="mean"),
    linewidth=1.3
  )+
  
  geom_line(
    aes(y=beta_q10,color="Q10"),
    linetype="22",
    linewidth=0.7
  )+
  
  geom_line(
    aes(y=beta_median,color="median"),
    linewidth=1
  )+
  
  geom_line(
    aes(y=beta_q90,color="Q90"),
    linetype="22",
    linewidth=0.7
  )+
  
  scale_color_manual(
    breaks = c("mean","Q10","median","Q90"),
    values = c(
      "mean"   = "dodgerblue4",
      "Q10"    = "black",
      "median" = "magenta4",
      "Q90"    = "orange"
    )
  )+
  
  scale_fill_manual(
    values=c("mean ± SD"="dodgerblue4")
  )+
  
  coord_cartesian(xlim=c(0,80))+
  
  labs(
    x="Time (days)",
    y="Transmission rate",
    color="",
    fill=""
  )+
  
  theme_minimal(base_size=14)

print(p)

# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-09-transmission-rates-fct-summary-fcts.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-09-transmission-rates-fct-summary-fcts.png")

ggsave(
  out_pdf,
  plot = p,
  width = 20,
  height = 8,
  units = "cm"
)

ggsave(
  out_png,
  plot = p,
  width = 20,
  height = 8,
  units = "cm"
)

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")
cat("✅ Summary matrix saved to OUTPUT/beta_overall_transmitters.qs\n")
