# ============================================================
# File: plot-Figure-C4-time-evolution-quantities.R
# Purpose: 6-panel plot showing time evolution of V, beta, I, Psi, FB, FU
# ============================================================

suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

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

# Ensure capitalized status
cohort_df$status <- factor(cohort_df$status, levels = c("Mild", "ICU", "Dead"))

# ------------------------------------------------------------
# 2️⃣ Parameters
# ------------------------------------------------------------
xi_c <- 0.001
xi_h <- 75
alpha <- 16.422
k_v <-  7.49
V_star <- k_v * (xi_c / (1 - xi_c))^(1/alpha)

# ------------------------------------------------------------
# 3️⃣ Compute beta_hat
# ------------------------------------------------------------
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

cohort_df$beta_hat <- compute_beta_hat(cohort_df$V)

# ------------------------------------------------------------
# 4️⃣ Sort + derivative
# ------------------------------------------------------------
cohort_df <- cohort_df %>%
  arrange(individual_id, time) %>%
  group_by(individual_id) %>%
  mutate(d_beta = beta_hat - lag(beta_hat)) %>%
  ungroup()

# ------------------------------------------------------------
# 5️⃣ Compute tau_c, tau_r, tau_h
# ------------------------------------------------------------
tau_df <- cohort_df %>%
  group_by(individual_id) %>%
  summarise(
    tau_c = {
      idx <- which(beta_hat >= xi_c)
      if (length(idx) > 0) min(time[idx]) else Inf
    },
    tau_r = {
      idx <- which(beta_hat >= xi_c)
      if (length(idx) == 0) Inf else {
        t_c <- min(time[idx])
        idx_r <- which(time >= t_c &
                         beta_hat <= xi_c &
                         d_beta < 0)
        if (length(idx_r) > 0) min(time[idx_r]) else Inf
      }
    },
    tau_h_start = {
      idx <- which(Psi >= xi_h)
      if (length(idx) > 0) min(time[idx]) else Inf
    },
    tau_h_end = {
      idx <- which(Psi >= xi_h)
      if (length(idx) > 0) max(time[idx]) else -Inf
    },
    .groups = "drop"
  )

tau_c_vec       <- setNames(tau_df$tau_c, tau_df$individual_id)
tau_r_vec       <- setNames(tau_df$tau_r, tau_df$individual_id)
tau_h_start_vec <- setNames(tau_df$tau_h_start, tau_df$individual_id)
tau_h_end_vec   <- setNames(tau_df$tau_h_end, tau_df$individual_id)

# ------------------------------------------------------------
# 6️⃣ Build effective beta_c(a)
# ------------------------------------------------------------
cohort_df <- cohort_df %>%
  mutate(
    tau_c       = tau_c_vec[as.character(individual_id)],
    tau_r       = tau_r_vec[as.character(individual_id)],
    tau_h_start = tau_h_start_vec[as.character(individual_id)],
    tau_h_end   = tau_h_end_vec[as.character(individual_id)],
    
    beta_hat_window = if_else(
      time >= tau_c &
        time <= tau_r &
        !(time >= tau_h_start & time <= tau_h_end),
      beta_hat,
      0
    )
  )

# ------------------------------------------------------------
# 7️⃣ Compute mean ± SD per status
# ------------------------------------------------------------
mean_df <- cohort_df %>%
  group_by(status, time) %>%
  summarise(
    Psi_mean = mean(Psi, na.rm = TRUE),
    Psi_sd   = sd(Psi,   na.rm = TRUE),
    V_mean   = mean(V,   na.rm = TRUE),
    V_sd     = sd(V,     na.rm = TRUE),
    I_mean   = mean(I,   na.rm = TRUE),
    I_sd     = sd(I,     na.rm = TRUE),
    F_B_mean = mean(F_B, na.rm = TRUE),
    F_B_sd   = sd(F_B,   na.rm = TRUE),
    F_U_mean = mean(F_U, na.rm = TRUE),
    F_U_sd   = sd(F_U,   na.rm = TRUE),
    beta_mean = mean(beta_hat_window, na.rm = TRUE),
    beta_sd   = sd(beta_hat_window,   na.rm = TRUE),
    .groups = "drop"
  )

status_cols <- c(
  Mild = "dodgerblue4",
  ICU  = "orange",
  Dead = "red"
)

# ------------------------------------------------------------
# 8️⃣ Helper function
# ------------------------------------------------------------
make_plot <- function(y_mean, y_sd, ylab, xlim=NULL, ylim=NULL, legend=FALSE) {
  
  ggplot(mean_df, aes(x=time,
                      y=!!sym(y_mean),
                      color=status,
                      fill=status)) +
    geom_ribbon(
      aes(ymin=pmax(!!sym(y_mean)-!!sym(y_sd),0),
          ymax=     !!sym(y_mean)+!!sym(y_sd)),
      alpha=0.22, color=NA
    ) +
    geom_line(linewidth=1.1) +
    scale_color_manual(values=status_cols) +
    scale_fill_manual(values=status_cols) +
    labs(x="Time (days)", y=ylab) +
    { if(!is.null(ylim)) coord_cartesian(ylim=ylim) else NULL } +
    {if(!is.null(xlim))  coord_cartesian(xlim=xlim)}+
    theme_minimal(base_size=13) +
    theme(
      legend.position = if (legend) "top" else "none",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ------------------------------------------------------------
# 9️⃣ Create the 6 panels
# ------------------------------------------------------------
p_V    <- make_plot("V_mean","V_sd",
                    expression("Viral load"~log[10]*"(copies/ml)"))

p_beta <- make_plot("beta_mean","beta_sd",
                    expression("Transmission rate"),
                    xlim=c(0,30))

p_I    <- make_plot("I_mean","I_sd",
                    expression("Infected cells"~10^9*" cells/ml"))

p_Psi  <- make_plot("Psi_mean","Psi_sd",
                    expression("Lung tissue damage"~Psi(t)~"%"),
                    ylim=c(0,100),
                    legend=TRUE)

p_FB   <- make_plot("F_B_mean","F_B_sd",
                    expression("Bound IFN"~10^{-5}~"pg/ml"))

p_FU   <- make_plot("F_U_mean","F_U_sd",
                    expression("Unbound IFN (pg/ml)"))

# ------------------------------------------------------------
# 🔟 Combine 2 rows × 3 columns
# ------------------------------------------------------------
combined_all <-
  (p_V | p_beta | p_I) /
  (p_Psi | p_FB | p_FU)

combined_all <- combined_all +
  plot_layout(guides="collect") &
  theme(legend.position="top")

print(combined_all)

# ------------------------------------------------------------
# Save single-page figure
# ------------------------------------------------------------

dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "FIGS/Figure-C4-time-evolution-quantities.pdf"
out_png <- "FIGS/Figure-C4-time-evolution-quantities.png"

ggsave(
  out_pdf,
  plot = combined_all,
  width = 26, height = 16, units = "cm", dpi = 300
)

ggsave(
  out_png,
  plot = combined_all,
  width = 26, height = 16, units = "cm", dpi = 300
)

cat("✅ Extracted 6-panel evolution plot saved to", out_pdf, "and", out_png, "\n")
