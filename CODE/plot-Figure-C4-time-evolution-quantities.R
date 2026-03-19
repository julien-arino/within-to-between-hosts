#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-C4-time-evolution-quantities.R
# Purpose: 6-panel plot showing time evolution of V, beta, I, Psi, FB, FU
# ============================================================

suppressWarnings(suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(here)
  # library(parallel)
  library(future.apply)
}))

# ------------------------------------------------------------
# 1. Automatically find the latest files
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-C4-time-evolution-quantities.R ...\n\n")

# Set project root automatically relative to the .git tracking directory
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

# Load sim_state instead of truncated state
files <- list.files(output_dir, pattern = "^cohort_sim_state_.*\\.qs$", full.names = TRUE)
if (length(files) == 0) stop("No cohort_sim_state file found in ", output_dir)
latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort trajectories:", basename(latest_file), "\n")
cohort_list <- qs_read(latest_file, nthreads = N_QS_THREADS)

# Load baseline statuses
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xih_75_xid_85\\.qs$|^cohort_status_P.*_xid_85\\.qs$", full.names = TRUE)
if (length(status_files) == 0) stop("No baseline cohort_status file found")
latest_status_file <- status_files[which.max(file.mtime(status_files))]
cat("Loading newest cohort baseline statuses:", basename(latest_status_file), "\n")
cohort_status <- qs_read(latest_status_file, nthreads = N_QS_THREADS)

status_map <- setNames(as.character(cohort_status$status), as.character(cohort_status$ID))

# ------------------------------------------------------------
# 2. Parallel interpolation of trajectories
# ------------------------------------------------------------
cat("\nStarting future multisession for parallel interpolation...\n")
n_cores <- parallel::detectCores()
workers_to_use <- max(2, n_cores - 2)
if (n_cores >= 64) workers_to_use <- max(2, round(n_cores * 2 / 3))

plan(multisession, workers = workers_to_use)

# We need a unified global max time to interpolate everyone out to
max_time <- max(sapply(cohort_list, function(lst) max(lst$time, na.rm = TRUE)))

interpolate_trajectory <- function(lst_element) {
  if (is.null(lst_element)) return(NULL)
  
  df <- as.data.frame(lst_element)
  if (nrow(df) == 0) return(NULL)
  
  # Determine death naturally from the trajectory
  xi_d <- 85
  is_dead <- FALSE
  t_death <- max_time
  
  if ("Psi" %in% names(df)) {
    psi_max <- max(df$Psi, na.rm = TRUE)
    if (psi_max >= xi_d) {
      is_dead <- TRUE
      idx_death <- which(df$Psi >= xi_d)[1]
      t_death <- df$time[idx_death]
    }
  }
  
  # Interpolate strictly to the global max_time on a 0.1 grid
  t_interp <- seq(0, max_time, by = 0.1)
  
  res <- data.frame(
    time = t_interp,
    stringsAsFactors = FALSE
  )
  
  vars_to_interp <- c("V", "I", "F_B", "F_U", "Psi")
  
  for (var in vars_to_interp) {
    if (var %in% names(df)) {
      active_df <- df[df$time <= t_death, ]
      
      interpolated_vals <- approx(active_df$time, active_df[[var]], xout = t_interp, rule = 2, ties = mean)$y
      if (is_dead) {
        interpolated_vals[t_interp > t_death] <- NA
      }
      
      res[[var]] <- interpolated_vals
    } else {
      res[[var]] <- NA
    }
  }
  
  return(res)
}

cohort_list_interp <- future_lapply(cohort_list, interpolate_trajectory, future.seed = TRUE)

# Explicitly ensure the resulting list names match the original cohort list layout
if (!is.null(names(cohort_list))) {
  names(cohort_list_interp) <- names(cohort_list)
} else {
  names(cohort_list_interp) <- as.character(seq_along(cohort_list_interp))
}

# Flush the (potentially massive) original trajectory list from memory
rm(cohort_list)
gc()

# Assemble the complete new cohort metric tensor
cohort_df <- dplyr::bind_rows(Filter(Negate(is.null), cohort_list_interp), .id = "ID")

# Stitch on their severity status purely from the official baseline output
cohort_df <- cohort_df %>%
  inner_join(cohort_status %>% select(ID, status) %>% mutate(ID = as.character(ID)), 
             by = "ID")

cohort_df$status <- factor(cohort_df$status, levels = c("Mild", "ICU", "Dead"))

# ------------------------------------------------------------
# 3. Parameters
# ------------------------------------------------------------
xi_c <- 4       # start of infectious period (viral load threshold)
xi_r <- 1.0     # end of infectious period (viral load threshold)
xi_h <- 75
alpha <- 16.422
k_v <-  7.49

# ------------------------------------------------------------
# 4. Compute beta_hat
# ------------------------------------------------------------
compute_beta_hat <- function(V, alpha = 16.422, k_v = 7.49) {
  (V^alpha) / (V^alpha + k_v^alpha)
}

cohort_df$beta_hat <- compute_beta_hat(cohort_df$V)

# ------------------------------------------------------------
# 5. Compute tau_c, tau_r, tau_h
# ------------------------------------------------------------
tau_df <- cohort_df %>%
  group_by(ID) %>%
  summarise(
    tau_c = {
      idx <- which(V >= xi_c)
      if (length(idx) > 0) min(time[idx]) else Inf
    },
    tau_r = {
      idx <- which(V >= xi_c)
      if (length(idx) == 0) Inf else {
        V_peak_idx <- which.max(V)
        t_peak <- time[V_peak_idx]
        idx_r <- which(time > t_peak & V <= xi_r)
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

tau_c_vec       <- setNames(tau_df$tau_c, tau_df$ID)
tau_r_vec       <- setNames(tau_df$tau_r, tau_df$ID)
tau_h_start_vec <- setNames(tau_df$tau_h_start, tau_df$ID)
tau_h_end_vec   <- setNames(tau_df$tau_h_end, tau_df$ID)

# ------------------------------------------------------------
# 7. Build effective beta_c(a)
# ------------------------------------------------------------
cohort_df <- cohort_df %>%
  mutate(
    tau_c       = tau_c_vec[as.character(ID)],
    tau_r       = tau_r_vec[as.character(ID)],
    tau_h_start = tau_h_start_vec[as.character(ID)],
    tau_h_end   = tau_h_end_vec[as.character(ID)],
    
    beta_hat_window = if_else(
      time >= tau_c &
        time <= tau_r &
        !(time >= tau_h_start & time <= tau_h_end),
      beta_hat,
      0
    )
  )

# ------------------------------------------------------------
# 8. Compute mean ± SD per status
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
# 9. Helper function
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
# 10. Create the 6 panels
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

dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-C4-time-evolution-quantities.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-C4-time-evolution-quantities.png")

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

cat("✅ 6-panel evolution plot saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
