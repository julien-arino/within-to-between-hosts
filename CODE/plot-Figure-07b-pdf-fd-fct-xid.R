#!/usr/bin/env Rscript
# ============================================================
# File: plot-Figure-07b-pdf-fd-fct-xid.R
# Description:
# ============================================================

# First things first: locate project directory and load helper functions
# and constants
project_dir <- here::here()
source(file.path(project_dir, "CODE", "functions-and-definitions.R"))

# Say what we are running and start the clock
start_time <- start_time_and_hello("plot-Figure-07b-pdf-fd-fct-xid.R")

# Load libraries
load_libraries(c("ggplot2", "dplyr", "qs2", "tidyr"))

# ------------------------------------------------------------
# 1. Setup paths and load dynamically matched cohort_status files
# ------------------------------------------------------------

# --- Load pre-computed distributions for xid=85 ---
dist_files <- list.files(output_dir, pattern = "^cohort_distribution_filters_P.*_xih_75_xid_85_xic_4_xir_1\\.qs$", full.names = TRUE)
if (length(dist_files) == 0) stop("No distributions files found for xid=85 in OUTPUT")

latest_dist_file <- dist_files[which.max(file.mtime(dist_files))]
cat("Loading newest distribution for xid=85:", basename(latest_dist_file), "\n")
df_85_raw <- qs_read(latest_dist_file, nthreads = N_QS_THREADS)

df_85 <- df_85_raw$rolling_average %>%
  mutate(
    a = time,
    f_d_val = f_d,
    xi_d = factor("85", levels = c("85", "95"))
  ) %>%
  select(a, f_d_val, xi_d)

# --- Compute distributions locally for xid=95 ---
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xih_75_xid_95_xic_4_xir_1\\.qs$", full.names = TRUE)
if (length(status_files) == 0) stop("No cohort_status file found for xid=95 in OUTPUT")
latest_status_file <- status_files[which.max(file.mtime(status_files))]
cat("Loading newest cohort status for xid=95:", basename(latest_status_file), "\n")
cohort_status_95 <- qs_read(latest_status_file, nthreads = N_QS_THREADS)

cat("Computing f_d for xid=95 locally...\n")
time_pts <- seq(0, 100.0, by = 0.1)

dead_cohort_95 <- cohort_status_95 %>% filter(status == "Dead")
dens_95 <- density(dead_cohort_95$tau_d, from = 0, to = max(time_pts), n = length(time_pts))
f_d_95 <- dens_95$y

df_95 <- data.frame(
  a = time_pts,
  f_d_val = f_d_95,
  xi_d = factor("95", levels = c("85", "95"))
)

# Combine both for plotting
mu_all <- dplyr::bind_rows(df_85, df_95)


# ------------------------------------------------------------
# 5. Plot all μ_P(a) densities together
# ------------------------------------------------------------
p_mu <- ggplot(mu_all, aes(x = a, y = f_d_val, color = xi_d, fill = xi_d)) +
  geom_area(alpha = 0.15, position = "identity") +
  geom_line(linewidth = 0.8) +
  labs(
    x = expression("Time since infection (days)"),
    y = expression(f[d](a)),
    color = expression(xi^d ~ "(%)"),
    fill = expression(xi^d ~ "(%)")
  ) +
  scale_color_manual(
    values = c("85" = "firebrick3", "95" = "mistyrose3"),
    labels = c(
      "85" = expression(xi^d == 85),
      "95" = expression(xi^d == 95)
    )
  ) +
  scale_fill_manual(
    values = c("85" = "firebrick3", "95" = "mistyrose3"),
    labels = c(
      "85" = expression(xi^d == 85),
      "95" = expression(xi^d == 95)
    )
  ) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, max(mu_all$f_d_val, na.rm = TRUE) + 0.01)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "italic"),
    legend.text = element_text(size = 11)
  )

print(p_mu)

# ------------------------------------------------------------
# 3. Save outputs
# ------------------------------------------------------------

out_pdf <- file.path(figs_dir, "Figure-07b-pdf-fd-fct-xid.pdf")
out_png <- file.path(figs_dir, "Figure-07b-pdf-fd-fct-xid.png")

ggsave(out_pdf, plot = p_mu, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(out_png, plot = p_mu, width = 25, height = 15, units = "cm", dpi = 300)

cat("Figures saved:\n -", out_pdf, "\n -", out_png, "\n")

print_end_time(start_time, "plot-Figure-07b-pdf-fd-fct-xid.R")
