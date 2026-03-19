#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

# ------------------------------------------------------------
# 1. Setup paths and load dynamically matched cohort_status files
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-07b-pdf-fd-fct-xid.R ...\n\n")

suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

# --- Find and load pre-computed distributions ---
dist_files <- list.files(output_dir, pattern = "^cohort_distributions_P.*_mu\\.qs$", full.names = TRUE)
if (length(dist_files) == 0) stop("No mu distributions file found in OUTPUT")
latest_dist <- dist_files[which.max(file.mtime(dist_files))]
cat("Loading newest mu distributions:", basename(latest_dist), "\n")
mu_overall <- qs_read(latest_dist, nthreads = N_QS_THREADS)

# ------------------------------------------------------------
# 2. Format to long (mu_all) for plotting
# ------------------------------------------------------------
mu_all <- mu_overall %>%
  tidyr::pivot_longer(
    cols = starts_with("mu_xid_"),
    names_to = "xi_label",
    values_to = "mu_a"
  ) %>%
  mutate(
    xi_d = sub("mu_xid_", "", xi_label),
    a = time
  ) %>%
  mutate(xi_d = factor(xi_d, levels = c("95", "85")))


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
    values = c("85" = "firebrick3", "95" = "mistyrose3"),
    labels = c(
      "85" = expression(xi^d == 85),
      "95" = expression(xi^d == 95)
    )) +
  scale_fill_manual(
    values = c("85" = "firebrick3", "95" = "mistyrose3"),
    labels = c(
      "85" = expression(xi^d == 85),
      "95" = expression(xi^d == 95)
    ))+
  coord_cartesian(xlim = c(0, 20), ylim = c(0, max(mu_all$mu_a, na.rm=TRUE) + 0.1)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "italic"),
    legend.text  = element_text(size = 11)
  )

print(p_mu)

# ------------------------------------------------------------
# 3. Save outputs
# ------------------------------------------------------------

dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-07b-pdf-fd-fct-xid.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-07b-pdf-fd-fct-xid.png")

ggsave(out_pdf, plot = p_mu, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(out_png, plot = p_mu, width = 25, height = 15, units = "cm", dpi = 300)

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")
