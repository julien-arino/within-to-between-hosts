suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

# ------------------------------------------------------------
# 1. Setup paths and load the full continuous state array
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-08-pdf-fr-fct-xic.R ...\n\n")

suppressWarnings(suppressPackageStartupMessages(library(here)))
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}
output_dir <- file.path(project_dir, "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-all.R")) else source("functions-all.R")
}

# --- Find and load pre-computed distributions ---
dist_files <- list.files(output_dir, pattern = "^cohort_distributions_P.*_gamma\\.qs$", full.names = TRUE)
if (length(dist_files) == 0) stop("No gamma distributions file found in OUTPUT")
latest_dist <- dist_files[which.max(file.mtime(dist_files))]
cat("Loading newest gamma distributions:", basename(latest_dist), "\n")
gamma_overall <- qs_read(latest_dist, nthreads = N_QS_THREADS)

# ------------------------------------------------------------
# 2. Format to long (gamma_all) for plotting
# ------------------------------------------------------------
gamma_all <- gamma_overall %>%
  tidyr::pivot_longer(
    cols = starts_with("gamma_xi_"),
    names_to = "xi_label",
    values_to = "gamma_a"
  ) %>%
  mutate(
    xi_c = sub("gamma_xi_", "", xi_label),
    a = time
  ) %>%
  mutate(xi_c = factor(xi_c, levels = c("4", "6")))

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
    x = "Age of infection (days)",
    y = "Probability density function",
    color = expression(xi^c),
    fill  = expression(xi^c)
  ) +
  
  coord_cartesian(xlim = c(0,35)) +
  
  scale_color_manual(
    values = c("#e7298a", "#1b9e77"),
    labels = c(expression(4), expression(6))
  ) +
  
  scale_fill_manual(
    values = c("#e7298a", "#1b9e77"),
    labels = c(expression(4), expression(6))
  ) +
  
  theme_minimal(base_size = 14)

print(p)

# ------------------------------------------------------------
# 3. Save results
# ------------------------------------------------------------
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-08-pdf-fr-fct-xic.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-08-pdf-fr-fct-xic.png")

ggsave(out_pdf, plot = p, width = 20, height = 8, units = "cm", dpi = 300)
ggsave(out_png, plot = p, width = 20, height = 8, units = "cm", dpi = 300)

cat("\n✅ Figures saved to", out_pdf, "and", out_png, "\n")
