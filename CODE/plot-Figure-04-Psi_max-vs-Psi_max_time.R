# ============================================================
# File: plot-Psi_max-vs-Psi_max_time-Figure-04.R
# Goal:
#   1) Load truncated cohort results to get Psi_max & status per patient
#   2) Load process-virtual-cohort-results.csv to get R0_within and tau_Psi_max
#   3) Join and make figure (Psi_max vs t_Psi_max, split by R0_within<1 vs >=1)
# Outputs:
#   - FIGS/plot_Psi_max_vs_Psi_max_time_Figure_04.png and .pdf
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(qs2)

dir.create("OUTPUT", showWarnings = FALSE, recursive = TRUE)
dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)

output_dir <- file.path(getwd(), "OUTPUT")

# 1a) Find latest processed times .qs file
times_files <- list.files(output_dir, pattern = "^cohort_times_P.*\\.qs$", full.names = TRUE)
if (length(times_files) == 0) stop("Cannot find times qs file in ", output_dir)
TIMES_FILE <- times_files[which.max(file.mtime(times_files))]

# ---------------------------
# 1) Load data
# ---------------------------
cat("Loading processed times:", basename(TIMES_FILE), "\n")
summary_df <- qs_read(TIMES_FILE)

# ---------------------------
# 2) Formatting for plotting (Create status, rename columns to match plot)
# ---------------------------
# The manuscript definitions:
# tau_d present = Dead
# tau_h present = ICU (Hospitalization)
# otherwise = Mild
summary_df <- summary_df %>%
  mutate(
    status = case_when(
      !is.na(tau_d) ~ "Dead",
      !is.na(tau_h_start) ~ "ICU",
      TRUE ~ "Mild"
    ),
    R0 = R0_within,
    t_Psi_max = tau_Psi_max,
    Psi_max = Psi_max
  ) %>%
  filter(!is.na(R0))

summary_df$status <- factor(summary_df$status, levels = c("Mild", "ICU", "Dead"))

status_colors <- c(
  Mild = "dodgerblue4",
  ICU  = "orange",
  Dead = "red"
)

# ---------------------------
# 5) Create plot
# ---------------------------
p2 <- ggplot(summary_df, aes(x = t_Psi_max, y = Psi_max)) +
  
  # Patients with R0 >= 1 (colored by status)
  geom_point(
    data = subset(summary_df, R0 >= 1),
    aes(color = status),
    alpha = 0.18,
    size = 0.6
  ) +
  
  # Patients with R0 < 1 (DARK GREEN)
  geom_point(
    data = subset(summary_df, R0 < 1),
    color = "darkgreen",
    alpha = 0.4,
    size = 0.8
  ) +
  
  coord_cartesian(xlim=c(0,100), ylim=c(0,100)) +
  
  scale_color_manual(
    values = status_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    x = expression("Time in days to maximum tissue damage"~tau[Psi[max]]),
    y = expression("Maximum tissue damage"~Psi[max]),
    parse = TRUE
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 2, alpha = 1)
    )
  )

p2 <- p2 +
  annotate("text",
           x = 1, y = 15,
           label = "R[0] < 1",
           size = 4,
           color = "darkgreen",
           parse = TRUE,
           fontface = "bold") +
  annotate("text",
           x = 45, y = 95,
           label = "R[0] >= 1",
           size = 4,
           color = "black",
           parse = TRUE,
           fontface = "bold")
print(p2)

# ---------------------------
# 6) Save figure
# ---------------------------
out_png <- "FIGS/plot_Psi_max_vs_Psi_max_time_Figure_04.png"
out_pdf <- "FIGS/plot_Psi_max_vs_Psi_max_time_Figure_04.pdf"

ggsave(filename = out_png, plot = p2, width = 25, height = 15, units = "cm", dpi = 300)
ggsave(filename = out_pdf, plot = p2, width = 25, height = 15, units = "cm", dpi = 300)

cat("\nFigure 4 saved to:\n  -", out_png, "\n  -", out_pdf, "\n")
