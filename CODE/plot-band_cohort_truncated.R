# ============================================================
# File: plot_cohort_truncated_bands.R
# Purpose: Plot mean ± SD bands per status (mild / ICU / dead)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(qs2)
  library(scales)
})

# ------------------------------------------------------------
# 1️⃣ Load truncated cohort
# ------------------------------------------------------------
cohort_df <- qs_read("OUTPUT/cohort_results_truncated.qs")
cat("Loaded:", nrow(cohort_df), "rows\n")


cohort_df$status <- recode(cohort_df$status,
                           rest = "Mild",
                           ICU  = "ICU",
                           dead = "Dead")
# ------------------------------------------------------------
# 2️⃣ Compute mean & SD per time *and* per status
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
    .groups = "drop"
  )

status_cols <- c(
  "Mild" = "dodgerblue3",
  "ICU"  = "orange",
  "Dead" = "red3"
)


status_fills <- status_cols  # fixed
# ------------------------------------------------------------
# 🔧 Helper: function to make ribbon plot
# ------------------------------------------------------------

make_band_plot <- function(df, y_mean, y_sd, title, ylab, ylim=NULL, show_legend = FALSE) {
  
  ggplot(df, aes(x = time, y = !!sym(y_mean), color = status, fill = status)) +
    
    geom_ribbon(
      aes(ymin = pmax(!!sym(y_mean) - !!sym(y_sd), 0),
          ymax =      !!sym(y_mean) + !!sym(y_sd)),
      alpha = 0.22, color = NA
    ) +
    
    geom_line(linewidth = 1.1) +
    
    scale_color_manual(values = status_cols) +
    scale_fill_manual(values = status_cols) +
    
    labs(title = "", x = "Time (days)", y = ylab) +
    
    { if(!is.null(ylim)) coord_cartesian(ylim = ylim) else NULL } +
    
    theme_bw(base_size = 12) +
    theme(
      #legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face="bold", hjust=0.5),
      legend.position = if (show_legend) "top" else "none"
    )
}

# ------------------------------------------------------------
# 3️⃣ Generate the 5 plots
# ------------------------------------------------------------

p1 <- make_band_plot(mean_df, "Psi_mean", "Psi_sd",
                     "Lung tissue damage Ψ(t)",
                     expression("Lung tissue damage"~Psi(t)~"%"),
                     ylim=c(0,100), show_legend = TRUE)

p2 <- make_band_plot(mean_df, "V_mean", "V_sd",
                     "Viral load",
                     expression("Viral load"~log[10]*"(copies/ml)"),
                     show_legend = FALSE)

p3 <- make_band_plot(mean_df, "I_mean", "I_sd",
                     "Infected cells",
                     expression("Infected cells"~10^9*" cells/ml"))

p4 <- make_band_plot(mean_df, "F_B_mean", "F_B_sd",
                     "Bound IFN", expression("Bound IFN"~10^{-5}~"pg/ml"),
                     show_legend = FALSE)+
  scale_y_continuous(
    breaks = c(0, 5e-05, 1e-04, 1.5e-04),
    labels = c("0", "0.5", "1", "1.5")
  )

p5 <- make_band_plot(mean_df, "F_U_mean", "F_U_sd",
                     "Unbound IFN",   expression("Unbound IFN pg/ml"),
                     show_legend = FALSE)

# ------------------------------------------------------------
# 4️⃣ Combine layout to match your previous panel
# ------------------------------------------------------------

combined_band <- p1 | (p2 / p4 | p3 / p5)

# ------------------------------------------------------------
# 5️⃣ Save
# ------------------------------------------------------------

print(combined_band)

# ------------------------------------------------------------
# 9️⃣ Save and print
# ------------------------------------------------------------
ggsave("/home/cdjuikem/overleaf-within-to-between-hosts-new/FIGS/band_cohort_truncated.pdf",
       plot = combined_band, width = 22, height = 15, units = "cm", dpi = 300)

