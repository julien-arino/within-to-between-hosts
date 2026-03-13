# ============================================================
# File: plot-Figure-B3-PCA.R
# Purpose: PCA biplot including ψ[max] (maximum Psi per patient)
# ============================================================

# 1️⃣ Load libraries
suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
})

# ------------------------------------------------------------
# 2. Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
output_dir <- file.path(getwd(), "OUTPUT")
files <- list.files(output_dir, pattern = "cohort_censored_.*\\.qs$|cohort-censored_.*\\.qs$|cohort_results_truncated\\.qs$", full.names = TRUE)

if (length(files) == 0) {
  stop("No truncated/censored cohort file found in ", output_dir)
}

latest_file <- files[which.max(file.mtime(files))]
cat("Loading newest cohort results:", basename(latest_file), "\n")

cohort_df <- qs_read(latest_file)
cat("✅ Dataset loaded:", nrow(cohort_df), "rows\n")

# 3️⃣ Summarize by individual
# Mean for V, F_B, F_U — but MAX for Psi
patient_summary <- cohort_df %>%
  group_by(individual_id, status) %>%
  summarise(
    V       = mean(V, na.rm = TRUE),
    F_B     = mean(F_B, na.rm = TRUE),
    F_U     = mean(F_U, na.rm = TRUE),
    psi_max = max(Psi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(status))

# Keep complete cases only
patient_summary <- patient_summary %>%
  filter(if_all(c(V, F_B, F_U, psi_max), ~ !is.na(.x)))

cat("✅ Individuals summarized:", nrow(patient_summary), "\n")

# 4️⃣ Prepare PCA data
data_pca_full <- patient_summary %>%
  mutate(status = factor(status, levels = c("Mild", "ICU", "Dead"))) %>%
  select(V, F_B, F_U, psi_max, status) %>%
  as.data.frame()

quali_pos <- ncol(data_pca_full)

# 5️⃣ PCA computation
res.pca <- tryCatch(
  PCA(data_pca_full, quali.sup = quali_pos, graph = FALSE),
  error = function(e) { message("❌ PCA failed: ", e$message); NULL }
)

# 6️⃣ PCA biplot visualization
if (!is.null(res.pca)) {
  
  status_colors <- c(Mild = "dodgerblue4", ICU = "orange", Dead = "red")
  status_labels <- c(Mild = "Mild", ICU = "ICU", Dead = "Dead")
  
  # Base biplot (no variable labels yet)
  p_biplot <- fviz_pca_biplot(
    res.pca,
    geom.ind      = "point",
    habillage     = "status",
    addEllipses   = TRUE,
    ellipse.level = 0.68,
    label         = "none",    # we will add our own labels
    col.var       = "black") +
    scale_color_manual(values = status_colors,
                       labels = status_labels,
                       name   = "Status") +
    guides(fill = "none", shape = "none") +
    theme_minimal(base_size = 14) +
    labs(title = NULL) +
    theme(
      plot.title = element_blank(), 
      legend.position  = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey85")
    )
  
  # Extract variable coordinates
  var_coord <- as.data.frame(res.pca$var$coord[, 1:2])
  var_coord$varname <- rownames(var_coord)
  
  # Map R names -> math labels
  label_map <- c(
    "psi_max" = "psi[max]",
    "V"       = "V",
    "F_U"     = "F[U]",
    "F_B"     = "F[B]"
  )
  var_coord$label <- label_map[var_coord$varname]
  
  stretch_factor <- 10.5
  var_coord$x_lab <- var_coord$Dim.1 * stretch_factor
  var_coord$y_lab <- var_coord$Dim.2 * stretch_factor
  
  # Add labels at vector ends
  p_biplot <- p_biplot +
    geom_text(
      data = var_coord,
      aes(x = x_lab, y = y_lab, label = label),
      parse = TRUE,
      size = 4.5,
      fontface = "bold",
      hjust = 0.5,
      vjust = 0.5,
      inherit.aes = FALSE
    )
  
  print(p_biplot)
  
  # ------------------------------------------------------------
  # Save
  # ------------------------------------------------------------
  dir.create("FIGS", showWarnings = FALSE, recursive = TRUE)
  out_pdf <- "FIGS/plot_Figure_B3_PCA.pdf"
  out_png <- "FIGS/plot_Figure_B3_PCA.png"
  
  ggsave(
    filename = out_pdf,
    plot = p_biplot,
    width = 25, height = 15, units = "cm", dpi = 300
  )
  
  ggsave(
    filename = out_png,
    plot = p_biplot,
    width = 25, height = 15, units = "cm", dpi = 300
  )
  
  cat("✅ PCA biplot saved to", out_pdf, "and", out_png, "\n")
  
} else {
  message("⚠️ PCA computation failed. No plot generated.")
}
