# ============================================================
# File: plot-Figure-C3-PCA.R
# Purpose: PCA biplot including ψ[max] (maximum Psi per patient)
# ============================================================

# 1. Load libraries
suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
})
suppressWarnings(suppressPackageStartupMessages(library(here)))

# ------------------------------------------------------------
# 2. Automatically find the latest cohort_truncated file
# ------------------------------------------------------------
cat("\n\n>>> Running plot-Figure-C3-PCA.R ...\n\n")

# Set project root automatically relative to the .git tracking directory
project_dir <- here()
if (basename(project_dir) == "CODE") {
  project_dir <- dirname(project_dir)
}

output_dir <- file.path(project_dir, "OUTPUT")

# 1a. Load simulation parameters (contains the max metrics)
param_files <- list.files(output_dir, pattern = "^cohort_sim_parameters_P.*\\.qs$", full.names = TRUE)
if (length(param_files) == 0) {
  stop("No cohort_sim_parameters file found in ", output_dir)
}
latest_param_file <- param_files[which.max(file.mtime(param_files))]
cat("Loading newest cohort parameters for MAX metrics:", basename(latest_param_file), "\n")
cohort_params <- qs_read(latest_param_file)

# 1b. Load corresponding patient status summary (classifier target: xi_h=75, xi_d=85)
status_files <- list.files(output_dir, pattern = "^cohort_status_P.*_xih_75_xid_85\\.qs$|^cohort_status_P.*_xid_85\\.qs$", full.names = TRUE)
if (length(status_files) == 0) {
  stop("No baseline cohort_status file found in ", output_dir)
}
latest_status_file <- status_files[which.max(file.mtime(status_files))]
cat("Loading newest cohort baseline statuses:", basename(latest_status_file), "\n")
cohort_status <- qs_read(latest_status_file)

# ------------------------------------------------------------
# 3. Summarize and merge by individual
# ------------------------------------------------------------
if("individual_id" %in% names(cohort_params) && !("ID" %in% names(cohort_params))) {
  cohort_params <- rename(cohort_params, ID = individual_id)
}

cohort_merged <- cohort_status %>%
  select(ID, status) %>%
  inner_join(cohort_params, by = "ID")

# 4. Prepare PCA data natively resolving Max Values
data_pca_full <- cohort_merged %>%
  filter(!is.na(status)) %>%
  filter(if_all(c(max_V, max_F_B, max_F_U, max_Psi), ~ !is.na(.x))) %>%
  mutate(status = factor(status, levels = c("Mild", "ICU", "Dead"))) %>%
  select(max_V, max_F_B, max_F_U, max_Psi, status) %>%
  as.data.frame()

cat("Individuals summarized for PCA:", nrow(data_pca_full), "\n")

quali_pos <- ncol(data_pca_full)

# 5. PCA computation
res.pca <- tryCatch(
  PCA(data_pca_full, quali.sup = quali_pos, graph = FALSE),
  error = function(e) { message("❌ PCA failed: ", e$message); NULL }
)

# 6. PCA biplot visualization
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
    "max_Psi" = "psi[max]",
    "max_V"   = "V[max]",
    "max_F_U" = "F[U]^max",
    "max_F_B" = "F[B]^max"
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
dir.create(file.path(project_dir, "FIGS"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(project_dir, "FIGS", "Figure-C3-PCA.pdf")
out_png <- file.path(project_dir, "FIGS", "Figure-C3-PCA.png")
  
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
  
  cat("✅ PCA plot saved to:\n  -", out_pdf, "\n  -", out_png, "\n")
  
} else {
  message("⚠️ PCA computation failed. No plot generated.")
}
