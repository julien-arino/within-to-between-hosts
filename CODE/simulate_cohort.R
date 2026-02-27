library(ggplot2)
library(dplyr)
library(furrr)
library(future)
library(patchwork)
library(purrr)
library(qs)
# ------------------------------------------------------------
# 1. Model
# ------------------------------------------------------------
source("functions-wb.R")  # contains simulate_model_ode() and compute_severity_ode()

# ------------------------------------------------------------
# 2. Setup parallel plan
# ------------------------------------------------------------
total_cores <- parallel::detectCores()
n_cores <- max(1L, floor(0.8 * total_cores))  # use ~80% of cores
plan(multisession, workers = n_cores)
opts <- furrr_options(seed = TRUE, scheduling = 1)

cat("Using", n_cores, "cores for parallel simulation\n")

# ------------------------------------------------------------
# 3. Parameters
# ------------------------------------------------------------
set.seed(123)
N <-300000  # number of patients

param_sd <- list(
  beta_V = 0.1994,
  d_V = 0.67,
  p = 158.65,
  p_FI = 1.8741,
  k_lin_f = 2.49,
  k_int_f = 2.54
)

patients <- data.frame(
  beta_V  = pmax(rnorm(N, default_params$beta_V, param_sd$beta_V), 0),
  d_V     = pmax(rnorm(N, default_params$d_V, param_sd$d_V), 0),
  p       = pmax(rnorm(N, default_params$p, param_sd$p), 0),
  p_FI    = pmax(rnorm(N, default_params$p_FI, param_sd$p_FI), 0),
  k_lin_f = pmax(rnorm(N, default_params$k_lin_f, param_sd$k_lin_f), 0),
  k_int_f = pmax(rnorm(N, default_params$k_int_f, param_sd$k_int_f), 0)
)

# ------------------------------------------------------------
# 4. Parallel simulation
# ------------------------------------------------------------
simulate_one_patient <- function(i, patient_df) {
  pars_i <- modifyList(default_params, as.list(patient_df[i, ]))
  out <- tryCatch(as.data.frame(simulate_model_ode(pars_i)), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  sev <- compute_severity_ode(out, S_max = pars_i$S_max)
  out$patient_id <- i
  out$status <- sev$status
  out$Psi <- sev$Psi
  out
}

t0 <- Sys.time()
results <- future_map(1:N, ~ simulate_one_patient(.x, patients),
                      .options = opts, .progress = TRUE)
res1 <- simulate_one_patient(1, patients)
object.size(res1) |> format(units = "auto")

t1 <- Sys.time()

cat("Simulation finished in", round(difftime(t1, t0, units = "mins"), 2), "minutes\n")


cat("Combining results...\n")
cohort_df <- list_rbind(results)  

#cat("Saving cohort results (gzip compression)...\n")
#saveRDS(cohort_df, "OUTPUT_clo/cohort_results.rds", compress = "gzip")

qsave(cohort_df, "OUTPUT_clo/cohort_results.qs", preset = "balanced")


#cat("✅ Data successfully saved.\n")


# ------------------------------------------------------------
# 5a. Compute mean trajectories for all key variables
# ------------------------------------------------------------
# mean_df <- cohort_df %>%
#   group_by(time) %>%
#   summarise(
#     Psi_mean = mean(Psi, na.rm = TRUE),
#     Psi_sd   = sd(Psi, na.rm = TRUE),
#     V_mean   = mean(V, na.rm = TRUE),
#     V_sd     = sd(V, na.rm = TRUE),
#     I_mean   = mean(I, na.rm = TRUE),
#     I_sd     = sd(I, na.rm = TRUE),
#     F_B_mean = mean(F_B, na.rm = TRUE),
#     F_B_sd   = sd(F_B, na.rm = TRUE),
#     F_U_mean = mean(F_U, na.rm = TRUE),
#     F_U_sd   = sd(F_U, na.rm = TRUE)
#   )
# 
# 
# x_center <- mean(range(cohort_df$time, na.rm = TRUE))
# 
# # 1️⃣ Lung tissue damage
# p1 <- ggplot(cohort_df, aes(x = time, y = Psi, group = patient_id)) +
#   geom_line(alpha = 0.05, color = "gray50") +
#   geom_line(data = mean_df, aes(x = time, y = Psi_mean), inherit.aes = FALSE,
#             color = "blue", linewidth = 1) +
#   annotate("rect", xmin = -Inf, xmax = Inf, ymin = 75, ymax = 85, fill = "gold", alpha = 0.25) +
#   annotate("rect", xmin = -Inf, xmax = Inf, ymin = 85, ymax = Inf, fill = "red", alpha = 0.25) +
#   annotate("text", x = x_center, y = 80, label = "ICU patients", size = 3.5) +
#   annotate("text", x = x_center, y = 92, label = "Dead patients", size = 3.5) +
#   labs(title = "Lung tissue damage", x = "Time (days)", y = "Percentage Ψ of tissue damage") +
#   theme_bw(base_size = 12)
# 
# # 2️⃣ Viral load
# p2 <- ggplot(cohort_df, aes(x = time, y = V, group = patient_id)) +
#   geom_line(alpha = 0.05, color = "gray60") +
#   geom_line(data = mean_df, aes(x = time, y = V_mean), inherit.aes = FALSE,
#             color = "blue", linewidth = 1) +
#   labs(title = "Viral load", x = "Time (days)", y = expression(log[10]*"(copies/ml)")) +
#   theme_bw(base_size = 12)
# 
# # 3️⃣ Infected cells
# p3 <- ggplot(cohort_df, aes(x = time, y = I, group = patient_id)) +
#   geom_line(alpha = 0.05, color = "gray60") +
#   geom_line(data = mean_df, aes(x = time, y = I_mean), inherit.aes = FALSE,
#             color = "blue", linewidth = 1) +
#   labs(title = "Infected cells", x = "Time (days)", y = expression(10^9*" cells/ml")) +
#   theme_bw(base_size = 12)
# 
# # 4️⃣ Bound IFN
# p4 <- ggplot(cohort_df, aes(x = time, y = F_B, group = patient_id)) +
#   geom_line(alpha = 0.05, color = "gray60") +
#   geom_line(data = mean_df, aes(x = time, y = F_B_mean), inherit.aes = FALSE,
#             color = "blue", linewidth = 1) +
#   labs(title = "Bound IFN", x = "Time (days)", y = "pg/ml") +
#   scale_y_continuous(
#     breaks = c(0, 2e-4, 4e-4, 6e-4),
#     labels = c("0", "2×10⁻⁴", "4×10⁻⁴", "6×10⁻⁴")
#   )+
#   theme_bw(base_size = 12)
# 
# # 5️⃣ Unbound IFN
# p5 <- ggplot(cohort_df, aes(x = time, y = F_U, group = patient_id)) +
#   geom_line(alpha = 0.05, color = "gray60") +
#   geom_line(data = mean_df, aes(x = time, y = F_U_mean), inherit.aes = FALSE,
#             color = "blue", linewidth = 1) +
#   labs(title = "Unbound IFN", x = "Time (days)", y = "pg/ml") +
#   theme_bw(base_size = 12)
# 
# # 6️⃣ Combine layout (like your example)
# right_side <- (p2 | p3) / (p4 | p5)
# combined <- p1 | right_side + plot_layout(widths = c(0.8, 2.2))
# 
# ggsave("FIGS_clo/simulate_cohort.png",
#        plot = combined, width = 20, height = 15, units = "cm", dpi = 300)
# 
# print(combined)
# 
# # ------------------------------------------------------------
# #6. Clean-up parallel workers
# # ------------------------------------------------------------
#  plan(sequential)
