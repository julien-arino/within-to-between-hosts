# ==============================================================================
# File: process-right-censor-dead-individuals.R
# Description:
#   This script processes the full virtual cohort simulation results.
#   It identifies individuals based on their maximum lung tissue damage (Psi_max):
#     - Mild: Psi_max < 75%
#     - ICU: 75% <= Psi_max < 85%
#     - Dead: Psi_max >= 85%
#
#   For individuals classified as "Dead", their simulation trajectories are
#   truncated exactly at the time they hit their maximum damage score (t_max),
#   removing any biologically irrelevant data points simulated after death.
#   The cleaned and truncated data is saved for further analysis.
# ==============================================================================

library(dplyr)
library(qs2)

# Load the non-truncated dataset
cohort_df <- qs_read("OUTPUT/cohort_results.qs")

# 1️⃣ Summary per individual: Psi_max, t_max, computed status
summary_df <- cohort_df %>%
  group_by(individual_id) %>%
  summarise(
    Psi_max = max(Psi, na.rm = TRUE),
    t_max   = time[which.max(Psi)],
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      Psi_max < 75 ~ "Mild",
      Psi_max < 85 ~ "ICU",
      TRUE ~ "Dead"
    )
  )

# 2️⃣ Merge and truncate trajectories
cohort_truncated <- cohort_df %>%
  inner_join(summary_df, by = "individual_id") %>%
  filter(Psi_max < 85 | time <= t_max)

# 3️⃣ Select Columns
cohort_truncated <- cohort_truncated %>%
  select(time, V, I, F_U, F_B,
    Psi,
    individual_id,
    status = status.y, # <-- you asked to KEEP THIS EXACTLY
    Psi_max, t_max
  )

# 5️⃣ Save truncated data
qs_save(cohort_truncated,
  "OUTPUT/cohort_results_truncated.qs"
)

cat("✅ Saved: cohort_results_truncated.qs\n")
