library(dplyr)
library(qs)

# Load the non-truncated dataset
cohort_df <- qread("OUTPUT_clo/cohort_results.qs")

# 1️⃣ Summary per patient: Psi_max, t_max, computed status
summary_df <- cohort_df %>%
  group_by(patient_id) %>%
  summarise(
    Psi_max = max(Psi, na.rm = TRUE),
    t_max   = time[which.max(Psi)],
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      Psi_max < 75 ~ "Mild",
      Psi_max < 85 ~ "ICU",
      TRUE         ~ "Dead"
    )
  )

# 2️⃣ Merge and truncate trajectories
cohort_truncated <- cohort_df %>%
  inner_join(summary_df, by = "patient_id") %>%
  filter(Psi_max < 85 | time <= t_max)

# 3️⃣ Select Columns
cohort_truncated <- cohort_truncated %>%
  select(time, V, I, F_U, F_B,
         Psi,
         patient_id,
         status = status.y,   # <-- you asked to KEEP THIS EXACTLY
         Psi_max, t_max)

# 5️⃣ Save truncated data
qsave(cohort_truncated,
      "OUTPUT_clo/cohort_results_truncated.qs",
      preset = "balanced")

cat("✅ Saved: cohort_results_truncated.qs\n")


