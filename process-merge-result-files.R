source("functions_all.R")

NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"

list_files = c(
  "sim_10000invididuals_2023_07_15-21_41_14",
  "sim_10000invididuals_2023_07_17-21_10_34"#,
  #"sim_50000invididuals_2023_07_18-07_44_19"
)

# Load the files. Let's do the latest one for now..
for (f in list_files) {
  if (f == list_files[1]) {
    cohort = readRDS(sprintf("%s/%s.Rds", NAS, f))
  } else {
    tmp = readRDS(sprintf("%s/%s.Rds", NAS, f))
    cohort$cohort = append(cohort$cohort, tmp$cohort)
    cohort$parameters = rbind(cohort$parameters, tmp$parameters)
  }
}

saveRDS(cohort, sprintf("%s/%s_merged.Rds", NAS, list_files[1]))
