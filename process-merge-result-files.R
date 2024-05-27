# Merge files generated in multiple parts

library(dplyr)

source("functions_all.R")

NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
OUTPUT_LOCAL = "/home/jarino/OUTPUT_local/"

list_files = 
  data.frame(name = list.files(
    path = sprintf("%s/within-to-between-hosts",
                   OUTPUT_LOCAL),
    pattern = glob2rx("*.Rds"))
  )
list_files = list_files %>%
  filter(grepl("part", name)) %>%
  filter(grepl("sim", name))

# Load the files. Let's do the latest one for now..
for (f in list_files$name) {
  if (f == list_files$name[1]) {
    cohort = readRDS(sprintf("%s/%s", NAS, f))
  } else {
    tmp = readRDS(sprintf("%s/%s", NAS, f))
    cohort$cohort = append(cohort$cohort, tmp$cohort)
    cohort$parameters = rbind(cohort$parameters, tmp$parameters)
  }
}

saveRDS(cohort, sprintf("%s/%s_merged.Rds", NAS, list_files$name[1]))
