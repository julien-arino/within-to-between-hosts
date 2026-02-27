# Merge files generated in multiple parts

library(dplyr)

source("functions-all.R")

NAS = "/home/jarino/OUTPUT_NAS_small/within-to-between-hosts"
OUTPUT_LOCAL = "/home/jarino/OUTPUT_local/within-to-between-hosts"

list_files = 
  data.frame(name = list.files(
    path = sprintf("%s",
                   NAS),
    pattern = glob2rx("*.Rds"))
  )
list_files = list_files %>%
  filter(grepl("part", name))
list_files$short_name = gsub(".Rds", "", list_files$name)
tmp = strsplit(list_files$short_name, "_")
list_files$type = unlist(lapply(tmp, function(x) x[1]))
list_files$nb_sims = unlist(lapply(tmp, function(x) x[2]))
list_files$date_time = unlist(lapply(tmp, function(x) x[3]))
tmp_tmp = strsplit(unlist(lapply(tmp, function(x) x[4])), "-")
list_files$part_curr = unlist(lapply(tmp_tmp, function(x) x[2]))
list_files$part_outof = unlist(lapply(tmp_tmp, function(x) x[4]))

# Check we have everything we are expecting

# Load and merge the files
for (dt in unique(list_files$date_time)) {
  writeLines(paste("Processing", dt))
  tmp_files = list_files %>%
    filter(date_time == dt)
  for (f in tmp_files$name) {
    writeLines(paste("  Processing", f))
    if (f == tmp_files$name[1]) {
      cohort = readRDS(sprintf("%s/%s", NAS, f))
    } else {
      tmp = readRDS(sprintf("%s/%s", NAS, f))
      cohort$cohort = append(cohort$cohort, tmp$cohort)
      cohort$parameters = rbind(cohort$parameters, tmp$parameters)
    }
  }
  saveRDS(cohort,
          file = sprintf("%s/%s_%s_%s_merged.Rds", 
                         OUTPUT_LOCAL, 
                         tmp_files$type[1],
                         tmp_files$nb_sims[1],
                         tmp_files$date_time[1]))
}


