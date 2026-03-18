library(qs2)
library(dplyr)
output_dir <- "/home/jarino/github/within-to-between-hosts/OUTPUT"
files <- list.files(output_dir, pattern = "^cohort_P.*\\.qs$", full.names = TRUE)
files <- files[!grepl("-times\\.qs$", files)]
latest_file <- files[which.max(file.mtime(files))]

cat("Loading latest cohort results:", basename(latest_file), "\n")
cohort_df <- qs_read(latest_file)

# Extract individual 3
ind3 <- cohort_df$cohort[[3]]

if (is.list(ind3) && "vars" %in% names(ind3)) {
    traj <- ind3$vars
} else {
    traj <- as.data.frame(ind3)
}

print("Head of trajectory:")
print(head(traj, 20))

print("Max V encountered:")
print(max(traj$V))

print("Did V ever exceed 4.0? (assuming xi_c = 4)")
print(any(traj$V >= 4.0))

# Print the points around the maximum
cat("Peak interval:\n")
max_idx <- argmax(traj$V)
print(traj[(max(1, max_idx - 5)):(min(nrow(traj), max_idx + 5)), c("time", "V", "Psi")])
