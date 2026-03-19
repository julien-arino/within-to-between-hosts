library(qs2)

# Find the original unfiltered cohort simulation files
output_dir <- file.path(getwd(), "OUTPUT")
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-all.R")) else source("functions-all.R")
}

# Match cohort_P*.qs but exclude cohort_times and cohort_censored
all_files <- list.files(output_dir, pattern = "^cohort_P.*\\.qs$", full.names = TRUE)
unfiltered_files <- all_files[!grepl("(_times_|_censored)", all_files)]

if (length(unfiltered_files) == 0) {
    stop("No unfiltered cohort_P*.qs files found in OUTPUT directory.")
}

latest_file <- unfiltered_files[which.max(file.mtime(unfiltered_files))]
cat("Loading newest UNFILTERED cohort results:", basename(latest_file), "\n")

# The original file is a nested list/dict containing :cohort and :parameters
save_data <- qs_read(latest_file, nthreads = N_QS_THREADS)
cohort <- save_data$cohort
parameters <- save_data$parameters

# Extract max_V from the nested :maxima elements
N <- length(cohort)
v_max_vals <- numeric(N)
for (i in seq_len(N)) {
  # Depending on Julia version, the dictionary keys may be strings or symbols
  if (!is.null(cohort[[i]]$maxima$max_V)) {
    v_max_vals[i] <- cohort[[i]]$maxima$max_V
  } else {
    v_max_vals[i] <- NA
  }
}

cat("\n==== Percentiles of Unfiltered V_max (0 to 100 by 10) ====\n")
quantiles_V <- quantile(v_max_vals, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
print(quantiles_V)
cat("============================================================\n")

cat("\n==== Percentiles of Unfiltered R0_within (0 to 100 by 10) ====\n")
quantiles_R0 <- quantile(parameters$R0_within, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
print(quantiles_R0)
cat("================================================================\n")
