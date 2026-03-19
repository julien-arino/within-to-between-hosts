#!/usr/bin/env Rscript
cat("\n\n>>> Running debug-read-write-big-data.R ...\n\n")
suppressWarnings(suppressPackageStartupMessages(library(arrow)))
suppressWarnings(suppressPackageStartupMessages(library(qs2)))
suppressWarnings(suppressPackageStartupMessages(library(tictoc)))
if(!exists("N_QS_THREADS")) {
  if(exists("project_dir")) source(file.path(project_dir, "CODE", "functions-and-definitions.R")) else source("functions-and-definitions.R")
}

cat("\n=======================================================\n")
cat("OPTION A: Legacy Pipeline (Direct List Load via qs2)\n")
cat("=======================================================\n")
cat("  >> Loading 1M DataFrames recursively from qs2...\n")
tic()
list_from_qs <- qs_read("shared_data_from_julia.qs", nthreads = N_QS_THREADS)
toc()

cat("\n=======================================================\n")
cat("OPTION B: Modern Pipeline (Arrow Load -> Reshape in Memory)\n")
cat("=======================================================\n")
cat("  >> Zero-copy parsing the 2.2GB Flattened DataFrame from Arrow...\n")
tic()
df_arrow <- read_feather("shared_data_from_julia.arrow")
toc()

cat("  >> Reshaping Arrow Flattened DataFrame back into a List of 1M DataFrames using split()...\n")
tic()
# Note: Since the ID column in Julia is named 'ID', we use $ID here.
list_from_arrow <- split(df_arrow, df_arrow$ID)
toc()

cat("\nBenchmark Complete!\n")

cat("\n\n>>> debug-read-write-big-data.R successfully finished running ✅\n\n")
