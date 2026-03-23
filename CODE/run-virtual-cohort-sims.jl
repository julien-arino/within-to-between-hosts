#!/usr/bin/env julia
# ============================================================
# File: run-virtual-cohort-sims.jl
# Description:
#   This script creates a virtual cohort of individuals and runs simulations
#   for each individual in parallel, saving the results in various formats.
#   The script saves results as Rds files for later exploitation in R.
#   Load required packages
# ============================================================
start_time_total = time()
println("\n\n>>> Starting run-virtual-cohort-sims.jl ...\n\n")
using Dates
using Distributed
using Serialization
using Printf  # Import Printf for @sprintf
using CSV     # Import the CSV package for reading/writing CSV files
using RCall   # Import RCall for interacting with R
using DataFrames # Import DataFrames

# Some of the output files can be quite large, so we include the option
# to save them on a large capacity disk (e.g., a NAS) if available.
OUTPUT_NAS = "/mnt/NAS-small-OUTPUT/within-to-between-hosts/"
# Use absolute paths so that RCall (which may have a different
# working directory context) can reliably find output files.
SCRIPT_DIR = @__DIR__  # Directory containing this script
OUTPUT_LOCAL = normpath(joinpath(SCRIPT_DIR, "..", "OUTPUT"))

# Select the output path based on knowledge of available disk space.
OUTPUT = OUTPUT_LOCAL

# Ensure the output directory exists (both Julia & R save operations rely on it).
mkpath(OUTPUT)

# Load external functions
include("functions-and-definitions.jl")

# Run parallel?
PARALLEL = true

# Save as jls?
SAVE_JLS = false

# Save as qs (using R's qs2/qs package)?
SAVE_QS = true

# Number of individuals in the virtual cohort
N = 1_000_000

# Set general parameters and (common) initial conditions
params = set_parameters()
IC = set_IC()

# Generate virtual cohort
individuals = generate_params_cohort(params, N)
individuals_idx = 1:N

# Run computation sequentially for all individuals
println("Starting computation for all $N individuals")
start_time = time()

if PARALLEL
    # Prepare parallel processing environment.

    # If there are pre-existing workers attached to this Julia session, remove
    # them first. Leftover workers from prior runs or other masters can try to
    # reconnect with a different cookie and trigger "Invalid connection
    # credentials sent by remote." Remove workers owned by this session so we
    # start with a clean slate.
    if nprocs() > 1
        println("Found existing workers (nprocs=$(nprocs())). Removing them to avoid stale/invalid connections...")
        try
            rmprocs(workers())
        catch e
            @warn "Failed to remove existing workers" exception = (e, catch_backtrace())
        end
    end

    # In case julia was not started with multiple processes, add some here. 
    # For large CPU counts, we add two thirds of the CPUs.
    # For smaller ones, we leave two free.
    if nprocs() < 2
        if Sys.CPU_THREADS >= 64
            addprocs(max(2, Int(round(Sys.CPU_THREADS * 2 / 3))))
        else
            addprocs(max(2, Sys.CPU_THREADS - 2))
        end
        println("Done setting up $(nprocs()) workers. Moving on to distribute variables...")
    end
    # Ensure all workers have the required functions and modules
    @everywhere using DifferentialEquations  # Import DifferentialEquations
    @everywhere using Serialization
    @everywhere using DataFrames
    @everywhere include("functions-and-definitions.jl")

    # Ensure all workers have the required variables
    @everywhere IC = $IC
    @everywhere params = $params
end

# Run computation
raw_results = if PARALLEL
    pmap(x -> run_one_individual(x, individuals, IC), individuals_idx)
else
    map(x -> run_one_individual(x, individuals, IC), individuals_idx)
end

println("Extracting maxima and R0 into the parameters table...")
individuals[!, :max_V] = [r[:maxima][:max_V] for r in raw_results]
individuals[!, :max_F_U] = [r[:maxima][:max_F_U] for r in raw_results]
individuals[!, :max_F_B] = [r[:maxima][:max_F_B] for r in raw_results]
individuals[!, :max_Psi] = [r[:maxima][:max_Psi] for r in raw_results]

individuals[!, :tau_max_V] = [r[:maxima][:tau_max_V] for r in raw_results]
individuals[!, :tau_max_F_U] = [r[:maxima][:tau_max_F_U] for r in raw_results]
individuals[!, :tau_max_F_B] = [r[:maxima][:tau_max_F_B] for r in raw_results]
individuals[!, :tau_max_Psi] = [r[:maxima][:tau_max_Psi] for r in raw_results]

individuals[!, :R0_within] = [r[:R0_within] for r in raw_results]

# Strip standard array for trajectory metadata
COHORT = [r[:vars] for r in raw_results]

# Print elapsed time
elapsed_computation = time() - start_time
println("Computation completed in $(elapsed_computation) seconds")

# Record date-time for unique file naming (capture current time for each run)
# Note: use UTC to avoid issues with compute nodes with time set incorrectly
date_time_start = Dates.format(now(UTC), "yyyymmdd-HHMMSS")

## Save the results as a JLS file
# Only save if SAVE_JLS is true
if SAVE_JLS
    println("Saving results as JLS")
    save_path_params = joinpath(OUTPUT, @sprintf("cohort_sim_parameters_P%07d_DT%s.jls", N, date_time_start))
    save_path_ic = joinpath(OUTPUT, @sprintf("cohort_sim_IC_P%07d_DT%s.jls", N, date_time_start))
    save_path_state = joinpath(OUTPUT, @sprintf("cohort_sim_state_P%07d_DT%s.jls", N, date_time_start))

    serialize(save_path_params, individuals)
    serialize(save_path_ic, IC)
    serialize(save_path_state, COHORT)
    println("Results saved to JLS files with DT$date_time_start")
end

## Save the results as a qs file (faster than Rds)
# Only save if SAVE_QS is true
# Beware:
# - R must be installed and available in the PATH, with libraries qs2 (preferably) or qs installed
# - RCall must be installed so julia can call R
# - This copies the SAVE variable to R, so if SAVE is large, this will be slow and likely to fail if RAM
#   is insufficient.
if SAVE_QS
    println("Saving results as QS (via RCall)")
    save_path_qs_params = joinpath(OUTPUT, @sprintf("cohort_sim_parameters_P%07d_DT%s.qs", N, date_time_start))
    save_path_qs_ic = joinpath(OUTPUT, @sprintf("cohort_sim_IC_P%07d_DT%s.qs", N, date_time_start))
    save_path_qs_state = joinpath(OUTPUT, @sprintf("cohort_sim_state_P%07d_DT%s.qs", N, date_time_start))

    @rput individuals
    @rput IC
    @rput COHORT
    @rput save_path_qs_params
    @rput save_path_qs_ic
    @rput save_path_qs_state

    R"""
    if(!requireNamespace("qs2", quietly=TRUE)) {
        warning("'qs2' package is not installed in R. Cannot save in qs format.")
    } else {
        qs2::qs_save(individuals, file = save_path_qs_params)
        qs2::qs_save(IC, file = save_path_qs_ic)
        qs2::qs_save(COHORT, file = save_path_qs_state)
    }
    """
    println("Results saved to QS files with DT$date_time_start")
end


# Close the cluster if parallel processing was used
if PARALLEL
    println("Shutting down workers...")
    redirect_stderr(devnull) do
        redirect_stdout(devnull) do
            try
                rmprocs(workers())  # Remove all worker processes
            catch
                # Ignore teardown warnings
            end
        end
    end
end

println("\nFiles saved in ", OUTPUT, ":")
if SAVE_JLS
    println("  - ", basename(save_path_params))
    println("  - ", basename(save_path_ic))
    println("  - ", basename(save_path_state))
end
if SAVE_QS
    println("  - ", basename(save_path_qs_params))
    println("  - ", basename(save_path_qs_ic))
    println("  - ", basename(save_path_qs_state))
end

elapsed_total = time() - start_time_total
mins = floor(Int, elapsed_total / 60)
secs = round(elapsed_total - mins * 60, digits=1)
println("\n>>> Finished run-virtual-cohort-sims.jl in $mins minutes and $(secs) seconds ✅\n")
