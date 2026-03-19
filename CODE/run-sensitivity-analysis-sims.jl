#!/usr/bin/env julia
# This script creates a set of parameter values for individuals, varying
# all parameters according to a Sobol sequence, and runs simulations
# for each individual in parallel, saving the results in various formats.
# The script calls on R to create the cohort parameters and saves results as 
# Rds or qs files for later exploitation in R.

# Load required packages
using Dates
using Distributed
using Serialization
using Printf  # Import Printf for @sprintf
using CSV # Import the CSV package for reading/writing CSV files
using RCall  # Import RCall for interacting with R

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

# Save as CSV?
SAVE_CSV = false

# Save as Rds?
SAVE_RDS = false

# Save as qs?
SAVE_QS = true

## Type of output
# "maxima" = save only the maxima and their time of occurrence
# "select_variables" = select variables to save
# "all" = save all variables
type_output = "maxima"

# Number of individuals in the virtual cohort
N = 1_000_000

# Build the absolute path to the R scripts to use FOR WITH_V0 PIPELINE
r_functions_path = joinpath(SCRIPT_DIR, "functions-and-definitions.R")
r_script_path = joinpath(SCRIPT_DIR, "create-sample-for-sensitivity.R")

## Generate the sample in R
# Send the Julia objects to R
@rput N r_functions_path r_script_path
# Source the R script to generate the sample
R"""
source(r_functions_path)
source(r_script_path)
"""
@rget individuals  # Get the sample from R

# Print the dimension of pars.sobol
println("individuals size (with V0): ", size(individuals))

# Establish standard IC vector
IC = set_IC()

# Set the individual indices
individuals_idx = 1:N

# Run computation sequentially for all individuals
println("Starting computation for all $N individuals")
start_time = time()

if PARALLEL
    # Prepare parallel processing environment. 
    # Clean up any previously attached workers to avoid mismatched connection cookies
    if nprocs() > 1
        println("Cleaning up pre-existing workers (nprocs=$(nprocs())) before spawning new ones...")
        try
            rmprocs(workers())
        catch e
            @warn "Failed to remove existing workers in run-sensitivity-analysis-sims" exception = (e, catch_backtrace())
        end
    end

    # Give the system a brief moment to clear the ports before spinning up new workers
    sleep(1.0)

    # In case julia was not started with multiple processes, add some here. 
    # For large CPU counts, we add two thirds of the CPUs.
    # For smaller ones, we leave two free.
    num_to_add = if Sys.CPU_THREADS >= 64
        max(2, Int(round(Sys.CPU_THREADS * 2 / 3)))
    else
        max(2, Sys.CPU_THREADS - 2)
    end

    println("Spawning $num_to_add new workers...")
    # Explicitly pass the project path to ensure credential domains align
    addprocs(num_to_add, exeflags="--project=@.")

    # Ensure all workers have the required functions and modules
    @everywhere using DifferentialEquations  # Import the DifferentialEquations package
    @everywhere using Serialization
    @everywhere include("functions-and-definitions.jl")
end

# INTERCEPT function to extract dynamic V0 from the individual dataframe before running the ODE wrapper
@everywhere function run_individual_with_dynamic_V0(idx, individuals_df, base_IC, type_output)
    # create a locally modified mutable copy of the IC array for this individual
    local_IC = copy(base_IC)
    # The sampled V0 is attached directly to the individuals dataframe now
    local_IC[1] = individuals_df[idx, :V0]

    return run_one_individual(idx, individuals_df, local_IC)
end

# Run computation via our intercept function
COHORT = if PARALLEL
    pmap(x -> run_individual_with_dynamic_V0(x, individuals, IC, type_output), individuals_idx)
else
    map(x -> run_individual_with_dynamic_V0(x, individuals, IC, type_output), individuals_idx)
end

# Close the cluster if parallel processing was used
if PARALLEL
    println("Shutting down workers...")
    rmprocs(workers())  # Remove all worker processes
end

# Print elapsed time
elapsed_time = time() - start_time
println("Computation completed in $(elapsed_time) seconds")

# Record date-time for unique file naming (capture current time for each run)
# Note: use UTC to avoid issues with compute nodes with time set incorrectly
date_time_start = Dates.format(now(UTC), "yyyymmdd-HHMMSS")

# Compute R0 for each individual and add as a column
println("Computing R0 for each individual in the cohort...")
compute_R0_cohort!(individuals)
println("R0 computation completed")

# Start preparing the save variable
SAVE = Dict()
SAVE[:parameters] = individuals
# Add IC and results to save variable
SAVE[:IC] = IC
SAVE[:cohort] = COHORT

## Save the results as a JLS file
if SAVE_JLS
    println("Saving results as JLS")
    save_path = joinpath(OUTPUT, @sprintf("sensitivity_P%07d_DT%s_%s.jls", N, date_time_start, type_output))
    serialize(save_path, SAVE)
    println("Results saved to $save_path")
end

## Save the results as an Rds file
if SAVE_RDS
    println("Saving results as Rds (via RCall)")
    save_path_rds = joinpath(OUTPUT, @sprintf("sensitivity_P%07d_DT%s_%s.Rds", N, date_time_start, type_output))
    @rput SAVE  # Send the Julia object to R
    @rput save_path_rds  # Send the absolute path to R
    R"""
    saveRDS(SAVE, file = save_path_rds)
    """
    println("Results saved to $save_path_rds")
end

## Save the results as a qs file (faster than Rds)
if SAVE_QS
    println("Saving results as QS (via RCall)")
    save_path_qs = joinpath(OUTPUT, @sprintf("sensitivity_P%07d_DT%s_%s.qs", N, date_time_start, type_output))
    @rput SAVE          # Send the Julia object to R
    @rput save_path_qs  # Send the absolute path to R
    R"""
    if(!requireNamespace("qs2", quietly=TRUE)) {
        warning("'qs2' package is not installed in R. Cannot save in qs format.")
    } else {
        qs2::qs_save(SAVE, file = save_path_qs)
    }
    """
    println("Results saved to $save_path_qs")
end

## Save the results as a CSV file
if SAVE_CSV
    println("Saving results as CSV")
    save_path_csv = joinpath(OUTPUT, @sprintf("sensitivity_P%07d_DT%s_%s.csv", N, date_time_start, type_output))

    # Prepare the appropriate DataFrame based on the type of output
    if type_output == "select_variables"
        # Prepare a long-format DataFrame for selected variables
        long_table = DataFrame(sim_nb=Int[], time=Float64[], F_B=Float64[], F_U=Float64[], I=Float64[], V=Float64[])

        for (sim_nb, result) in enumerate(COHORT)
            times = result[:time]
            F_B = result[:F_B]
            F_U = result[:F_U]
            I = result[:I]
            V = result[:V]

            # Append rows for this simulation
            append!(long_table, DataFrame(sim_nb=fill(sim_nb, length(times)), time=times, F_B=F_B, F_U=F_U, I=I, V=V))
        end

    elseif type_output == "maxima"
        # Prepare a DataFrame for maxima
        long_table = DataFrame(sim_nb=Int[], variable=String[], value=Float64[])

        for (sim_nb, result) in enumerate(COHORT)
            for (var, value) in result
                append!(long_table, DataFrame(sim_nb=[sim_nb], variable=[string(var)], value=[value]))
            end
        end

    else
        error("Unknown type_output: $type_output")
    end

    # Save the DataFrame as a CSV
    CSV.write(save_path_csv, long_table)

    println("Results saved to $save_path_csv")
end
