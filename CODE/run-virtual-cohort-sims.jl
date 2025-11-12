# This script creates a virtual cohort of individuals and runs simulations
# for each individual in parallel, saving the results in various formats.
# The script saves results as Rds files for later exploitation in R.

# Load required packages
using Dates
using Distributed
using Serialization
using Printf  # Import Printf for @sprintf
using CSV     # Import the CSV package for reading/writing CSV files
using RCall   # Import RCall for interacting with R

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
include("functions-all.jl")

# Run parallel?
PARALLEL = true

# Save as jls?
SAVE_JLS = false

# Save as CSV?
SAVE_CSV = false

## Type of output
# "maxima" = save only the maxima and their time of occurrence
# "select_variables" = select variables to save
# "all" = save all variables
type_output = "select_variables"

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
            @warn "Failed to remove existing workers" exception=(e, catch_backtrace())
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
    @everywhere include("functions-all.jl")

    # Ensure all workers have the required variables
    @everywhere IC = $IC
    @everywhere params = $params
end

# Run computation
COHORT = if PARALLEL
    pmap(x -> run_one_individual(x, individuals, IC; type_output = type_output), individuals_idx)
else
    map(x -> run_one_individual(x, individuals, IC; type_output = type_output), individuals_idx)
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
# Only save if SAVE_JLS is true
if SAVE_JLS
    println("Saving results as JLS")
    save_path = joinpath(OUTPUT, @sprintf("cohort_P%07d_DT%s_%s.jls", N, date_time_start, type_output))
    serialize(save_path, SAVE)
    println("Results saved to $save_path")
end

## Save the results as an Rds file
save_path_rds = joinpath(OUTPUT, @sprintf("cohort_P%07d_DT%s_%s.Rds", N, date_time_start, type_output))
@rput SAVE          # Send the Julia object to R
@rput save_path_rds # Send the absolute path to R (now guaranteed absolute)
R"""
saveRDS(SAVE, file = save_path_rds)
"""
println("Results saved to $save_path_rds")

## Save the results as a CSV file
# Only save if SAVE_CSV is true
if SAVE_CSV
    println("Saving results as CSV")
    save_path_csv = joinpath(OUTPUT, @sprintf("cohort_P%07d_DT%s_%s.csv", N, date_time_start, type_output))

    # Prepare the appropriate DataFrame based on the type of output
    if type_output == "select_variables"
        # Prepare a long-format DataFrame for selected variables
        long_table = DataFrame(sim_nb = Int[], time = Float64[], F_B = Float64[], F_U = Float64[], I = Float64[], V = Float64[])

        for (sim_nb, result) in enumerate(COHORT)
            times = result[:time]
            F_B = result[:F_B]
            F_U = result[:F_U]
            I = result[:I]
            V = result[:V]

            # Append rows for this simulation
            append!(long_table, DataFrame(sim_nb = fill(sim_nb, length(times)), time = times, F_B = F_B, F_U = F_U, I = I, V = V))
        end

    elseif type_output == "maxima"
        # Prepare a DataFrame for maxima
        long_table = DataFrame(sim_nb = Int[], variable = String[], value = Float64[])

        for (sim_nb, result) in enumerate(COHORT)
            for (var, value) in result
                append!(long_table, DataFrame(sim_nb = [sim_nb], variable = [string(var)], value = [value]))
            end
        end

    else
        error("Unknown type_output: $type_output")
    end

    # Save the DataFrame as a CSV
    CSV.write(save_path_csv, long_table)

    println("Results saved to $save_path_csv")
end