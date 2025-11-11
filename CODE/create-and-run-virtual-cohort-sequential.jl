using Dates
using Serialization
using Printf  # Import Printf for @sprintf
using DifferentialEquations  # Import the DifferentialEquations package
using CSV # Import the CSV package for reading/writing CSV files

# Load external functions
include("functions-all.jl")

OUTPUT_NAS = "/mnt/NAS-small-OUTPUT/within-to-between-hosts/"
OUTPUT_LOCAL = "OUTPUT/"
OUTPUT = OUTPUT_LOCAL

# Run parallel? (Disabled for sequential testing)
PARALLEL = false

# Number of patients in the virtual cohort
N = 1_000

# Set general parameters and (common) initial conditions
params = set_parameters()
IC = set_IC()

# Generate virtual cohort
patients = generate_params_patients(params, N)
patients_idx = 1:N

println("Starting computations")

# Record date-time at start to have common file name
date_time_start = Dates.format(now(UTC), "yyyyMMdd-HHmmss")

# Run computation sequentially for all patients
println("Starting computation for all $N patients")
start_time = time()

# COHORT = map(x -> run_one_patient(x, patients, IC), patients_idx)

COHORT = map(x -> run_one_patient(x, patients, IC; type_output = "select_variables"), patients_idx)


elapsed_time = time() - start_time
println("Computation completed in $(elapsed_time) seconds")

# Start preparing the save variable
SAVE = Dict()
SAVE[:parameters] = patients
# Add IC and results to save variable
SAVE[:IC] = IC
SAVE[:cohort] = COHORT

println("Saving results")
save_path = joinpath(OUTPUT_LOCAL, @sprintf("sim_P%07d_DT%s.jls", N, date_time_start))
serialize(save_path, SAVE)

println("Results saved to $save_path")

println("Saving results as CSV")
save_path_csv = joinpath(OUTPUT_LOCAL, @sprintf("sim_P%07d_DT%s.csv", N, date_time_start))

# Prepare a long-format DataFrame
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

# Save the long-format table as a CSV
CSV.write(save_path_csv, long_table)

println("Results saved to $save_path_csv")