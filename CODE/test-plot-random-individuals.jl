# plot-random-individuals.jl

using RCall
using Glob
using Random
using Plots

# Default inputs
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = normpath(joinpath(SCRIPT_DIR, "..", "OUTPUT"))

# Find the latest Rds file 
rds_files = glob("*.Rds", OUTPUT_DIR)
if isempty(rds_files)
    error("No .Rds result files found in $OUTPUT_DIR")
end

# Get the latest Rds file by modification time
latest_rds = sort(rds_files, by=mtime, rev=true)[1]
println("Loading latest result file: $latest_rds")

# Load the Rds file via RCall into a Julia object
R"""
save_data <- readRDS($latest_rds)
"""
SAVE = rcopy(R"save_data")

# Extract the cohort simulation results
cohort = SAVE[:cohort]
N = length(cohort)

println("Loaded results for cohort of size: $N")

num_to_plot = min(500, N)
random_indices = randperm(N)[1:num_to_plot]

println("Plotting $num_to_plot random individuals...")

# Initialize plots
plt_V = plot(title="Viral Load (V) for $num_to_plot Random Individuals",
    xlabel="Time (days)", ylabel="V (log10 copies/ml)",
    legend=false, grid=true)

plt_Psi = plot(title="Tissue Damage (Psi) for $num_to_plot Random Individuals",
    xlabel="Time (days)", ylabel="Psi (%)",
    legend=false, grid=true)

for i in random_indices
    # Handle structure difference depending on type_output used
    if haskey(cohort[i], :vars)
        data = cohort[i][:vars]
    else
        data = cohort[i]
    end

    time_pts = data[:time]
    V = data[:V]

    # Check if Psi is natively saved or if we are reading an old file where it needs to be computed
    if haskey(data, :Psi)
        Psi = data[:Psi]
    else
        S_max = SAVE[:parameters][i, :S_max]
        S = data[:S]
        R = data[:R]
        Psi = 100 .* (S_max .- (S .+ R)) ./ S_max
    end

    # Plot
    plot!(plt_V, time_pts, V, color=:blue, alpha=0.3)
    plot!(plt_Psi, time_pts, Psi, color=:red, alpha=0.3)
end

# Save the plots
out_v = joinpath(OUTPUT_DIR, "plot_V_random100.png")
out_psi = joinpath(OUTPUT_DIR, "plot_Psi_random100.png")

savefig(plt_V, out_v)
savefig(plt_Psi, out_psi)

println("Saved V plot to: $out_v")
println("Saved Psi plot to: $out_psi")
