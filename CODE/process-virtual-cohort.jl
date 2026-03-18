## process-virtual-cohort.jl
# Processes output files from run-virtual-cohort-sims.jl and computes the
# various tau times defined in the manuscript.

using RCall
using DataFrames
using Glob
using Statistics
using Printf

# Default inputs
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = normpath(joinpath(SCRIPT_DIR, "..", "OUTPUT"))

# Find the latest results file (either qs or Rds)
result_files = vcat(glob("cohort_P*.qs", OUTPUT_DIR), glob("cohort_P*.Rds", OUTPUT_DIR))
result_files = filter(f -> !occursin("_times", f) && !occursin("_censored", f), result_files)

if isempty(result_files)
    error("No .qs or .Rds result files found in $OUTPUT_DIR")
end

# Get the latest file by modification time
latest_file = sort(result_files, by=mtime, rev=true)[1]
println("Loading latest result file: $latest_file")

# Load the file via RCall into a Julia object based on extension
if endswith(latest_file, ".qs")
    R"""
    if(!requireNamespace("qs2", quietly=TRUE)) {
        stop("'qs2' package is not installed in R.")
    } else {
        save_data <- qs2::qs_read($latest_file)
    }
    """
else
    R"""
    save_data <- readRDS($latest_file)
    """
end

SAVE = rcopy(R"save_data")

# Extract the cohort simulation results and parameters
cohort = SAVE[:cohort]
params_df = SAVE[:parameters]
N = length(cohort)

println("Loaded results for cohort of size: $N")

# Thresholds defined in the manuscript
# Disease severity & Death
const xi_h = 75.0  # Hospitalization threshold (%)
const xi_d = 85.0  # Death threshold (%)

# Infectiousness
const xi_c = 9.0    # Start of infectious period (assumed base on log10 init condition 4.5 * 2)
const xi_r = 1.0    # End of infectious period (assumed threshold)

# Arrays to store computed tau values
tau_Psi_max = fill(NaN, N)
tau_h_start = fill(NaN, N)
tau_h_end = fill(NaN, N)
tau_d = fill(NaN, N)
tau_c = fill(NaN, N)
tau_r = fill(NaN, N)
tau_V_max = fill(NaN, N)
R0_P2P = fill(0.0, N)
Psi_max = fill(NaN, N)

println("Computing tau values for each individual...")

# Constants for transmission beta function
const alpha_i = 16.422
const k_i = 7.49

function beta_i(V::Float64)
    if V <= 0.0
        return 0.0
    end
    # Using log-transform to handle large powers safely
    # V^alpha / (V^alpha + k^alpha) = 1 / (1 + (k/V)^alpha)
    return 1.0 / (1.0 + (k_i / V)^alpha_i)
end

function compute_tau_for_individual_idx(i, data, S_P_0=2000.0)
    time_pts = data[:time]
    Psi = data[:Psi]
    V = data[:V]

    # Preallocate returns: [tau_Psi_max, tau_h_start, tau_h_end, tau_d, tau_c, tau_r, tau_V_max, R0_P2P, Psi_max_val]
    res = fill(NaN, 9)
    res[8] = 0.0

    # ------------------- #
    # Eq 2.2 / Eq 5.2     #
    # ------------------- #

    Psi_max_val = maximum(Psi)
    Psi_max_idx = argmax(Psi)
    res[1] = time_pts[Psi_max_idx] # tau_Psi_max
    res[9] = Psi_max_val           # actual max value

    # ------------------- #
    # Eq 5.3 (Hospital)   #
    # ------------------- #
    h_indices = findall(>=(xi_h), Psi)
    if !isempty(h_indices)
        res[2] = time_pts[h_indices[1]]   # tau_h_start
        res[3] = time_pts[h_indices[end]] # tau_h_end
    end

    # ------------------- #
    # Eq 5.4 (Death)      #
    # ------------------- #
    if Psi_max_val >= xi_d
        d_idx = findfirst(x -> (x >= xi_d && x ≈ Psi_max_val), Psi)
        if d_idx !== nothing
            res[4] = time_pts[d_idx] # tau_d
            res[3] = NaN             # tau_h_end = NaN if dead
        end
    end

    # ------------------- #
    # Eq 5.5 (V_max)      #
    # ------------------- #
    V_max_val = maximum(V)
    V_max_idx = argmax(V)
    res[7] = time_pts[V_max_idx] # tau_V_max

    # ------------------- #
    # Eq 5.7 (Start Inf)  #
    # ------------------- #
    c_idx = findfirst(>=(xi_c), V)
    if c_idx !== nothing
        res[5] = time_pts[c_idx] # tau_c
    end

    # ------------------- #
    # Eq 5.8 (End Inf)    #
    # ------------------- #
    if isnan(res[4])
        post_peak_indices = (V_max_idx+1):length(V)
        r_offset = findfirst(<=(xi_r), @view V[post_peak_indices])
        if r_offset !== nothing
            r_idx = post_peak_indices[1] + r_offset - 1
            res[6] = time_pts[r_idx] # tau_r
        end
    end

    # ------------------- #
    # Eq R0_P             #
    # ------------------- #
    tau_end = min(isnan(res[6]) ? Inf : res[6], isnan(res[4]) ? Inf : res[4])

    if isinf(tau_end)
        end_idx = length(time_pts)
    else
        end_idx = findfirst(>=(tau_end), time_pts)
        if end_idx === nothing
            end_idx = length(time_pts)
        end
    end

    # Trapezoidal numerical integration up to tau_end
    integ_val = 0.0
    for j in 1:(end_idx-1)
        dt = time_pts[j+1] - time_pts[j]
        integ_val += 0.5 * dt * (beta_i(V[j]) + beta_i(V[j+1]))
    end

    res[8] = S_P_0 * integ_val

    return res
end

using Distributed

if N >= 50000
    println("N >= 50,000. Configuring distributed processing...")

    # Remove existing workers if any to prevent stale connections
    if nprocs() > 1
        println("Found existing workers (nprocs=$(nprocs())). Removing them...")
        try
            rmprocs(workers())
        catch e
            @warn "Failed to remove existing workers" exception = (e, catch_backtrace())
        end
    end

    # Add workers based on available CPUs
    if Sys.CPU_THREADS >= 64
        addprocs(max(2, Int(round(Sys.CPU_THREADS * 2 / 3))))
    else
        addprocs(max(2, Sys.CPU_THREADS - 2))
    end
    println("Done setting up $(nprocs()) workers. Moving on to distribute parameters...")

    @everywhere begin
        using OrderedCollections

        # Distribute the computation function and thresholds to all workers
        const xi_h_w = $xi_h
        const xi_d_w = $xi_d
        const xi_c_w = $xi_c
        const xi_r_w = $xi_r

        const alpha_i_w = 16.422
        const k_i_w = 7.49

        function beta_i_w(V::Float64)
            if V <= 0.0
                return 0.0
            end
            return 1.0 / (1.0 + (k_i_w / V)^alpha_i_w)
        end

        function compute_tau_for_individual_idx(data, S_P_0=2000.0)
            time_pts = data[:time]
            Psi = data[:Psi]
            V = data[:V]

            res = fill(NaN, 9)
            res[8] = 0.0

            Psi_max_val = maximum(Psi)
            Psi_max_idx = argmax(Psi)
            res[1] = time_pts[Psi_max_idx]
            res[9] = Psi_max_val

            h_indices = findall(>=(xi_h_w), Psi)
            if !isempty(h_indices)
                res[2] = time_pts[h_indices[1]]
                res[3] = time_pts[h_indices[end]]
            end

            if Psi_max_val >= xi_d_w
                d_idx = findfirst(x -> (x >= xi_d_w && x ≈ Psi_max_val), Psi)
                if d_idx !== nothing
                    res[4] = time_pts[d_idx]
                    res[3] = NaN
                end
            end

            V_max_val = maximum(V)
            V_max_idx = argmax(V)
            res[7] = time_pts[V_max_idx]

            c_idx = findfirst(>=(xi_c_w), V)
            if c_idx !== nothing
                res[5] = time_pts[c_idx]
            end

            if isnan(res[4])
                post_peak_indices = (V_max_idx+1):length(V)
                r_offset = findfirst(<=(xi_r_w), @view V[post_peak_indices])
                if r_offset !== nothing
                    r_idx = post_peak_indices[1] + r_offset - 1
                    res[6] = time_pts[r_idx]
                end
            end

            tau_end = min(isnan(res[6]) ? Inf : res[6], isnan(res[4]) ? Inf : res[4])

            if isinf(tau_end)
                end_idx = length(time_pts)
            else
                end_idx = findfirst(>=(tau_end), time_pts)
                if end_idx === nothing
                    end_idx = length(time_pts)
                end
            end

            integ_val = 0.0
            for j in 1:(end_idx-1)
                dt = time_pts[j+1] - time_pts[j]
                integ_val += 0.5 * dt * (beta_i_w(V[j]) + beta_i_w(V[j+1]))
            end

            res[8] = S_P_0 * integ_val

            return res
        end
    end

    println("Processing $(N) individuals in parallel...")

    # We build an array of pairs/tuples containing just the pieces needed 
    # for each worker to run smoothly without sending the entire `cohort` memory.
    inputs = [(cohort[i][:vars]) for i in 1:N]

    # Run the computation in parallel
    results = pmap(inputs) do data
        compute_tau_for_individual_idx(data)
    end

    # Unpack the results
    for i in 1:N
        res = results[i]
        tau_Psi_max[i] = res[1]
        tau_h_start[i] = res[2]
        tau_h_end[i] = res[3]
        tau_d[i] = res[4]
        tau_c[i] = res[5]
        tau_r[i] = res[6]
        tau_V_max[i] = res[7]
        R0_P2P[i] = res[8]
        Psi_max[i] = res[9]
    end

    # Shut down workers
    println("Shutting down workers...")
    rmprocs(workers())

else
    println("Running sequentially...")
    for i in 1:N
        data = cohort[i][:vars]

        res = compute_tau_for_individual_idx(i, data)

        tau_Psi_max[i] = res[1]
        tau_h_start[i] = res[2]
        tau_h_end[i] = res[3]
        tau_d[i] = res[4]
        tau_c[i] = res[5]
        tau_r[i] = res[6]
        tau_V_max[i] = res[7]
        R0_P2P[i] = res[8]
        Psi_max[i] = res[9]
    end
end

# Create a much cleaner output DataFrame retaining only exactly what the user requested
out_df = DataFrame()
out_df[!, :ID] = params_df[!, :ID]
out_df[!, :R0_within] = params_df[!, :R0_within]
out_df[!, :R0_P2P] = R0_P2P
out_df[!, :Psi_max] = Psi_max

# Note: The manuscript specifies xi_h=75% and xi_d=85% for severity, and V>=4.5 for transmission
out_df[!, :tau_Psi_max] = tau_Psi_max
out_df[!, :tau_h_start] = tau_h_start
out_df[!, :tau_h_end] = tau_h_end
out_df[!, :tau_d] = tau_d
out_df[!, :tau_c] = tau_c
out_df[!, :tau_r] = tau_r
out_df[!, :tau_V_max] = tau_V_max

# Extract the original maximums natively tracked during the ODE solver run back in `run-virtual-cohort-sims.jl`
out_df[!, :max_V] = [cohort[i][:maxima][:max_V] for i in 1:N]
out_df[!, :max_F_U] = [cohort[i][:maxima][:max_F_U] for i in 1:N]
out_df[!, :max_F_B] = [cohort[i][:maxima][:max_F_B] for i in 1:N]
out_df[!, :tau_max_V] = [cohort[i][:maxima][:tau_max_V] for i in 1:N]
out_df[!, :tau_max_F_U] = [cohort[i][:maxima][:tau_max_F_U] for i in 1:N]
out_df[!, :tau_max_F_B] = [cohort[i][:maxima][:tau_max_F_B] for i in 1:N]

# Display a quick summary of the results
println("\n--- Summary of Computed Taus ---")
survived = isnan.(tau_d)
died = .!survived
hosp = .!isnan.(tau_h_start)
inf_start = .!isnan.(tau_c)

println("Total Individuals: ", N)
println(@sprintf("Survived: %d (%.2f%%)", sum(survived), 100 * sum(survived) / N))
println(@sprintf("Died (tau_d present): %d (%.2f%%)", sum(died), 100 * sum(died) / N))
println(@sprintf("Hospitalized (tau_h_start present): %d (%.2f%%)", sum(hosp), 100 * sum(hosp) / N))
println(@sprintf("Started Infectious Period (tau_c present): %d (%.2f%%)", sum(inf_start), 100 * sum(inf_start) / N))
println("Mean time to death:     ", round(mean(filter(!isnan, tau_d)); digits=2), " days")
println("Mean time to hospital:  ", round(mean(filter(!isnan, tau_h_start)); digits=2), " days")
println("Mean time to infectious:", round(mean(filter(!isnan, tau_c)); digits=2), " days")

println("\nSaving updated DataFrame out as process-virtual-cohort-results.csv ...")
r_csv_lib = RCall.rcopy(R"library(readr)") # Quick hack if CSV.jl isn't around, but let's just use CSV package
# We assume the user has DataFrames + CSV installed. We will use CSV.jl
using CSV
out_path = joinpath(OUTPUT_DIR, "process-virtual-cohort-results.csv")
CSV.write(out_path, out_df)
println("Saved successfully to $out_path")

# Save as qs file with "cohort_times_" prefix
base_name = splitext(basename(latest_file))[1]
new_base = replace(base_name, "cohort_" => "cohort_times_")
out_path_qs = joinpath(OUTPUT_DIR, new_base * ".qs")
println("\nSaving updated DataFrame out as $out_path_qs ...")
@rput out_df
@rput out_path_qs
R"""
if(!requireNamespace("qs2", quietly=TRUE)) {
    warning("'qs2' package is not installed in R. Cannot save in qs format.")
} else {
    qs2::qs_save(out_df, file = out_path_qs)
}
"""
println("Saved successfully to $out_path_qs")

# ------------------------------------------------------------
# Format and Save Truncated Cohort Data
# ------------------------------------------------------------
println("\nConstructing and saving the right-censored cohort dataset (truncated at death)...")

# We want roughly these columns: time, V, I, F_U, F_B, Psi, individual_id, status, Psi_max, t_max
n_total_rows = sum(length(cohort[i][:vars][:time]) for i in 1:N)
# Preallocate slightly larger than we'll need because we'll crop some out
all_time = Vector{Float64}(undef, n_total_rows)
all_V = Vector{Float64}(undef, n_total_rows)
all_I = Vector{Float64}(undef, n_total_rows)
all_F_U = Vector{Float64}(undef, n_total_rows)
all_F_B = Vector{Float64}(undef, n_total_rows)
all_Psi = Vector{Float64}(undef, n_total_rows)

all_individual_id = Vector{Int}(undef, n_total_rows)
all_status = Vector{String}(undef, n_total_rows)
all_Psi_max = Vector{Float64}(undef, n_total_rows)
all_t_max = Vector{Float64}(undef, n_total_rows)

row_idx = 1
for i in 1:N
    global row_idx
    data = cohort[i][:vars]
    pts = length(data[:time])

    # Calculate patient status natively
    c_Psi_max = maximum(data[:Psi])
    c_t_max = data[:time][argmax(data[:Psi])]

    c_status = "Mild"
    if c_Psi_max >= xi_d
        c_status = "Dead"
    elseif c_Psi_max >= xi_h
        c_status = "ICU"
    end

    c_tau_d = tau_d[i]

    for j in 1:pts
        # Truncate strictly at time of death
        if !isnan(c_tau_d) && data[:time][j] > c_tau_d
            continue
        end

        all_time[row_idx] = data[:time][j]
        all_V[row_idx] = data[:V][j]
        all_I[row_idx] = data[:I][j]
        all_F_U[row_idx] = data[:F_U][j]
        all_F_B[row_idx] = data[:F_B][j]
        all_Psi[row_idx] = data[:Psi][j]

        all_individual_id[row_idx] = i
        all_status[row_idx] = c_status
        all_Psi_max[row_idx] = c_Psi_max
        all_t_max[row_idx] = c_t_max

        row_idx += 1
    end
end

# Crop down to exact length filled
idx_end = row_idx - 1
cohort_censored_df = DataFrame(
    time=all_time[1:idx_end],
    V=all_V[1:idx_end],
    I=all_I[1:idx_end],
    F_U=all_F_U[1:idx_end],
    F_B=all_F_B[1:idx_end],
    Psi=all_Psi[1:idx_end],
    individual_id=all_individual_id[1:idx_end],
    status=all_status[1:idx_end],
    Psi_max=all_Psi_max[1:idx_end],
    t_max=all_t_max[1:idx_end]
)

out_truncated_qs = joinpath(OUTPUT_DIR, replace(base_name, "cohort_" => "cohort_censored_") * ".qs")

@rput cohort_censored_df
@rput out_truncated_qs
R"""
if(!requireNamespace("qs2", quietly=TRUE)) {
    warning("'qs2' package is not installed in R. Cannot save censored data.")
} else {
    qs2::qs_save(cohort_censored_df, file = out_truncated_qs)
}
"""
println("Saved right-censored cohort successfully to $out_truncated_qs")
