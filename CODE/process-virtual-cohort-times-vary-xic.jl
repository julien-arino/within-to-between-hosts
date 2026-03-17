## process-virtual-cohort-times-vary-xic.jl
# Processes output files from run-virtual-cohort-sims.jl and computes the
# various tau times defined in the manuscript across a range of infectiousness
# onset thresholds (xi_c), holding hospitalisation (xi_h) and death (xi_d) fixed.

using RCall
using DataFrames
using Glob
using Statistics
using Printf
using Distributed

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

# Load the file via RCall
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

# Fixed thresholds for hospitalisation, deaths and end of infectious period
const xi_h = 75.0
const xi_d = 85.0
const xi_r = 1.0    # End of infectious period

# Grid of infectious start thresholds to iterate over (down to the V(0)=1.0 initial condition)
xi_c_vals = collect(1.0:1.0:10.0)

# Constants for transmission beta function
const alpha_i = 16.422
const k_i = 7.49

function beta_i(V::Float64)
    if V <= 0.0
        return 0.0
    end
    return 1.0 / (1.0 + (k_i / V)^alpha_i)
end

function compute_tau_for_individual_idx(data, xi_c, S_P_0=2000.0)
    time_pts = data[:time]
    Psi = data[:Psi]
    V = data[:V]

    res = fill(NaN, 9)
    res[8] = 0.0

    Psi_max_val = maximum(Psi)
    Psi_max_idx = argmax(Psi)
    res[1] = time_pts[Psi_max_idx]
    res[9] = Psi_max_val

    # Uses the globally fixed xi_h
    h_indices = findall(>=(xi_h), Psi)
    if !isempty(h_indices)
        res[2] = time_pts[h_indices[1]]
        res[3] = time_pts[h_indices[end]]
    end

    # Uses the globally fixed xi_d
    if Psi_max_val >= xi_d
        d_idx = findfirst(x -> (x >= xi_d && x ≈ Psi_max_val), Psi)
        if d_idx !== nothing
            res[4] = time_pts[d_idx]
            res[3] = NaN
        end
    end

    V_max_val = maximum(V)
    V_max_idx = argmax(V)
    res[7] = time_pts[V_max_idx]

    # Uses the passed xi_c value
    c_idx = findfirst(>=(xi_c), V)
    if c_idx !== nothing
        res[5] = time_pts[c_idx]
    end

    if isnan(res[4])
        post_peak_indices = (V_max_idx+1):length(V)
        r_offset = findfirst(<=(xi_r), @view V[post_peak_indices])
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
        integ_val += 0.5 * dt * (beta_i(V[j]) + beta_i(V[j+1]))
    end

    res[8] = S_P_0 * integ_val

    return res
end

if N >= 50000
    println("N >= 50,000. Configuring distributed processing...")
    if nprocs() > 1
        try
            rmprocs(workers())
        catch e
        end
    end
    if Sys.CPU_THREADS >= 64
        addprocs(max(2, Int(round(Sys.CPU_THREADS * 2 / 3))))
    else
        addprocs(max(2, Sys.CPU_THREADS - 2))
    end
    println("Done setting up $(nprocs()) workers.")

    @everywhere begin
        const xi_h_w = 75.0
        const xi_d_w = 85.0
        const xi_r_w = 1.0
        const alpha_i_w = 16.422
        const k_i_w = 7.49

        function beta_i_w(V::Float64)
            if V <= 0.0
                return 0.0
            end
            return 1.0 / (1.0 + (k_i_w / V)^alpha_i_w)
        end

        function compute_tau_for_individual_idx_w(data, xi_c, S_P_0=2000.0)
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

            c_idx = findfirst(>=(xi_c), V)
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
end

inputs = [Dict(cohort[i][:vars]) for i in 1:N]
base_name = splitext(basename(latest_file))[1]

for xi_c in xi_c_vals
    println("\n========================================")
    println("Processing for fixed xi_h = $xi_h, xi_d = $xi_d, varying xi_c = $xi_c")
    println("========================================")

    local tau_Psi_max = fill(NaN, N)
    local tau_h_start = fill(NaN, N)
    local tau_h_end = fill(NaN, N)
    local tau_d = fill(NaN, N)
    local tau_c = fill(NaN, N)
    local tau_r = fill(NaN, N)
    local tau_V_max = fill(NaN, N)
    local R0_P2P = fill(0.0, N)
    local Psi_max = fill(NaN, N)

    if N >= 50000
        # Run the computation in parallel
        local results = pmap(inputs) do data
            compute_tau_for_individual_idx_w(data, xi_c)
        end

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
    else
        for i in 1:N
            res = compute_tau_for_individual_idx(inputs[i], xi_c)
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

    # Create output DataFrame
    local out_df = DataFrame()
    out_df[!, :ID] = params_df[!, :ID]
    out_df[!, :R0_within] = params_df[!, :R0_within]
    out_df[!, :R0_P2P] = R0_P2P
    out_df[!, :Psi_max] = Psi_max
    out_df[!, :tau_Psi_max] = tau_Psi_max
    out_df[!, :tau_h_start] = tau_h_start
    out_df[!, :tau_h_end] = tau_h_end
    out_df[!, :tau_d] = tau_d
    out_df[!, :tau_c] = tau_c
    out_df[!, :tau_r] = tau_r
    out_df[!, :tau_V_max] = tau_V_max

    out_df[!, :max_V] = [cohort[i][:maxima][:max_V] for i in 1:N]
    out_df[!, :max_F_U] = [cohort[i][:maxima][:max_F_U] for i in 1:N]
    out_df[!, :max_F_B] = [cohort[i][:maxima][:max_F_B] for i in 1:N]
    out_df[!, :max_V_t] = [cohort[i][:maxima][:max_V_t] for i in 1:N]
    out_df[!, :max_F_U_t] = [cohort[i][:maxima][:max_F_U_t] for i in 1:N]
    out_df[!, :max_F_B_t] = [cohort[i][:maxima][:max_F_B_t] for i in 1:N]

    c_str = isinteger(xi_c) ? string(Int(xi_c)) : string(xi_c)

    local new_base = replace(base_name, "cohort_" => "cohort_times_")
    new_base = new_base * "_xih_75_xid_85_xic_" * c_str

    local out_path_qs = joinpath(OUTPUT_DIR, new_base * ".qs")

    @rput out_df
    @rput out_path_qs
    R"""
    if(!requireNamespace("qs2", quietly=TRUE)) {
        warning("'qs2' package is not installed in R.")
    } else {
        qs2::qs_save(out_df, file = out_path_qs)
    }
    """
    println("Saved successfully to $out_path_qs")
end

if nprocs() > 1 && N >= 50000
    println("Shutting down workers...")
    rmprocs(workers())
end

println("\nAll thresholds processed successfully.")
