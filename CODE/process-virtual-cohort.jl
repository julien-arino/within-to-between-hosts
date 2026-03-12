# process-virtual-cohort.jl
# Processes output files from run-virtual-cohort-sims.jl and computes the
# various tau times defined in the manuscript.

using RCall
using DataFrames
using Glob
using Statistics

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
    
    # Preallocate returns: [tau_Psi_max, tau_h_start, tau_h_end, tau_d, tau_c, tau_r, tau_V_max, R0_P2P]
    res = fill(NaN, 8)
    res[8] = 0.0
    
    # ------------------- #
    # Eq 2.2 / Eq 5.2     #
    # ------------------- #
    
    Psi_max_val = maximum(Psi)
    Psi_max_idx = argmax(Psi)
    res[1] = time_pts[Psi_max_idx] # tau_Psi_max
    
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
            
            res = fill(NaN, 8)
            res[8] = 0.0
            
            Psi_max_val = maximum(Psi)
            Psi_max_idx = argmax(Psi)
            res[1] = time_pts[Psi_max_idx]
            
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
        tau_h_end[i]   = res[3]
        tau_d[i]       = res[4]
        tau_c[i]       = res[5]
        tau_r[i]       = res[6]
        tau_V_max[i]   = res[7]
        R0_P2P[i]      = res[8]
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
        tau_h_end[i]   = res[3]
        tau_d[i]       = res[4]
        tau_c[i]       = res[5]
        tau_r[i]       = res[6]
        tau_V_max[i]   = res[7]
        R0_P2P[i]      = res[8]
    end
end

# Create a much cleaner output DataFrame retaining only exactly what the user requested
out_df = DataFrame()
out_df[!, :ID] = params_df[!, :ID]
out_df[!, :R0_within] = params_df[!, :R0_within]
out_df[!, :R0_P2P] = R0_P2P

# Note: The manuscript specifies xi_h=75% and xi_d=85% for severity, and V>=4.5 for transmission
out_df[!, :tau_Psi_max] = tau_Psi_max
out_df[!, :tau_h_start] = tau_h_start
out_df[!, :tau_h_end]   = tau_h_end
out_df[!, :tau_d]       = tau_d
out_df[!, :tau_c]       = tau_c
out_df[!, :tau_r]       = tau_r
out_df[!, :tau_V_max]   = tau_V_max

# Extract the original maximums natively tracked during the ODE solver run back in `run-virtual-cohort-sims.jl`
out_df[!, :max_V] = [cohort[i][:maxima][:max_V] for i in 1:N]
out_df[!, :max_F_U] = [cohort[i][:maxima][:max_F_U] for i in 1:N]
out_df[!, :max_F_B] = [cohort[i][:maxima][:max_F_B] for i in 1:N]
out_df[!, :max_V_t] = [cohort[i][:maxima][:max_V_t] for i in 1:N]
out_df[!, :max_F_U_t] = [cohort[i][:maxima][:max_F_U_t] for i in 1:N]
out_df[!, :max_F_B_t] = [cohort[i][:maxima][:max_F_B_t] for i in 1:N]

# Display a quick summary of the results
println("\n--- Summary of Computed Taus ---")
survived = isnan.(tau_d)
died = .!survived
hosp = .!isnan.(tau_h_start)
inf_start = .!isnan.(tau_c)

println("Total Individuals: ", N)
println("Survived: ", sum(survived))
println("Died (tau_d present): ", sum(died))
println("Hospitalized (tau_h_start present): ", sum(hosp))
println("Started Infectious Period (tau_c present): ", sum(inf_start))
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
