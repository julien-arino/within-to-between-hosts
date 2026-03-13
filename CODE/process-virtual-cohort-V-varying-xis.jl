## process-virtual-cohort-V-varying-xis.jl
# Processes output files from run-virtual-cohort-sims.jl and times output files 
# to compute the V values when individuals are transmitting.

using RCall
using DataFrames
using Glob
using Distributed
using Printf

# Default inputs
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = normpath(joinpath(SCRIPT_DIR, "..", "OUTPUT"))

# 1. Find the latest base cohort result file
base_files = vcat(glob("cohort_P*.qs", OUTPUT_DIR), glob("cohort_P*.Rds", OUTPUT_DIR))
base_files = filter(f -> !occursin("_times", f) && !occursin("_censored", f) && !occursin("_V_", f), base_files)

if isempty(base_files)
    error("No raw cohort .qs or .Rds result files found in $OUTPUT_DIR")
end

latest_base_file = sort(base_files, by=mtime, rev=true)[1]
println("Loading base result file: $latest_base_file")

# Load the base file
if endswith(latest_base_file, ".qs")
    R"""
    if(!requireNamespace("qs2", quietly=TRUE)) {
        stop("'qs2' package is not installed in R.")
    } else {
        base_data <- qs2::qs_read($latest_base_file)
    }
    """
else
    R"""
    base_data <- readRDS($latest_base_file)
    """
end

BASE_SAVE = rcopy(R"base_data")
cohort = BASE_SAVE[:cohort]
N = length(cohort)
println("Loaded base cohort size: $N")

# Set up raw data arrays to avoid passing the whole cohort down
times_raw = [cohort[i][:vars][:time] for i in 1:N]
V_raw = [cohort[i][:vars][:V] for i in 1:N]

# 2. Find all matching times files
base_name_pure = splitext(basename(latest_base_file))[1]
times_files = glob(replace(base_name_pure, "cohort_" => "cohort_times_") * "*.qs", OUTPUT_DIR)

if isempty(times_files)
    error("No corresponding cohort_times_P*.qs files found for $base_name_pure in $OUTPUT_DIR")
end

# Linear interpolation helper
function interpolate_v(t_target::Float64, times::Vector{Float64}, V::Vector{Float64})
    # If exactly on a point
    if t_target in times
        idx = findfirst(==(t_target), times)
        return V[idx]
    end
    
    # Boundary conditions
    if t_target <= times[1]
        return V[1]
    end
    if t_target >= times[end]
        return V[end]
    end
    
    # Find surrounding indices
    idx_upper = findfirst(>(t_target), times)
    idx_lower = idx_upper - 1
    
    t0 = times[idx_lower]
    t1 = times[idx_upper]
    v0 = V[idx_lower]
    v1 = V[idx_upper]
    
    # Linear interpolation
    fraction = (t_target - t0) / (t1 - t0)
    return v0 + fraction * (v1 - v0)
end

function process_individual_v!(output_t::Vector{Float64}, output_v::Vector{Float64}, 
                              times::Vector{Float64}, V::Vector{Float64},
                              tau_c::Float64, tau_h_start::Float64, tau_h_end::Float64, 
                              tau_d::Float64, tau_r::Float64; dt::Float64=0.1)
    
    # If the individual never reached infectious period, skip
    if isnan(tau_c)
        return
    end

    # Define segment 1: [tau_c, min(tau_h_start, tau_d, tau_r)]
    # tau_h_start, tau_d, tau_r may be NaN. Treat NaN as Inf for finding the end bound.
    end1 = min(isnan(tau_h_start) ? Inf : tau_h_start, 
               isnan(tau_d) ? Inf : tau_d, 
               isnan(tau_r) ? Inf : tau_r)
    
    # Safeguard
    if isinf(end1)
        end1 = times[end]
    end

    if end1 > tau_c
        t_grid1 = collect(tau_c:dt:end1)
        if t_grid1[end] < end1
            push!(t_grid1, end1)
        end
        
        for t in t_grid1
            push!(output_t, t)
            push!(output_v, interpolate_v(t, times, V))
        end
    end

    # Define segment 2: If hospitalized AND did not die, transmit from tau_h_end to tau_r
    if !isnan(tau_h_end) && isnan(tau_d)
        end2 = isnan(tau_r) ? Inf : tau_r
        if isinf(end2)
            end2 = times[end]
        end
        
        if end2 > tau_h_end
            t_grid2 = collect(tau_h_end:dt:end2)
            if t_grid2[end] < end2
                push!(t_grid2, end2)
            end
            
            for t in t_grid2
                push!(output_t, t)
                push!(output_v, interpolate_v(t, times, V))
            end
        end
    end
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
        addprocs(max(2, Sys.CPU_THREADS * 2 ÷ 3))
    else
        addprocs(max(2, Sys.CPU_THREADS - 2))
    end
    println("Done setting up $(nprocs()) workers.")

    @everywhere begin
        function interpolate_v_w(t_target::Float64, times::Vector{Float64}, V::Vector{Float64})
            if t_target in times
                idx = findfirst(==(t_target), times)
                return V[idx]
            end
            if t_target <= times[1] return V[1] end
            if t_target >= times[end] return V[end] end
            
            idx_upper = findfirst(>(t_target), times)
            idx_lower = idx_upper - 1
            t0 = times[idx_lower]; t1 = times[idx_upper]
            v0 = V[idx_lower]; v1 = V[idx_upper]
            fraction = (t_target - t0) / (t1 - t0)
            return v0 + fraction * (v1 - v0)
        end

        function process_individual_v_w(times::Vector{Float64}, V::Vector{Float64},
                                        tau_c::Float64, tau_h_start::Float64, tau_h_end::Float64, 
                                        tau_d::Float64, tau_r::Float64; dt::Float64=0.1)
            output_t = Float64[]
            output_v = Float64[]
            
            if isnan(tau_c)
                return output_t, output_v
            end

            end1 = min(isnan(tau_h_start) ? Inf : tau_h_start, 
                       isnan(tau_d) ? Inf : tau_d, 
                       isnan(tau_r) ? Inf : tau_r)
            if isinf(end1) end1 = times[end] end

            if end1 > tau_c
                t_grid1 = collect(tau_c:dt:end1)
                if t_grid1[end] < end1 push!(t_grid1, end1) end
                for t in t_grid1
                    push!(output_t, t)
                    push!(output_v, interpolate_v_w(t, times, V))
                end
            end

            if !isnan(tau_h_end) && isnan(tau_d)
                end2 = isnan(tau_r) ? Inf : tau_r
                if isinf(end2) end2 = times[end] end
                
                if end2 > tau_h_end
                    t_grid2 = collect(tau_h_end:dt:end2)
                    if t_grid2[end] < end2 push!(t_grid2, end2) end
                    for t in t_grid2
                        push!(output_t, t)
                        push!(output_v, interpolate_v_w(t, times, V))
                    end
                end
            end
            return output_t, output_v
        end
    end
end

for times_file in times_files
    println("\n========================================")
    println("Processing $times_file")
    println("========================================")

    # Load the times df
    R"""
    times_data <- qs2::qs_read($times_file)
    """
    times_df = rcopy(R"times_data")

    # Double check it matches the base data
    if nrow(times_df) != N
        @warn "Times file rows ($(nrow(times_df))) do not match cohort size ($N). Skipping."
        continue
    end

    if N >= 50000
        # Package inputs for pmap
        # Each input tuple: (times, V, tau_c, tau_h_start, tau_h_end, tau_d, tau_r)
        inputs_p = Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64}[]
        for i in 1:N
            push!(inputs_p, (times_raw[i], V_raw[i], 
                             times_df[i, :tau_c], times_df[i, :tau_h_start], 
                             times_df[i, :tau_h_end], times_df[i, :tau_d], times_df[i, :tau_r]))
        end

        results = pmap(inputs_p) do args
            process_individual_v_w(args...)
        end
        
        # Package into DataFrame
        out_df = DataFrame(ID = times_df[!, :ID],
                           time = [Vector{Float64}() for _ in 1:N], 
                           V = [Vector{Float64}() for _ in 1:N])
                           
        for i in 1:N
            out_df[i, :time] = results[i][1]
            out_df[i, :V] = results[i][2]
        end
    else
        # Serial processing
        out_df = DataFrame(ID = times_df[!, :ID],
                           time = [Vector{Float64}() for _ in 1:N], 
                           V = [Vector{Float64}() for _ in 1:N])

        for i in 1:N
            tau_c = times_df[i, :tau_c]
            tau_h_start = times_df[i, :tau_h_start]
            tau_h_end = times_df[i, :tau_h_end]
            tau_d = times_df[i, :tau_d]
            tau_r = times_df[i, :tau_r]
            
            output_t = Float64[]
            output_v = Float64[]
            process_individual_v!(output_t, output_v, times_raw[i], V_raw[i], tau_c, tau_h_start, tau_h_end, tau_d, tau_r)
            
            out_df[i, :time] = output_t
            out_df[i, :V] = output_v
        end
    end

    # Extract varying xis from the input file name to construct the new output file name
    # cohort_times_P1000000_DT20260312-210843_vars_and_max_xih_50_xid_75.qs
    time_base = splitext(basename(times_file))[1]
    new_base = replace(time_base, "cohort_times_" => "cohort_V_")
    out_path_qs = joinpath(OUTPUT_DIR, new_base * ".qs")

    @rput out_df
    @rput out_path_qs
    R"""
    qs2::qs_save(out_df, file = out_path_qs)
    """
    println("Saved successfully to $out_path_qs")
end

if nprocs() > 1 && N >= 50000
    println("Shutting down workers...")
    rmprocs(workers())
end

println("\nAll times files processed successfully.")
