# process-envelopes.jl
# Interpolates virtual cohort solutions onto a fixed time grid and computes
# cross-sectional descriptive statistics (envelopes) including mean and percentiles.
# Usage: julia -t auto process-envelopes.jl (using threads reduces RAM usage significantly compared to Distributed)

using RCall
using DataFrames
using Glob
using Statistics

# Default inputs
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = normpath(joinpath(SCRIPT_DIR, "..", "OUTPUT"))

# Find the latest results file (either qs or Rds)
result_files = vcat(glob("*.qs", OUTPUT_DIR), glob("*.Rds", OUTPUT_DIR))
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
        if(!requireNamespace("qs", quietly=TRUE)) {
            stop("Neither 'qs2' nor 'qs' packages are installed in R.")
        } else {
            save_data <- qs::qread($latest_file)
        }
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

# Extract the cohort simulation results
cohort = SAVE[:cohort]
N = length(cohort)

println("Loaded results for cohort of size: $N")

# Time grid Configuration
const t_start = 0.0
const t_end = 100.0
const t_step = 0.1
t_grid = collect(t_start:t_step:t_end)
num_timepoints = length(t_grid)

# Define the variables to interpolate
vars_to_process = [:V, :I, :Psi, :F_B, :F_U]

function interp1d(x, y, xi)
    yi = zeros(length(xi))
    for k in 1:length(xi)
        xk = xi[k]
        if xk <= x[1]
            yi[k] = y[1]
        elseif xk >= x[end]
            yi[k] = y[end]
        else
            idx = searchsortedlast(x, xk)
            if idx == 0
                yi[k] = y[1]
            elseif idx == length(x)
                yi[k] = y[end]
            else
                t = (xk - x[idx]) / (x[idx+1] - x[idx])
                yi[k] = y[idx] + t * (y[idx+1] - y[idx])
            end
        end
    end
    return yi
end

println("Computing envelopes...")
println("Active threads: ", Threads.nthreads())

out_df = DataFrame(
    time = Float64[],
    variable = String[],
    mean = Float64[],
    p05 = Float64[],
    p10 = Float64[],
    p50 = Float64[],
    p90 = Float64[],
    p95 = Float64[]
)

for var in vars_to_process
    println("Processing $var ...")
    mat = zeros(Float64, N, num_timepoints)
    
    # Fill the matrix safely with shared memory threading
    Threads.@threads for i in 1:N
        data = haskey(cohort[i], :vars) ? cohort[i][:vars] : cohort[i]
        
        if !haskey(data, var)
            mat[i, :] .= NaN
        else
            time_pts = data[:time]
            vals = data[var]
            
            # Filter NaNs or identical points
            if length(time_pts) > 0
                mat[i, :] .= interp1d(time_pts, vals, t_grid)
            else
                mat[i, :] .= NaN
            end
        end
    end
    
    # If all NaNs, skip
    if all(isnan.(mat))
        println("Skipping $var (Data contains NaN natively or not tracked)")
        continue
    end
    
    # Compute cross-sectional metrics
    println("  Computing statistics for $var ...")
    v_mean = zeros(num_timepoints)
    v_p05 = zeros(num_timepoints)
    v_p10 = zeros(num_timepoints)
    v_p50 = zeros(num_timepoints)
    v_p90 = zeros(num_timepoints)
    v_p95 = zeros(num_timepoints)
    
    # Using threads for the 1000 time column percentiles computation
    Threads.@threads for t_idx in 1:num_timepoints
        slice = filter(!isnan, mat[:, t_idx]) # handle potential NaNs
        if length(slice) > 0
            v_mean[t_idx] = mean(slice)
            v_p05[t_idx] = quantile(slice, 0.05)
            v_p10[t_idx] = quantile(slice, 0.10)
            v_p50[t_idx] = quantile(slice, 0.50)
            v_p90[t_idx] = quantile(slice, 0.90)
            v_p95[t_idx] = quantile(slice, 0.95)
        else
            v_mean[t_idx] = NaN
            v_p05[t_idx] = NaN
            v_p10[t_idx] = NaN
            v_p50[t_idx] = NaN
            v_p90[t_idx] = NaN
            v_p95[t_idx] = NaN
        end
    end
    
    sub_df = DataFrame(
        time = t_grid,
        variable = fill(string(var), num_timepoints),
        mean = v_mean,
        p05 = v_p05,
        p10 = v_p10,
        p50 = v_p50,
        p90 = v_p90,
        p95 = v_p95
    )
    
    append!(out_df, sub_df)
    
    # Free memory
    mat = nothing
    GC.gc()
end

using CSV
out_path = joinpath(OUTPUT_DIR, "process-envelopes-results.csv")
CSV.write(out_path, out_df)
println("Successfully saved envelopes data to $out_path")
