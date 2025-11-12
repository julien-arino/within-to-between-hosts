###################
## FUNCTIONS_ALL.jl
###################
#
# This file contains all functions used in the code.
#

using DifferentialEquations
using Random
using Statistics
using DataFrames
using Distributions
using Distributed
using Roots

@everywhere using Distributed

##
## rhs_within_host_ODE
##
# The right-hand side of the within-host ODE
function rhs_within_host_ODE!(du, u, pp, t)
    V, S, I, R, D, F_U, F_B, A_I, A_R = u
    p, τ_I, V0, S0, λ_S, S_max, β, d_V, d_I, d_D, ψ_F_prod, p_FI, η_FI, k_lin_f, k_int_f, k_B_F, T_star, k_U_F, δ, ε_FI, A_F = pp

    du[1] = p * I - d_V * V
    du[2] = λ_S * (1 - (S + I + D + R) / S_max) * S - β * S * V
    du[3] = β * S * V * (1 - F_B / (ε_FI + F_B)) * A_I - d_I * I
    du[4] = λ_S * (1 - (S + I + D + R) / S_max) * R +
            β * S * V * (F_B / (ε_FI + F_B)) * A_R
    du[5] = d_I * I - d_D * D
    du[6] = ψ_F_prod + p_FI * I / (I + η_FI) -
            k_lin_f * F_U - k_B_F * ((T_star + I) * A_F - F_B) * F_U + k_U_F * F_B
    du[7] = -k_int_f * F_B + k_B_F * ((T_star + I) * A_F - F_B) * F_U - k_U_F * F_B
    du[8] = δ * A_I * (I - I)
    du[9] = δ * A_R * (R - R)
end

##
## rhs_within_host_DDE
## (The original model used DDEs; this function is kept for reference)
##
# The right-hand side of the within-host DDE
function rhs_within_host_DDE!(du, u, h, pp, t)
    V, S, I, R, D, F_U, F_B, A_I, A_R = u
    p, τ_I, V0, S0, λ_S, S_max, β, d_V, d_I, d_D, ψ_F_prod, p_FI, η_FI, k_lin_f, k_int_f, k_B_F, T_star, k_U_F, δ, ε_FI, A_F = pp

    lagged_values = h(pp, t-τ_I)
    V_t, S_t, I_t, R_t = lagged_values[1], lagged_values[2], lagged_values[3], lagged_values[4]

    du[1] = p * I - d_V * V
    du[2] = λ_S * (1 - (S + I + D + R) / S_max) * S - β * S * V
    du[3] = β * S_t * V_t * (1 - F_B / (ε_FI + F_B)) * A_I - d_I * I
    du[4] = λ_S * (1 - (S + I + D + R) / S_max) * R +
            β * S_t * V_t * (F_B / (ε_FI + F_B)) * A_R
    du[5] = d_I * I - d_D * D
    du[6] = ψ_F_prod + p_FI * I / (I + η_FI) -
            k_lin_f * F_U - k_B_F * ((T_star + I) * A_F - F_B) * F_U + k_U_F * F_B
    du[7] = -k_int_f * F_B + k_B_F * ((T_star + I) * A_F - F_B) * F_U - k_U_F * F_B
    du[8] = δ * A_I * (I_t - I)
    du[9] = δ * A_R * (R_t - R)
end

##
## set_IC
##
# Set initial conditions
function set_IC()
    return [1.0, 0.16, 0.0, 0.0, 0.0, 0.015, 1.1e-8, 1.0, 1.0]
end

##
## set_parameters
##
# Set parameters
function set_parameters()
    params = Dict(
        :λ_S => 0.74,
        :S_max => 0.16,
        :d_I => 0.1,
        :d_D => 8.0,
        :τ_I => 0.17,
        :β => 0.3,
        :β_stddev => 0.1994,
        :d_V => 8.4,
        :d_V_stddev => 0.67,
        :p => 394.0,
        :p_stddev => 158.65,
        :k_U_F => 6.072,
        :p_FI => 2.8235,
        :p_FI_stddev => 1.8741,
        :ψ_F_prod => 0.25,
        :k_B_F => 0.0107,
        :k_B_F_stddev => 0.01,  # Estimated based on parameter range [0.001, 0.05] pattern; Jenner et al. 2021 ref [54] Sheahan et al. 2020
        :k_lin_f => 16.635,
        :k_lin_f_stddev => 2.49,
        :k_int_f => 16.968,
        :ε_FI => 2e-4,
        :η_FI => 0.022328,
        :T_star => 1.104e-4,
        :δ => 0.1,
        :avo => 6.02214e23,
        :MM_F => 19000.0,
        :R_F_T => 1000.0,
        :R_F_I => 1300.0
    )
    params[:A_F] = (params[:MM_F] / params[:avo]) *
                   (params[:R_F_I] + params[:R_F_T]) *
                   (1 / 5000) * (10^9 * 1e12)

    # Add initial condition keys
    IC = set_IC()
    params[:V0] = IC[1]
    params[:S0] = IC[2]
    params[:I0] = IC[3]
    params[:R0] = IC[4]

    return params
end

##
## add_IC_to_params
##
# Add initial conditions to parameters
function add_IC_to_params(params, IC)
    params[:V0] = IC[1]
    params[:S0] = IC[2]
    params[:I0] = IC[3]
    params[:R0] = IC[4]
    return params
end

##
## generate_params_cohort
##
# Generate parameters for individuals in the virtual cohort
function generate_params_cohort(params, n = 1000)
    # Step 1: Extract stddev parameters and create params_varying_stddev
    all_keys = collect(keys(params))
    stddev_keys = filter(x -> occursin("_stddev", string(x)), all_keys)
    
    params_varying_stddev = Dict()
    for key in stddev_keys
        # Remove "_stddev" suffix from the key name
        new_key = Symbol(replace(string(key), "_stddev" => ""))
        params_varying_stddev[new_key] = params[key]
    end
    
    # Step 2: Create params_varying_mean (corresponding mean values)
    params_varying_mean = Dict()
    for key in keys(params_varying_stddev)
        params_varying_mean[key] = params[key]
    end
    
    # Step 3: Create params_fixed (everything else, excluding _stddev entries)
    varying_param_names = keys(params_varying_mean)
    params_fixed = Dict()
    for key in all_keys
        # Include if it's not a varying parameter and not a _stddev entry
        if !(key in varying_param_names) && !occursin("_stddev", string(key))
            params_fixed[key] = params[key]
        end
    end
    
    # Step 4: Build the output DataFrame
    OUT = DataFrame()
    
    # Add fixed parameters (same value for all individuals)
    for (key, val) in params_fixed
        OUT[!, key] = fill(val, n)
    end
    
    # Add varying parameters (sampled from Normal distribution)
    for key in keys(params_varying_mean)
        mean_val = params_varying_mean[key]
        stddev_val = params_varying_stddev[key]
        
        # Sample with symmetric variation and a scaling factor of 1 on the provided stddev
        OUT[!, key] = rand(Normal(mean_val, 1 * stddev_val), n)
        
        # Enforce non-negativity with a small minimum to avoid numerical issues
        min_val = max(1e-10, mean_val * 0.001)  # At least 0.1% of mean value
        OUT[!, key] = map(x -> max(x, min_val), OUT[!, key])
    end
    
    # Add ID column
    OUT[!, :ID] = 1:n
    
    return OUT
end

##
## history_function
##
# Define a standalone history function
function history_function(t, p)
    return IC  # Always return the initial conditions
end

##
## run_one_individual
##
# Simulate the within-host model for one individual
function run_one_individual(idx, individuals, IC; type_output = "solution")
    # Extract individual-specific parameters
    params_tmp = individuals[idx, :]
    params_tmp = add_IC_to_params(params_tmp, IC)

    # Convert DataFrameRow to a vector of parameter values
    params_vector = collect(values(params_tmp))

    # Define the time span
    tspan = (0.0, 200.0)

    ## The original model used DDEs; here we use ODEs for simplicity, but the DDE code
    ## is kept for reference.
    # # Define the lag
    # lag = params_tmp[:τ_I]
    # # Set the integration method
    # integrator = MethodOfSteps(Rodas5())
    # # Define the DDE problem
    # prob = DDEProblem(rhs_within_host!, history_function, tspan, params_vector; constant_lags = [lag])

    # Set the integration method
    integrator = Rodas5()
    # Define the ODE problem
    prob = ODEProblem(rhs_within_host_ODE!, IC, tspan, params_vector)

    # --- Nonnegativity callback defined locally ---
    condition_nonneg(u, t, integrator) = any(x -> x < 0, u)
    affect_nonneg!(integrator) = (integrator.u .= max.(integrator.u, 0.0))
    nonneg_callback = DiscreteCallback(condition_nonneg, affect_nonneg!)
    # ---------------------------------------------

    # Solve the DDE problem with the nonnegativity callback
    sol = solve(prob, integrator; callback=nonneg_callback)

    # Return based on type_output
    if type_output == "solution"
        return sol
    elseif type_output == "maxima"
        return find_maxima(sol)
    elseif type_output == "select_variables"
        # Extract time and selected variables
        selected_data = Dict(
            :time => sol.t,
            :V => sol[1, :],  # Viral load
            :I => sol[3, :],  # Infected cells
            :F_B => sol[7, :],  # Bound IFN
            :F_U => sol[6, :]   # Unbound IFN
        )
        return selected_data
    else
        error("Invalid type_output: must be 'solution', 'maxima', or 'select_variables'")
    end
end

##
## find_maxima
##
# Extract maxima and their corresponding times from the solution
function find_maxima(sol)
    # Ensure the solution has enough data points
    if length(sol.t) == 0
        error("Solution has no time points")
    end

    # Extract the desired maxima and their corresponding times
    OUT = Dict()
    OUT[:max_V] = maximum(sol[1, :])  # Maximum viral load
    OUT[:max_F_U] = maximum(sol[6, :])  # Maximum unbound IFN
    OUT[:max_F_B] = maximum(sol[7, :])  # Maximum bound IFN
    OUT[:max_A_I] = maximum(sol[8, :])  # Maximum active infected cells
    OUT[:max_A_R] = maximum(sol[9, :])  # Maximum active recovered cells

    # Ensure indices are valid before accessing
    OUT[:max_V_t] = sol.t[argmax(sol[1, :])]  # Time of maximum viral load
    OUT[:max_F_U_t] = sol.t[argmax(sol[6, :])]  # Time of maximum unbound IFN
    OUT[:max_F_B_t] = sol.t[argmax(sol[7, :])]  # Time of maximum bound IFN
    OUT[:max_A_I_t] = sol.t[argmax(sol[8, :])]  # Time of maximum active infected cells
    OUT[:max_A_R_t] = sol.t[argmax(sol[9, :])]  # Time of maximum active recovered cells

    return OUT
end

##
## format_df
##
# Format data for plotting
function format_df(time, data; line_plotted = "mean", lower = "2.5%", upper = "97.5%")
    df = DataFrame(
        time = time,
        lower = data[:, lower],
        line = data[:, line_plotted],
        upper = data[:, upper]
    )
    return df
end

##
## value_indicators
##
# Compute indicators for all parameter sets
function value_indicators(params_change, params_fixed; t_f = 200, ncpus = 60, parallel = true)
    println("Finalizing parameters data frame")
    col_names = names(params_change)[2:end]
    col_names_fixed = setdiff(names(params_fixed), col_names)

    individuals_idx = 1:size(params_change, 1)
    if parallel
        addprocs(ncpus)
        results = pmap(idx -> run_one_individual_indicators(idx, params_change, params_fixed), individuals_idx)
    else
        results = map(idx -> run_one_individual_indicators(idx, params_change, params_fixed), individuals_idx)
    end
    return results
end


# Compute virus-free equilibrium value of F_U (unbound IFN)
# Solves: ψ_F_prod - k_lin_f*F_U = k_int_f * (k_B_F * T_star * A_F * F_U) / (k_U_F + F_U)
function equilibrium_FU(params)
    ψ_F_prod = params[:ψ_F_prod]
    k_lin_f  = params[:k_lin_f]
    k_int_f  = params[:k_int_f]
    k_B_F    = params[:k_B_F]
    T_star   = params[:T_star]
    A_F      = params[:A_F]
    k_U_F    = params[:k_U_F]

    f(FU) = ψ_F_prod - k_lin_f*FU - k_int_f * (k_B_F * T_star * A_F * FU) / (k_U_F + FU)
    FU_root = find_zero(f, (0.0, 1e3), Bisection())  # search in positive range
    return FU_root
end

# Compute F_U factor in R0 formula
function F_U_factor(FU, params)
    ε_FI   = params[:ε_FI]
    k_U_F  = params[:k_U_F]
    A_F    = params[:A_F]
    T_star = params[:T_star]
    k_B_F  = params[:k_B_F]
    numerator = ε_FI * (FU + k_U_F)
    denominator = FU * A_F * T_star * k_B_F + ε_FI * (FU + k_U_F)
    return numerator / denominator
end

# Compute R0
function reproduction_number(params)
    FU = equilibrium_FU(params)
    Ffac = F_U_factor(FU, params)
    p    = params[:p]
    β    = params[:β]
    Smax = params[:S_max]
    dI   = params[:d_I]
    dV   = params[:d_V]
    R0 = p * β * Smax * Ffac / (dI * dV)
    return R0
end

##
## compute_R0_cohort
##
# Compute R0 for each individual in the cohort DataFrame
# and add it as a new column
function compute_R0_cohort!(individuals::DataFrame)
    n = nrow(individuals)
    R0_values = zeros(n)
    
    for i in 1:n
        # Extract parameters for individual i as a Dict-like object
        params_i = individuals[i, :]
        R0_values[i] = reproduction_number(params_i)
    end
    
    # Add R0 column to the DataFrame
    individuals[!, :R0_within] = R0_values
    return individuals
end

# Example usage of R0 (uncomment to test):
# include("functions-all.jl")
# params = set_parameters()
# println("Equilibrium F_U: ", equilibrium_FU(params))
# println("R0: ", reproduction_number(params))
