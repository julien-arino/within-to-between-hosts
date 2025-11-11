##################
## FUNCTIONS_ALL.jl
##################
#
# This file contains all functions used in the code.
#

using DifferentialEquations
using Random
using Statistics
using DataFrames
using Distributions
using Distributed

@everywhere using Distributed

##
## rhs_within_host
##
# The right-hand side of the within-host DDE
function rhs_within_host!(du, u, h, pp, t)
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
        :d_V => 8.4,
        :p => 394.0,
        :k_U_F => 6.072,
        :p_FI => 2.8235,
        :ψ_F_prod => 0.25,
        :k_B_F => 0.0107,
        :k_lin_f => 16.635,
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
## generate_params_patients
##
# Generate parameters for patients
function generate_params_patients(params, n = 1000)
    names_params = keys(params)
    idx_stddev = filter(x -> occursin("_stddev", string(x)), names_params)  # Convert Symbol to String
    params_with_stddev = map(x -> replace(string(x), "_stddev" => ""), collect(idx_stddev))  # Convert Set to Array

    OUT = DataFrame()
    for curr_col in names_params
        if !(curr_col in params_with_stddev)
            OUT[!, curr_col] = fill(params[curr_col], n)  # Use [!, column] syntax
        else
            mean_val = params[curr_col]
            stddev_val = params[Symbol(string(curr_col) * "_stddev")]  # Convert back to Symbol
            OUT[!, curr_col] = rand(Normal(mean_val, 3 * stddev_val), n)  # Use [!, column] syntax
            OUT[!, curr_col] = map(x -> max(x, mean_val), OUT[!, curr_col])  # Use [!, column] syntax
        end
    end

    # Add initial condition columns
    OUT[!, :V0] = fill(params[:V0], n)
    OUT[!, :S0] = fill(params[:S0], n)
    OUT[!, :I0] = fill(params[:I0], n)
    OUT[!, :R0] = fill(params[:R0], n)

    OUT[!, :ID] = 1:n  # Use [!, column] syntax
    return OUT
end

##
## history_function
##
# Define a standalone history function
function history_function(t, p)
    return IC  # Always return the initial conditions
end

# # Define a custom callable struct for the history function
# struct HistoryFunction
#     IC::Vector{Float64}
# end

# # Make the struct callable
# function (hf::HistoryFunction)(t::Float64)
#     # Return the initial conditions for all state variables
#     return hf.IC
# end

##
## run_one_patient
##
# Simulate the within-host model for one patient
function run_one_patient(idx, patients, IC; type_output = "solution")
    # Extract patient-specific parameters
    params_tmp = patients[idx, :]
    params_tmp = add_IC_to_params(params_tmp, IC)

    # Convert DataFrameRow to a vector of parameter values
    params_vector = collect(values(params_tmp))

    # Define the lag
    lag = params_tmp[:τ_I]

    # Define the time span
    tspan = (0.0, 200.0)

    # Set the integration method
    integrator = MethodOfSteps(Rodas5())

    # Define the DDE problem
    prob = DDEProblem(rhs_within_host!, history_function, tspan, params_vector; constant_lags = [lag])

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

    patients_idx = 1:size(params_change, 1)
    if parallel
        addprocs(ncpus)
        results = pmap(idx -> run_one_patient_indicators(idx, params_change, params_fixed), patients_idx)
    else
        results = map(idx -> run_one_patient_indicators(idx, params_change, params_fixed), patients_idx)
    end
    return results
end