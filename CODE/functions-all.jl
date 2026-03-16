###################
## FUNCTIONS_ALL.jl
###################
#
# This file contains most functions used in the code.
#

using DifferentialEquations
using DiffEqCallbacks
using Random
using Statistics
using DataFrames
using Distributions
using Distributed
using Roots

##
## rhs_within_host_ODE
##
# The right-hand side of the within-host ODE
function rhs_within_host_ODE!(du, u, pp, t)
    V, S, I, R, D, F_U, F_B = u

    # Extract ODE parameters by name to guarantee correct mapping regardless of dataframe column order
    p = pp.p
    λ_S = pp.λ_S
    S_max = pp.S_max
    β_V = pp.β_V
    d_V = pp.d_V
    d_I = pp.d_I
    d_D = pp.d_D
    ψ_F_prod = pp.ψ_F_prod
    p_FI = pp.p_FI
    η_FI = pp.η_FI
    k_lin_f = pp.k_lin_f
    k_int_f = pp.k_int_f
    k_B_F = pp.k_B_F
    c_star = pp.c_star
    k_U_F = pp.k_U_F
    ε_FI = pp.ε_FI
    a_F = pp.a_F

    du[1] = p * I - d_V * V
    du[2] = λ_S * (1 - (S + I + D + R) / S_max) * S - β_V * S * V
    du[3] = β_V * S * V * (1 - F_B / (ε_FI + F_B)) - d_I * I
    du[4] = λ_S * (1 - (S + I + D + R) / S_max) * R + β_V * S * V * (F_B / (ε_FI + F_B))
    du[5] = d_I * I - d_D * D
    du[6] = ψ_F_prod + p_FI * I / (I + η_FI) -
            k_lin_f * F_U - k_B_F * ((c_star + I) * a_F - F_B) * F_U + k_U_F * F_B
    du[7] = -k_int_f * F_B + k_B_F * ((c_star + I) * a_F - F_B) * F_U - k_U_F * F_B
end

##
## set_IC
##
# Set initial conditions
function set_IC()
    return [1.0, 0.16, 0.0, 0.0, 0.0, 0.015, 1.1e-8]
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
        :β_V => 0.3,
        :β_V_stddev => 0.1994,
        :d_V => 8.4,
        :d_V_stddev => 0.67,
        :p => 394.0,
        :p_stddev => 158.65,
        :k_U_F => 6.072,
        :p_FI => 2.8235,
        :p_FI_stddev => 1.8741,
        :ψ_F_prod => 0.25,
        :k_B_F => 0.0107,
        :k_B_F_stddev => 0.01,  # Estimated based on parameter range [0.001, 0.05] pattern; Jenner et al. 2021, eahan et al. 2020
        :k_lin_f => 16.635,
        :k_lin_f_stddev => 2.49,
        :k_int_f => 16.968,
        :ε_FI => 2e-4,
        :η_FI => 0.022328,
        :c_star => 1.104e-4,
        :avo => 6.02214e23,
        :MM_F => 19000.0,
        :R_F_T => 1000.0,
        :R_F_I => 1300.0
    )
    params[:a_F] = (params[:MM_F] / params[:avo]) *
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
function generate_params_cohort(params, n=1000)
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

    # Add varying parameters (sampled from LogNormal distribution)
    for key in keys(params_varying_mean)
        mean_val = params_varying_mean[key]
        stddev_val = params_varying_stddev[key]

        # Calculate parameters for the underlying Normal distribution
        # so that the resulting LogNormal has the desired mean and stddev
        var_val = stddev_val^2
        sigma2 = log(1.0 + var_val / mean_val^2)
        mu = log(mean_val) - sigma2 / 2.0
        sigma = sqrt(sigma2)

        # Sample from LogNormal distribution
        OUT[!, key] = rand(LogNormal(mu, sigma), n)
    end

    # Add ID column
    OUT[!, :ID] = 1:n

    return OUT
end


##
## run_one_individual
##
# Simulate the within-host model for one individual
function run_one_individual(idx, individuals, IC; type_output="solution")
    # Extract individual-specific parameters
    params_tmp = individuals[idx, :]
    params_tmp = add_IC_to_params(params_tmp, IC)

    # Convert DataFrameRow to a correctly named tuple to ensure parameters are explicitly mapped
    params_nt = NamedTuple(params_tmp)

    # Define the time span
    tspan = (0.0, 100.0)

    # Set the integration method
    integrator = Rodas5()
    # Define the ODE problem
    prob = ODEProblem(rhs_within_host_ODE!, IC, tspan, params_nt)

    # Nonnegativity callback
    nonneg_callback = PositiveDomain()

    # Solve the ODE problem with the nonnegativity callback
    sol = solve(prob, integrator; callback=nonneg_callback)

    # Return based on type_output
    if type_output == "solution"
        return sol
    elseif type_output == "maxima"
        return find_maxima(sol)
    elseif type_output == "select_variables"
        # Extract time and selected variables
        S_max = params_tmp[:S_max]
        S = sol[2, :]
        R = sol[4, :]
        Psi = 100 .* (S_max .- (S .+ R)) ./ S_max
        selected_data = Dict(
            :time => sol.t,
            :Psi => Psi,      # Tissue damage
            :V => sol[1, :],  # Viral load
            :I => sol[3, :],  # Infected cells
            :F_B => sol[7, :],  # Bound IFN
            :F_U => sol[6, :]   # Unbound IFN
        )
        return selected_data
    elseif type_output == "vars_and_max"
        S_max = params_tmp[:S_max]
        S = sol[2, :]
        R = sol[4, :]
        Psi = 100 .* (S_max .- (S .+ R)) ./ S_max
        selected_data = Dict(
            :time => sol.t,
            :Psi => Psi,      # Tissue damage
            :I => sol[3, :],  # Infected cells
            :V => sol[1, :],  # Viral load
            :F_B => sol[7, :],  # Bound IFN
            :F_U => sol[6, :]   # Unbound IFN
        )
        maxima_data = find_maxima(sol)
        return Dict(:vars => selected_data, :maxima => maxima_data)
    else
        error("Invalid type_output: must be 'solution', 'maxima', 'select_variables', or 'vars_and_max'")
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

    # Ensure indices are valid before accessing
    OUT[:max_V_t] = sol.t[argmax(sol[1, :])]  # Time of maximum viral load
    OUT[:max_F_U_t] = sol.t[argmax(sol[6, :])]  # Time of maximum unbound IFN
    OUT[:max_F_B_t] = sol.t[argmax(sol[7, :])]  # Time of maximum bound IFN

    return OUT
end

##
## format_df
##
# Format data for plotting
function format_df(time, data; line_plotted="mean", lower="2.5%", upper="97.5%")
    df = DataFrame(
        time=time,
        lower=data[:, lower],
        line=data[:, line_plotted],
        upper=data[:, upper]
    )
    return df
end

##
## value_indicators
##
# Compute indicators for all parameter sets
function value_indicators(params_change, params_fixed; t_f=200, ncpus=60, parallel=true)
    println("Finalizing parameters data frame")
    col_names = names(params_change)[2:end]
    col_names_fixed = setdiff(names(params_fixed), col_names)

    individuals_idx = 1:size(params_change, 1)
    if parallel
        # Remove any pre-existing workers to avoid credential/cookie mismatches
        if nprocs() > 1
            println("Removing existing workers before addprocs (nprocs=$(nprocs()))...")
            try
                rmprocs(workers())
            catch e
                @warn "Failed to remove existing workers in value_indicators" exception = (e, catch_backtrace())
            end
        end
        addprocs(ncpus)
        results = pmap(idx -> run_one_individual_indicators(idx, params_change, params_fixed), individuals_idx)
    else
        results = map(idx -> run_one_individual_indicators(idx, params_change, params_fixed), individuals_idx)
    end
    return results
end


# Compute virus-free equilibrium value of F_U (unbound IFN) and F_B (bound IFN)
function equilibrium_F(params)
    ψ_F_prod = params[:ψ_F_prod]
    k_lin_f = params[:k_lin_f]
    k_int_f = params[:k_int_f]
    k_B_F = params[:k_B_F]
    c_star = params[:c_star]
    a_F = params[:a_F]
    k_U_F = params[:k_U_F]

    a_0 = -ψ_F_prod - k_U_F * ψ_F_prod / k_int_f
    a_1 = k_lin_f + k_B_F * (c_star * a_F - ψ_F_prod / k_int_f) + k_U_F * k_lin_f / k_int_f
    a_2 = k_B_F * k_lin_f / k_int_f

    FU = (-a_1 + sqrt(a_1^2 - 4 * a_2 * a_0)) / (2 * a_2)
    FB = (ψ_F_prod - k_lin_f * FU) / k_int_f
    return (FU, FB)
end

# Compute theta^0 in R0 formula
function theta_zero(FB, params)
    ε_FI = params[:ε_FI]
    return ε_FI / (ε_FI + FB)
end

# Compute R0
function reproduction_number(params)
    FU, FB = equilibrium_F(params)
    θ0 = theta_zero(FB, params)
    p = params[:p]
    β_V = params[:β_V]
    Smax = params[:S_max]
    dI = params[:d_I]
    dV = params[:d_V]
    R0 = p * β_V * Smax * θ0 / (dI * dV)
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

# ---- Compute disease severity (Ψ_i) and classification ----
function compute_severity(out_dict::Dict, S_max::Float64)
    Psi = 100 .* (S_max .- (out_dict[:S] .+ out_dict[:R])) ./ S_max
    Psi_max = maximum(Psi)

    status = if Psi_max >= 85
        "dead"
    elseif Psi_max >= 75
        "ICU"
    else
        "rest"
    end

    idx_ICU = findfirst(>=(75), Psi)
    t_ICU = isnothing(idx_ICU) ? NaN : out_dict[:time][idx_ICU]

    idx_death = findfirst(>=(85), Psi)
    t_death = isnothing(idx_death) ? NaN : out_dict[:time][idx_death]

    return Dict(
        :Psi => Psi,
        :Psi_max => Psi_max,
        :status => status,
        :t_ICU => t_ICU,
        :t_death => t_death
    )
end

##
## rhs_between_host_DID!
##
# The right-hand side of the between-host integro-differential equation
# Here, we use a simple discretization of the history to compute the integral.
# Since U_P(t) is defined by an integral over U_P(t-a), we can either use Delays
# or since U_P(t) acts as an algebraic equation, solve it implicitly or manually 
# integrate it at each timestep using an explicit Euler/RK method stepping algorithm.

function solve_between_host_DID(params, t_end; dt=0.1)
    # Extract parameters
    bP = params[:b_P]
    dP = params[:d_P]

    # These must be provided as functions of infection-age `a`
    β_P = params[:β_P]
    γ_P = params[:γ_P]
    μ_P = params[:μ_P]

    # Maximum age of infection to consider in the integral (to avoid integrating to infinity)
    a_max = get(params, :a_max, 100.0)

    # Initial conditions
    S0 = params[:S_P0]
    U0 = params[:U_P0]

    # Time grid
    t_span = 0.0:dt:t_end
    N_t = length(t_span)

    # State history arrays
    S_P = zeros(N_t)
    U_P = zeros(N_t)

    # Initialize
    S_P[1] = S0
    U_P[1] = U0

    # Precompute the survival probability kernel: exp[-∫_0^a (d_P + γ_P(ξ) + μ_P(ξ)) dξ]
    # We compute this on the same dt grid to optimize the integral computation later
    age_grid = 0.0:dt:a_max
    N_a = length(age_grid)
    survival_prob = zeros(N_a)
    integral_val = 0.0
    for j in 1:N_a
        age = age_grid[j]
        survival_prob[j] = exp(-integral_val)
        if j < N_a
            # Trapezoidal rule for the exponent integral
            integral_val += 0.5 * dt * ((dP + γ_P(age) + μ_P(age)) + (dP + γ_P(age + dt) + μ_P(age + dt)))
        end
    end

    # Precompute β_P on the age grid
    beta_vals = β_P.(age_grid)

    # Time stepping (Euler method for simplicity and robustness in IDEs, matching the integration grid)
    for i in 1:(N_t-1)
        t = t_span[i]

        # 1. Compute incidence U_P(t) by evaluating the integral
        # U_P(t) = S_P(t) ∫ β_P(a) U_P(t-a) exp(...) da
        # We approximate the integral using the history of U_P.
        # Note: For t < a, U_P(t-a) is historically 0 unless we assume an initial endemic history.
        # By default we assume U_P(t) = 0 for t < 0, but there's an initial pulse U0 at t=0.
        integral_U = 0.0

        # Number of historical points available
        max_j = min(i, N_a)

        for j in 1:max_j
            a = age_grid[j]
            idx_history = i - j + 1
            U_hist = U_P[idx_history]

            # Rectangular integration
            integral_U += beta_vals[j] * U_hist * survival_prob[j] * dt
        end

        # Update current U_P(t) based on the computed integral and current S_P
        # If this is the first step, U_P is already set to the initial condition constraint.
        if i > 1
            U_P[i] = S_P[i] * integral_U
        end

        # 2. Update Susceptibles S_P
        # dS_P/dt = b_P - d_P S_P - U_P
        dS = bP - dP * S_P[i] - U_P[i]
        S_P[i+1] = S_P[i] + dt * dS
    end

    # Update the final U_P
    integral_U = 0.0
    max_j = min(N_t, N_a)
    for j in 1:max_j
        idx_history = N_t - j + 1
        U_hist = U_P[idx_history]
        integral_U += beta_vals[j] * U_hist * survival_prob[j] * dt
    end
    U_P[N_t] = S_P[N_t] * integral_U

    return Dict(:t => t_span, :S_P => S_P, :U_P => U_P)
end
