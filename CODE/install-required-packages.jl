#!/usr/bin/env julia
## install-required-packages.jl
# Ensure that all Julia packages required by the project are installed.

using Pkg

# List of required packages based on the project's source code
required_packages = [
    "CSV",
    "DataFrames",
    "Dates",
    "DiffEqCallbacks",
    "DifferentialEquations",
    "Distributed",
    "Distributions",
    "Glob",
    "OrderedCollections",
    "Plots",
    "Printf",
    "RCall",
    "Random",
    "Roots",
    "Serialization",
    "Statistics"
]

println("Checking and installing required Julia packages...")

# Iterate through the list and install if missing
for pkg in required_packages
    try
        # Using Pkg.add directly will check if it's already installed and 
        # resolve it. It does not reinstall if everything is up to date, 
        # so calling Pkg.add on an already installed package is safe.
        println("Ensuring package is installed: ", pkg)
        Pkg.add(pkg)
    catch e
        println("Error installing $pkg: ", e)
    end
end

println("All required Julia packages are installed.")
