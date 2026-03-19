#!/usr/bin/env julia
println("\n\n>>> Running debug-read-write-big-data.jl ...\n\n")
using RCall, Arrow, DataFrames, Serialization, Random, Printf

println("Generating 1 Million DataFrames inside an Array... (This may take a minute)")
N_individuals = 1_000_000
N_rows_per = 40

list_of_dfs = Vector{DataFrame}(undef, N_individuals)
for i in 1:N_individuals
    list_of_dfs[i] = DataFrame(
        ID = fill(i, N_rows_per),
        time = 1.0:1.0:N_rows_per,
        v1 = rand(N_rows_per),
        v2 = rand(N_rows_per),
        v3 = rand(N_rows_per),
        v4 = rand(N_rows_per),
        v5 = rand(N_rows_per)
    )
end

println("-> DataFrames Array generated. Memory roughly ~2.2 GB.")

println("\n=======================================================")
println("OPTION A: Legacy Pipeline (RCall memory transfer -> qs2)")
println("=======================================================")
println("  >> Timing the transfer of 1M DataFrames to R memory via @rput...")
@time @rput list_of_dfs

println("  >> Timing the save to disk via R's `qs2::qs_save`...")
@time R"qs2::qs_save(list_of_dfs, 'shared_data_from_julia.qs', nthreads = N_QS_THREADS)"


println("\n=======================================================")
println("OPTION B: Modern Pipeline (Julia Flatten -> Arrow native)")
println("=======================================================")
println("  >> Timing the flattening of 1M DataFrames into a single continuous table...")
@time df_flat = reduce(vcat, list_of_dfs)

println("  >> Timing the save to disk natively via `Arrow.write`...")
@time Arrow.write("shared_data_from_julia.arrow", df_flat)

println("\nJulia sequence complete! Run the R script next to benchmark the loading.")

println("\n\n>>> debug-read-write-big-data.jl successfully finished running ✅\n\n")
