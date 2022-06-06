using CairoMakie
using Statistics
using DelimitedFiles

include("../src/pottsmc.jl")
include("../src/observables.jl")
include("../src/statutils.jl")

#=
Perform simulation
=#

L = 16
q = 3
cold_start = true
err_nblocks = 100
esteps = 1000  # Number of steps for equilibration
nconfigs = 5000  # Number of steps for measurements


data_path = "D:\\Projects\\Potts-QCD\\potts-model-monte-carlo\\data\\corrtime_data.txt"
data = readdlm(data_path, ',', Float64)
Temps = data[1, :]
autocorr_steps_wolff = convert.(Int64, data[3, :])

println(".==================================")
println("| Lattice Size: $(L) x $(L)        ")
println(".==================================")
println("|  ")

for i = 1:length(Temps)
    T = Temps[i]
    println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")

    location="ising_uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"
    potts = initialize_model_2d(L, q; cold_start=cold_start)
    uncorrelated_spins = zeros(Float64, (L, L, nconfigs))
    for step=1:esteps
        wolff_cluster_update!(potts, T)
    end

    twice_τ = 2*autocorr_steps_wolff[i]
    nsteps = twice_τ*nconfigs
    for j in 1:nsteps
        wolff_cluster_update!(potts, T)
        if j%twice_τ == 0
            uncorrelated_spins[:, :, j÷twice_τ] = copy(potts.lattice)
        end
    end

    println("| Process complete on thread #$(Threads.threadid()): T = $T")
end

println("| Program Finished!")
println(".==================================")