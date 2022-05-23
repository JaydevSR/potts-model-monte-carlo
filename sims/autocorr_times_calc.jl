include("../src/potts.jl")
include("../src/statutils.jl")

using CairoMakie
using Statistics

#=
Perform simulation
=#

L = 16
q = 3
println("================================\n")
println("    Lattice Size: $(L) x $(L)")
println("================================\n")

# Temps = [1.0]
Temps = [i for i=0.4:0.1:1.6]
nsteps = 10000  # Number of steps for measurements

auto_corr_times_metro = zeros(Int64, length(Temps))
auto_corr_times_wolff = zeros(Int64, length(Temps))

potts = initialize_model_2d(L, q)
lattice_0 = copy(potts.lattice)

for i = 1:length(Temps)
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    E_arr = zeros(Float64, nsteps)

    # For Metropolis
    potts.lattice = lattice_0
    println("   | For Metropolis algorithm ...")
    E_arr[1] = potts_hamiltonian(potts) / potts.L^2
    for step=2:nsteps
        delE = metropolis_batch_update!(potts, T)
        E_arr[step] = E_arr[step-1] + (delE / potts.L^2)
    end
    corrfn = autocorrelation_fn(E_arr)
    τ = sum(corrfn[1:500])
    # println(sum(isnan.(corrfn[1:200])))
    auto_corr_times_metro[i] = convert(Int64, ceil(τ))

    # For Wolff
    potts.lattice = lattice_0
    println("   | For Wolff algorithm ...")
    E_arr[1] = potts_hamiltonian(potts) / potts.L^2
    for step=2:nsteps
        wolff_cluster_update!(potts, T)
        E_arr[step] = potts_hamiltonian(potts) / potts.L^2
    end
    corrfn = autocorrelation_fn(E_arr)
    τ = sum(corrfn[1:500])
    auto_corr_times_wolff[i] = convert(Int64, ceil(τ))
    println("   |          ")
    println("   +-> Done.\n")
end

println("Using Metropolis: $(auto_corr_times_metro)")
println("Using Wolff: $(auto_corr_times_wolff)")

#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "integrated autocorrelation time",
    title = "PottsModel2D $(L)x$(L):  Metropolis Algorithm")

ax2 = Axis(f[2, 1], xlabel = "temperature, T", ylabel = "integrated autocorrelation time",
    title = "PottsModel2D $(L)x$(L): Wolff Cluster Algorithm")

scatterlines!(
    ax1, Temps, auto_corr_times_metro,
    markersize = 7
)

scatterlines!(
    ax2, Temps, auto_corr_times_wolff,
    markersize = 7
)

save("plots/2Dmodel/integrated_auto_corr_time_comparison_$(L).pdf", f)

println("Program Finished!")
println("===========================\n")