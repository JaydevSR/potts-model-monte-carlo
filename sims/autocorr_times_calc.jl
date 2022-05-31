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
println("================================\n")
println("    Lattice Size: $(L) x $(L)")
println("================================\n")

# Temps = [1.0]
Temps = [0.5, 0.7, 0.9, 0.98, 1.0, 1.02, 1.1, 1.3, 1.5]
nsteps = 20000  # Number of steps for measurements
wsteps = 500

auto_corr_times_metro = zeros(Int64, length(Temps))
auto_corr_times_wolff = zeros(Int64, length(Temps))

potts = initialize_model_2d(L, q; cold_start=true)
lattice_0 = copy(potts.lattice)

data_mat = zeros(Float64, (3, length(Temps)))
data_mat[1, :] = copy(Temps)

for i = 1:length(Temps)
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    m_arr = zeros(Float64, nsteps)

    # For Metropolis
    potts.lattice = copy(lattice_0)
    println("   | For Metropolis algorithm ...")

    print("   |   | Progress: [")
    for step=1:nsteps
        if step%1000==0
            print("*")
        end
        m_arr[step] = magnetisation(potts) / potts.L^2
        delE = metropolis_batch_update!(potts, T)
    end
    print("]\n")

    corrfn = autocorrelation_fn(m_arr)
    τ = sum(corrfn[1:wsteps])
    # println(sum(isnan.(corrfn[1:200])))
    auto_corr_times_metro[i] = convert(Int64, ceil(τ))
    println("   |   | τ = $(auto_corr_times_metro[i])")

    # For Wolff
    potts.lattice = copy(lattice_0)
    println("   | For Wolff algorithm ...")

    print("   |   | Progress: [")
    for step=1:nsteps
        if step%1000==0
            print("*")
        end
        m_arr[step] = magnetisation(potts) / potts.L^2
        wolff_cluster_update!(potts, T)
    end
    print("]\n")

    corrfn = autocorrelation_fn(m_arr)
    τ = sum(corrfn[1:wsteps])
    # scale value for comparison with Metropolis algorithm
    if T > 1.0
        τ = τ * (succeptibility(m_arr, T, potts.L^potts.d) * T / (potts.L^potts.d))
    end
    auto_corr_times_wolff[i] = convert(Int64, ceil(τ))
    println("   |   | τ = $(auto_corr_times_wolff[i])")
    println("   +----> Done.\n")
end

println("Using Metropolis: $(auto_corr_times_metro)")
println("Using Wolff: $(auto_corr_times_wolff)")
data_mat[2, :] = auto_corr_times_metro
data_mat[3, :] = auto_corr_times_wolff

open("D:\\Projects\\Potts-QCD\\potts-model-monte-carlo\\data\\corr_data.txt", "w") do io
    writedlm(io, data_mat, ',')
end;

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