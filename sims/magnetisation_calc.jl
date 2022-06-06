using CairoMakie
using Statistics

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

Temps = [0.5, 0.7, 0.9, 0.98, 1.0, 1.02, 1.1, 1.3, 1.5]
esteps = 1000  # Number of steps for equilibration
nsteps = 20000  # Number of steps for measurements

m_T = zeros(Float64, length(Temps))  # Array of magnetisation per site
err_m_T = zeros(Float64, length(Temps))

χ_T = zeros(Float64, length(Temps))  # Array of specific heat
err_χ_T = zeros(Float64, length(Temps))

potts = initialize_model_2d(L, q; cold_start=true)

for i = 1:length(Temps)
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    m_arr = zeros(Float64, nsteps)

    println("   | Making measurements ...")
    for step=1:esteps
        metropolis_batch_update!(potts, T)
        # wolff_cluster_update!(potts, T)
    end

    m_arr[1] = magnetisation(potts) / potts.L^2
    for step=2:nsteps
        metropolis_batch_update!(potts, T)
        # wolff_cluster_update!(potts, T)

        m_arr[step] = magnetisation(potts) / potts.L^2
    end

    println("   | Calculating observables ...")
    m_T[i] = mean(m_arr) 
    err_m_T[i] = blocking_err(m_arr, A -> mean(A); blocks=50)

    χ_T[i] = succeptibility(m_arr, T, potts.L^potts.d)
    err_χ_T[i] = blocking_err(m_arr, specific_heat, T, potts.L^potts.d; blocks=50)
    println("   |          ")
    println("   +-> Done.\n")
end


#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "magnetisation, u",
    title = "PottsModel2D $(L)x$(L): Magnetisation (per site) v/s temperature")

ax2 = Axis(f[2, 1], xlabel = "temperature, T", ylabel = "succeptibility, χ",
    title = "PottsModel2D $(L)x$(L): Succeptibility (per site) v/s temperature")

errorbars!(
    ax1, Temps, m_T, err_m_T,
    whiskerwidth = 10
)
scatter!(
    ax1, Temps, m_T,
    markersize = 7
)

errorbars!(
    ax2, Temps, χ_T, err_χ_T,
    whiskerwidth = 10)
scatter!(
    ax2, Temps, χ_T,
    markersize = 7
)

save("plots/2Dmodel/mag_and_suzz_vs_Temp_$(L)_metro.pdf", f)

println("Program Finished!")
println("===========================\n")