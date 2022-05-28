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

Temps = [0.5, 0.98, 1.02, 1.5]
nsteps = 10000  # Number of steps for measurements

for i=1:length(Temps)
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    f = Figure(resolution = (1200, 800), title="3 state potts 2D, Lattice size = $(L)x$(L), temperature=$(T)")
    ax1 = Axis(f[1, 1:2], xlabel="step", ylabel="magnetisation", title="metropolis", xticks=LinearTicks(10))
    ax2 = Axis(f[1, 3], xlabel="step", ylabel="magnetisation", title="metropolis", xticks=LinearTicks(10))
    ax3 = Axis(f[2, 1:2], xlabel="step", ylabel="magnetisation", title="wolff", xticks=LinearTicks(10))
    ax4 = Axis(f[2, 3], xlabel="step", ylabel="magnetisation", title="wolff", xticks=LinearTicks(10))

    potts2d = initialize_model_2d(L, q; cold_start=false)
    hot_save = copy(potts2d.lattice)
    m_arr = zeros(Float64, nsteps)

    # metropolis_batch_update
    println("   | Metropolis algorithm ...")
    println("   |   | Hot start ...")
    for step=1:nsteps
        m_arr[step] = magnetisation(potts2d) / potts2d.L^2
        metropolis_batch_update!(potts2d, T)
    end
    scatterlines!(ax1, m_arr, color=:red, label="hot start", markersize=3, markercolor=:red)
    scatterlines!(ax2, m_arr[1:1000], color=:red, label="hot start", markersize=3, markercolor=:red)

    println("   |   | Cold start ...")
    potts2d.lattice = zero(potts2d.lattice)
    for step=1:nsteps
        m_arr[step] = magnetisation(potts2d) / potts2d.L^2
        metropolis_batch_update!(potts2d, T)
    end
    scatterlines!(ax1, m_arr, color=:blue, label="cold start", markersize=3, markercolor=:blue)
    scatterlines!(ax2, m_arr[1:1000], color=:blue, label="cold start", markersize=3, markercolor=:blue)
    println("   | Done.\n   |")

    # wolff_cluster_update
    println("   | Wolff algorithm ...")
    println("   |   | Hot start ...")
    potts2d.lattice = copy(hot_save)
    for step=1:nsteps
        m_arr[step] = magnetisation(potts2d) / potts2d.L^2
        wolff_cluster_update!(potts2d, T)
    end
    scatterlines!(ax3, m_arr, color=:red, label="hot start", markersize=3, markercolor=:red)
    scatterlines!(ax4, m_arr[1:1000], color=:red, label="hot start", markersize=3, markercolor=:red)

    println("   |   | Cold start ...")
    potts2d.lattice = zero(potts2d.lattice)
    for step=1:nsteps
        m_arr[step] = magnetisation(potts2d) / potts2d.L^2
        wolff_cluster_update!(potts2d, T)
    end
    scatterlines!(ax3, m_arr, color=:blue, label="cold start", markersize=3, markercolor=:blue)
    scatterlines!(ax4, m_arr[1:1000], color=:blue, label="cold start", markersize=3, markercolor=:blue)
    println("   | Done.\n   |")

    println("   | Saving plots ...")
    axislegend(ax1)
    axislegend(ax2)
    axislegend(ax3)
    axislegend(ax4)

    save("plots/2Dmodel/equilibration_and_start_temp$(T)_size$(L).png", f)
end