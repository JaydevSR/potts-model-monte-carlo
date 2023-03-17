include("../../src/pottsmc.jl")
using CairoMakie

# load data
lattice_sizes = [48, 56, 64, 72, 80]
cols = Dict([(32, :blue), (48, :red), (56, :pink), (64, :green), (72, :orange), (80, :purple), (96, :cyan), (128, :deepskyblue)])

temps = reshape(readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), ',', Float64), :)
max_order = 7

figures = [Figure(resolution=(800, 600)) for i in 1:max_order]
axes = [Axis(figures[1][1, 1], xlabel="T", ylabel="m(T)", title="Mean Magnetization Per Site")]
append!(axes, [
    Axis(
        figures[i][1, 1],
        xlabel="T", ylabel="χₙ(T)",
        title="Order $(i-1) Susceptibility Per Site"
    ) for i in 2:max_order]
)

for i in 2:max_order
    xlims!(axes[i], (0.98, 1.02))
end

# plot susceptibilities
for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    suzzs = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_L$(L).txt"), '\t', Float64)
    errors = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_errors_L$(L).txt"), '\t', Float64)

    for order_k in 1:max_order
        ax = axes[order_k]
        errorbars!(ax, temps, suzzs[order_k, :], errors[order_k, :], color=:black, whiskerwidth=12, label="L = $L")
        scatterlines!(ax, temps, suzzs[order_k, :], color = cols[L], markersize=15, label="L = $L")
    end
end

axislegend(axes[1], position=:rt, merge=true)
for idx in 2:max_order
    axislegend(axes[idx], position=:lt, merge=true)
end

save(joinpath("plots", "2DModel", "final_plots", "potts_magnetisation.svg"), figures[1])
for idx in 2:max_order
    save(joinpath("plots", "2DModel", "final_plots", "potts_suzzs_order$(idx-1).svg"), figures[idx])
end