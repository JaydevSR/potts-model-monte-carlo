include("../../src/pottsmc.jl")
using CairoMakie

# load data
lattice_sizes = [48, 56, 64, 72, 80]
cols = Dict([(32, :blue), (48, :red), (56, :pink), (64, :green), (72, :orange), (80, :purple), (96, :cyan), (128, :deepskyblue)])

temps = reshape(readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), ',', Float64), :)
max_order = 6

figures = [Figure(resolution = (1000, 800), fontsize = 28) for i in 1:max_order]
axes = [
    Axis(
        figures[1][1, 1],
        xlabel=L"T", ylabel=L"\langle m\rangle",
        title=L"\langle m \rangle = \frac{1}{L^2\beta} \left[\frac{\partial \ln Z}{\partial h}\right]_{h \rightarrow 0}",
        xticks = temps[1:2:end],
        xlabelsize = 46, ylabelsize = 46,
        xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
    ),

    Axis(
        figures[2][1, 1],
        xlabel=L"T", ylabel=L"\chi",
        title=L"\chi = \frac{1}{L^2\beta} \left[\frac{\partial^{2} \ln Z}{\partial h^2}\right]_{h \rightarrow 0}",
        xticks = temps[1:2:end],
        xlabelsize = 46, ylabelsize = 46,
        xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
    )
]

append!(axes, [
    Axis(
        figures[i][1, 1],
        xlabel = L"T", ylabel = L"\chi_{%$i}",
        title = L"\chi_{%$i} = \frac{1}{L^2\beta} \left[\frac{\partial^{%$i} \ln Z}{\partial h^{%$i}}\right]_{h \rightarrow 0}",
        xticks = temps[1:2:end],
        xlabelsize = 46, ylabelsize = 46,
        xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
    ) for i in 3:max_order
])

for i in 1:max_order
    xlims!(axes[i], (0.982, 1.022))
end

# plot susceptibilities
for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    suzzs = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_L$(L).txt"), '\t', Float64)
    errors = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_errors_L$(L).txt"), '\t', Float64)

    for order_k in 1:max_order
        ax = axes[order_k]
        errorbars!(ax, temps, suzzs[order_k, :], errors[order_k, :], color=:black, whiskerwidth=12, label="L = $L")
        scatterlines!(ax, temps, suzzs[order_k, :],
            # color = cols[L],
            markersize=25, linewidth=5, 
            label="L = $L")
    end
end

for idx in 1:max_order
    Legend(figures[idx][1, 2], axes[idx], merge = true)
end

save(joinpath("plots", "2DModel", "final_plots", "potts_magnetisation.svg"), figures[1])
save(joinpath("plots", "2DModel", "final_plots", "potts_susceptibility.svg"), figures[2])
for idx in 3:max_order
    save(joinpath("plots", "2DModel", "final_plots", "potts_susceptibility_order$(idx).svg"), figures[idx])
end