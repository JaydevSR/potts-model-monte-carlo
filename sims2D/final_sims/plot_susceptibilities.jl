include("../../src/pottsmc.jl")
using CairoMakie

# load data
lattice_sizes = [48, 56, 64, 80, 72, 96, 128]
cols = Dict([(128, :steelblue3), (96, :mediumseagreen), (80, :orange), (72, :darkcyan), (64, :orchid3), (56, :green), (48, :purple), (32, :pink)])

temps = reshape(readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), ',', Float64), :)
max_order = 6

figures = [Figure(resolution = (1000, 800), fontsize = 32) for i in 1:max_order]
axes = [
    Axis(
        figures[1][1, 1],
        xlabel=L"T (J/k_B)", ylabel=L"\mathcal{O}",
        # title=L"\langle m \rangle = \frac{1}{L^2\beta} \left[\frac{\partial \ln Z}{\partial h}\right]_{h \rightarrow 0}",
        xticks = temps[1:2:end],
        xlabelsize = 46, ylabelsize = 46,
        xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
    ),

    Axis(
        figures[2][1, 1],
        xlabel=L"T (J/k_B)", ylabel=L"\chi",
        # title=L"\chi = \frac{1}{L^2\beta} \left[\frac{\partial^{2} \ln Z}{\partial h^2}\right]_{h \rightarrow 0}",
        xticks = temps[1:2:end],
        xlabelsize = 46, ylabelsize = 46,
        xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
    )
]

append!(axes, [
    Axis(
        figures[i][1, 1],
        xlabel = L"T (J/k_B)", ylabel = L"\chi_{%$i}",
        # title = L"\chi_{%$i} = \frac{1}{L^2\beta} \left[\frac{\partial^{%$i} \ln Z}{\partial h^{%$i}}\right]_{h \rightarrow 0}",
        xticks = temps[1:2:end],
        xlabelsize = 46, ylabelsize = 46,
        xgridstyle = :dot, xgridwidth = 1.1, xgridcolor = :gray23,
        ygridstyle = :dot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
    ) for i in 3:max_order
])

for i in 1:max_order
    xlims!(axes[i], (0.982, 1.022))
    lines!(axes[i], temps, zero.(temps), linestyle=:dash, color=:black,
                linewidth=4)
end

# plot susceptibilities
for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    suzzs = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_L$(L).txt"), '\t', Float64)
    errors = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_errors_L$(L).txt"), '\t', Float64)

    for order_k in 1:max_order
        if order_k <= 4 || L != 128 # to not plot 128 on higher orders
            ax = axes[order_k]
            errorbars!(ax, temps, suzzs[order_k, :], errors[order_k, :],
                color=:black, whiskerwidth=16, linewidth=3, label="L = $L")
            scatterlines!(ax, temps, suzzs[order_k, :],
                color = cols[L],
                markersize=25, linewidth=5, 
                label="L = $L")
        end
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