include("../../src/pottsmc.jl")
using CairoMakie

# C3/C1 > C4/C2 > C5/C1 > C6/C2

# load data
lattice_sizes = [48, 64, 80, 96, 128]
xranges = Dict(
    48 => (0.98, 1.025),
    64 => (0.98, 1.025),
    80 => (0.98, 1.025),
    96 => (0.98, 1.025),
    128 => (0.98, 1.025)
)

temps = reshape(readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), ',', Float64), :)
max_order = 6

order_colors = Dict(:C3C1 => :dodgerblue2, :C4C2 => :orange, :C5C1 => :mediumseagreen, :C6C2 => :orchid3)
markersize = 18
linewidth = 4
marker = :diamond

psuedoTc = Dict(
            48 => 1.004098703447542,
            64 => 1.0009103256327856,
            80 => 0.9997345330528584,
            96 => 0.9988572671581376,
            128 => 0.997941002359664)

# Plot Ratios
for stepL in eachindex(lattice_sizes)    
    L = lattice_sizes[stepL]

    println("Calculation for size: $L ...")

    suzzs = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_L$(L).txt"), '\t', Float64)
    errors = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_errors_L$(L).txt"), '\t', Float64)

    suzzs_L = suzzs .± errors

    # psuedoTc[L] = temps[findmax(suzzs[2, :])[2]]
    temps_scaled = temps ./ psuedoTc[L]

    temps_scaled = round.(temps_scaled, digits=3)

    fig = Figure(resolution=(900, 600), fontsize = 19)
    axes = [
        Axis(fig[1, 1],
            title=L"\text{Comparison between } \chi_4/\chi_2 \text{ and } \chi_3/\chi_1 \text{ for L=%$L}",
            xticks = temps_scaled[1:2:end], yticks = WilkinsonTicks(5),
            xlabelsize = 24, ylabelsize = 24,
            xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
            ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23),

        Axis(fig[2, 1], ylabel="Ratio",
            title=L"\text{Comparison between } \chi_5/\chi_1 \text{ and } \chi_4/\chi_2 \text{ for L=%$L}",
            xticks = temps_scaled[1:2:end], yticks = WilkinsonTicks(5),
            xlabelsize = 24, ylabelsize = 24,   
            xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
            ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23),

        Axis(fig[3, 1], xlabel=L"T / T_c(L)",
            title=L"\text{Comparison between } \chi_6/\chi_2 \text{ and } \chi_5/\chi_1 \text{ for L=%$L}",
            xticks = temps_scaled[1:2:end], yticks = WilkinsonTicks(5),
            xlabelsize = 24, ylabelsize = 24,
            xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
            ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23)
    ]

    for ax in axes
        hlines!(ax, [0], linewidth = 1.5, color = :red)
    end

    # χ_3 / χ_1    
    r31 = suzzs_L[3, :] ./ suzzs_L[1, :]
    r31_values = Measurements.value.(r31)
    r31_errors = Measurements.uncertainty.(r31)

    # band!(axes[1], temps_scaled, r31_values .- r31_errors, r31_values .+ r31_errors,
    #     color = (order_colors[:C3C1], 0.25), label=L"\chi_3 / \chi_1")

    errorbars!(axes[1], temps_scaled, r31_values, r31_errors, label=L"\chi_3 / \chi_1", whiskerwidth=12, linewidth=2)

    scatterlines!(axes[1], temps_scaled, r31_values, label=L"\chi_3 / \chi_1",
        color=order_colors[:C3C1], marker=marker, markersize=markersize, linewidth=linewidth)

    # χ_4 / χ_2
    r42 = suzzs_L[4, :] ./ suzzs_L[2, :]
    r42_values = Measurements.value.(r42)
    r42_errors = Measurements.uncertainty.(r42)

    # band!(axes[1], temps_scaled, r42_values .- r42_errors, r42_values .+ r42_errors,
    #     color = (order_colors[:C4C2], 0.25), label=L"\chi_4 / \chi_2")
    # band!(axes[2], temps_scaled, r42_values .- r42_errors, r42_values .+ r42_errors,
    #     color = (order_colors[:C4C2], 0.25), label=L"\chi_4 / \chi_2")

    errorbars!(axes[1], temps_scaled, r42_values, r42_errors, label=L"\chi_4 / \chi_2", whiskerwidth=12, linewidth=2)
    errorbars!(axes[2], temps_scaled, r42_values, r42_errors, label=L"\chi_4 / \chi_2", whiskerwidth=12, linewidth=2)

    scatterlines!(axes[1], temps_scaled, r42_values, label=L"\chi_4 / \chi_2",
        color=order_colors[:C4C2], marker=marker, markersize=markersize, linewidth=linewidth)
    scatterlines!(axes[2], temps_scaled, r42_values, label=L"\chi_4 / \chi_2",
        color=order_colors[:C4C2], marker=marker, markersize=markersize, linewidth=linewidth)

    # χ_5 / χ_1
    r51 = suzzs_L[5, :] ./ suzzs_L[1, :]
    r51_values = Measurements.value.(r51)
    r51_errors = Measurements.uncertainty.(r51)

    # band!(axes[2], temps_scaled, r51_values .- r51_errors, r51_values .+ r51_errors,
    #     color = (order_colors[:C5C1], 0.25), label=L"\chi_5 / \chi_1")
    # band!(axes[3], temps_scaled, r51_values .- r51_errors, r51_values .+ r51_errors,
    #     color = (order_colors[:C5C1], 0.25), label=L"\chi_5 / \chi_1")

    errorbars!(axes[2], temps_scaled, r51_values, r51_errors, label=L"\chi_5 / \chi_1", whiskerwidth=12, linewidth=2)
    errorbars!(axes[3], temps_scaled, r51_values, r51_errors, label=L"\chi_5 / \chi_1", whiskerwidth=12, linewidth=2)

    scatterlines!(axes[2], temps_scaled, r51_values,
        label=L"\chi_5 / \chi_1", color=order_colors[:C5C1], marker=marker, markersize=markersize, linewidth=linewidth)
    scatterlines!(axes[3], temps_scaled, r51_values,
        label=L"\chi_5 / \chi_1", color=order_colors[:C5C1], marker=marker, markersize=markersize, linewidth=linewidth)

    # χ_6 / χ_2
    r62 = suzzs_L[6, :] ./ suzzs_L[2, :]
    r62_values = Measurements.value.(r62)
    r62_errors = Measurements.uncertainty.(r62)

    # band!(axes[3], temps_scaled, r62_values .- r62_errors, r62_values .+ r62_errors,
    #     color = (order_colors[:C6C2], 0.25), label=L"\chi_6 / \chi_2")

    errorbars!(axes[3], temps_scaled, r62_values, r62_errors, label=L"\chi_6 / \chi_2", whiskerwidth=12, linewidth=2)

    scatterlines!(axes[3], temps_scaled, r62_values,
        label=L"\chi_6 / \chi_2", color=order_colors[:C6C2], marker=marker, markersize=markersize, linewidth=linewidth)

    plots_in_fig = AbstractPlot[]
    labels_in_fig = AbstractString[]
    for ax in axes
        pl, lb = Makie.get_labeled_plots(ax, merge=false, unique=false)
        append!(plots_in_fig, pl)
        append!(labels_in_fig, lb)
        xlims!(ax, xranges[L])
    end

    ulabels = Base.unique(labels_in_fig)
    mergedplots = [[lp for (i, lp) in enumerate(plots_in_fig) if labels_in_fig[i] == ul]
            for ul in ulabels]

    Legend(fig[:, 2], mergedplots, ulabels)

    save(joinpath("plots", "2DModel", "final_plots", "higher_order_susceptibility_ratios_size$(L).svg"), fig)
end

println(psuedoTc)

# f = Figure(resolution=(1000, 600));
# axes = [
#     Axis(f[1, 1], xlabel=L"T", ylabel="Ratio", title="Ratio of odd orders for L=$L"),
#     Axis(f[1, 2], xlabel=L"T", ylabel="Ratio", title="Ratio of even orders for L=$L")
# ]

# is_odd_pair(x) = x[1] > x[2] && isodd(x[1]) && isodd(x[2])
# is_even_pair(x) = x[1] > x[2] && iseven(x[1]) && iseven(x[2])

# pair_range = [(i, j) for i=1:max_order-1, j=1:max_order]

# all_odd_pairs = filter(is_odd_pair, pair_range)
# all_even_pairs = filter(is_even_pair, pair_range)

# for (oh, ol) in all_odd_pairs
#     ratio = suzzs_L[oh, :] ./ suzzs_L[ol, :]
#     ratio_values = Measurements.value.(ratio)
#     ratio_errors = Measurements.uncertainty.(ratio)

#     band!(axes[1], temps, ratio_values .- ratio_errors, ratio_values .+ ratio_errors,
#         color = (order_colors[oh-ol], 0.25), label=L"\chi^{%$oh} / \chi^{%$ol}")

#     errorbars!(axes[1], temps, ratio_values, ratio_errors, label=L"\chi^{%$oh} / \chi^{%$ol}", whiskerwidth=12)

#     scatterlines!(axes[1], temps, ratio_values, label=L"\chi^{%$oh} / \chi^{%$ol}", color=order_colors[oh-ol])
# end

# for (oh, ol) in all_even_pairs
#     ratio = suzzs_L[oh, :] ./ suzzs_L[ol, :]
#     ratio_values = Measurements.value.(ratio)
#     ratio_errors = Measurements.uncertainty.(ratio)

#     band!(axes[2], temps, ratio_values .- ratio_errors, ratio_values .+ ratio_errors,
#         color = (order_colors[oh-ol], 0.25), label=L"\chi^{%$oh} / \chi^{%$ol}")

#     errorbars!(axes[2], temps, ratio_values, ratio_errors, label=L"\chi^{%$oh} / \chi^{%$ol}", whiskerwidth=12)

#     scatterlines!(axes[2], temps, ratio_values, label=L"\chi^{%$oh} / \chi^{%$ol}", color=order_colors[oh-ol])
# end

# axislegend(axes[1], position=:rt, merge=true)
# axislegend(axes[2], position=:rt, merge=true)

# # display(f)
# save(joinpath("plots", "2DModel", "final_plots", "higher_order_susceptibility_ratios_size$(L).svg"), f)