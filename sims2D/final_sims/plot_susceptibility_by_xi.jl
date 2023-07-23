include("../../src/pottsmc.jl")
using CairoMakie

# load data
lattice_sizes = [32, 48, 64, 128]

γ = 13//9
ν = 5//6

temps = reshape(readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), ',', Float64), :)

f_L_combined = Figure(resolution = (1000, 800), fontsize = 32)
f_ξ_combined = Figure(resolution = (1000, 800), fontsize = 32)

ax_L_combined = Axis(
    f_L_combined[1, 1],
    xlabel=L"T", ylabel=L"\frac{\chi}{L^{\gamma/\nu}}",
    xticks = temps[1:2:end],
    xlabelsize = 46, ylabelsize = 46,
    xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
    ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
)
ax_ξ_combined = Axis(
    f_ξ_combined[1, 1],
    xlabel=L"T", ylabel=L"\frac{\chi}{ξ^\gamma}",
    xticks = temps[1:2:end],
    xlabelsize = 46, ylabelsize = 46,
    xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
    ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
)

for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    suzzs = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_L$(L).txt"), '\t', Float64)
    errors = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_errors_L$(L).txt"), '\t', Float64)

    suzz_2 = suzzs[2, :]
    error_2 = errors[2, :]

    f_L = Figure(resolution = (1000, 800), fontsize = 32)
    f_ξ = Figure(resolution = (1000, 800), fontsize = 32)

    ax_L = Axis(
        f_L[1, 1],
        xlabel=L"T", ylabel=L"\frac{\chi}{L^{\gamma/\nu}}",
        title=L"$%$L \times %$L$ lattice",
        xticks = temps[1:2:end],
        xlabelsize = 46, ylabelsize = 46,
        xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
    )
    ax_ξ = Axis(
        f_ξ[1, 1],
        xlabel=L"T", ylabel=L"\frac{\chi}{ξ^\gamma}",
        title=L"$%$L \times %$L$ lattice",
        xticks = temps[1:2:end],
        xlabelsize = 46, ylabelsize = 46,
        xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23, xticklabelrotation = pi/4
    )

    # xlims!(ax_L, (0.982, 1.022))
    # xlims!(ax_ξ, (0.982, 1.022))

    ratio_L = (suzz_2 .± error_2) / L^(γ/ν)
    ratio_L_vals = Measurements.value.(ratio_L)
    ratio_L_errs = Measurements.uncertainty.(ratio_L)
    errorbars!(ax_L, temps, ratio_L_vals, ratio_L_errs, color=:black, whiskerwidth=20, linewidth=2, label="ratio")
    scatterlines!(ax_L, temps, ratio_L_vals, markersize=25, linewidth=3.5, linestyle=:dot, label="ratio")

    errorbars!(ax_L_combined, temps, ratio_L_vals, ratio_L_errs, color=:black, whiskerwidth=20, linewidth=3, label="L = $L")
    scatterlines!(ax_L_combined, temps, ratio_L_vals, markersize=25, linewidth=3.5, linestyle=:dashdot, label="L = $L")

    ξ_vals = reshape(readdlm(joinpath("data", "2DModel", "correlations", "potts_corrlen_values_L$(L).txt"), '\t', Float64), :)
    ξ_errs = reshape(readdlm(joinpath("data", "2DModel", "correlations", "potts_corrlen_errors_L$(L).txt"), '\t', Float64), :)

    ratio_ξ = (suzz_2 .± error_2) ./ ((ξ_vals .± ξ_errs).^γ)

    ratio_ξ_vals = Measurements.value.(ratio_ξ)
    ratio_ξ_errs = Measurements.uncertainty.(ratio_ξ)
    errorbars!(ax_ξ, temps, ratio_ξ_vals, ratio_ξ_errs, color=:black, whiskerwidth=12, linewidth=3, label="ratio")
    scatterlines!(ax_ξ, temps, ratio_ξ_vals, markersize=25, linewidth=3.5, linestyle=:dot, label="ratio")

    errorbars!(ax_ξ_combined, temps, ratio_ξ_vals, ratio_ξ_errs, color=:black, whiskerwidth=12, linewidth=3, label="L = $L")
    scatterlines!(ax_ξ_combined, temps, ratio_ξ_vals, markersize=25, linewidth=3.5, linestyle=:dashdot, label="L = $L")

    # axislegend(ax_L; position=:rt, merge=true)
    # axislegend(ax_ξ; position=:rt, merge=true)

    save("plots/2Dmodel/correlations/potts_susceptibility_L_gamma_nu_L$L.svg", f_L)
    save("plots/2Dmodel/correlations/potts_susceptibility_ξ_gamma_L$L.svg", f_ξ)
end

axislegend(ax_L_combined; position=:rt, merge=true)
axislegend(ax_ξ_combined; position=:rt, merge=true)

save("plots/2Dmodel/correlations/potts_susceptibility_L_gamma_nu_combined.svg", f_L_combined)
save("plots/2Dmodel/correlations/potts_susceptibility_ξ_gamma_combined.svg", f_ξ_combined)