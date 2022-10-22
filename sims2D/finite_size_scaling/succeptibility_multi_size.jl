include("../../src/pottsmc.jl")
using CairoMakie

Lvals = [32, 48, 64, 80]
cols = Dict([(32, :blue), (48, :red), (64, :green), (80, :purple)])
temps = [0.900, 0.920, 0.940, 0.960, 0.980, 0.984, 0.988, 0.992, 0.996, 1.000, 1.004,
         1.008, 1.012, 1.016, 1.020, 1.024, 1.028, 1.032, 1.036, 1.040, 1.060, 1.080, 1.100]

mags_def = 1
suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

mean_mags = zeros(Float64, (length(Lvals), length(temps)))
err_mean_mags = zeros(Float64, (length(Lvals), length(temps)))
suzz = zeros(Float64, (length(Lvals), length(temps)))
err_suzz = zeros(Float64, (length(Lvals), length(temps)))

for Lidx in eachindex(Lvals)
    L = Lvals[Lidx]
    for tidx in eachindex(temps)
        T = temps[tidx]
        mags = readdlm(joinpath("data", "2DModel", "Size$(L)", "mags", "potts_mags_temp$(T)_L$(L).txt"), ',', Float64)
        mags ./= L^2

        @views mean_mags[Lidx, tidx] = mean(mags[mags_def, :])
        @views err_mean_mags[Lidx, tidx] = bootstrap_err(mags[mags_def, :], A -> mean(A); r=100)
        @views suzz[Lidx, tidx] = suzz_kth(mags[mags_def, :], T, L^2, 2)
        @views err_suzz[Lidx, tidx] = bootstrap_err(mags[mags_def, :], A -> suzz_kth(A, T, L^2, 2); r=100)
    end
end

f1 = Figure();
f2 = Figure();

ax1 = Axis(f1[1, 1], xlabel = "temperature, T", ylabel = "magnetization (per site), m",
    title = "PottsModel2D: magnetization (per site) v/s temperature");

ax2 = Axis(f2[1, 1], xlabel = "temperature, T", ylabel = "susceptibility, χ",
    title = "PottsModel2D: susceptibility (per site) v/s temperature");

for Lidx in eachindex(Lvals)
    L = Lvals[Lidx]

    errorbars!(ax1, temps, mean_mags[Lidx, :], err_mean_mags[Lidx, :]; color=:black, whiskerwidth=14)
    errorbars!(ax2, temps, suzz[Lidx, :], err_suzz[Lidx, :]; color=:black, whiskerwidth=14)

    scatterlines!(ax1, temps, mean_mags[Lidx, :]; color=cols[L], label="L = $L")
    scatterlines!(ax2, temps, suzz[Lidx, :]; color=cols[L], label="L = $L")
end

axislegend(ax1)
axislegend(ax2)
display(f1)
display(f2)

save("plots/2DModel/finite_size_scaling/mag_full_plot.svg", f1)
save("plots/2DModel/finite_size_scaling/suzz_full_plot.svg", f2)