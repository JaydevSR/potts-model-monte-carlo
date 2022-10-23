include("../../src/pottsmc.jl")
using CairoMakie

lattice_sizes = [32, 48, 64, 80]
cols = Dict([(32, :blue), (48, :red), (64, :green), (80, :purple)])
temps = [0.900, 0.920, 0.940, 0.960, 0.980, 0.982, 0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996, 0.998,
        1.000, 1.002, 1.004, 1.006, 1.008, 1.010, 1.012, 1.014, 1.016, 1.018, 1.020, 1.022, 1.024, 1.026,
        1.028, 1.030, 1.032, 1.034, 1.036, 1.038, 1.040, 1.060, 1.080, 1.100]

mags_def = 1

binder = zeros(Float64, (length(lattice_sizes), length(temps)))
err_binder = zeros(Float64, (length(lattice_sizes), length(temps)))

for Lidx in eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]
    Threads.@threads for tidx in eachindex(temps)
        T = temps[tidx]
        mags = readdlm(joinpath("data", "2DModel", "Size$(L)", "mags", "potts_mags_temp$(T)_L$(L).txt"), ',', Float64)
        mags ./= L^2

        @views binder[Lidx, tidx] = binders_cumulant(mags[mags_def, :]) 
        @views err_binder[Lidx, tidx] = bootstrap_err(mags[mags_def, :], A -> binders_cumulant(A); r=100)
    end
end

# plots
trange = 10:14
f = Figure(resolution=(800, 600));

ax1 = Axis(f[1, 1], xlabel = "T", ylabel = "U",
    title = "PottsModel2D: Binder's cumulant v/s temperature")
ax2 = Axis(f[2, 1], xlabel = "T", ylabel = "U")

for Lidx in eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]
    errorbars!(ax1, temps[:], binder[Lidx, :], err_binder[Lidx, :], whiskerwidth = 12)
    scatterlines!(ax1, temps[:], binder[Lidx, :], markersize = 10, label="Size=$(L)x$(L)", linestyle = :dot, color=cols[L])
    errorbars!(ax2, temps[trange], binder[Lidx, trange], err_binder[Lidx, trange], whiskerwidth = 12)
    scatterlines!(ax2, temps[trange], binder[Lidx, trange], markersize = 10, label="Size=$(L)x$(L)", linestyle = :dot, color=cols[L])
end

axislegend(ax1)
display(f)
save("plots/2Dmodel/cumulants/binders_cumulant.svg", f)