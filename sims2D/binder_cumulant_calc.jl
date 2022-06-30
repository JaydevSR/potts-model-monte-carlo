include("../src/pottsmc.jl")
using CairoMakie

Lvals = [16, 24, 32, 40, 48]
cols = Dict([(16, :blue), (24, :red), (32, :green), (40, :purple), (48, :black)])
temps = [0.98, 0.984, 0.990, 0.996, 1.000]
basepath = joinpath(["data", "magdata"])
err_nblocks = 20

f1 = Figure()

ax1 = Axis(f1[1, 1], xlabel = "T", ylabel = "U",
    title = "PottsModel2D: Binder's cumulant v/s temperature")

for L in Lvals
    binder = zeros(Float64, length(temps))  # Array of magnetisation per site

    for i in eachindex(temps)
        T = temps[i]
        loc = joinpath([basepath, "Size$L", "potts_mags_temp$(T)_L$(L).txt"])
        m_arr = readdlm(loc, ',', Float64) ./ L^2

        binder[i] = binders_cumulant(m_arr) 
        # err_mags[i] = blocking_err(m_arr, A -> binders_cumulant(A); blocks=err_nblocks)
    end

    # errorbars!(ax1, temps, mags, err_mags, whiskerwidth = 10)
    scatterlines!(ax1, temps, binder ,markersize = 7, label="Size=$(L)x$(L)", linestyle = :dot, color=cols[L])
end

axislegend(ax1)
display(f1)
