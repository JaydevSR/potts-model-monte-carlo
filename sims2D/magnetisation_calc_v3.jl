include("../src/pottsmc.jl")
using CairoMakie

Lvals = [16, 24, 32, 40]#, 48]
cols = Dict([(16, :blue), (24, :red), (32, :green), (40, :purple), (48, :black)])
temps = [0.5, 0.7, 0.9, 0.94, 0.98, 0.984, 0.988, 0.992, 0.996, 1.0, 1.004, 1.008, 1.012, 1.016, 1.02, 1.06,  1.1, 1.3, 1.5]
basepath = joinpath(["data", "uncorr_configs"])
bootstrap_samples = 100

f1 = Figure()
f2 = Figure()

ax1 = Axis(f1[1, 1], xlabel = "temperature, T", ylabel = "magnetization (per site), m",
    title = "PottsModel2D: magnetization (per site) v/s temperature")

ax2 = Axis(f2[1, 1], xlabel = "temperature, T", ylabel = "succeptibility, χ",
    title = "PottsModel2D: Succeptibility (per site) v/s temperature")

for L in Lvals
    println("L = ", L)
    mags = zeros(Float64, length(temps))  # Array of magnetization per site
    err_mags = zeros(Float64, length(temps))

    suzzs = zeros(Float64, length(temps))  # Array of specific heat
    err_suzzs = zeros(Float64, length(temps))

    for i in eachindex(temps)
        T = temps[i]
        println("T = ", T)
        loc = joinpath([basepath, "Size$L", "potts_uncorr_configs_temp$(T)_L$(L).txt"])
        uncorr_configs = readdlm(loc, ',', Int64)
        m_arr = zeros(size(uncorr_configs)[2])
        Threads.@threads for col in 1:size(uncorr_configs)[2]
            m_arr[col] = magnetization(uncorr_configs[:, col], L, 3, 2) / L^2
        end

        mags[i] = mean(m_arr) 
        err_mags[i] = bootstrap_err(m_arr, A -> mean(A); r=bootstrap_samples)

        suzz_kth(m_arr, T, nsites, k) = (1/T) * nsites * cumulant(m_arr, k)
        order_k = 2
        suzzs[i] = suzz_kth(m_arr, T, L*L, order_k)
        err_suzzs[i] = bootstrap_err(m_arr, suzz_kth, T, L*L, order_k; r=bootstrap_samples)
    end

    errorbars!(ax1, temps, mags, err_mags, whiskerwidth = 10)
    scatterlines!(ax1, temps, mags ,markersize = 7, label="Size=$(L)x$(L)", linestyle = :dot, color=cols[L])
    
    errorbars!(ax2, temps, suzzs, err_suzzs, whiskerwidth = 10)
    scatterlines!(ax2, temps, suzzs, markersize = 7, label="Size=$(L)x$(L)", linestyle = :dot, color=cols[L])
end

axislegend(ax1)
axislegend(ax2)
display(f1)
display(f2)