include("../src/pottsmc.jl")

Lvals = [16, 24, 32, 40, 48]
cols = Dict([(16, :blue), (24, :red), (32, :green), (40, :purple), (48, :black)])
temps = [0.90, 0.94, 0.98, 0.984, 0.990, 0.996, 1.000, 1.002, 1.004, 1.006, 1.008, 1.01, 1.016, 1.02, 1.024, 1.028, 1.04, 1.06, 1.10]
basepath = joinpath(["data", "magdata"])
err_nblocks = 20

f1 = Figure()
f2 = Figure()

ax1 = Axis(f1[1, 1], xlabel = "temperature, T", ylabel = "magnetisation (per site), m",
    title = "PottsModel2D: Magnetisation (per site) v/s temperature")

ax2 = Axis(f2[1, 1], xlabel = "temperature, T", ylabel = "succeptibility, χ",
    title = "PottsModel2D: Succeptibility (per site) v/s temperature")

T_star_arr = []

for L in Lvals
    mags = zeros(Float64, length(temps))  # Array of magnetisation per site
    err_mags = zeros(Float64, length(temps))

    suzzs = zeros(Float64, length(temps))  # Array of specific heat
    err_suzzs = zeros(Float64, length(temps))

    for i in eachindex(temps)
        T = temps[i]
        loc = joinpath([basepath, "Size$L", "potts_mags_temp$(T)_L$(L).txt"])
        m_arr = readdlm(loc, ',', Float64) ./ L^2

        mags[i] = mean(m_arr) 
        err_mags[i] = blocking_err(m_arr, A -> mean(A); blocks=err_nblocks)

        suzz_kth(m_arr, T, nsites, k) = (1/T) * nsites * cumulant(m_arr, k)
        suzzs[i] = suzz_kth(m_arr, T, L*L, 2)
        err_suzzs[i] = blocking_err(m_arr, suzz_kth, T, L*L, 2; blocks=err_nblocks)
    end

    errorbars!(ax1, temps, mags, err_mags, whiskerwidth = 10)
    scatterlines!(ax1, temps, mags ,markersize = 7, label="Size=$(L)x$(L)", linestyle = :dot, color=cols[L])
    
    # errorbars!(ax2, temps, suzzs, err_suzzs, whiskerwidth = 10)
    scatterlines!(ax2, temps, suzzs, markersize = 7, label="Size=$(L)x$(L)", linestyle = :dot, color=cols[L])
    push!(T_star_arr, temps[argmax(suzzs)])
end

open("data/max_arg_suzz.txt", "w") do io
    writedlm(io, [Lvals T_star_arr], ',')
end

axislegend(ax1)
axislegend(ax2)
display(f1)
display(f2)

