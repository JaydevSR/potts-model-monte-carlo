include("../../src/pottsmc.jl")
using CairoMakie
using LsqFit

lattice_sizes = [32, 48, 56, 64, 72, 80, 96]
cols = Dict([(32, :blue), (48, :red), (56, :pink), (64, :green), (72, :orange), (80, :purple), (96, :cyan), (128, :deepskyblue)])
temps = [0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996, 0.998,
         1.000, 1.002, 1.004, 1.006, 1.008, 1.010, 1.012, 1.014]

mags_def = 1
suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

mag = zeros(Float64, (length(lattice_sizes), length(temps)))
err_mag = zeros(Float64, (length(lattice_sizes), length(temps)))
suzz = zeros(Float64, (length(lattice_sizes), length(temps)))
err_suzz = zeros(Float64, (length(lattice_sizes), length(temps)))

println("Calculating susceptibilities ...")
for Lidx in eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]
    Threads.@threads for tidx in eachindex(temps)
        T = temps[tidx]
        mags = readdlm(joinpath("data", "2DModel", "Size$(L)", "mags", "potts_mags_temp$(T)_L$(L).txt"), ',', Float64)
        mags ./= L^2

        @views mag[Lidx, tidx] = mean(mags[mags_def, :])
        @views err_mag[Lidx, tidx] = bootstrap_err(mags[mags_def, :], A -> mean(A); r=200)
        @views suzz[Lidx, tidx] = suzz_kth(mags[mags_def, :], T, L^2, 2)
        @views err_suzz[Lidx, tidx] = bootstrap_err(mags[mags_def, :], A -> suzz_kth(A, T, L^2, 2); r=100)
    end
end

println("Making fits ...")
fs = Figure();
axs = Axis(fs[1, 1], xlabel="T", ylabel="χ")
dfl = open("data/susceptibility_fit_params.txt", "w")
for Lidx=eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]

    # susceptibility to gaussian fit
    suzz_model(x, p) = p[1]*exp(1).^(- p[2] * (x .- p[3]).^2)
    p0s = [0.5, 1.0, 1.0]
    wt = 1 ./ err_suzz[Lidx, :].^2

    errorbars!(axs, temps, suzz[Lidx, :], err_suzz[Lidx, :], color=cols[L], label="L = $L", whiskerwidth=12)
    scatter!(axs, temps, suzz[Lidx, :], color=cols[L], label="L = $L")

    gaussian_fit = curve_fit(suzz_model, temps, suzz[Lidx, :], p0s)
    # residual_sum = sum((suzz[Lidx, :] .- suzz_model(temps, gaussian_fit.param)).^2 ./ err_suzz[Lidx, :].^2)
    # chi_sq = residual_sum / (length(temps) - length(p0s))
    println("L=$L: ", sum(gaussian_fit.resid.^2))
    # println("L=$L: χ^2 = $chi_sq")
    par = gaussian_fit.param
    parerr = stderror(gaussian_fit)
    println(dfl, "Size $L :")
    println(dfl, " * scaling parameter = $(par[1]) ± $(parerr[1])")
    println(dfl, " * width parameter = $(par[2]) ± $(parerr[2])")
    println(dfl, " * center parameter = $(par[3]) ± $(parerr[3])")
    println(dfl, "\n")
    lines!(axs, temps[1]:0.0001:temps[end], suzz_model(temps[1]:0.0001:temps[end], gaussian_fit.param), color=cols[L], label="L = $L")
end
axislegend(axs, position=:rt, merge=true)
display(fs)
close(dfl)
save("plots/susceptibility_gaussian_fits.svg", fs)
