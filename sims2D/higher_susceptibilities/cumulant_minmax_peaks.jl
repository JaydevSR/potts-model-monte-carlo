include("../../src/pottsmc.jl")
using CairoMakie
using LsqFit

Tc_inf = 0.99628906
lattice_sizes = [48, 56, 64, 72, 80, 96, 128]
cols = Dict([(32, :blue), (48, :red), (56, :pink), (64, :green), (72, :orange), (80, :purple), (96, :cyan), (128, :deepskyblue)])
temps = [0.982, 0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996, 0.998,
         1.000, 1.002, 1.004, 1.006, 1.008, 1.010, 1.012, 1.014, 1.016, 1.018, 1.020]

mags_def = 1
suzz_order = 3
suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

mag = zeros(Float64, (length(lattice_sizes), length(temps)))
err_mag = zeros(Float64, (length(lattice_sizes), length(temps)))
suzz = zeros(Float64, (length(lattice_sizes), length(temps)))
err_suzz = zeros(Float64, (length(lattice_sizes), length(temps)))

parabola(x, p) = p[1] .+ p[2]*x .+ p[3]*x.^2

parabola_extrema(p) = -p[2] / (2 * p[3])
p0 = [0.0, 0.0, 1.0]

println("Calculating susceptibilities ...")
for Lidx in eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]
    Threads.@threads for tidx in eachindex(temps)
        T = temps[tidx]
        mags = readdlm(joinpath("data", "2DModel", "Size$(L)", "mags", "potts_mags_temp$(T)_L$(L).txt"), ',', Float64)
        mags ./= L^2

        @views mag[Lidx, tidx] = mean(mags[mags_def, :])
        @views err_mag[Lidx, tidx] = bootstrap_err(mags[mags_def, :], A -> mean(A); r=200)
        @views suzz[Lidx, tidx] = suzz_kth(mags[mags_def, :], T, L^2, suzz_order)
        @views err_suzz[Lidx, tidx] = bootstrap_err(mags[mags_def, :], A -> suzz_kth(A, T, L^2, suzz_order); r=100)
    end
end

min_locations = fill(0.0 ± 0.0, length(lattice_sizes))
max_locations = fill(0.0 ± 0.0, length(lattice_sizes))

println("Plotting ...")
fs = Figure();
axs = Axis(fs[1, 1], xlabel="T", ylabel="χ^($suzz_order)")
for Lidx=eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]
    max_peak = findmax(suzz[Lidx, :])
    min_peak = findmin(suzz[Lidx, :])
    min_range = (min_peak[2]-2):(min_peak[2]+2)
    max_range = (max_peak[2]-2):(max_peak[2]+2)

    fit_min = curve_fit(parabola, temps[min_range], suzz[Lidx, min_range], p0)
    fit_max = curve_fit(parabola, temps[max_range], suzz[Lidx, max_range], p0)

    min_locations[Lidx] = parabola_extrema(fit_min.param .± stderror(fit_min))
    max_locations[Lidx] = parabola_extrema(fit_max.param .± stderror(fit_max))

    lines!(axs, temps[min_range[1]]:0.0001:temps[min_range[end]], parabola(temps[min_range[1]]:0.0001:temps[min_range[end]], fit_min.param), color=:grey7, label="L = $L")
    lines!(axs, temps[max_range[1]]:0.0001:temps[max_range[end]], parabola(temps[max_range[1]]:0.0001:temps[max_range[end]], fit_max.param), color=:grey7, label="L = $L")
    errorbars!(axs, temps, suzz[Lidx, :], err_suzz[Lidx, :], color=cols[L], label="L = $L", whiskerwidth=12)
    scatterlines!(axs, temps, suzz[Lidx, :], color=cols[L], label="L = $L")
end
axislegend(axs, position=:rb, merge=true)
display(fs)

# log(T - Tc_inf) ∝ log(L)

println("Plotting ...")
f = Figure();
ax = Axis(f[1, 1], xlabel="1/L", ylabel="T - Tc_inf")
scatterlines!(ax, inv.(lattice_sizes), Measurements.value.(max_locations) .- Tc_inf, color=:red, label="Maxima")
scatterlines!(ax, inv.(lattice_sizes), Measurements.value.(min_locations) .- Tc_inf, color=:blue, label="Minima")
axislegend(ax, merge=true)
display(f)

amax, bmax = linear_least_squares(inv.(lattice_sizes), max_locations .- Tc_inf)
amin, bmin = linear_least_squares(inv.(lattice_sizes), min_locations .- Tc_inf)

# TODO: THE ERROR ESTIMATES SUCK ASS