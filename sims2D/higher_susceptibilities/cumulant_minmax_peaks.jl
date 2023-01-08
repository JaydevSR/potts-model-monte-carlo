include("../../src/pottsmc.jl")
using CairoMakie
using LsqFit
using PyCall

Tc_inf = 0.99691951 ± 7.4795e-04
lattice_sizes = [48, 56, 64, 72, 80, 96, 128]
cols = Dict([
            (32, :blue), (48, :red), (56, :pink),
            (64, :green), (72, :orange), (80, :purple),
            (96, :cyan), (128, :deepskyblue)
])

temps = [0.982, 0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996, 0.998,
         1.000, 1.002, 1.004, 1.006, 1.008, 1.010, 1.012, 1.014, 1.016, 1.018, 1.020]

mags_def = 1
suzz_order = 6
suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

mag = zeros(Float64, (length(lattice_sizes), length(temps)))
err_mag = zeros(Float64, (length(lattice_sizes), length(temps)))
suzz = zeros(Float64, (length(lattice_sizes), length(temps)))
err_suzz = zeros(Float64, (length(lattice_sizes), length(temps)))

parabola(x, p) = p[1] .+ p[2]*x .+ p[3]*x.^2

parabola_extrema(p) = -p[2] / (2 * p[3])
p0 = [1.0, 1.0, 1.0]

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

println("Calculating peaks ...")
for Lidx=eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]
    max_peak = findmax(suzz[Lidx, :])
    min_peak = findmin(suzz[Lidx, :])
    min_range = (min_peak[2]-2):(min_peak[2]+1)
    max_range = (max_peak[2]-1):(max_peak[2]+2)

    fit_min = curve_fit(parabola, temps[min_range], suzz[Lidx, min_range], p0)
    fit_max = curve_fit(parabola, temps[max_range], suzz[Lidx, max_range], p0)

    min_locations[Lidx] = parabola_extrema(fit_min.param .± stderror(fit_min))
    max_locations[Lidx] = parabola_extrema(fit_max.param .± stderror(fit_max))
end

min_loc_values = Measurements.value.(min_locations)
min_loc_errs = Measurements.uncertainty.(min_locations)
max_loc_values = Measurements.value.(max_locations)
max_loc_errs = Measurements.uncertainty.(max_locations)

println("Fitting peak location to find Tc ...")
py"""
import numpy as np
import lmfit

def residuals(params, xvals, yvals=None, eps=None):
    tinf = params['a']
    scaling = params['b']
    expnu = params['c']

    model = tinf + scaling * (xvals ** (- 1 / expnu))
    if yvals is None:
        return model
    if eps is None:
        return model - yvals
    return (model - yvals) / eps

pars_min = lmfit.Parameters()
pars_min.add_many(('a', 1.0), ('b', 1.0), ('c', 1.0))
pars_max = lmfit.Parameters()
pars_max.add_many(('a', 1.0), ('b', 1.0), ('c', 1.0))

out_min = lmfit.minimize(residuals, pars_min, args=($lattice_sizes, $min_loc_values, $min_loc_errs), method='leastsq')
out_max = lmfit.minimize(residuals, pars_max, args=($lattice_sizes, $max_loc_values, $max_loc_errs), method='leastsq')

tc_min = out_min.params['a'].value
tc_max = out_max.params['a'].value
tc_min_err = out_min.params['a'].stderr
tc_max_err = out_max.params['a'].stderr
"""

tc_min = py"tc_min" ± py"tc_min_err"
tc_max = py"tc_max" ± py"tc_max_err"

println("Tc_min = $tc_min")
println("Tc_max = $tc_max")
println("Tc_inf = $Tc_inf")

# write results to a file
open("data/cumulant_minmax_peak_collapse.txt", "a") do f
    println(f, "Order = $suzz_order ---------------------------------------------------+")
    println(f, "|   Tc_min = $(Measurements.value.(tc_min)) ± $(Measurements.uncertainty.(tc_min))")
    println(f, "|   Tc_max = $(Measurements.value.(tc_max)) ± $(Measurements.uncertainty.(tc_max))")
    println(f, "|   Tc_inf = $(Measurements.value.(Tc_inf)) ± $(Measurements.uncertainty.(Tc_inf))")
    println(f, "+----------------------------------------------------------------------+\n")
end