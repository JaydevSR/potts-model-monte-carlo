include("../../src/pottsmc.jl")
using PyCall
using CairoMakie

lattice_sizes = [48, 64, 128]
temps = reshape(readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), ',', Float64), :)

py"""
import lmfit
import numpy as np

def residuals(params, N, xvals, yvals=None, eps=None):
    a = params['a']
    b = params['b']
    c = params['c']

    model = a * (np.exp(xvals / -b) + np.exp((N - xvals) / -b)) + c
    if yvals is None:
        return model
    if eps is None:
        return model - yvals
    return (model - yvals) / eps

"""

for lidx in eachindex(lattice_sizes)
    L = lattice_sizes[lidx]
    r_vals = 2:L
    r_vals = collect(r_vals) .- 1
    ss_corr_data = readdlm(joinpath("data", "2DModel", "correlations", "potts_corrfun_L$L.txt"), '\t', Float64)
    ξ_vals = zeros(Float64, length(temps))
    ξ_errs = zeros(Float64, length(temps))
    for tidx in eachindex(temps)
        T = temps[tidx]
        corrfun_T = ss_corr_data[r_vals, tidx] .- ss_corr_data[L ÷ 2 - 1, tidx]

        py"""
        pars = lmfit.Parameters()
        pars.add_many(('a', 1.0), ('b', 1.0), ('c', 1.0))
        
        parabola_fit = lmfit.minimize(residuals, pars, args=($L, $r_vals, $corrfun_T), method='leastsq')
        chi_red = parabola_fit.redchi
        
        fitted_pars = {}
        for name, param in parabola_fit.params.items():
            fitted_pars[name] = (param.value, param.stderr)
        
        """
        ξ_vals[tidx] = py"fitted_pars['b'][0]"
        ξ_errs[tidx] = py"fitted_pars['b'][1]"
    end

    f = Figure()
    ax = Axis(f[1, 1], xlabel="T", ylabel="ξ(T)")
    errorbars!(ax, temps, ξ_vals, ξ_errs, color=:black, whiskerwidth=12)
    scatter!(ax, temps, ξ_vals)

    display(f)

    # save data
    open(joinpath("data", "2DModel", "correlations", "potts_corrlen_values_L$(L).txt"), "w") do f
        writedlm(f, ξ_vals)
    end

    open(joinpath("data", "2DModel", "correlations", "potts_corrlen_errors_L$(L).txt"), "w") do f
        writedlm(f, ξ_errs)
    end
end