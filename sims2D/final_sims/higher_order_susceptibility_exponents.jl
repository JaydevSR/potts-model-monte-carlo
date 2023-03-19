include("../../src/pottsmc.jl")
using CairoMakie
using LsqFit
using PyCall

lattice_sizes = [48, 56, 64, 72, 80, 96, 128]
temps = reshape(readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), ',', Float64), :)
max_order = 6

max_values = zeros(Float64, (max_order-1, length(lattice_sizes)))
max_errors = zeros(Float64, (max_order-1, length(lattice_sizes)))

for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]

    suzzs = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_L$(L).txt"), '\t', Float64)
    errors = readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_errors_L$(L).txt"), '\t', Float64)

    @assert size(suzzs) == size(errors) == (max_order, length(temps))

    for r in 2:max_order
        (mx, imx) = findmax(suzzs[r, :])
        max_values[r-1, stepL] = mx
        max_errors[r-1, stepL] = errors[r, imx]
    end
end

py"""
import numpy as np
import lmfit

def residuals_exp(params, xvals, yvals=None, eps=None):
    exp = params['exp']
    scale = params['scale']

    model = scale * xvals**exp
    if yvals is None:
        return model
    if eps is None:
        return model - yvals
    return (model - yvals) / eps

def get_scaling_parameters(lattice_sizes, values, errors):
    pars = lmfit.Parameters()
    pars.add_many(('scale', 0.0), ('exp', 1.0))
    out = lmfit.minimize(residuals_exp, pars, args=(lattice_sizes, values, errors), method='leastsq')
    return {
            'scale': (out.params['scale'].value, out.params['scale'].stderr),
            'exp': (out.params['exp'].value, out.params['exp'].stderr)
        }
"""

scaling_exps = fill(0.0 ± 0.0, max_order-1)
for r in 1:max_order-1
    order = r + 1

    vals = max_values[r, :] .± max_errors[r, :]
    vals ./= max_values[r, 1]

    max_params = py"""get_scaling_parameters($lattice_sizes, $(Measurements.value.(vals)), $(Measurements.uncertainty.(vals)))"""

    scaling_exps[r] = measurement(max_params["exp"]...)
end

println(scaling_exps)