include("../../src/pottsmc.jl")
using PyCall
using CairoMakie

lattice_size = 64
eqsteps = 5_000
n_steps = 1_00_000

temperature = 0.9950
potts = PottsModel2D(lattice_size, 3, :cold)
stack=LazyStack(Int)
cluster=falses(size(potts.lattice))

for i=1:eqsteps  # equilibration
    wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
end

println("Calculating the correlation function ...")
ss_corr = ss_correlation_fn(potts.lattice, lattice_size)
for i in 1:n_steps
    i%1000 == 0 ? print("=") : nothing
    i%10000 == 0 ? print("\n") : nothing
    for j in 1:20
        wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
    end
    ss_corr .+= ss_correlation_fn(potts.lattice, lattice_size)
end
ss_corr ./= n_steps
ss_corr ./= maximum(ss_corr)

# Fit corr_model to the tail of ss_corr

fcorr = Figure();
axcorr = Axis(fcorr[1, 1], xlabel=L"r", ylabel=L"C(r) = \langle \sigma(0) \sigma(r) \rangle", title="Site-Site Correlation Function for L=$lattice_size, T=$temperature")

scatterlines!(axcorr, 0:lattice_size, [ss_corr; 1.0], linestyle=:dashdot)
save(joinpath("plots", "2DModel", "final_plots", "spin_correlation_function_L$lattice_size.svg"), fcorr)

r_vals = lattice_size÷4 : 3lattice_size÷4 + 2
c_r = ss_corr[r_vals]
r_vals = collect(r_vals) .- 1

println("Making fits ...")
py"""
import lmfit
import numpy as np

def residuals(params, xvals, yvals=None, eps=None):
    a = params['a']
    b = params['b']
    N = $lattice_size

    model = a * (np.exp(xvals / -b) + np.exp((N - xvals) / -b))
    if yvals is None:
        return model
    if eps is None:
        return model - yvals
    return (model - yvals) / eps

pars = lmfit.Parameters()
pars.add_many(('a', 1.0), ('b', 1.0))

parabola_fit = lmfit.minimize(residuals, pars, args=($r_vals, $c_r), method='leastsq')

fitted_pars = {}
for name, param in parabola_fit.params.items():
    fitted_pars[name] = (param.value, param.stderr)

min_series = residuals(parabola_fit.params, $r_vals)
"""
println(py"fitted_pars")

a_fit = py"fitted_pars['a'][0]" ± py"fitted_pars['a'][1]"
ξ_fit = py"fitted_pars['b'][0]" ± py"fitted_pars['b'][1]"

println("Plotting ...")
ffit = Figure();

fit_title = "Fitting Site-Site Correlation Function for L=$lattice_size, T=$temperature \n Fit Parameters: A=$a_fit, ξ=$ξ_fit"
axfit = Axis(ffit[1, 1], xlabel=L"r", ylabel=L"C(r) = \langle \sigma(0) \sigma(r) \rangle", title=fit_title)

scatter!(axfit, r_vals, c_r, label="numerical data")

fit_eqn = L"A \left[ \exp\left( \frac{-r}{\xi} \right) + \exp\left( \frac{L - r}{\xi} \right) \right]"
lines!(axfit, r_vals, py"min_series", label=fit_eqn)
axislegend(axfit, position=:ct)
# display(ffit)
save(joinpath("plots", "2DModel", "final_plots", "fit_correlation_length_L$lattice_size.svg"), ffit)