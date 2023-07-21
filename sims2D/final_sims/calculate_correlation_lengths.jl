include("../../src/pottsmc.jl")
using PyCall
using CairoMakie

krondelta(x, y) = Int(x == y)

function potts_correlation_fn(sites::Matrix, L::Int, q::Int)
    corrfn = zeros(Float64, L)
    nsamples = zeros(Float64, L)
    for α=0:q-1
        @inbounds for i=1:L
            for j=1:L
                r = abs(i-j)
                corrfn[r+1] += krondelta(sites[i,i], α) * krondelta(sites[i,j], α) 
                corrfn[r+1] += krondelta(sites[i,i], α) * krondelta(sites[i,j], α)
                nsamples[r+1] += 2
            end
        end
    end
    return (corrfn ./ nsamples) .- inv(q^2)
end

lattice_size = 48
eqsteps = 20_000
n_steps = 50_000

psuedoTc = Dict(
            48 => 1.004098703447542,
            64 => 1.0009103256327856,
            80 => 0.9997345330528584,
            96 => 0.9988572671581376,
            128 => 0.997941002359664)

temperature_rel = 1.01
temperature = temperature_rel * psuedoTc[lattice_size]
potts = PottsModel2D(lattice_size, 3, :cold)
stack=LazyStack(Int)
cluster=falses(size(potts.lattice))

for i=1:eqsteps  # equilibration
    wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
end

println("Calculating the correlation function ...")
ss_corr = potts_correlation_fn(potts.lattice, lattice_size, 3)
for i in 1:n_steps
    i%1000 == 0 ? print("=") : nothing
    i%10000 == 0 ? print("\n") : nothing
    for j in 1:15
        wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
    end
    ss_corr .+= potts_correlation_fn(potts.lattice, lattice_size, 3)
end
ss_corr ./= n_steps
ss_corr[2:end] .= (ss_corr[2:end] .+ ss_corr[end:-1:2]) ./ 2

# Fit corr_model to the tail of ss_corr

fcorr = Figure(resolution = (1000, 800), fontsize = 32);
axcorr = Axis(fcorr[1, 1],
    xlabel=L"r",
    ylabel=L"C_{\text{Potts}}(r)",
    title="Site-Site Correlation Function for L=$lattice_size, T / T_c(L) = $(temperature_rel)",
    xlabelsize = 46, ylabelsize = 46,
    xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
    ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23)

scatterlines!(axcorr, 0:lattice_size, [ss_corr; ss_corr[1]], linestyle=:dashdot, linewidth=5, markersize=25)
save(joinpath("plots", "2DModel", "final_plots", "spin_correlation_function_Tr$(temperature_rel)_L$lattice_size.svg"), fcorr)
display(fcorr)

r_vals = 2 : lattice_size
c_r = ss_corr[r_vals]
r_vals = collect(r_vals) .- 1

println("Making fits ...")
py"""
import lmfit
import numpy as np

def residuals(params, xvals, yvals=None, eps=None):
    a = params['a']
    b = params['b']
    c = params['c']
    N = $lattice_size

    model = a * (np.exp(xvals / -b) + np.exp((N - xvals) / -b)) + c
    if yvals is None:
        return model
    if eps is None:
        return model - yvals
    return (model - yvals) / eps

pars = lmfit.Parameters()
pars.add_many(('a', 1.0), ('b', 1.0), ('c', 1.0))

parabola_fit = lmfit.minimize(residuals, pars, args=($r_vals, $c_r), method='leastsq')
chi_red = parabola_fit.redchi

fitted_pars = {}
for name, param in parabola_fit.params.items():
    fitted_pars[name] = (param.value, param.stderr)

min_series = residuals(parabola_fit.params, $r_vals)
"""
println(py"fitted_pars")

a_fit = py"fitted_pars['a'][0]" ± py"fitted_pars['a'][1]"
ξ_fit = py"fitted_pars['b'][0]" ± py"fitted_pars['b'][1]"
chi_red_fit = round(py"chi_red", sigdigits=3)

println("Plotting ...")
ffit = Figure(resolution = (1000, 800), fontsize = 32);
fit_title = "L=$lattice_size, T/T_c(L) = $temperature_rel \n A=$a_fit, ξ=$ξ_fit"

axfit = Axis(ffit[1, 1],
    xlabel=L"r",
    ylabel=L"C_{\text{Potts}}(r)",
    title=fit_title,
    xlabelsize = 46, ylabelsize = 46,
    xgridstyle = :dashdot, xgridwidth = 1.1, xgridcolor = :gray23,
    ygridstyle = :dashdot, ygridwidth = 1.1, ygridcolor = :gray23)

scatter!(axfit, r_vals, c_r, label="numerical data", markersize=25)

fit_eqn = L"A \left[ \exp\left(-\frac{r}{\xi} \right) + \exp\left(-\frac{L - r}{\xi} \right) \right]"
lines!(axfit, r_vals, py"min_series", label=fit_eqn, linewidth=5, linestyle=:dash, color=:black)

axislegend(axfit, position=:ct)
display(ffit)
save(joinpath("plots", "2DModel", "final_plots", "fit_correlation_length_Tr$(temperature_rel)_L$lattice_size.svg"), ffit)