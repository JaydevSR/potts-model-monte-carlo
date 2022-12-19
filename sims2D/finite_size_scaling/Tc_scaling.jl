include("../../src/pottsmc.jl")
using CairoMakie
using LsqFit

lattice_sizes = [48, 56, 64, 72, 80, 96]
cols = Dict([(32, :blue), (48, :red), (64, :green), (80, :purple)])
Tc_L = [
    1.004098558580445 ± 0.00012838372820415597,
    1.002880459951727 ± 0.00016614142151693173,
    1.0019690678472521 ± 0.00020775901663022396,
    1.0011767926962136 ± 0.00024041071306175725,
    1.0005831960780498 ± 0.0002808961827846406,
    0.9995852248655276 ± 0.0003212285341357723
]
Tc_vals = Measurements.value.(Tc_L)
Tc_err = Measurements.uncertainty.(Tc_L)

fit_model(x, p) = p[1] .- p[2] * x.^-(1/p[3])
p0 = [1.0, 1.0, 1.0]
# wt = inv.(Tc_err)

fit = curve_fit(fit_model, lattice_sizes, Tc_vals, p0)
par = fit.param
par_err = stderror(fit)

println(par .± par_err)

# Plotting
f = Figure();
ax = Axis(f[1, 1], xlabel="L", ylabel="Tc(L)", title="Peak location v/s lattice size")
lines!(ax, lattice_sizes[1]:lattice_sizes[end], 
    fit_model(lattice_sizes[1]:lattice_sizes[end], par), label="Fit", color=:red)
errorbars!(ax, lattice_sizes, Tc_vals, Tc_err, label="Data", color=:black)
scatter!(ax, lattice_sizes, Tc_vals, label="Data", color=:black)
axislegend(ax, merge=true)
display(f)

residual_sum = sum((fit_model(lattice_sizes, par) .- Tc_vals).^2 ./ Tc_err.^2)

#= NOTE ON RESULTS

These fits are without considering the errors on data points. I am getting huge parameter errors if I use this information. I am not sure why maybe I should try the python package for chi^2 fitting.
Measurement{Float64}[0.99251 ± 0.0006, -0.181 ± 0.023, 1.408 ± 0.092]

My first guess is that the gaussian fits are not able to capture the peak location, nu should be close to 0.84 but it is coming out to be 1.4. Mostly the higher the value of nu the slower the peak moves towards Tc. It can be seen in the fit plots that gaussian peak moving behind the actual susceptibility peak. Still the Tc value should be reliable as as L -> \infty the peak location should be at Tc for both the cases.

=#
