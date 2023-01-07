include("../../src/pottsmc.jl")
using CairoMakie
using LsqFit
using PyCall

lattice_sizes = [32, 48, 56, 64, 72, 80, 96, 128]
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

## Fit to parabola
println("Making fits ...")
fs = Figure();
axs = Axis(fs[1, 1], xlabel="T", ylabel="χ")
dfl = open("data/susceptibility_fit_params_peak_parabola_pycall.txt", "w")

suzz_model(xvals, params) = params["a"] .+ params["b"] * xvals .+ params["c"] * xvals.^2

for Lidx=eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]
    center = findmax(suzz[Lidx, :])[2]
    trange = center-2:center+2

    wt = 1 ./ err_suzz[Lidx, trange].^2

    py"""
    import lmfit

    def residuals(params, xvals, yvals=None, eps=None):
        a = params['a']
        b = params['b']
        c = params['c']

        model = a + b*xvals + c*xvals**2
        if yvals is None:
            return model
        if eps is None:
            return model - yvals
        return (model - yvals) / eps
    
    pars = lmfit.Parameters()
    pars.add_many(('a', 1.0), ('b', 1.0), ('c', 1.0))

    parabola_fit = lmfit.minimize($residuals, pars, args=($(temps[trange]), $(suzz[Lidx, trange]), $wt), method='leastsq')
    
    fitted_pars = {}
    for name, param in parabola_fit.params.items():
        fitted_pars[name] = (param.value, param.stderr)
    """

    println(py"fitted_pars")

    pars = Dict()
    pars_errs = Dict()
    for name in ["a", "b", "c"]
        pars[name] = py"fitted_pars[$name]"[1]
        pars_errs[name] = py"fitted_pars[$name]"[2]
    end

    b, c = pars["b"] ± pars_errs["b"], pars["c"] ± pars_errs["c"]
    peak = -b / (2c)

    println(peak)

    println(dfl, "Size $L :")
    println(dfl, " * peak location = $(Measurements.value(peak)) ± $(Measurements.uncertainty(peak))")
    println(dfl, "\n")

    errorbars!(axs, temps, suzz[Lidx, :], err_suzz[Lidx, :], color=cols[L], label="L = $L", whiskerwidth=12)
    scatter!(axs, temps, suzz[Lidx, :], color=cols[L], label="L = $L")
    linesrange = temps[trange[1]]:0.0005:temps[trange[end]]
    lines!(axs, linesrange, suzz_model(linesrange, pars), color=cols[L], label="L = $L")
end
axislegend(axs, position=:rt, merge=true)
display(fs)
close(dfl)
save("plots/2Dmodel/finite_size_scaling/susceptibility_parabolic_fits_pycall.svg", fs)
