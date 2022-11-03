include("../../src/pottsmc.jl")
using CairoMakie
using LsqFit

lattice_sizes = [32, 48, 64, 80]
cols = Dict([(32, :blue), (48, :red), (64, :green), (80, :purple)])
temps = [0.900, 0.920, 0.940, 0.960, 0.980, 0.982, 0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996, 0.997, 0.998,
         0.999, 1.000, 1.001, 1.002, 1.003, 1.004, 1.006, 1.008, 1.010, 1.012, 1.014, 1.016, 1.018, 1.020, 1.022, 
         1.024, 1.026, 1.028, 1.030, 1.032, 1.034, 1.036, 1.038, 1.040, 1.060, 1.080, 1.100]

mags_def = 1
suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

suzz = zeros(Float64, (length(lattice_sizes), length(temps)))
err_suzz = zeros(Float64, (length(lattice_sizes), length(temps)))

for Lidx in eachindex(lattice_sizes)
    L = lattice_sizes[Lidx]
    Threads.@threads for tidx in eachindex(temps)
        T = temps[tidx]
        mags = readdlm(joinpath("data", "2DModel", "Size$(L)", "mags", "potts_mags_temp$(T)_L$(L).txt"), ',', Float64)
        mags ./= L^2

        @views suzz[Lidx, tidx] = suzz_kth(mags[mags_def, :], T, L^2, 2)
        @views err_suzz[Lidx, tidx] = bootstrap_err(mags[mags_def, :], A -> suzz_kth(A, T, L^2, 2); r=100)
    end
end

suzz_star = zeros(Float64, length(lattice_sizes))
error_suzz_star = zeros(Float64, length(lattice_sizes))
T_star = zeros(Float64, length(lattice_sizes))
T_star_err = zeros(Float64, length(lattice_sizes))

for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    ss, ts = findmax(suzz[stepL, :])
    suzz_star[stepL] = ss
    error_suzz_star[stepL] = err_suzz[stepL, ts]
    T_star[stepL] = temps[ts]
    T_star_err[stepL] = 0.5*(temps[ts+1] - temps[ts-1])
end

## Least squares fit to find the critical temperature
model(x, p) = p[1] .+ p[2] * x.^(-inv(p[3]))
p0 = [1.0, 1.0, 1.0]
# wt = 1 ./ T_star_err

fit = curve_fit(model, lattice_sizes, T_star, p0)  # weighted LsqFit (χ² minimization)
Tc, A, ν = fit.param .± standard_errors(fit)


## Plots
# fig = Figure(resolution = (800, 600));
# ax = Axis(fig[1, 1]; xlabel = "L^(-1/ν)", ylabel = "T*(L)", title="Scaling of peak location of susceptibility with size")

# # Plot the fit
# x = 
# lines!(ax, x, model.(x, Ref(p0)), color = :black, linestyle=:dot,
#     label="Fit (slope=$(a1_with_err), intercept=$(a0_with_err))")
# band!(ax, x_log, y_log_fit_val .- y_log_fit_val_err, y_log_fit_val .+ y_log_fit_val_err,
#     color = (:black, 0.2), label="Fit (slope=$(a1_with_err), intercept=$(a0_with_err))")

# errorbars!(ax, x_log, Measurements.value.(y_log), Measurements.uncertainty.(y_log),
#     color = :black, whiskerwidth = 15, label="Data")
# scatter!(ax, x_log, Measurements.value.(y_log), color = :red, label = "Data")

# axislegend(ax, merge=true, position=:lt)
# display(fig)
# save("plots/2Dmodel/finite_size_scaling/succeptibility_Tc_fss.svg", fig)