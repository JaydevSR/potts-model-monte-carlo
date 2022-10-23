include("../../src/pottsmc.jl")
using CairoMakie

suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

lattice_sizes = [32, 48, 64, 80]
temps = [0.980, 0.982, 0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996, 0.998,
        1.000, 1.002, 1.004, 1.006, 1.008, 1.010, 1.012, 1.014, 1.016, 1.018, 1.020, 1.022, 1.024, 1.026,
        1.028, 1.030, 1.032, 1.034, 1.036, 1.038, 1.040]
max_order = 6
errors = zeros(Float64, (length(lattice_sizes), max_order, length(temps)))
suzzs = zeros(Float64, (length(lattice_sizes), max_order, length(temps)))

for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    Threads.@threads for idx in eachindex(temps)
        T = temps[idx]
        mags = readdlm(joinpath("data", "2DModel", "Size$(L)", "mags", "potts_mags_temp$(T)_L$(L).txt"), ',', Float64)
        mags ./= L^2
        mags_def = 1
        suzzs[stepL, 1, idx] = mean(mags[mags_def, :])
        errors[stepL, 1, idx] = bootstrap_err(mags[mags_def, :], mean)
        for order_k in 2:max_order
            suzzs[stepL, order_k, idx] = suzz_kth(mags[mags_def, :], T, L*L, order_k)
            errors[stepL, order_k, idx] = bootstrap_err(mags[mags_def, :], suzz_kth, T, L*L, order_k; r=200)
        end
    end
end

# Plot the data
for order in 1:max_order
    let f = Figure(resolution=(800, 600))
        ax = Axis(f[1, 1], xlabel="T", ylabel="Order $(order) value", title="Order $(order) derivative of ln(Z) with respect to h as h → 0")
        for stepL in eachindex(lattice_sizes)
            L = lattice_sizes[stepL]
            errorbars!(ax, temps, suzzs[stepL, order, :], errors[stepL, order, :], whiskerwidth = 12)
            scatterlines!(ax, temps, suzzs[stepL, order, :], markersize = 15, label="size=$(L)×$(L)", linestyle = :dashdot,
            marker=:diamond, linewidth=3)
        end
        axislegend(ax, position=:rt)
        save("plots/2DModel/cumulants/order_$(order).svg", f)
    end
end

order_colors = [:red, :blue, :dodgerblue, :dodgerblue, :orange, :orange]
# Plot Ratios
for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    suzzs_L = suzzs[stepL, :, :] .± errors[stepL, :, :]
    f = Figure(resolution=(1000, 600));
    axes = [Axis(f[1, 1], xlabel="T", ylabel="Ratio", title="Ratio of odd orders for L=$L"),
        Axis(f[1, 2], xlabel="T", ylabel="Ratio", title="Ratio of even orders for L=$L")]

    for order in 3:2:max_order
        ratio = suzzs_L[order, :] ./ suzzs_L[order-2, :]
        ratio_values = Measurements.value.(ratio)
        ratio_errors = Measurements.uncertainty.(ratio)

        band!(axes[1], temps, ratio_values .- ratio_errors, ratio_values .+ ratio_errors,
            color = (order_colors[order], 0.25), label="Order $(order) / Order 1")
        errorbars!(axes[1], temps, ratio_values, ratio_errors, label="Order $(order) / Order 1", whiskerwidth=12)
        scatterlines!(axes[1], temps, ratio_values, label="Order $(order) / Order 1", color=order_colors[order])
    end

    for order in 4:2:max_order
        ratio = suzzs_L[order, :] ./ suzzs_L[2, :]
        ratio_values = Measurements.value.(ratio)
        ratio_errors = Measurements.uncertainty.(ratio)

        band!(axes[2], temps, ratio_values .- ratio_errors, ratio_values .+ ratio_errors,
            color = (order_colors[order], 0.25), label="Order $(order) / Order 2")
        errorbars!(axes[2], temps, ratio_values, ratio_errors, label="Order $(order) / Order 2", whiskerwidth=12)
        scatterlines!(axes[2], temps, ratio_values, label="Order $(order) / Order 2", color=order_colors[order])
    end

    axislegend(axes[1], position=:lt, merge=true)
    axislegend(axes[2], position=:rt, merge=true)

    # display(f)
    save(joinpath("plots", "2DModel", "cumulants", "cumulant_ratios_size$(L).svg"), f)
end

# save values
# for stepL in eachindex(lattice_sizes)
#     L = lattice_sizes[stepL]
#     suzzs_L = suzzs[stepL, :, :]
#     errors_L = errors[stepL, :, :]
#     mkpath(joinpath("data", "2DModel", "Size$(L)", "cumulants"))
#     open(joinpath("data", "2DModel", "Size$(L)", "cumulants", "cumulant_values.txt"), "w") do io
#         writedlm(io, suzzs_L, ',')
#     end

#     open(joinpath("data", "2DModel", "Size$(L)", "cumulants", "cumulant_errors.txt"), "w") do io
#         writedlm(io, errors_L, ',')
#     end
# end