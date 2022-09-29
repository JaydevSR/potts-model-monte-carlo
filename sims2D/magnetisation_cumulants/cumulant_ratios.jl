include("../../src/pottsmc.jl")
using CairoMakie

suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

lattice_sizes=[16, 32, 48, 64]
temps = [0.980, 0.984, 0.988, 0.992, 0.996, 1.000, 1.004,
     1.008, 1.012, 1.016, 1.020, 1.024, 1.028, 1.032, 1.036, 1.040]
max_order = 6
errors = zeros(Float64, (length(lattice_sizes), max_order, length(temps)))
suzzs = zeros(Float64, (length(lattice_sizes), max_order, length(temps)))
for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    for idx in eachindex(temps)
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

for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    suzzs_L = suzzs[stepL, :, :]
    T_star = temps[argmax(suzzs_L[2, :])]
    f = Figure(resolution=(1000, 600));
    axes = [Axis(f[1, 1], xlabel="T", ylabel="Ratio", title="Ratio: Order 3 by Order 1 for L=$L"),
        Axis(f[2, 1], xlabel="T", ylabel="Ratio", title="Ratio: Order 4 by Order 2 for L=$L"),
        Axis(f[3, 1], xlabel="T", ylabel="Ratio", title="Ratio: Order 5 by Order 3 for L=$L"),
        ]

    for idx in eachindex(axes)
        ax = axes[idx]
        scatterlines!(ax, temps, suzzs_L[idx+3, :] ./ suzzs_L[idx+1, :])
        vlines!(ax, T_star, color=:blue, linestyle=:dash)
        hlines!(ax, 0, color=:blue, linestyle=:dash)
    end

    save(joinpath("plots", "2DModel", "cumulants", "cumulant_ratios_size$(L).svg"), f)
end