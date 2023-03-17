include("../../src/pottsmc.jl")

suzz_kth(M_arr, T, nsites, k) = cumulant(M_arr, k) / (T * nsites)

lattice_sizes = [32, 48, 56, 64, 72, 80, 96, 128]

temps = [
    0.980, 0.982, 0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996, 0.998,
    1.000, 1.002, 1.004, 1.006, 1.008, 1.010, 1.012, 1.014, 1.016, 1.018, 1.020, 1.022, 1.024, 1.026,
    1.028, 1.030, 1.032, 1.034, 1.036, 1.038, 1.040
]

mags_def = 1 # maximum definition

# save temps
open(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), "w") do f
    writedlm(f, temps)
end

max_order = 7
errors = zeros(Float64, (max_order, length(temps)))
suzzs = zeros(Float64, (max_order, length(temps)))

for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    println("L = $L")
    Threads.@threads for idx in eachindex(temps)
        T = temps[idx]
        mags = readdlm(joinpath("data", "2DModel", "Size$(L)", "mags", "potts_mags_temp$(T)_L$(L).txt"), ',', Float64)

        suzzs[1, idx] = mean(mags[mags_def, :]) / L^2
        errors[1, idx] = bootstrap_err(mags[mags_def, :] ./ L^2, mean)
        for order_k in 2:max_order
            suzzs[order_k, idx] = suzz_kth(mags[mags_def, :], T, L*L, order_k)
            errors[order_k, idx] = bootstrap_err(mags[mags_def, :], suzz_kth, T, L*L, order_k; r=200)
        end
    end

    # save data
    open(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_L$(L).txt"), "w") do f
        writedlm(f, suzzs)
    end

    open(joinpath("data", "2DModel", "susceptibilities", "potts_suzzs_errors_L$(L).txt"), "w") do f
        writedlm(f, errors)
    end
end