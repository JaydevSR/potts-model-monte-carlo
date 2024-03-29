include("../../src/pottsmc.jl")
using CairoMakie
using Interpolations

suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

lattice_sizes=[16, 32, 48, 64]
temps = [0.500, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.920,
     0.940, 0.960, 0.980, 0.984, 0.988, 0.992, 0.996, 1.000, 1.004,
     1.008, 1.012, 1.016, 1.020, 1.024, 1.028, 1.032, 1.036, 1.040,
     1.060, 1.080, 1.100, 1.150, 1.200, 1.250, 1.300, 1.350, 1.400, 1.500]
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

# Determine widths
widths = zeros(Float64, length(lattice_sizes))
starts = zeros(Float64, length(lattice_sizes))
ends = zeros(Float64, length(lattice_sizes))
for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    lerp = LinearInterpolation(temps, suzzs[stepL, 2, :])
    # calculate the width non zero part of the curve
    mx = maximum(abs.(suzzs[stepL, 2, :]))
    delta = 0.005
    startp, endp = 0, 0
    for t in first(temps):delta:last(temps)
        if abs(lerp(t)) > 0.5mx
            startp = t - delta
            break
        end
    end
    for t in last(temps):-delta:first(temps)
        if abs(lerp(t)) > 0.5mx
            endp = t + delta
            break
        end
    end

    starts[stepL] = startp
    ends[stepL] = endp
    widths[stepL] = endp - startp
end

# starts, ends, widths = ([0.98, 0.985, 0.99, 0.99], [1.090, 1.035, 1.020, 1.015], [0.11, 0.050, 0.030, 0.025])