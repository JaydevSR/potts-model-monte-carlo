include("../../src/pottsmc.jl")
using CairoMakie
using StatsKit

suzz_kth(m_arr, T, nsites, k) = (1/T) * (nsites) * cumulant(m_arr, k)

lattice_sizes=[16, 32, 48, 64]
temps = [0.500, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.920,
        0.940, 0.960, 0.980, 0.984, 0.988, 0.992, 0.996, 1.000, 1.004,
        1.008, 1.012, 1.016, 1.020, 1.024, 1.028, 1.032, 1.036, 1.040,
        1.060, 1.080, 1.100, 1.150, 1.200, 1.250, 1.300, 1.350, 1.400, 1.500]

order_k = 2
errors = zeros(Float64, (length(lattice_sizes), length(temps)))
suzzs = zeros(Float64, (length(lattice_sizes), length(temps)))
for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    for idx in eachindex(temps)
        T = temps[idx]
        mags = readdlm(joinpath("data", "2DModel", "Size$(L)", "mags", "potts_mags_temp$(T)_L$(L).txt"), ',', Float64)
        mags ./= L^2
        mags_def = 1
        suzzs[stepL, idx] = suzz_kth(mags[mags_def, :], T, L*L, order_k)
        errors[stepL, idx] = bootstrap_err(mags[mags_def, :], suzz_kth, T, L*L, order_k; r=200)
    end
end

suzz_star = zeros(Float64, length(lattice_sizes))
T_star = zeros(Float64, length(lattice_sizes))

for stepL in eachindex(lattice_sizes)
    L = lattice_sizes[stepL]
    ss, ts = findmax(suzzs[stepL, :])
    suzz_star[stepL] = ss
    T_star[stepL] = temps[ts]
end

# linear fit
data_all = DataFrame(X=log.(lattice_sizes), Y=log.(suzz_star))
ols = lm(@formula(Y~X), data_all)
c_reg, m_reg = coef(ols)
x_reg = log.(lattice_sizes)
y_reg = m_reg .* x_reg .+ c_reg

##
f = Figure();
axes = [
    Axis(f[1,1], xlabel="log(L)", ylabel="log(χ*)", title="Maxima of Susceptibility: χ*(L) (γ=13/5, ν=5/6, γ/ν = 1.733)"),
]
lines!(axes[1], x_reg, y_reg, label="Linear Fit: log(χ*(L)) = $(round(m_reg, digits=3)) * log(L) + $(round(c_reg, digits=3))",
         color=:red, linewidth=2)
scatter!(axes[1], log.(lattice_sizes), log.(suzz_star), label="log(χ*)", color=:black, markersize=10)
axislegend(axes[1], position=:lt)

save("plots/2Dmodel/finite_size_scaling/susceptibility_peak_fss.png", f)
display(f)
