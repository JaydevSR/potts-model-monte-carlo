# Caclulate correlation length using distance dependent mass
include("../../src/pottsmc.jl")
using CairoMakie
using PyCall

lattice_size = 48
eqsteps = 5_000
n_steps = 1_00_000

psuedoTc = Dict(
            48 => 1.004098703447542,
            64 => 1.0009103256327856,
            80 => 0.9997345330528584,
            96 => 0.9988572671581376,
            128 => 0.997941002359664)

temperature_rel = 1.011
temperature = temperature_rel * psuedoTc[lattice_size]
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
    for j in 1:40
        wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
    end
    ss_corr .+= ss_correlation_fn(potts.lattice, lattice_size)
end
ss_corr ./= n_steps
ss_corr ./= maximum(ss_corr)

py"""
import scipy.optimize as opt
import math

def corrfn_ansatz(r, m, L, i):
    return math.exp(-m * r) / r**i + math.exp(m * (L - r)) / (L - r)**i

def corrfn_ratio_ansatz(r, m, L, i):
    return corrfn_ansatz(r + 1, m, L, i) / corrfn_ansatz(r, m, L, i)
"""

rng = 1:lattice_size÷2-1
# mvals0 = zeros(length(range))
mvals1 = zeros(length(rng))
for r in rng
    ratio = ss_corr[r+1] / ss_corr[r]

    py"""
    def f(x):
        return corrfn_ratio_ansatz($r, x, $lattice_size, 1) - $ratio

    sol = opt.root_scalar(f, bracket=[-1.0, 1.0], method='brentq')

    root, fnroot = sol.root, f(sol.root)
    """

    println("\n<== r = $r ==========>")
    println(py"root", py"fnroot")
end

figtitle = L"Distance dependent mass: $m_1(r) = \frac{1}{\xi_1(r)}$ for $L = %$lattice_size, T_r = %$temperature_rel$"

fig = Figure(fontsize=22);
ax = Axis(
        fig[1,1],
        xlabel=L"r",
        ylabel=L"m_1(r)",
        title=figtitle,
        xticks = range[1:2:end],
        xlabelsize = 22, ylabelsize = 22,
        xgridstyle = :dashdot, xgridwidth = 1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1, ygridcolor = :gray23,
        );
scatter!(ax, range, mvals1, marker=:xcross, markersize=16)
display(fig)
save(joinpath("plots", "2DModel", "final_plots", "ddm_lattice$(lattice_size)_T$temperature_rel.svg"), fig)