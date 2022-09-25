include("../../src/pottsmc.jl")
using CairoMakie
using StatsKit

function spin_prod(s1, s2; q=3)
    theta1 = 2π/q * s1
    theta2 = 2π/q * s2
    return real(exp(im*(theta1 + theta2)))
end

lattice_sizes = [32, 48, 64]
eqsteps = 10000

T_c = 1/log1p(sqrt(3))
potts = PottsModel2D(lattice_sizes[end], 3, :cold)
stack=LazyStack(Int)
cluster=falses(size(potts.lattice))

for i=1:eqsteps  # equilibration
    wolff_cluster_update!(potts, T_c, fix_vacuum=false, stack=stack, cluster=cluster)
end

ss_corr = ss_correlation_fn(potts.lattice, lattice_sizes[end])
for i in 1:10_000
    for j in 1:12
        wolff_cluster_update!(potts, T_c, fix_vacuum=false, stack=stack, cluster=cluster)
    end
    ss_corr .+= ss_correlation_fn(potts.lattice, lattice_sizes[end])
end
ss_corr ./= 10_000
ss_corr ./= maximum(ss_corr)

ll = length(ss_corr)
lh = ll ÷ 2

f1, f2 = Figure(), Figure();
ax1 = Axis(f1[1, 1], xlabel="r", ylabel="C(r)",
            title="Two point correlation function (periodic boundary, $L×$L)",
            xticks=0:2:L, yticks=0:0.02:1);
ax2 = Axis(f2[1, 1], xlabel="r", ylabel="C(r)",
            title="Two point correlation function (periodic boundary, $L×$L)",
            xticks=0:2:ll, yticks=0:0.02:1);

scatterlines!(ax1, 0:ll-1, ss_corr)
scatterlines!(ax2, 0:lh, ss_corr[1:lh+1])
save("plots/2DModel/correlations/two_point_correlations_complete.png", f1)
save("plots/2DModel/correlations/two_point_correlations_half.png", f2)

# EXpoenents

# linear fit
XX = log.(1:5)
YY = log.(ss_corr[2:6])
data_all = DataFrame(X=XX, Y=YY)
ols = lm(@formula(Y~X), data_all)
c_reg, m_reg = coef(ols)
y_reg = m_reg .* XX .+ c_reg

f3 = Figure();
ax3 = Axis(f3[1, 1], xlabel="log(r)", ylabel="log(C(r))",
            title="Two point correlation function (periodic boundary, $L×$L)");
scatterlines!(ax3, XX, YY, label="log(C(r))")
scatterlines!(ax3, XX, y_reg, color=:red, 
                label="linear fit: slope=$(round(m_reg, digits=4)), intercept=$(round(c_reg, digits=4))")
axislegend(ax3, position=:rt)
save("plots/2DModel/correlations/two_point_correlations_fit.png", f3)
