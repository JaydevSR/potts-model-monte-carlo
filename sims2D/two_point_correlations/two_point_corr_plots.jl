using CairoMakie
using DelimitedFiles
using StatsKit

lattice_size = 128
ss_corr = readdlm("data/ss_correlation_fn_size_$lattice_size.txt", ',', Float64)[:, 1]

ll = length(ss_corr)
lh = ll ÷ 2

ss_corr_half = ss_corr[1:lh] + ss_corr[end:-1:lh+1] ./ 2

f1, f2 = Figure(), Figure();
ax1 = Axis(f1[1, 1], xlabel="r", ylabel="C(r)",
            title="Two point correlation function (periodic boundary, $lattice_size×$lattice_size)",
            xticks=0:2:lattice_size, yticks=0:0.02:1);
ax2 = Axis(f2[1, 1], xlabel="r", ylabel="C(r)",
            title="Two point correlation function (periodic boundary, $lattice_size×$lattice_size)",
            xticks=0:2:ll, yticks=0:0.02:1);

scatterlines!(ax1, 0:ll-1, ss_corr)
scatterlines!(ax2, 0:lh, ss_corr[1:lh+1])
# display(f1)
# display(f2)
save("plots/2DModel/correlations/two_point_correlations_complete_$lattice_size.svg", f1)
save("plots/2DModel/correlations/two_point_correlations_half_$lattice_size.svg", f2)

# EXpoenents

# linear fit
XX = log.(1:lh-1)
YY = log.(ss_corr[2:lh])
data_all = DataFrame(X=XX, Y=YY)
ols = lm(@formula(Y~X), data_all)
c_reg, m_reg = coef(ols)
y_reg = m_reg .* XX .+ c_reg

f3 = Figure();
ax3 = Axis(f3[1, 1], xlabel="log(r)", ylabel="log(C(r))",
            title="Two point correlation function (periodic boundary, $lattice_size×$lattice_size)");
scatterlines!(ax3, XX, YY, label="log(C(r))")
scatterlines!(ax3, XX, y_reg, color=:red, 
                label="linear fit: slope=$(round(m_reg, digits=4)), intercept=$(round(c_reg, digits=4))")
axislegend(ax3, position=:rt)
# display(f3)
save("plots/2DModel/correlations/two_point_correlations_fit_$lattice_size.svg", f3)
