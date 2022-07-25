include("../../src/pottsmc.jl")
using CairoMakie

data = readdlm("data/max_arg_suzz.txt", ',', Float64)

ν = 5/6
Lvals = data[:, 1]
L_ν = Lvals.^-(1/ν)
T_star = data[:, 2]

f = Figure()
ax = Axis(f[1, 1], xlabel = "L^-1/ν", ylabel = "T_c(L)",
title = "PottsModel2D: Scaling of critical point using succeptibility")
xlims!(-0.0001, 0.05)
ylims!(0.98, 1.04)

# linear fit
data_all = DataFrame(X=L_ν, Y=T_star)
ols = lm(@formula(Y~X), data_all)
c_reg, m_reg = coef(ols) 
x_reg = append!([0.0], L_ν, [0.05])
y_reg = m_reg .* x_reg .+ c_reg


# scatter!(1 ./ (log.(Lvals).^2), T_star)
lines!(x_reg, y_reg, linestyle=:dot, color=:black, label="linear fit")
scatter!(L_ν, T_star, label="datapoints")
scatter!([0], [c_reg], color=:red, marker=:xcross, label="intercept=$(round(c_reg, digits=3))", markersize=16)

axislegend(ax, position=:rb)

display(f)