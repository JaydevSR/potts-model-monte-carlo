include("../src/pottsmc.jl")

data = readdlm("data/max_arg_suzz.txt", ',', Float64)

Lvals = data[:, 1]
T_star = data[:, 2]

f = Figure()
ax = Axis(f[1, 1], xlabel = "1/[log(L)]^2", ylabel = "T_c(L)",
title = "PottsModel2D: Scaling of critical point using succeptibility")
# xlims!(0, 0.4)
xlims!(0, 0.01)
ylims!(0.98, 1.04)

# scatter!(1 ./ (log.(Lvals).^2), T_star)
scatter!(Lvals.^-1.83, T_star)

display(f)