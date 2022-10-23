using CairoMakie

exp_vals = [1.749, 1.599, 1.461, 1.075, 1.077]
exp_errs = [0.011, 0.033, 0.062, 0.044, 0.056]

f = Figure();
ax = Axis(f[1, 1], xlabel="Order", ylabel="Exponent", title="Scaling exponents of Order $order derivative of ln(Z) with size");

col = :dodgerblue
band!(ax, 2:6, exp_vals .- exp_errs, exp_vals .+ exp_errs, color = (col, 0.2))
errorbars!(ax, 2:6, exp_vals, exp_errs, color = col, whiskerwidth = 15)
scatterlines!(ax, 2:6, exp_vals, color = col)

display(f)
save("plots/2Dmodel/finite_size_scaling/cumulant_exponents_with_order.svg", f)