using CairoMakie

orders = 2:7
exp_vals = [1.743, 1.608, 1.459, 1.408, 1.204, 1.125]
exp_errs = [0.011, 0.031, 0.049, 0.035, 0.053, 0.073]

f = Figure();
ax = Axis(f[1, 1], xlabel="Order", ylabel="Exponent", title="Scaling exponents of Order $order derivative of ln(Z) with size");

col = :dodgerblue
band!(ax, orders, exp_vals .- exp_errs, exp_vals .+ exp_errs, color = (col, 0.2))
errorbars!(ax, orders, exp_vals, exp_errs, color = col, whiskerwidth = 15)
scatterlines!(ax, orders, exp_vals, color = col)

display(f)
save("plots/2Dmodel/finite_size_scaling/cumulant_exponents_with_order.svg", f)