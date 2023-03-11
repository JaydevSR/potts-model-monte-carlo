include("../../src/pottsmc.jl")
using CairoMakie

# Simulation Parameters
L = 64
q = 3
temp = 1.0
nsteps = 1_00_000
wsteps = 1_000
esteps = 10_000
verbose = false
max_order = 6

println(".==================================")
println("| Lattice Size: $(L) x $(L)        ")
println(".==================================")
println("|                                  ")

model_wolff = PottsModel2D(L, q, :cold)
model_metro = PottsModel2D(L, q, :cold)

wolff_Padd = -(expm1(-1/temp))
spin_stack = LazyStack(Int)
spin_cluster = falses(L, L)

metro_paccept = Dict(append!([(i, 1.0) for i = -4:0], [(i, exp(-i / temp)) for i = 1:4]))

for step in 1:esteps
    wolff_cluster_update!(model_wolff, temp, P_add = wolff_Padd, stack = spin_stack, cluster = spin_cluster)
    metropolis_batch_update!(model_metro, temp; accept_probs = metro_paccept)
end

mag_wolff = zeros(Float64, (max_order, nsteps))
mag_metro = zeros(Float64, (max_order, nsteps))

for step in 1:nsteps
    wolff_cluster_update!(model_wolff, temp, P_add = wolff_Padd, stack = spin_stack, cluster = spin_cluster)
    metropolis_batch_update!(model_metro, temp; accept_probs = metro_paccept)
    mag_wolff[:, step] = susceptibility_kth(model_wolff.lattice, temp, L^2, 1:max_order)
    mag_metro[:, step] = susceptibility_kth(model_metro.lattice, temp, L^2, 1:max_order)
end

fig = Figure(resolution = (800, 600));
ax = Axis(fig[1, 1], xlabel = "Step", ylabel = "Magnetization")
ylims!(ax, (-1, 1))
lines!(ax, 1:nsteps, mag_wolff[1, :] ./ (L^2), label = "Wolff")
lines!(ax, 1:nsteps, mag_metro[1, :] ./ (L^2), label = "Metropolis")
axislegend(ax)
display(fig)

# Autocorrelation
acf_wolff = autocorrelation_fn(mag_wolff[1, :])
acf_metro = autocorrelation_fn(mag_metro[1, :])

fig_acf = Figure(resolution = (800, 600));
ax_acf = Axis(fig_acf[1, 1], xlabel = "Step", ylabel = "Autocorrelation")
ylims!(ax_acf, (-0.5, 1))
lines!(ax_acf, 1:10000, acf_wolff[1:10000], label = "Wolff")
lines!(ax_acf, 1:10000, acf_metro[1:10000], label = "Metropolis")
axislegend(ax_acf)
display(fig_acf)

sum(acf_metro[1:5000])
sum(acf_wolff[1:5000])