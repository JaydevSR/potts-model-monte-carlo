include("../src/potts.jl")
include("../src/statutils.jl")

using CairoMakie

#=
Perform simulation
=#

L = 16
q = 3
println("================================\n")
println("    Lattice Size: $(L) x $(L)")
println("================================\n")

T = 1.0
nsteps = 10000

println("Calculating for T = $(T) ...")

E_arr = zeros(Float64, nsteps)

E_arr[1] = hamiltonian(potts) / potts.L^2
for step=2:nsteps
    # delE = metropolis_batch_update!(potts, T)
    wolff_cluster_update!(potts, T)
    E_arr[step] = hamiltonian(potts) / potts.L^2 # E_arr[step-1] + (delE / potts.L^2)
end

corrfn = autocorrelation_fn(E_arr)

s = scatter(corrfn[1:250])
display(s)