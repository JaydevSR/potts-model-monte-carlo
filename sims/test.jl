include("../src/pottsmc.jl")

#=
Perform simulation
=#

L = 16
q = 3
println("================================\n")
println("    Lattice Size: $(L) x $(L)")
println("================================\n")

T = 1.0
nsteps = 50000

println("Calculating for T = $(T) ...")

m_arr = zeros(Float64, nsteps)
potts = PottsModel2D(L, q, :cold)

for step=1:nsteps
    if step%5000==0
        println("   |   | $(step / nsteps * 100) %")
    end
    m_arr[step] = magnetisation(potts) / potts.L^2
    # delE = metropolis_batch_update!(potts, T)
    wolff_cluster_update!(potts, T)
end

corrfn1 = autocorrelation_fn(m_arr[1:nsteps÷2])
corrfn2 = autocorrelation_fn(m_arr[(nsteps÷2+1):nsteps])
corrfn3 = autocorrelation_fn(m_arr)
s = lines(corrfn1[1:500], color=:blue, label="first half")
lines!(corrfn2[1:500], color=:red, label="second half")
lines!(corrfn3[1:500], color=:black, label="complete ($nsteps)")
axislegend()
display(s)