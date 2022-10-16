include("../../src/pottsmc.jl")
using StatsKit
using DelimitedFiles

function spin_prod(s1, s2; q=3)
    theta1 = 2π/q * s1
    theta2 = 2π/q * s2
    return real(exp(im*(theta1 + theta2)))
end

# lattice_size = 64
# eqsteps = 10000
# n_steps = 50000
lattice_size = 128
eqsteps = 2000
n_steps = 10000

T_c = 1/log1p(sqrt(3))
potts = PottsModel2D(lattice_size, 3, :cold)
stack=LazyStack(Int)
cluster=falses(size(potts.lattice))

for i=1:eqsteps  # equilibration
    wolff_cluster_update!(potts, T_c, fix_vacuum=false, stack=stack, cluster=cluster)
end

ss_corr = ss_correlation_fn(potts.lattice, lattice_size)
for i in 1:n_steps
    for j in 1:12
        wolff_cluster_update!(potts, T_c, fix_vacuum=false, stack=stack, cluster=cluster)
    end
    ss_corr .+= ss_correlation_fn(potts.lattice, lattice_size)
end
ss_corr ./= n_steps
ss_corr ./= maximum(ss_corr)

# Write ss_corr function
open("data/ss_correlation_fn_size_$(lattice_size).txt", "w") do f
    writedlm(f, ss_corr, ',')
end
