include("../../src/pottsmc.jl")
using CairoMakie

temperature = 1.5
external_field = 0.05
lattice_size = 32
eqsteps = 1000
msteps = 1_00_000
mag_definition = :vac
accepted = 0

potts = PottsModel2D(lattice_size, 3, :cold)

stack=LazyStack(Int)
cluster=falses(size(potts.lattice))

mags = zeros(Float64, msteps)

for i=1:eqsteps  # equilibration
    wolff_cluster_update!(potts, temperature, external_field; fix_vacuum=true, stack=stack, cluster=cluster)
end

for j in 1:msteps
    accepted += wolff_cluster_update!(potts, temperature, external_field; fix_vacuum=true, stack=stack, cluster=cluster)
    mags[j] = magnetization(potts, use_definition=mag_definition)
end

@info "Acceptance ratio: $(accepted/msteps)"

f = Figure();
ax = Axis(f[1,1], xlabel="Magnetization", ylabel="density", title="Magnetization distribution: T = $temperature, L = $lattice_size")
density!(mags, npoints=500, color=("blue", 0.5), strokewidth = 1.5,
strokecolor = :grey20)
display(f)