include("../../src/pottsmc.jl")
using CairoMakie

temperatures = [0.9, 1.0, 1.1]
external_fields = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
lattice_size = 32
eqsteps = 10_000
msteps = 1_00_000
mag_definition = :vac

potts = PottsModel2D(lattice_size, 3, :cold)

stack=LazyStack(Int)
cluster=falses(size(potts.lattice))

mags = zeros(Float64, msteps)
acceptance_rates = zeros(Float64, length(temperatures), length(external_fields))

for tidx in eachindex(temperatures)
    temp = temperatures[tidx]
    for hidx in eachindex(external_fields)
        h = external_fields[hidx]

        println("T = $temp, h = $h")
        accepted = 0
        for i=1:eqsteps  # equilibration
            wolff_cluster_update!(potts, temp, h; fix_vacuum=true, stack=stack, cluster=cluster)
        end
        
        for j in 1:msteps
            accepted += wolff_cluster_update!(potts, temp, h; fix_vacuum=true, stack=stack, cluster=cluster)
            mags[j] = magnetization(potts, use_definition=mag_definition)
        end

        acceptance_rates[tidx, hidx] = accepted/msteps
    end
end


f = Figure();
ax = Axis(f[1,1], xscale=log10,
    xlabel="external field", ylabel="acceptance ratio", title="Acceptance rations for L = $lattice_size")
for tidx in eachindex(temperatures)
    temp = temperatures[tidx]
    scatterlines!(external_fields[2:end], acceptance_rates[tidx, 2:end], label="T = $temp")
end

axislegend(ax)
display(f)
save("plots/2Dmodel/external_field/acceptance_rates_with_field_L$lattice_size.svg", f)