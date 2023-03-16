include("../../src/pottsmc.jl")
using CairoMakie
using BenchmarkTools

lattice_sizes = [32, 48, 64, 80]

temperature = 0.9950 #J
metro_accept = Dict(append!([(i, 1.0) for i = -4:0], [(i, exp(-i / temperature)) for i = 1:4]))
wolff_P_add = -(expm1(-1 / temperature))
wolff_stack = LazyStack(Int)

# compute auto-correlation times
nmags = 50_000
sum_window_metro = 3_000
sum_window_wolff = 500
nequib = 5_000

mags = zeros(Float64, nmags) # a shared array, Don't use multithreading or move to for loop

τ_values_metro = zero(lattice_sizes)
τ_values_wolff = zero(lattice_sizes)
stepping_times_metro = zeros(Float64, length(lattice_sizes))
stepping_times_wolff = zeros(Float64, length(lattice_sizes))

println("Calculating ac times ...")
for idx in eachindex(lattice_sizes)
    println("Lattice size = $(lattice_sizes[idx])")
    model = PottsModel2D(lattice_sizes[idx], 3, :cold)

    # Metropolis Algorithm
    for i in 1:nequib
        metropolis_batch_update!(model, temperature; accept_probs=metro_accept)
    end

    for i in 1:nmags
        mags[i] = magnetization(model; use_definition = :max)
        metropolis_batch_update!(model, temperature; accept_probs=metro_accept)
    end

    acf_metro = autocorrelation_fn(mags)
    τ_values_metro[idx] = view(acf_metro, 1:sum_window_metro) |> sum |> ceil |> x -> convert(Int64, x)
    stepping_times_metro[idx] = @belapsed metropolis_batch_update!(
                                            $model,
                                            $temperature;
                                            accept_probs=$metro_accept)

    # Wolff Cluster Algorithm
    cluster = falses(size(model.lattice))
    for i in 1:nequib
        wolff_cluster_update!(model, temperature; P_add=wolff_P_add, stack=wolff_stack, cluster=cluster)
    end

    for i in 1:nmags
        mags[i] = magnetization(model; use_definition = :max)
        wolff_cluster_update!(model, temperature; P_add=wolff_P_add, stack=wolff_stack, cluster=cluster)
    end

    acf_wolff = autocorrelation_fn(mags)
    τ_values_wolff[idx] = view(acf_wolff, 1:sum_window_wolff) |> sum |> ceil |> x -> convert(Int64, x)
    stepping_times_wolff[idx] = @belapsed wolff_cluster_update!(
                                                $model,
                                                $temperature;
                                                P_add=$wolff_P_add,
                                                stack=$wolff_stack,
                                                cluster=$cluster)
end

println("Metropolis: ", τ_values_metro)
println("Wolff: ", τ_values_wolff)

f_time = Figure();
ax_time = Axis(f_time[1,1], xlabel = "lattice size", ylabel = "time (in seconds)", title = "Computational time required for an uncorrelated sample")
ylims!(ax_time, (-0.05, 0.8))

ax_time_min = Axis(f_time, bbox = BBox(200, 400, 300, 500), xticklabelsize = 12, yticklabelsize = 12, title = "wolff only")

scatterlines!(ax_time, lattice_sizes, τ_values_wolff .* stepping_times_wolff, label="Wolff")
scatterlines!(ax_time_min, lattice_sizes, τ_values_wolff .* stepping_times_wolff, label="Wolff")

scatterlines!(ax_time, lattice_sizes, τ_values_metro .* stepping_times_metro, label="Metropolis")

axislegend(ax_time)

save("plots/2Dmodel/final_plots/computational_time_comparison.svg", f_time)