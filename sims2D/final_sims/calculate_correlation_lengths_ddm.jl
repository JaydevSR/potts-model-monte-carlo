# Caclulate correlation length using distance dependent mass
include("../../src/pottsmc.jl")
using CairoMakie
using NLsolve

krondelta(x, y) = Int(x == y)

function potts_correlation_fn(sites::Matrix, L::Int, q::Int)
    corrfn = zeros(Float64, L)
    nsamples = zeros(Float64, L)
    for α=0:q-1
        @inbounds for i=1:L
            for j=1:L
                r = abs(i-j)
                corrfn[r+1] += krondelta(sites[i,i], α) * krondelta(sites[i,j], α) 
                corrfn[r+1] += krondelta(sites[i,i], α) * krondelta(sites[i,j], α)
                nsamples[r+1] += 2
            end
        end
    end
    return (corrfn ./ nsamples) .- inv(q^2)
end

lattice_size = 48
eqsteps = 5_000
n_steps = 1_00_000

psuedoTc = Dict(
            48 => 1.004098703447542,
            64 => 1.0009103256327856,
            80 => 0.9997345330528584,
            96 => 0.9988572671581376,
            128 => 0.997941002359664)

temperature_rel = 1.011
temperature = temperature_rel * psuedoTc[lattice_size]
potts = PottsModel2D(lattice_size, 3, :cold)
stack=LazyStack(Int)
cluster=falses(size(potts.lattice))

for i=1:eqsteps  # equilibration
    wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
end

println("Calculating the correlation function ...")
ss_corr = ss_correlation_fn(potts.lattice, lattice_size)
for i in 1:n_steps
    i%1000 == 0 ? print("=") : nothing
    i%10000 == 0 ? print("\n") : nothing
    for j in 1:40
        wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
    end
    ss_corr .+= ss_correlation_fn(potts.lattice, lattice_size)
end
ss_corr ./= n_steps
ss_corr[2:end] .= (ss_corr[2:end] .+ ss_corr[end:-1:2]) ./ 2

#= '_' =#

lattice_size = 48
ss_corr = copy(eval(Symbol("ss_corr_$lattice_size")))

corrfn_ansatz(r, m, L, i) = exp(-m * r) / r^i + exp(m * (L - r)) / (L - r)^i

corrfn_ratio_ansatz(r, m; L::Int64, i::Int64=1) = corrfn_ansatz(r + 1, m, L, i) / corrfn_ansatz(r, m, L, i)

function residual!(F, x; ratio, r, L, i)
    F[1] = ratio - corrfn_ratio_ansatz(r, x[1]; L=L, i=i)
end

rng = 2:lattice_size÷2-2
mvals1 = zeros(length(rng))

power_i = 0
for r in rng
    ratio = ss_corr[r+1] / ss_corr[r]

    f1!(F, x) = residual!(F, x; ratio=ratio, r=r-1, L=lattice_size, i=power_i)
    initial_x = [1.0]

    soln1 = nlsolve(f1!, initial_x, autodiff=:forward, ftol=1e-10, iterations=10_000)

    mvals1[r-1] = soln1.zero[1]

    # println("\n<== r = $r ==========>")
    # println(soln1.zero)
end

figtitle = L"Distance dependent mass: $m_%$power_i(r) = \frac{1}{\xi_%$power_i(r)}$ for $L = %$lattice_size, T_r = %$temperature_rel$"

fig = Figure(fontsize=22);
ax = Axis(
        fig[1,1],
        xlabel=L"r",
        ylabel=L"m_%$power_i(r)",
        title=figtitle,
        xticks = rng[1:2:end],
        xlabelsize = 22, ylabelsize = 22,
        xgridstyle = :dashdot, xgridwidth = 1, xgridcolor = :gray23,
        ygridstyle = :dashdot, ygridwidth = 1, ygridcolor = :gray23,
        );
scatter!(ax, rng.-1, mvals1, marker=:xcross, markersize=16)
display(fig)
save(joinpath("plots", "2DModel", "final_plots", "ddm_lattice$(lattice_size)_T$(temperature_rel)_i$(power_i)_sign.svg"), fig)