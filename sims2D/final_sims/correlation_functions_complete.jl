include("../../src/pottsmc.jl")
using ProgressMeter

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

lattice_sizes = [32, 48, 64]
eqsteps = 20_000
n_steps = 50_000

temps = reshape(readdlm(joinpath("data", "2DModel", "susceptibilities", "potts_temps.txt"), ',', Float64), :)

stack = LazyStack(Int)

for lidx in eachindex(lattice_sizes)
    L = lattice_sizes[lidx]
    Pmeter = Progress(length(temps) + 1, "Calculating for L = $L ...")
    next!(Pmeter)
    potts = PottsModel2D(L, 3, :cold)
    cluster = falses(size(potts.lattice))
    ss_corr = zeros(Float64, (L, length(temps)))
    
    @inbounds for tidx in eachindex(temps)
        temperature = temps[tidx]
        for i=1:eqsteps  # equilibration
            wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
        end
        
        ss_corr[:, tidx] = potts_correlation_fn(potts.lattice, L, 3)
        for i in 1:n_steps
            for j in 1:15
                wolff_cluster_update!(potts, temperature, stack=stack, cluster=cluster)
            end
            ss_corr[:, tidx] += potts_correlation_fn(potts.lattice, L, 3)
        end
        ss_corr[:, tidx] /= n_steps
        ss_corr[2:end, tidx] = (ss_corr[2:end, tidx] + ss_corr[end:-1:2, tidx]) / 2
        next!(Pmeter)
    end
    finish!(Pmeter)

    # save data
    open(joinpath("data", "2DModel", "correlations", "potts_corrfun_L$(L).txt"), "w") do f
        writedlm(f, ss_corr)
    end
end