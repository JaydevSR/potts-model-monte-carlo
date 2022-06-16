# METROPOLIS SINGLE FLIP ALGORITHM

function metropolis_batch_update!(model::AbstractPottsModel, temp::Float64)
    accept_probs = Dict(
        append!(
            [(i, 1.0) for i=-2*model.d:0], 
            [(i, exp(-i / temp)) for i=1:2*model.d]
        ))
    ΔE = 0
    for i = 1:model.L^model.d
        k = CartesianIndex(Tuple(rand(1:model.L, model.d)))
        kval = model.lattice[k]
        flip_val = mod(rand((kval + 1):(kval + model.q - 1)), model.q)
        δE = δE_single_flip(model, k, flip_val)
        if rand() < accept_probs[δE]
            model.lattice[k] = flip_val
            ΔE += δE 
        end
    end
    return ΔE
end

function δE_single_flip(model::AbstractPottsModel, flip_site::CartesianIndex, flip_val::Int64)
    nnbrs = get_nearest_neighbors(model, flip_site)
    nnbrs_vals = [model.lattice[i] for i in nnbrs]
    s1 = sum(nnbrs_vals .== model.lattice[flip_site])
    s2 = sum(nnbrs_vals .== flip_val)
    return convert(Int64, -(s2 - s1))
end

# WOLFF CLUSTER ALGORITHM

function wolff_cluster_update!(model::AbstractPottsModel, temp::Float64)
    P_add = 1 - exp(-1/temp)
    cluster = falses(size(model.lattice))
    seed = rand(1:model.L, model.d)

    # initialize static stack
    ptrcol=[1]
    stack = zeros(Int64, (model.d, model.L^model.d))
    stack[:, ptrcol] = seed
    cluster[seed...] = true

    # choose a random spin out of other values
    sval = model.lattice[seed...]
    new_val = mod(rand((sval + 1):(sval + model.q - 1)), model.q)
    while ptrcol[] > 0
        k, ptrcol[] = stack[:, ptrcol], ptrcol[]-1  # pop
        kval = model.lattice[k...]
        model.lattice[k...] = new_val  # set new value
        nnbrs = get_nearest_neighbors(model, k)
        for nn ∈ nnbrs
            nnval = model.lattice[nn...]
            if kval == nnval && !cluster[nn...] && rand() < P_add
                stack[:, ptrcol[]+1], ptrcol[] = nn, ptrcol[]+1  # push
                cluster[nn...] = true
            end
        end
    end
end
