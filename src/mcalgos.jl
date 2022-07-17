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
        @inbounds kval = model.lattice[k]
        flip_val = mod(rand((kval + 1):(kval + model.q - 1)), model.q)
        δE = δE_single_flip(model, k, flip_val)
        if rand() < accept_probs[δE]
            @inbounds model.lattice[k] = flip_val
            model.counts[flip_val+1] += 1
            model.counts[kval+1] -= 1
            ΔE += δE 
        end
    end
    return ΔE
end

function δE_single_flip(model::AbstractPottsModel, flip_site::CartesianIndex, flip_val::Int64)
    nnbrs = get_nearest_neighbors(model, flip_site)
    @inbounds nnbrs_vals = model.lattice[nnbrs]
    s1 = sum(nnbrs_vals .== model.lattice[flip_site])
    s2 = sum(nnbrs_vals .== flip_val)
    return convert(Int64, -(s2 - s1))
end

# WOLFF CLUSTER ALGORITHM
function wolff_cluster_update!(model::AbstractPottsModel, temp::Float64; fix_vacuum::Bool=true, stack=nothing)
    P_add = 1 - exp(-1/temp)
    if isnothing(stack)
        stack = []
        sizehint!(stack, model.L^model.d)
    else
        empty!(stack)
    end
    cluster = falses(size(model.lattice))
    seed = CartesianIndex(Tuple(rand(1:model.L, model.d)))
    @inbounds sval = model.lattice[seed]
    push!(stack, seed)
    
    # choose a random spin out of other values
    new_val = mod(rand((sval + 1):(sval + model.q - 1)), model.q)
    @inbounds cluster[seed] = true
    while !isempty(stack)
        k = pop!(stack)
        @inbounds kval = model.lattice[k]
        @inbounds model.lattice[k] = new_val  # set new value
        model.counts[new_val+1] += 1
        model.counts[kval+1] -= 1
        nnbrs = get_nearest_neighbors(model, k)
        for nn ∈ nnbrs
            @inbounds nnval = model.lattice[nn]
            if kval == nnval && !cluster[nn] && rand() < P_add
                push!(stack, nn)
                @inbounds cluster[nn] = true
            end
        end
    end

    if fix_vacuum
        current_vacuum = argmax(model.counts) - 1
        rotation = Tuple(mod(s - current_vacuum, model.q) for s=0:model.q-1)
        map!(s -> rotation[s+1], model.lattice, model.lattice)
        unrot_counts = copy(model.counts)
        for i in eachindex(model.counts)
            model.counts[rotation[i]+1] = unrot_counts[i]
        end
    end
    nothing
end