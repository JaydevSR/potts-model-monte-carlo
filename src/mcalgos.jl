# METROPOLIS SINGLE FLIP ALGORITHM

@inbounds function metropolis_batch_update!(
                        model::PottsModel2D,
                        temp::Float64;
                        fix_vacuum::Bool=true,
                        accept_probs::Dict=Dict(append!([(i, 1.0) for i=-4:0], [(i, exp(-i/temp)) for i=1:4])))
    for _ in 1:model.L^2
        kx, ky = rand(1:model.L), rand(1:model.L)
        kval = model.lattice[kx, ky]
        flip_val = mod(rand((kval + 1):(kval + model.q - 1)), model.q)

        # compute δE
        s1, s2 = 0, 0
        for δ in model.shifts
            nnx, nny = mod1(kx + δ[1], model.L), mod1(ky + δ[2], model.L)
            nnval = model.lattice[nnx, nny]
            s1 += (nnval == kval)
            s2 += (nnval == flip_val)
        end
        δE = convert(Int64, (s1 - s2))
        if rand() < accept_probs[δE]
            model.lattice[kx, ky] = flip_val
            model.counts[flip_val+1] += 1
            model.counts[kval+1] -= 1
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

function δE_single_flip(model::AbstractPottsModel, flip_site::CartesianIndex, flip_val::Int64)
    nnbrs = get_nearest_neighbors(model, flip_site)
    @inbounds nnbrs_vals = model.lattice[nnbrs]
    s1 = sum(nnbrs_vals .== model.lattice[flip_site])
    s2 = sum(nnbrs_vals .== flip_val)
    return convert(Int64, -(s2 - s1))
end

# WOLFF CLUSTER ALGORITHM
function wolff_cluster_update!(
                    model::PottsModel2D,
                    temp::Float64;
                    P_add::Float64=-expm1(-1/temp), fix_vacuum::Bool=true,
                    stack::LazyStack{Int}=LazyStack(Int),
                    cluster::BitMatrix=falses(model.L, model.L))
    empty!(stack)

    if size(cluster) != size(model.lattice)
        error("Cluster must be of size $(size(model.lattice))")
    end
    cluster .= false
    
    seedx, seedy = rand(1:model.L), rand(1:model.L)
    push!(stack, seedx)
    push!(stack, seedy)
    sval = model.lattice[seedx, seedy]
    new_val = mod(rand((sval + 1):(sval + model.q - 1)), model.q)
    cluster[seedx, seedy] = true

    while !isempty(stack)
        ky = pop!(stack)
        kx = pop!(stack)
        @inbounds kval = model.lattice[kx, ky]
        @inbounds model.lattice[kx, ky] = new_val  # set new value
        model.counts[new_val+1] += 1
        model.counts[kval+1] -= 1
        for δ ∈ model.shifts
            nnx, nny = mod1(kx + δ[1], model.L), mod1(ky + δ[2], model.L)
            @inbounds nnval = model.lattice[nnx, nny]
            if kval == nnval && !cluster[nnx, nny] && rand() < P_add
                push!(stack, nnx)
                push!(stack, nny)
                @inbounds cluster[nnx, nny] = true
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

function wolff_cluster_update!(
                    model::PottsModel3D,
                    temp::Float64;
                    P_add::Float64=-expm1(-1/temp), fix_vacuum::Bool=true,
                    stack::LazyStack{Int}=LazyStack(Int),
                    cluster::BitArray{3}=falses(model.L, model.L, model.L)) 
    empty!(stack)

    if size(cluster) != size(model.lattice)
        error("Cluster must be of size $(size(model.lattice))")
    end
    cluster .= false
    
    seedx, seedy, seedz = rand(1:model.L), rand(1:model.L), rand(1:model.L)
    push!(stack, seedx)
    push!(stack, seedy)
    push!(stack, seedz)
    sval = model.lattice[seedx, seedy, seedz]
    new_val = mod(rand((sval + 1):(sval + model.q - 1)), model.q)
    cluster[seedx, seedy, seedz] = true

    while !isempty(stack)
        kz = pop!(stack)
        ky = pop!(stack)
        kx = pop!(stack)
        @inbounds kval = model.lattice[kx, ky, kz]
        @inbounds model.lattice[kx, ky, kz] = new_val  # set new value
        model.counts[new_val+1] += 1
        model.counts[kval+1] -= 1
        for δ ∈ model.shifts
            nnx, nny, nnz = mod1(kx + δ[1], model.L), mod1(ky + δ[2], model.L), mod1(kz + δ[3], model.L)
            @inbounds nnval = model.lattice[nnx, nny, nnz]
            if kval == nnval && !cluster[nnx, nny, nnz] && rand() < P_add
                push!(stack, nnx)
                push!(stack, nny)
                push!(stack, nnz)
                @inbounds cluster[nnx, nny, nnz] = true
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