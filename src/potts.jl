mutable struct PottsModel2D
    L::Int64
    q::Int64
    lattice::AbstractArray
end

mutable struct PottsModel3D
    L::Int64
    q::Int64
    lattice::AbstractArray
end

function initialize_model_2d(L, q, cold_start::Bool=false)
    if cold_start
        lattice=fill(1, (L, L)) # Third dimension is flattened along the column
    else
        lattice=rand(0:q-1, (L, L))
    end

    return PottsModel2D(L, q, lattice)
end

function initialize_model_3d(L, q, cold_start::Bool=false)
    if cold_start
        lattice=fill(1, (L, L*L)) # Third dimension is flattened along the column
    else
        lattice=rand(0:q-1, (L, L*L))
    end

    return PottsModel3D(L, q, lattice)
end

function wolff_cluster_update(model::PottsModel2D, temp::Float64)
    P_add = 1 - exp(-1/temp)
    cluster = falses(size(model.lattice))
    seed = CartesianIndex(rand(1:model.L, 2))
    stack = [seed]
    cluster[seed] = true
    while !isempty(stack)
        k = pop!(stack)
        kval = model.lattice[k]
        # choose a random spin out of other values
        model.lattice[k] = mod1(rand((kval + 1):(kval + model.q - 1)), model.q) ## FLAG
        nnbrs = get_nearest_neighbors(model, k)
        @inbounds for nn ∈ nnbrs
            nnval = model.lattice[nn]
            if kval == nnval && !cluster[nn] && rand() < P_add
                push!(stack, nn)
                cluster[nn] = true
            end
        end
    end
end

function get_nearest_neighbors(model::PottsModel2D, k::CartesianIndex)
    shifts = CartesianIndex.([(1, 0), (model.L - 1, 0), (0, 1), (0, model.L - 1)])
    nnbrs = [(k + δ) for δ in shifts]
    nnbrs = [CartesianIndex(mod1.(Tuple(nn), model.L)) for nn in nnbrs]  # Apply periodic boundary conditions
    return nnbrs
end

function get_nearest_neighbors(model::PottsModel3D, k)
    ## TODO 
end