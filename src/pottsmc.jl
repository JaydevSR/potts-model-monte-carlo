# TYPE DEFINITIONS AND HELPER FUNCTIONS

abstract type AbstractPottsModel{T} end

mutable struct PottsModel2D{T} <: AbstractPottsModel{T} where T <: Real
    L::Int64
    q::Int64
    d::Int64
    lattice::AbstractArray{T, 2}
end

mutable struct PottsModel3D{T} <: AbstractPottsModel{T} where T <: Real
    L::Int64
    q::Int64
    d::Int64
    lattice::AbstractArray{T, 3}
end

function initialize_model_2d(L, q, cold_start::Bool=false)
    if cold_start
        lattice=fill(1, (L, L)) # Third dimension is flattened along the column
    else
        lattice=rand(0:q-1, (L, L))
    end
    PottsModel2D{Int64}(L, q, 2, lattice)
end

function initialize_model_3d(L, q, cold_start::Bool=false)
    if cold_start
        lattice=fill(1, (L, L, L))
    else
        lattice=rand(0:q-1, (L, L, L))
    end
    PottsModel3D{Int64}(L, q, 3, lattice)
end

function get_nearest_neighbors(model::PottsModel2D, k::CartesianIndex)
    shifts = CartesianIndex.([
        (1, 0), (model.L - 1, 0), 
        (0, 1), (0, model.L - 1)
        ])
    nnbrs = [(k + δ) for δ in shifts]
    nnbrs = [CartesianIndex(mod1.(Tuple(nn), model.L)) for nn in nnbrs]  # Periodic boundary conditions
    return nnbrs
end

function get_nearest_neighbors(model::PottsModel3D, k::CartesianIndex)
    shifts = CartesianIndex.([
        (1, 0, 0), (model.L - 1, 0, 0),
        (0, 1, 0), (0, model.L - 1, 0),
        (0, 0, 1), (0, 0, model.L - 1)
        ])
    nnbrs = [(k + δ) for δ in shifts]
    nnbrs = [CartesianIndex(mod1.(Tuple(nn), model.L)) for nn in nnbrs]  # Periodic boundary conditions
    return nnbrs
end

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
        flip_val = mod1(rand((kval + 1):(kval + model.q - 1)), model.q)
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
    seed = CartesianIndex(Tuple(rand(1:model.L, model.d)))
    stack = [seed]
    sval = model.lattice[seed]
    # choose a random spin out of other values
    new_val = mod1(rand((sval + 1):(sval + model.q - 1)), model.q)
    cluster[seed] = true
    while !isempty(stack)
        k = pop!(stack)
        kval = model.lattice[k]
        model.lattice[k] = new_val  # set new value
        nnbrs = get_nearest_neighbors(model, k)
        for nn ∈ nnbrs
            nnval = model.lattice[nn]
            if kval == nnval && !cluster[nn] && rand() < P_add
                push!(stack, nn)
                cluster[nn] = true
            end
        end
    end
end
