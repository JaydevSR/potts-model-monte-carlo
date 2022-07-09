# TYPE DEFINITIONS AND HELPER FUNCTIONS

abstract type AbstractPottsModel end

mutable struct PottsModel2D{T<:Integer, N, C} <: AbstractPottsModel
    L::T
    q::T
    d::T
    lattice::Array{T, 2}
    shifts::SVector{N, C}
end

function PottsModel2D(L::T, q::T, start::Symbol=:cold) where T <: Integer
    if start==:cold
        lattice=fill(zero(T), (L, L))
    elseif start==:hot
        lattice=rand(zero(T):zero(T)+q-1, (L, L))
    else
        error("Start state can be one of symbols :$(:cold) or :$(:hot)")
    end
    shifts = SA[
        CartesianIndex(1, 0), CartesianIndex(L - 1, 0), 
        CartesianIndex(0, 1), CartesianIndex(0, L - 1)
        ]
    C = eltype(shifts)
    N = length(shifts)
    return PottsModel2D{T, N, C}(L, q, 2, lattice, shifts)
end

mutable struct PottsModel3D{T<:Integer, N, C} <: AbstractPottsModel
    L::T
    q::T
    d::T
    lattice::Array{T, 2}
    shifts::SVector{N, C}
end

function PottsModel3D(L::T, q::T, start::Symbol=:cold) where T <: Integer
    if start==:cold
        lattice=fill(zero(T), (L, L, L))
    elseif start==:hot
        lattice=rand(zero(T):zero(T)+q-1, (L, L, L))
    else
        error("Start state can be one of symbols :$(:cold) or :$(:hot)")
    end
    shifts = SA[
        CartesianIndex(1, 0, 0), CartesianIndex(L - 1, 0, 0), 
        CartesianIndex(0, 1, 0), CartesianIndex(0, L - 1, 0),
        CartesianIndex(0, 0, 1), CartesianIndex(0, 0, L - 1)
        ]
    C = eltype(shifts)
    N = length(shifts)
    return PottsModel3D{T, N, C}(L, q, 2, lattice, shifts)
end

function get_nearest_neighbors(model::AbstractPottsModel, k::CartesianIndex)
    SA[[CartesianIndex(mod1.((k+δ).I, model.L)) for δ in model.shifts]...]
end

function get_projection(model::AbstractPottsModel, proj_dir::Int=0)
    cos_dict = Tuple(cos(2 * π * mod(i - proj_dir, model.q) / model.q) for i=0:model.q-1)
    return Folds.map(s -> cos_dict[s+1], model.lattice)
end

function get_projection(lattice::AbstractArray{Int64}, q::Int64, proj_dir::Int=0)
    cos_vals = Tuple(cos(2 * π * mod(i - proj_dir, q) / q) for i=0:q-1)
    return Folds.map(s -> cos_vals[s+1], lattice)
end