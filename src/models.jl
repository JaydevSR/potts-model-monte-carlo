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
    nnbrs = SA[[CartesianIndex(mod1.((k+δ).I, model.L)) for δ in model.shifts]...]
    return nnbrs
end