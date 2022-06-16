# TYPE DEFINITIONS AND HELPER FUNCTIONS

abstract type AbstractPottsModel end

mutable struct PottsModel2D{T<:Integer} <: AbstractPottsModel
    L::T
    q::T
    d::T
    lattice::AbstractArray{T, 2}
    function PottsModel2D{T}(L::T, q::T, start::Symbol) where T<:Integer
        if start==:cold
            lattice=fill(zero(T), (L, L))
        elseif start==:hot
            lattice=rand(zero(T):zero(T)+q-1, (L, L))
        else
            error("Start state can be one of symbols :$(:cold) or :$(:hot)")
        end
        new{T}(L, q, 2, lattice)
    end
end

PottsModel2D(L::T, q::T, start::Symbol=:cold) where {T<:Integer} = PottsModel2D{T}(L, q, start)

mutable struct PottsModel3D{T<:Integer} <: AbstractPottsModel
    L::T
    q::T
    d::T
    lattice::AbstractArray{T, 3}
    function PottsModel3D{T}(L::T, q::T, start::Symbol) where T<:Integer
        if start==:cold
            lattice=fill(zero(T), (L, L, L))
        elseif start==:hot
            lattice=rand(zero(T):zero(T)+q-1, (L, L, L))
        else
            error("Start state can be one of symbols :$(:cold) or :$(:hot)")
        end
        new{T}(L, q, 3, lattice)
    end
end

PottsModel3D(L::T, q::T, start::Symbol=:cold) where {T<:Integer} = PottsModel3D{T}(L, q, start)

function get_nearest_neighbors(model::PottsModel2D, k::Array)
    shifts = [
        [1, 0], [model.L - 1, 0], 
        [0, 1], [0, model.L - 1]
        ]
    nnbrs = [k + δ for δ in shifts]
    nnbrs = [mod1.(nn, model.L) for nn in nnbrs]  # Periodic boundary conditions
    return nnbrs
end

function get_nearest_neighbors(model::PottsModel3D, k::CartesianIndex)
    shifts = [
        [1, 0, 0], [model.L - 1, 0, 0],
        [0, 1, 0], [0, model.L - 1, 0],
        [0, 0, 1], [0, 0, model.L - 1]
        ]
    nnbrs = [k + δ for δ in shifts]
    nnbrs = [mod1.(nn, model.L) for nn in nnbrs]  # Periodic boundary conditions
    return nnbrs
end

# function get_nearest_neighbors(model::PottsModel2D, k::CartesianIndex)
#     shifts = CartesianIndex.([
#         (1, 0), (model.L - 1, 0), 
#         (0, 1), (0, model.L - 1)
#         ])
#     nnbrs = [(k + δ) for δ in shifts]
#     nnbrs = [CartesianIndex(mod1.(Tuple(nn), model.L)) for nn in nnbrs]  # Periodic boundary conditions
#     return nnbrs
# end

# function get_nearest_neighbors(model::PottsModel3D, k::CartesianIndex)
#     shifts = CartesianIndex.([
#         (1, 0, 0), (model.L - 1, 0, 0),
#         (0, 1, 0), (0, model.L - 1, 0),
#         (0, 0, 1), (0, 0, model.L - 1)
#         ])
#     nnbrs = [(k + δ) for δ in shifts]
#     nnbrs = [CartesianIndex(mod1.(Tuple(nn), model.L)) for nn in nnbrs]  # Periodic boundary conditions
#     return nnbrs
# end