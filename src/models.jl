# TYPE DEFINITIONS AND HELPER FUNCTIONS

abstract type AbstractPottsModel{D} end

mutable struct PottsModel2D{T, LL, Q, D} <: AbstractPottsModel{D}
    L::Int
    q::Int
    d::Int
    lattice::Array{T, 2}
    counts::Vector{Int}
    shifts::NTuple{4, CartesianIndex{2}}
end

function PottsModel2D(L::Int, q::Int, start::Symbol=:cold)
    if start==:cold
        lattice=fill(0, (L, L))
        counts=zeros(Int, q)
        counts[1]=L*L
    elseif start==:hot
        lattice=rand(0:q-1, (L, L))
        counts=zeros(Int, q)
        for i in eachindex(lattice)
            @inbounds counts[lattice[i]+1] += 1
        end
    else
        error("Start state can be one of symbols :$(:cold) or :$(:hot)")
    end
    shifts = (
        CartesianIndex(1, 0), CartesianIndex(L - 1, 0), 
        CartesianIndex(0, 1), CartesianIndex(0, L - 1)
    )
    return PottsModel2D{Int, L, q, 2}(L, q, 2, lattice, counts, shifts)
end

mutable struct PottsModel3D{T, LL, Q, D} <: AbstractPottsModel{D}
    L::Int
    q::Int
    d::Int
    lattice::Array{T, 3}
    counts::Vector{Int}
    shifts::NTuple{6, CartesianIndex{3}}
end

function PottsModel3D(L::Int, q::Int, start::Symbol=:cold)
    if start==:cold
        lattice=fill(0, (L, L, L))
        counts=zeros(Int, q)
        counts[1]=L*L*L
    elseif start==:hot
        lattice=rand(0:q-1, (L, L, L))
        counts=zeros(Int, q)
        for i in eachindex(lattice)
            @inbounds counts[lattice[i]+1] += 1
        end
    else
        error("Start state can be one of symbols :$(:cold) or :$(:hot)")
    end
    shifts = (
        CartesianIndex(1, 0, 0), CartesianIndex(L - 1, 0, 0), 
        CartesianIndex(0, 1, 0), CartesianIndex(0, L - 1, 0),
        CartesianIndex(0, 0, 1), CartesianIndex(0, 0, L - 1)
    )
    return PottsModel3D{Int, L, q, 3}(L, q, 2, lattice, counts, shifts)
end

function get_nearest_neighbors(model::AbstractPottsModel, k::CartesianIndex)
    ns = length(model.shifts)
    return ntuple(i -> CartesianIndex(mod1(k[1] + model.shifts[i][1], model.L), mod1(k[2]+model.shifts[i][2], model.L)), ns)
end

function get_projection(model::AbstractPottsModel, proj_dir::Int=0)
    cos_dict = Tuple(cos(2 * π * mod(i - proj_dir, model.q) / model.q) for i=0:model.q-1)
    return Folds.map(s -> cos_dict[s+1], model.lattice)
end

function get_projection(lattice::AbstractArray{Int}, q::Int, proj_dir::Int=0)
    cos_vals = Tuple(cos(2 * π * mod(i - proj_dir, q) / q) for i=0:q-1)
    return Folds.map(s -> cos_vals[s+1], lattice)
end