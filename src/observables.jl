"""
    hamiltonian(model::AbstractPottsModel)

Calculate the Hamiltonian for q-state potts model (J=1).
"""
function hamiltonian(model::AbstractPottsModel)
    H=0
    for site in CartesianIndices(model.lattice)
        nnbrs = get_nearest_neighbors(model, site)
        for nn in nnbrs
            if model.lattice[nn] == model.lattice[site]
                H -= 1
            end
        end
    end
    return H / 2  # because of double counting
end

"""
    magnetization(model::AbstractPottsModel)

Calculate the total magnetization for q-state potts model.
"""

function magnetization(model::AbstractPottsModel; use_definition::Symbol=:max)
    if use_definition == :max
        return Float64((model.q // (model.q - 1)) * maximum(model.counts .- model.L^model.d // model.q))
    elseif use_definition == :vac
        return (2/3)*(sum(model.counts[p+1] * cos(2*pi*p/model.q) for p=0:model.q-1) + (model.L^model.d/2))
    elseif use_definition == :both
        M_max = Float64((model.q // (model.q - 1)) * maximum(model.counts .- model.L^model.d // model.q))
        M_vac = (2/3)*(sum(model.counts[p+1] * cos(2*pi*p/model.q) for p=0:model.q-1) + (model.L^model.d/2))
        return M_max, M_vac
    else
        error("Unknown magnetization definition")
    end
end

function magnetization(lattice::AbstractArray, L::Int, q::Int, d::Int; use_definition::Symbol=:max)
    counts = zeros(Int64, q)
    for site in CartesianIndices(lattice)
        counts[lattice[site] + 1] += 1
    end
    if use_definition == :max
        return Float64((q // (q - 1)) * maximum(counts .- L^d // q))
    elseif use_definition == :vac
        return (2/3)*(sum(counts[p+1] * cos(2*pi*p/q) for p=0:q-1) + (L^d/2))
    elseif use_definition == :both
        M_max = Float64((q // (q - 1)) * maximum(counts .- L^d // q))
        M_vac = (2/3)*(sum(counts[p+1] * cos(2*pi*p/q) for p=0:q-1) + (L^d/2))
        return M_max, M_vac
    else
        error("Unknown magnetization definition")
    end
end