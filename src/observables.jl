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

function magnetization(model::AbstractPottsModel; max_definition=true)
    counts = zeros(Int64, model.q)
    for site in CartesianIndices(model.lattice)
        counts[model.lattice[site] + 1] += 1
    end
    if max_definition
        return Float64((model.q // (model.q - 1)) * maximum(counts .- model.L^model.d // model.q))
    else
        return real(sum(counts[p+1]*exp(2*pi*im*p/model.q) for p in 0:model.q-1))
    end
end

function magnetization(lattice::AbstractArray, L::Int, q::Int, d::Int; max_definition=true)
    counts = zeros(Int64, q)
    for site in CartesianIndices(lattice)
        counts[lattice[site] + 1] += 1
    end
    if max_definition
        return Float64((q // (q - 1)) * maximum(counts .- L^d // q))
    else
        return real(sum(counts[p+1]*exp(2*pi*im*p//q) for p in 0:q-1))
    end
end