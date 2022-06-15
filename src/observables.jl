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
    magnetisation(model::AbstractPottsModel)

Calculate the total magnetisation for q-state potts model.
"""

function magnetisation(model::AbstractPottsModel)
    counts = zeros(Int64, model.q)
    for site in CartesianIndices(model.lattice)
        counts[model.lattice[site] + 1] += 1
    end
    return (model.q/(model.q - 1)) * maximum(counts .- model.L^model.d / model.q)
end