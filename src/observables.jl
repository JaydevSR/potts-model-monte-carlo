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
    M=0
    counts = zeros(Int64, model.q)
    for site in CartesianIndices(model.lattice)
        counts[model.lattice[site] + 1] += 1
    end
    return (model.q/(model.q - 1)) * maximum(counts .- model.L^model.d / model.q)
end

"""
    specific_heat(u_vals, T, ns)

Calculate the specific heat from given array of internal energy per site (total `ns` sites) at temperature `T`.
"""
function specific_heat(u_vals, T, ns)
    return (T^-2) * ns * var(u_vals, corrected = false)
end


"""
    succeptibility(m_vals, T, N)

Calculate the succeptibility from given array of mean magnetization per site (total `ns` sites) at temperature `T`.
"""
function succeptibility(m_vals, T, ns)
    return (T^-2) * ns * var(m_vals, corrected = false)
end