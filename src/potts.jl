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

function wolff_cluster_update!(model::PottsModel2D, temp::Float64)
    P_add = 1 - exp(-1/temp)
    cluster = falses(size(model.lattice))
    seed = CartesianIndex(Tuple(rand(1:model.L, 2)))
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

# Correct
function metropolis_batch_update!(model::PottsModel2D, temp::Float64)
    accept_probs = Dict(append!([(i, 1.0) for i=-4:0], [(i, exp(-i / temp)) for i=1:4]))
    change_in_energy = 0
    for i = 1:model.L^2
        k = CartesianIndex(Tuple(rand(1:model.L, 2)))
        kval = model.lattice[k]
        flip_val = mod1(rand((kval + 1):(kval + model.q - 1)), model.q)
        ΔE = deltaE_single_flip(model, k, flip_val)  # Flag 1
        if rand() < accept_probs[ΔE]  # Flag 2
            model.lattice[k] = flip_val
            change_in_energy += ΔE
        end
    end
    return change_in_energy
end

# Correct
function get_nearest_neighbors(model::PottsModel2D, k::CartesianIndex)
    shifts = CartesianIndex.([(1, 0), (model.L - 1, 0), (0, 1), (0, model.L - 1)])
    nnbrs = [(k + δ) for δ in shifts]
    nnbrs = [CartesianIndex(mod1.(Tuple(nn), model.L)) for nn in nnbrs]  # Apply periodic boundary conditions
    return nnbrs
end

function get_nearest_neighbors(model::PottsModel3D, k::CartesianIndex)
    ## TODO
end

# Flag: Check correctness
function deltaE_single_flip(model::PottsModel2D, flip_site::CartesianIndex, flip_val::Int64)
    nnbrs = get_nearest_neighbors(model, flip_site)
    nnbrs_vals = [model.lattice[i] for i in nnbrs]
    s1 = sum(nnbrs_vals .== model.lattice[flip_site])
    s2 = sum(nnbrs_vals .== flip_val)
    return convert(Int64, -(s2 - s1))
end

# Correct
function potts_hamiltonian(model::PottsModel2D)
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
    specific_heat(u_vals, T, N)

Calculate the specific heat from given array of internal energy per site (`N²` sites) at temperature `T`.
"""
function specific_heat(u_vals, T, N)
    return (T^-2) * N^2 * var(u_vals, corrected = false)
end


"""
    succeptibility(m_vals, T, N)

Calculate the succeptibility from given array of mean magnetization per site (`N²` sites) at temperature `T`.
"""
function succeptibility(m_vals, T, N)
    return (T^-2) * N^2 * var(m_vals, corrected = false)
end