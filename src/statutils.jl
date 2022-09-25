"""
    autocorrelation_fn(mags, N)

Calculate the autocorrelation function (normalized) of the given time series array.
"""
function autocorrelation_fn(series)
    tmax = length(series)
    autocorr = zeros(Float64, tmax)
    for t ∈ 1:tmax-1
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for tk ∈ 1:tmax-t
            sum1 += series[tk] * series[tk+t]
            sum2 += series[tk]
            sum3 += series[tk+t]
        end
        autocorr[t] = sum1 / (tmax - t) - (sum2 * sum3) / (tmax - t)^2
    end
    @. autocorr /= autocorr[1]
    return autocorr
end


"""
    bootstrap_err(samples, calc_qty; r=100)

Estimate the error in the given samples by bootstrap method.
Here, `calc_qty` is the function to calculate the quantity in which error has to be calculated.
And, `r` is a keyword arguments giving number of resamples.
"""
function bootstrap_err(samples, calc_qty, args...; r = 200)
    nob = length(samples)
    resample_arr = zeros(Float64, r)
    resample = zeros(eltype(samples), nob)
    for i = 1:r
        resample .= rand(samples, nob)
        resample_arr[i] = calc_qty(resample, args...)
    end
    err = std(resample_arr, corrected = false)
    return err
end


"""
    blocking_err(samples, calc_qty; blocks=20)

Estimate the error in the given samples by blocking method.
Here, `calc_qty` is the function to calculate the quantity in which error has to be calculated.
And, `blocks` is a keyword arguments giving number of blocks.
"""
function blocking_err(samples, calc_qty, args...; blocks = 20)
    block_array = zeros(Float64, blocks)
    blocklength = length(samples) ÷ blocks
    for i = 1:blocks
        sample_block = samples[(i-1)*blocklength+1:i*blocklength]
        block_array[i] = calc_qty(sample_block, args...)
    end
    err = std(block_array)
    return err
end

"""
Calculate the k'th cumulant of given samples.
"""
function cumulant(samples::Array{<:Real}, k::UnitRange{Int}, m::Real)
    cmoms = zeros(Float64, last(k))
    cumls = zeros(Float64, last(k))
    cmoms[1] = 0
    cumls[1] = m
    for i=2:last(k)
        @inbounds cmoms[i] =  moment(samples, i, m)
        @inbounds kn = cmoms[i]
        for j=2:i-2
            @inbounds kn -= binomial(i-1, j)*cmoms[j]*cumls[i-j]
        end
        @inbounds cumls[i] = kn
    end
    return cumls[k]
end

cumulant(samples::Array{<:Real}, k::Int, m::Real) = cumulant(samples, 1:k, m)[k]
cumulant(samples::Array{<:Real}, k::Union{Int, UnitRange{Int}}) = cumulant(samples, k, mean(samples))

"""
Calculate Binder's cumulant for given array of magnetization (per spin).
"""
function binders_cumulant(marr)
    m4 = m2 = 0
    for m in marr
        m4, m2 = m4 + m^4, m2 + m^2
    end
    return 1 - (m4 / (3*m2*m2))
end

"""
Calculate two-point correlation function for given lattice of spins.
"""
function ss_correlation_fn(sites::Matrix, L=Int64; metric=*)
    ss_corrs = zeros(Float64, L)
    nsamples = zeros(Float64, L)
    for i=1:L
        for j=1:L
            r = abs(i-j)
            ss_corrs[r+1] += metric(sites[i,i], sites[i,j]) 
            ss_corrs[r+1] += metric(sites[i,i], sites[j,i])
            nsamples[r+1] += 2
        end
    end
    return ss_corrs ./ nsamples
end
