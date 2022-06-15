function potts_getconfigdata_to_txt(
    lattice_sizes::AbstractArray{Int64},
    temps::AbstractArray{Float64},
    q::Int64,
    d::Int64,
    nconfigs::Int64,
    eqsteps::Int64,
    autocorr_times::AbstractArray{Int64}=ones(Int64, length(Temps));
    store_at::AbstractString="",
    start::Symbol=:cold,
    verbose=true,
    ntau=2
    )
    for L in lattice_sizes
        verbose && println(".==================================")
        verbose && println("| Lattice Size: $(L) x $(L)        ")
        verbose && println(".==================================")
        verbose && println("|  ")

        ispath(store_at*"Size$L") ? 1 : mkpath(store_at*"Size$L")
        Threads.@threads for stepT in 1:length(temps)
            filename="Size$L/potts_uncorr_configs_temp$(temps[stepT])_L$(L).txt"
            T = temps[stepT]
            verbose && println("| Process strarted on thread #$(Threads.threadid()) (T = $(T)).")
            
            if d==2
                potts = PottsModel2D(L, q, start)
            elseif d==3
                potts = PottsModel3D(L, q, start)
            else
                error("No model available for d=$d.")
            end

            for i=1:eqsteps  # equilibration
                wolff_cluster_update!(potts, T)
            end

            uncorrelated_spins = zeros(Int64, (potts.L^potts.d, nconfigs))

            nτ = ntau*autocorr_times[stepT]
            numsteps = nτ*nconfigs
            el = @elapsed for j in 1:numsteps
                wolff_cluster_update!(potts, T)
                if j%nτ == 0
                    uncorrelated_spins[:, j÷nτ] = reshape(potts.lattice, potts.L^potts.d)
                end
            end
            open(store_at*filename, "w") do io
                writedlm(io, uncorrelated_spins, ',')
            end;
            verbose && println("| Process complete on thread #$(Threads.threadid()) (T = $T) in $el seconds.")
        end
        verbose && println("| Done.")
        verbose && println(".==================================")
    end
    nothing
end

function potts_getmagdata_to_txt(
    lattice_sizes::AbstractArray{Int64},
    temps::AbstractArray{Float64},
    q::Int64,
    d::Int64,
    nconfigs::Int64,
    eqsteps::Int64,
    autocorr_times::AbstractArray{Int64}=ones(Int64, length(Temps));
    store_at::AbstractString="",
    start::Symbol=:cold,
    verbose=true,
    ntau=2
    )
    for L in lattice_sizes
        verbose && println(".==================================")
        verbose && println("| Lattice Size: $(L) x $(L)        ")
        verbose && println(".==================================")
        verbose && println("|  ")

        ispath(store_at*"Size$L") ? 1 : mkpath(store_at*"Size$L")
        Threads.@threads for stepT in 1:length(temps)
            filename="Size$L/potts_mags_temp$(temps[stepT])_L$(L).txt"
            T = temps[stepT]
            verbose && println("| Process strarted on thread #$(Threads.threadid()) (T = $(T))")
            
            if d==2
                potts = PottsModel2D(L, q, start)
            elseif d==3
                potts = PottsModel3D(L, q, start)
            else
                error("No model available for d=$d.")
            end

            for i=1:eqsteps  # equilibration
                wolff_cluster_update!(potts, T)
            end

            mags = zeros(Float64, nconfigs)

            nτ = ntau*autocorr_times[stepT]
            numsteps = nτ*nconfigs
            el = @elapsed for j in 1:numsteps
                wolff_cluster_update!(potts, T)
                if j%nτ == 0
                    mags[j÷nτ] = magnetisation(potts)
                end
            end
            open(store_at*filename, "w") do io
                writedlm(io, mags, ',')
            end;
            verbose && println("| Process complete on thread #$(Threads.threadid()) (T = $T) in $el seconds.")
        end
        verbose && println("| Done.")
        verbose && println(".==================================")
    end
    nothing
end