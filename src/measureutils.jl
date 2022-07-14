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
    ntau=2,
    mode="w",
    fix_vacuum=true
    )
    for L in lattice_sizes
        verbose && println(".==================================")
        verbose && println("| Lattice Size: $(L) x $(L)        ")
        verbose && println(".==================================")
        verbose && println("|  ")
        @sync for stepT in 1:length(temps)
            Threads.@spawn potts_getconfigdata_to_txt(
                L, temps[stepT], q, d, 
                nconfigs, eqsteps, autocorr_times[stepT];
                store_at=store_at, start=start, ntau=ntau,
                mode=mode, fix_vacuum=fix_vacuum
                )
        end
        verbose && println("| Done.")
        verbose && println(".==================================")
    end
    nothing
end

function potts_getconfigdata_to_txt(
    L::Int,
    T::Float64,
    q::Int64,
    d::Int64,
    nconfigs::Int64,
    eqsteps::Int64,
    τ::Int64=1;
    store_at::AbstractString="",
    start::Symbol=:cold,
    ntau=2,
    mode="w",
    verbose=true,
    fix_vacuum=true
    )
    verbose && println("| Process strarted on thread #$(Threads.threadid()) (T = $(T)).")
    el = @elapsed begin
        szpath = joinpath([store_at, "Size$L"])
        ispath(szpath) ? 1 : mkpath(szpath)
        
        if d==2
            potts = PottsModel2D(L, q, start)
        elseif d==3
            potts = PottsModel3D(L, q, start)
        else
            error("No model available for d=$d.")
        end
    
        for i=1:eqsteps  # equilibration
            wolff_cluster_update!(potts, T, fix_vacuum=fix_vacuum)
        end
    
        uncorrelated_spins = zeros(Int64, (potts.L^potts.d, nconfigs))
    
        nτ = ntau*τ
        numsteps = nτ*nconfigs
        el = @elapsed for j in 1:numsteps
            wolff_cluster_update!(potts, T, fix_vacuum=fix_vacuum)
            if j%nτ == 0
                uncorrelated_spins[:, j÷nτ] = reshape(potts.lattice, potts.L^potts.d)
            end
        end
    
        filename="potts_uncorr_configs_temp$(T)_L$(L).txt"
        open(joinpath([szpath, filename]), mode) do io
            writedlm(io, uncorrelated_spins, ',')
        end; 
    end
    verbose && println("| Process complete on thread #$(Threads.threadid()) (T = $T) in $el seconds.")
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
    ntau=2,
    mode="w",
    fix_vacuum=true
    )
    for L in lattice_sizes
        verbose && println(".==================================")
        verbose && println("| Lattice Size: $(L) x $(L)        ")
        verbose && println(".==================================")
        verbose && println("|  ")

        szpath = joinpath([store_at, "Size$L"])
        ispath(szpath) ? 1 : mkpath(szpath)
        @sync for stepT in 1:length(temps)
            Threads.@spawn potts_getmagdata_to_txt(
                L, temps[stepT], q, d, 
                nconfigs, eqsteps, autocorr_times[stepT];
                store_at=store_at, start=start, ntau=ntau,
                mode=mode, fix_vacuum=fix_vacuum
                )
        end
        verbose && println("| Done.")
        verbose && println(".==================================")
    end
    nothing
end

function potts_getmagdata_to_txt(
    L::Int,
    T::Float64,
    q::Int64,
    d::Int64,
    nconfigs::Int64,
    eqsteps::Int64,
    τ::Int64=1;
    store_at::AbstractString="",
    start::Symbol=:cold,
    ntau=2,
    mode="w",
    verbose=true,
    fix_vacuum=true
    )
    verbose && println("| Process strarted on thread #$(Threads.threadid()) (T = $(T)).")
    el = @elapsed begin
        szpath = joinpath([store_at, "Size$L"])
        ispath(szpath) ? 1 : mkpath(szpath)

        stk = LazyStack(typeof(CartesianIndex(Tuple(rand(1:L, d)))))
        
        if d==2
            potts = PottsModel2D(L, q, start)
        elseif d==3
            potts = PottsModel3D(L, q, start)
        else
            error("No model available for d=$d.")
        end
    
        for i=1:eqsteps  # equilibration
            wolff_cluster_update!(potts, T, fix_vacuum=fix_vacuum)
        end
    
        mags = zeros(Float64, nconfigs)
    
        nτ = ntau*τ
        numsteps = nτ*nconfigs
        el = @elapsed for j in 1:numsteps
            wolff_cluster_update!(potts, T, fix_vacuum=fix_vacuum)
            if j%nτ == 0
                mags[j÷nτ] = magnetization(potts)
            end
        end
    
        filename="potts_mags_temp$(T)_L$(L).txt"
        open(joinpath([szpath, filename]), mode) do io
            writedlm(io, mags, ',')
        end; 
    end
    verbose && println("| Process complete on thread #$(Threads.threadid()) (T = $T) in $el seconds.")
    nothing 
end