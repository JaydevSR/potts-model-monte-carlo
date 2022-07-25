include("../../src/pottsmc.jl")

##
Threads.nthreads()
lattice_sizes = [40]
temps = [0.5, 0.7, 0.9, 0.94, 0.98, 0.984, 0.988, 0.992, 0.996, 1.0, 1.004, 1.008, 1.012, 1.016, 1.02, 1.06,  1.1, 1.3, 1.5]
τvals = [  3,   5,  10,   25,   50,    60,    80,    90,   110, 115,   105,    95,    80,    60,   40,   20,   10,   5,   3]
##
q=3
d=2
nconfigs=5000
eqsteps=5000
store_at="data/metro/"
ntau=2
start=:cold
mag_definition=:both
get_configs=false
fix_vacuum=true
verbose=true
##

@time for L in lattice_sizes
    verbose && println(".==================================")
    verbose && println("| Lattice Size: $(L) x $(L)        ")
    verbose && println(".==================================")
    verbose && println("|  ")
    @sync for stepT in 1:length(temps)
        Threads.@spawn begin
            T = temps[stepT]
            τ = τvals[stepT]
            accept_probs=Dict(append!([(i, 1.0) for i=-4:0], [(i, exp(-i/T)) for i=1:4]))
            verbose && println("| Process strarted for T = $(T).")
            if d==2
                potts = PottsModel2D(L, q, start)
            elseif d==3
                potts = PottsModel3D(L, q, start)
            else
                error("No model available for d=$d.")
            end
    
            szpath = joinpath([store_at, "$(d)DModel", "Size$L"])
            mkpath(szpath)
        
            for i=1:eqsteps  # equilibration
                metropolis_batch_update!(potts, T, fix_vacuum=fix_vacuum, accept_probs=accept_probs)
            end
        
            if get_configs
                uncorrelated_spins = zeros(Int64, (potts.L^potts.d, nconfigs))
            end
    
            if mag_definition == :both
                mags = zeros(Float64, (2, nconfigs))
            else
                mags = zeros(Float64, (1, nconfigs))
            end
    
            nτ = ntau*τ
            numsteps = nτ*nconfigs
            el = @elapsed for j in 1:numsteps
                metropolis_batch_update!(potts, T, fix_vacuum=fix_vacuum, accept_probs=accept_probs)
                if j%nτ == 0
                    mags[:, j÷nτ] .= magnetization(potts, use_definition=mag_definition)
                    if get_configs
                        uncorrelated_spins[:, j÷nτ] = reshape(potts.lattice, potts.L^potts.d)
                    end
                end
            end
        
            if get_configs
                filename_conf="potts_uncorr_configs_temp$(T)_L$(L).txt"
                mkpath(joinpath([szpath, "uncorr_configs"]))
                open(joinpath([szpath, "uncorr_configs", filename_conf]), "w") do io
                    writedlm(io, uncorrelated_spins, ',')
                end; 
            end
    
            filename_mag="potts_mags_temp$(T)_L$(L).txt"
            mkpath(joinpath([szpath, "mags"]))
            open(joinpath([szpath, "mags", filename_mag]), "w") do io
                writedlm(io, mags, ',')
            end;
            verbose && println("| Process complete for T = $T in $(round(el, digits=0)) seconds.")
        end
    end
    verbose && println("| Done.")
    verbose && println(".==================================")
end