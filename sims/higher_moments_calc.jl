include("../src/pottsmc.jl")

#=
Perform simulation
=#

L = 16
q = 3
start = :cold
err_nblocks = 25
esteps = 1000  # Number of steps for equilibration
nconfigs = 20000  # Number of steps for measurements


data_path = "D:\\Projects\\Potts-QCD\\potts-model-monte-carlo\\data\\corrtimes_data.txt"
data = readdlm(data_path, ',', Float64)
Temps = data[1, :]
autocorr_steps_wolff = convert.(Int64, data[3, :])

println(".==================================")
println("| Lattice Size: $(L) x $(L)        ")
println(".==================================")
println("|  ")

m_T = zeros(Float64, length(Temps))  # Array of magnetisation per site
err_m_T = zeros(Float64, length(Temps))

suzz_T = zeros(Float64, length(Temps))  # Array of specific heat
err_suzz_T = zeros(Float64, length(Temps))

Threads.@threads for i = 1:length(Temps)
    T = Temps[i]
    println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
    m_arr = zeros(Float64, nconfigs)

    local potts = PottsModel2D(L, q, start)
    for step=1:esteps
        wolff_cluster_update!(potts, T)
    end

    twice_τ = 2*autocorr_steps_wolff[i]
    numsteps = twice_τ*nconfigs
    for j in 1:numsteps
        wolff_cluster_update!(potts, T)
        if j%twice_τ == 0
            m_arr[j÷twice_τ] = magnetisation(potts) / potts.L^2
        end
    end

    m_T[i] = mean(m_arr) 
    err_m_T[i] = blocking_err(m_arr, A -> mean(A); blocks=err_nblocks)

    suzz_T[i] = succeptibility(m_arr, T, potts.L^potts.d)
    err_suzz_T[i] = blocking_err(m_arr, specific_heat, T, potts.L^potts.d; blocks=err_nblocks)

    println("| Process complete on thread #$(Threads.threadid()): T = $T")
end


#=
Plots
=#
println("| Generating Plots ...")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "magnetisation, u",
    title = "PottsModel2D $(L)x$(L): Magnetisation (per site) v/s temperature")

ax2 = Axis(f[2, 1], xlabel = "temperature, T", ylabel = "succeptibility, χ",
    title = "PottsModel2D $(L)x$(L): Succeptibility (per site) v/s temperature")

errorbars!(
    ax1, Temps, m_T, err_m_T,
    whiskerwidth = 10
)
scatter!(
    ax1, Temps, m_T,
    markersize = 7
)

errorbars!(
    ax2, Temps, suzz_T, err_suzz_T,
    whiskerwidth = 10)
scatter!(
    ax2, Temps, suzz_T,
    markersize = 7
)

save("plots/2Dmodel/mag_and_suzz_vs_Temp_$(L)_uncorr_wolff.png", f)

println("| Program Finished!")
println(".==================================")