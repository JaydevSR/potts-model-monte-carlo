include("../src/pottsmc.jl")

#=
Perform simulation
=#

L = 16
q = 3
println(".==================================")
println("| Lattice Size: $(L) x $(L)        ")
println(".==================================")
println("|                                  ")

# Temps = [1.0]
Temps = [0.5, 0.7, 0.9, 0.94, 0.98, 1.02, 1.06, 1.1, 1.3, 1.5]
nsteps = 40000  # Number of steps for measurements
wsteps = 250
esteps = 2000
verbose = false

auto_corr_times_metro = zeros(Float64, length(Temps))
auto_corr_times_wolff = zeros(Float64, length(Temps))
auto_corr_times_wolff_scaled = zeros(Float64, length(Temps))

data_mat = zeros(Float64, (4, length(Temps)))
data_mat[1, :] = copy(Temps)

Threads.@threads for i = 1:length(Temps)
    T = Temps[i]
    println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
    m_arr = zeros(Float64, nsteps)

    # For Metropolis
    potts = PottsModel2D(L, q, :cold)
    for step=1:esteps
        metropolis_batch_update!(potts, T)
    end

    for step=1:nsteps
        m_arr[step] = magnetisation(potts) / potts.L^2
        delE = metropolis_batch_update!(potts, T)
    end

    corrfn = autocorrelation_fn(m_arr)
    auto_corr_times_metro[i] = sum(corrfn[1:wsteps])

    # For Wolff
    potts = PottsModel2D(L, q, :cold)
    for step=1:esteps
        wolff_cluster_update!(potts, T)
    end

    for step=1:nsteps
        m_arr[step] = magnetisation(potts) / potts.L^2
        wolff_cluster_update!(potts, T)
    end

    corrfn = autocorrelation_fn(m_arr)
    τ = sum(corrfn[1:wsteps])
    auto_corr_times_wolff[i] = τ
    # scale value for comparison with Metropolis algorithm
    if T >= 1.06
        τ = τ * (succeptibility(m_arr, T, potts.L^potts.d) * T / (potts.L^potts.d))
    end
    auto_corr_times_wolff_scaled[i] = τ
    println("| Process complete on thread #$(Threads.threadid()) (T = $T): τ = $(auto_corr_times_metro[i]), $(auto_corr_times_wolff_scaled[i]) ($(auto_corr_times_wolff[i]))")
end

println("| ")
println("| Using Metropolis: $(auto_corr_times_metro)")
println("| Using Wolff: $(auto_corr_times_wolff)")
println("| Using Wolff (Scaled): $(auto_corr_times_wolff_scaled)")
println("| ")
data_mat[2, :] = convert.(Int64, ceil.(auto_corr_times_metro))
data_mat[3, :] = convert.(Int64, ceil.(auto_corr_times_wolff))
data_mat[4, :] = convert.(Int64, ceil.(auto_corr_times_wolff_scaled))

open("D:\\Projects\\Potts-QCD\\potts-model-monte-carlo\\data\\corrtimes_data.txt", "w") do io
    writedlm(io, data_mat, ',')
end;

#=
Plots
=#
println("| Generating Plots ...")
println("| ")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "integrated autocorrelation time",
    title = "PottsModel2D $(L)x$(L):  Metropolis Algorithm")

ax2 = Axis(f[2, 1], xlabel = "temperature, T", ylabel = "integrated autocorrelation time",
    title = "PottsModel2D $(L)x$(L): Wolff Cluster Algorithm")

scatterlines!(
    ax1, Temps, auto_corr_times_metro,
    markersize = 7
)

scatterlines!(
    ax2, Temps, auto_corr_times_wolff_scaled,
    markersize = 7, label="Scaled for T > T_c"
)

scatterlines!(
    ax2, Temps, auto_corr_times_wolff, 
    markersize = 7, label="unscaled", linestyle=:dash
)

axislegend(ax2, position=:lt)

save("plots/2Dmodel/integrated_auto_corr_time_comparison_$(L)_v2.png", f)

println("| Program Finished!")
println(".==================================")