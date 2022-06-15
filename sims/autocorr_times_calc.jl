include("../src/pottsmc.jl")

#=
Perform simulation
=#

L = 48
q = 3
println(".==================================")
println("| Lattice Size: $(L) x $(L)        ")
println(".==================================")
println("|                                  ")

# Temps = [1.0]
Temps = [0.5, 0.7, 0.9, 0.94, 0.98, 1.02, 1.06, 1.1, 1.3, 1.5]
nsteps = 50000  # Number of steps for measurements
wsteps = 250
esteps = 10000
verbose = false

auto_corr_times_wolff = zeros(Float64, length(Temps))

data_mat = zeros(Float64, (2, length(Temps)))
data_mat[1, :] = copy(Temps)

Threads.@threads for i = 1:length(Temps)
    T = Temps[i]
    println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
    m_arr = zeros(Float64, nsteps)
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
    println("| Process complete on thread #$(Threads.threadid()) (T = $T): τ = $(auto_corr_times_wolff[i])")
end

println("| ")
println("| Using Wolff, τ=$(auto_corr_times_wolff)")
println("| ")
data_mat[2, :] = convert.(Int64, ceil.(auto_corr_times_wolff))

open("data/corrtimes_data_L$(L).txt", "w") do io
    writedlm(io, data_mat, ',')
end;