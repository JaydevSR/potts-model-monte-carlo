include("../src/pottsmc.jl")

#=
Perform simulation
=#

L = 16
q = 3
start = :cold
esteps = 10_000  # Number of steps for equilibration
nconfigs = 50_000  # Number of steps for measurements


data_path = "data/corrtimes_data.txt"
data = readdlm(data_path, ',', Float64)
Temps = data[1, :]
autocorr_steps_wolff = convert.(Int64, data[3, :])

function test_routine(T, L, q, start, τ, nconfigs, esteps)
    local potts = PottsModel2D(L, q, start)

    for step=1:esteps
        wolff_cluster_update!(potts, T)
    end

    twice_τ = 2*τ
    numsteps = twice_τ*nconfigs
    for j in 1:numsteps
        wolff_cluster_update!(potts, T)
    end
end

println(".==================================")
println("| Lattice Size: $(L) x $(L)        ")
println(".==================================")
println("|  ")

elapsed_times = zeros(Float64, length(Temps))

Threads.@threads for i = 1:length(Temps)
    T = Temps[i]
    τ = autocorr_steps_wolff[i]
    println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
    elapsed_times[i] = @elapsed test_routine(T, L, q, start, τ, nconfigs, esteps)
    println("| Process complete on thread #$(Threads.threadid()): T = $T (elapsed $(elapsed_times[i]) seconds)")
end

if Sys.iswindows()
    store_at="data/benchmark_times_laptop.txt"
else
    store_at="data/benchmark_times_pc_linux.txt"
end

open(store_at, "w") do io
    writedlm(io, [Temps elapsed_times], ',')
end

println("| ")
println("| Program Finished!")
println(".==================================")