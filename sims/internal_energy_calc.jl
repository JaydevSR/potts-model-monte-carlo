include("../src/pottsmc.jl")

#=
Perform simulation
=#

L = 16
q = 3
println("================================\n")
println("    Lattice Size: $(L) x $(L)")
println("================================\n")

Temps = [i for i=0.1:0.2:2.1]
esteps = 1000  # Number of steps for equilibration
nsteps = 5000  # Number of steps for measurements

u_T = zeros(Float64, length(Temps))  # Array of mean internal energy per site
err_u_T = zeros(Float64, length(Temps))

c_T = zeros(Float64, length(Temps))  # Array of specific heat
err_c_T = zeros(Float64, length(Temps))

potts = PottsModel2D(L, q, :cold)

for i = 1:length(Temps)
    T = Temps[i]
    println("Calculating for T = $(T) ...")

    E_arr = zeros(Float64, nsteps)

    for step=1:esteps
        metropolis_batch_update!(potts, T)
    end

    E_arr[1] = hamiltonian(potts) / potts.L^2
    for step=2:nsteps
        # delE = metropolis_batch_update!(potts, T)
        # E_arr[step] = E_arr[step-1] + (delE / potts.L^2)

        wolff_cluster_update!(potts, T)
        E_arr[step] = hamiltonian(potts) / potts.L^2
    end

    u_T[i] = mean(E_arr) 
    err_u_T[i] = blocking_err(E_arr, A -> mean(A))

    c_T[i] = specific_heat(E_arr, T, potts.L^potts.d)
    err_c_T[i] = blocking_err(E_arr, specific_heat, T, potts.L^potts.d)
    println("   |          ")
    println("   +-> Done.\n")
end


#=
Plots
=#
println("Generating Plots ...")
f = Figure()

ax1 = Axis(f[1, 1], xlabel = "temperature, T", ylabel = "internal energy, u",
    title = "PottsModel2D $(L)x$(L): Internal energy (per site) v/s temperature")

ax2 = Axis(f[2, 1], xlabel = "temperature, T", ylabel = "specific heat, c",
    title = "PottsModel2D $(L)x$(L): Specific heat (per site) v/s temperature")

errorbars!(
    ax1, Temps, u_T, err_u_T,
    whiskerwidth = 10
)
scatter!(
    ax1, Temps, u_T,
    markersize = 7
)

errorbars!(
    ax2, Temps, c_T, err_c_T,
    whiskerwidth = 10)
scatter!(
    ax2, Temps, c_T,
    markersize = 7
)

save("plots/2Dmodel/u_and_c_vs_Temp_$(L)_wolff.pdf", f)

println("Program Finished!")
println("===========================\n")