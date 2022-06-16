include("../src/pottsmc.jl")

#=
Perform simulation
=#

L = 16
q = 3
println("================================\n")
println("    Lattice Size: $(L) x $(L)")
println("================================\n")

T = 1.0
# thsteps = [1000, 1000, 1000, 1000]
# nsteps = [10000, 20000, 30000, 40000]
thsteps = [10000]
nsteps = [100000]
cols = [:blue, :red, :green, :black]

f=Figure(resolution=(800, 600))
ax=Axis(f[1,1])

for i=1:length(thsteps)
    col=cols[i]
    st=nsteps[i]
    th=thsteps[i]

    println("   | Th=$th, St=$st ...")
    potts = PottsModel2D(L, q, :cold)
    m_arr = zeros(Float64, st)
    for step=1:th
        wolff_cluster_update!(potts, T)
    end

    for step=1:st
        if step%5000==0
            println("   |   | $(step / st * 100) %")
        end
        m_arr[step] = magnetisation(potts) / potts.L^2
        wolff_cluster_update!(potts, T)
    end

    corrfn = autocorrelation_fn(m_arr)

    lines!(corrfn[1:1000] .+ (i-1), color=col, label="thermalization steps=$th, nsteps=$st")
end

axislegend(ax)
display(f)