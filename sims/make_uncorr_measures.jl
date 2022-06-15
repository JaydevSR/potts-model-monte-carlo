include("../src/pottsmc.jl")

Lvals = [16, 24, 32, 40, 48]
temps = [0.50, 0.70, 0.90, 0.94, 0.98, 1.02, 1.06, 1.10, 1.30, 1.50]
τvals = [   2,    3,    4,    5,   13,   23,   14,   14,   40,   61]
q=3
d=2
nconfigs=10000
eqsteps=10000
store_at="data/magdata/"

potts_getmagdata_to_txt(Lvals, temps, q, d, nconfigs, eqsteps, τvals; store_at=store_at, ntau=10)
