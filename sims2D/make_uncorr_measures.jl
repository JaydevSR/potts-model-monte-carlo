include("../src/pottsmc.jl")

Lvals = [16, 24, 32, 40, 48]
# temps = [0.50, 0.70, 0.90, 0.94, 0.98, 1.02, 1.06, 1.10, 1.30, 1.50]
# τvals = [   2,    3,    4,    5,   13,   23,   14,   14,   40,   61]
# temps = [0.984, 0.990, 0.996, 1.000, 1.004, 1.01, 1.016]
# τvals = [   15,    17,    20,    30,    26,   25,    23]
# temps = [1.002, 1.006, 1.008, 1.024, 1.028, 1.04]
# τvals = [  30,     26,    25,    20,    18,   16]
temps = [1.012, 1.014, 1.018, 1.022, 1.026, 1.03]
τvals = [   25,    24,    23,    22,    19,   16]
q=3
d=2
nconfigs=5000
eqsteps=5000
store_at="data/magdata/"

potts_getmagdata_to_txt(Lvals, temps, q, d, nconfigs, eqsteps, τvals; store_at=store_at, ntau=3, mode="w")
