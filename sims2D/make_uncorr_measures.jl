include("../src/pottsmc.jl")

Lvals = [16]#, 24, 32, 40]
temps = [0.5, 0.7, 0.9, 0.94, 0.98, 0.984, 0.988, 0.992, 0.996, 1.0, 1.004, 1.008, 1.012, 1.016, 1.02, 1.06,  1.1, 1.3, 1.5]
τvals = [  2,   3,   4,    5,   13,    15,    17,    19,   20,   30,    26,    25,    25,    23,   21,   14,   14,  40,  61]
q=3
d=2
nconfigs=2000
eqsteps=5000
store_at="data/uncorr_configs/"
fix_vacuum=true

# potts_getmagdata_to_txt(Lvals, temps, q, d, nconfigs, eqsteps, τvals; store_at=store_at, ntau=4, mode="w", fix_vacuum=fix_vacuum)
potts_getconfigdata_to_txt(Lvals, temps, q, d, nconfigs, eqsteps, τvals; store_at=store_at, ntau=2, mode="w", fix_vacuum=fix_vacuum)