include("../src/pottsmc.jl")

Lvals = [16, 32, 48, 64]
temps = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9]
append!(temps, [0.920, 0.940, 0.960, 0.980, 0.984, 0.988, 0.992, 0.996, 1.000])
append!(temps, [1.004, 1.008, 1.012, 1.016, 1.020, 1.040, 1.060, 1.080, 1.100])
append!(temps, [1.2, 1.25, 1.3, 1.4, 1.5])
τvals = [2, 2, 3, 4, 4, 4]
append!(τvals, [ 5,  5,  8, 13, 15, 17, 19, 20, 30])
append!(τvals, [26, 25, 25, 23, 21, 18, 16, 15, 14])
append!(τvals, [30, 35, 40, 50, 61])
q=3
d=2
nconfigs=5000
eqsteps=10000
store_at="data/"
ntau=5
start=:cold
mag_definition=:both
get_configs=true
fix_vacuum=true
verbose=true


potts_get_measurements_to_txt(
    Lvals, temps, q, d, nconfigs, eqsteps, τvals;
    ntau=ntau, start=start,
    mag_definition=mag_definition,
    get_configs=get_configs,
    fix_vacuum=fix_vacuum,
    store_at=store_at,
    verbose=verbose,
    mode="w"
    )

# potts_getmagdata_to_txt(Lvals, temps, q, d, nconfigs, eqsteps, τvals; store_at=store_at, ntau=4, mode="w", fix_vacuum=fix_vacuum)
# potts_getconfigdata_to_txt(Lvals, temps, q, d, nconfigs, eqsteps, τvals; store_at=store_at, ntau=2, mode="w", fix_vacuum=fix_vacuum)