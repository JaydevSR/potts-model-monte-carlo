include("../src/pottsmc.jl")

Lvals = [16, 32, 48, 64]

# temps = [0.500, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.920,
#          0.940, 0.960, 0.980, 0.984, 0.988, 0.992, 0.996, 1.000, 1.004,
#          1.008, 1.012, 1.016, 1.020, 1.040, 1.060, 1.080, 1.100, 1.150,
#          1.200, 1.250, 1.300, 1.350, 1.400, 1.500]

temps = [1.204, 1.208, 1.320, 1.360]
# τvals = [02, 02, 03, 03, 04, 04, 04, 04, 05,
#          05, 08, 13, 15, 17, 19, 20, 30, 26,
#          25, 25, 23, 21, 18, 16, 15, 14, 14,
#          30, 35, 40, 45, 50, 60]
τvals = [21, 20, 20, 19]

q, d = 3, 2
nconfigs=5000
eqsteps=10000
store_at="data/"
ntau=4
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