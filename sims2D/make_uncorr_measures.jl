include("../src/pottsmc.jl")

# lattice_sizes = [32, 48, 64, 80]
lattice_sizes = [56, 72]

temps = [0.500, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.920,
         0.940, 0.960, 0.980, 0.982, 0.984, 0.986, 0.988, 0.990, 0.992,
         0.994, 0.996, 0.998, 1.000, 1.002, 1.004, 1.006, 1.008, 1.010,
         1.012, 1.014, 1.016, 1.018, 1.020, 1.022, 1.024, 1.026, 1.028,
         1.030, 1.032, 1.034, 1.036, 1.038, 1.040, 1.060, 1.080, 1.100,
         1.150, 1.200, 1.250, 1.300, 1.350, 1.400, 1.500]

τvals = [02, 02, 03, 03, 04, 04, 04, 04, 05,
         05, 08, 13, 14, 15, 16, 17, 18, 19,
         20, 30, 28, 26, 25, 25, 25, 24, 23,
         22, 21, 21, 21, 20, 20, 20, 19, 19,
         18, 18, 16, 16, 15, 15, 14, 14, 14,
         20, 30, 35, 40, 45, 50, 60]

q, d = 3, 2
nconfigs=20000
eqsteps=10000
store_at="data/"
ntau=3
start=:cold
mag_definition=:both
get_configs=false
fix_vacuum=true
verbose=true

potts_get_measurements_to_txt(
    lattice_sizes, temps, q, d, nconfigs, eqsteps, τvals;
    ntau=ntau, start=start,
    mag_definition=mag_definition,
    get_configs=get_configs,
    fix_vacuum=fix_vacuum,
    store_at=store_at,
    verbose=verbose,
    mode="w"
    )