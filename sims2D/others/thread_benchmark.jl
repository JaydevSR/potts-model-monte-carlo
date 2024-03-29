include("../../src/utils.jl")
include("../../src/models.jl")
include("../../src/mcalgos.jl")
include("../../src/observables.jl")
include("../../src/statutils.jl")
include("../../src/measureutils.jl")

using DelimitedFiles

lattice_sizes = [32]
temps = [0.500, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.920,
         0.940, 0.960, 0.980, 0.984, 0.988, 0.992, 0.996, 1.000, 1.004,
         1.008, 1.012, 1.016, 1.020, 1.024, 1.028, 1.032, 1.036, 1.040,
         1.060, 1.080, 1.100, 1.150, 1.200, 1.250, 1.300, 1.350, 1.400, 1.500]

τvals = [02, 02, 03, 03, 04, 04, 04, 04, 05,
         05, 08, 13, 15, 17, 19, 20, 30, 26,
         25, 25, 23, 21, 21, 20, 20, 19, 18,
         16, 15, 14, 14, 30, 35, 40, 45, 50, 60]

q, d = 3, 2
nconfigs=4000
eqsteps=10000
store_at="dump/"
ntau=2
start=:cold
mag_definition=:both
get_configs=false
fix_vacuum=true
verbose=false

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
