import lmfit
import numpy as np
from matplotlib import pyplot as plt

"""
Size 32 :
 * peak location = 1.0088133938977524 ± 0.315729537030768


Size 48 :
 * peak location = 1.0037126712055784 ± 0.09902317928004663


Size 56 :
 * peak location = 1.0021600385069998 ± 0.18151555991193455


Size 64 :
 * peak location = 1.0009103256327856 ± 0.15560977544242188


Size 72 :
 * peak location = 1.0003331493025998 ± 0.06417313931737789


Size 80 :
 * peak location = 0.9997345330528584 ± 0.11894946740275754


Size 96 :
 * peak location = 0.9988572671581376 ± 0.17022948837868232


Size 128 :
 * peak location = 0.997941002359664 ± 0.17442585458631407
"""

lattice_sizes = np.array([48, 56, 64, 72, 80, 96, 128])
tc_vals = np.array([
    1.0037126712055784,
    1.0021600385069998,
    1.0009103256327856,
    1.0003331493025998,
    0.9997345330528584,
    0.9988572671581376,
    0.997941002359664
])

tc_err = np.array([
    0.09902317928004663,
    0.18151555991193455,
    0.15560977544242188,
    0.06417313931737789,
    0.11894946740275754,
    0.17022948837868232,
    0.17442585458631407
])

def residuals(params, xvals, yvals=None, eps=None):
    """Residuals function for finite size scaling of Tc."""
    tinf = params['a']
    scaling = params['b']
    expnu = params['c']

    model = tinf + scaling * (xvals ** (- 1 / expnu))
    if yvals is None:
        return model
    if eps is None:
        return model - yvals
    return (model - yvals) / eps

pars = lmfit.Parameters()
pars.add_many(('a', 1.0), ('b', 1.0), ('c', 1.0))

out = lmfit.minimize(residuals, pars, args=(lattice_sizes, tc_vals, tc_err), method='leastsq')

# print the results to a file
with open('data/Tc_scaling_fit_parabola_results_pycall.txt', 'w') as f:
    f.write(lmfit.fit_report(out))

print(lmfit.fit_report(out))
a, b, c = out.params['a'].value, out.params['b'].value, out.params['c'].value

# plot the results
plt.scatter(lattice_sizes, tc_vals, label='data')
lrange = np.linspace(lattice_sizes[0], lattice_sizes[-1], 100)
plt.plot(lrange, residuals(out.params, lrange), label='fit')
plt.xlabel('Lattice size')
plt.ylabel('Tc')
plt.title(f'Finite size scaling of Tc: T_c(L) = {a:.4f} + {b:.4f} * L^(-1 / {c:.4f})')
plt.legend()
plt.savefig('plots/2Dmodel/finite_size_scaling/Tc_scaling_fit_parabola_pycall.svg')
