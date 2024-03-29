import lmfit
import numpy as np
from matplotlib import pyplot as plt

"""
Size 32 :
 * peak location = 1.0087550077965644 ± 0.3897079104394817


Size 48 :
 * peak location = 1.0037076586533769 ± 0.10648230645344149


Size 56 :
 * peak location = 1.002158707554152 ± 0.10193685496743785


Size 64 :
 * peak location = 1.0008862078388885 ± 0.07007491365312664


Size 72 :
 * peak location = 1.000333823075783 ± 0.05558539643410161


Size 80 :
 * peak location = 0.9997234675238178 ± 0.05026909135561001


Size 96 :
 * peak location = 0.9988293864556604 ± 0.03661099661084414


Size 128 :
 * peak location = 0.9980082580114912 ± 0.02188981537738993
"""

lattice_sizes = np.array([48, 56, 64, 72, 80, 96, 128])
tc_vals = np.array([
    1.0037076586533769,
    1.002158707554152,
    1.0008862078388885,
    1.000333823075783,
    0.9997234675238178,
    0.9988293864556604,
    0.9980082580114912
])

tc_err = np.array([
    0.10648230645344149,
    0.10193685496743785,
    0.07007491365312664,
    0.05558539643410161,
    0.05026909135561001,
    0.03661099661084414,
    0.02188981537738993
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
with open('data/Tc_scaling_fit_parabola_results.txt', 'w') as f:
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
plt.savefig('plots/2Dmodel/finite_size_scaling/Tc_scaling_fit_parabola.svg')
