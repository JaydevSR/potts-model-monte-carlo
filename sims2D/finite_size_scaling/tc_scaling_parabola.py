import lmfit
import numpy as np
from matplotlib import pyplot as plt

lattice_sizes = np.array([32, 48, 56, 64, 72, 80, 96, 128])
tc_vals = np.array([
    1.008788264284399,
    1.0037077804743122,
    1.0021463679815212,
    1.000886609848195,
    1.000334873385511,
    0.9997025205774708,
    0.9988191144571779,
    0.9980114447110828
])

tc_err = np.array([
    0.3230155328751795,
    0.11565719794367375,
    0.09523867104088478,
    0.06996979656817764,
    0.06248921955640057,
    0.05469383449564052,
    0.034335474900586346,
    0.021541391255513305
])

def residuals(params, xvals, yvals=None, eps=None):
    """Residuals function for finite size scaling of Tc."""
    tinf = params['a']
    scaling = params['b']
    expnu = params['c']

    model = tinf - scaling * (xvals ** (- 1 / expnu))
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

# plot the results
plt.scatter(lattice_sizes, tc_vals, label='data')
lrange = np.linspace(lattice_sizes[0], lattice_sizes[-1], 100)
plt.plot(lrange, residuals(out.params, lrange), label='fit')
plt.xlabel('Lattice size')
plt.ylabel('Tc')
plt.legend()
plt.savefig('plots/2Dmodel/finite_size_scaling/Tc_scaling_fit_parabola.svg')
