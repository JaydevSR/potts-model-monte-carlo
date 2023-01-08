import lmfit
import numpy as np
import matplotlib.pyplot as plt

"""
Size 32 :
 * peak location = 1.009347198871323 ± 0.00017271331290456279


Size 48 :
 * peak location = 1.004098703447542 ± 4.558877215489571e-5


Size 56 :
 * peak location = 1.00293349803508 ± 3.540467136196057e-5


Size 64 :
 * peak location = 1.0020574408340042 ± 2.9553098374471333e-5


Size 72 :
 * peak location = 1.0013518959209922 ± 2.548910629662237e-5


Size 80 :
 * peak location = 1.0008723162246127 ± 2.316395387744226e-5
"""

lattice_sizes = np.array([48, 56, 64, 72, 80])
tc_vals = np.array([
    1.004098703447542,
    1.00293349803508,
    1.0020574408340042,
    1.0013518959209922,
    1.0008723162246127
])

tc_err = np.array([
    4.558877215489571e-5,
    3.540467136196057e-5,
    2.9553098374471333e-5,
    2.548910629662237e-5,
    2.316395387744226e-5
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
with open('data/Tc_scaling_fit_gaussian_results.txt', 'w') as f:
    f.write(lmfit.fit_report(out))

print(lmfit.fit_report(out))

a, b, c = out.params['a'].value, out.params['b'].value, out.params['c'].value

# plot the results
plt.scatter(lattice_sizes, tc_vals, label='data')
lrange = np.linspace(lattice_sizes[0], lattice_sizes[-1], 100)
plt.plot(lrange, residuals(out.params, lrange), label=f'fit')
plt.xlabel('Lattice size')
plt.ylabel('Tc')
plt.title(f'Finite size scaling of Tc: T_c(L) = {a:.4f} + {b:.4f} * L^(-1 / {c:.4f})')
plt.legend()
plt.savefig('plots/2Dmodel/finite_size_scaling/Tc_scaling_fit_gaussian.svg')
