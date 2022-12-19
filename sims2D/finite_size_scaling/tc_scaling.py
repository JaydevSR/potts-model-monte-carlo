import lmfit
import numpy as np

lattice_sizes = np.array([48, 56, 64, 72, 80, 96])
tc_vals = np.array([
    1.004098558580445,
    1.002880459951727,
    1.0019690678472521,
    1.0011767926962136,
    1.0005831960780498,
    0.9995852248655276
])

tc_err = np.array([
    0.00012838372820415597,
    0.00016614142151693173,
    0.00020775901663022396,
    0.00024041071306175725,
    0.0002808961827846406,
    0.0003212285341357723
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
with open('data/Tc_scaling_fit_results.txt', 'w') as f:
    f.write(lmfit.fit_report(out))

print(lmfit.fit_report(out))
