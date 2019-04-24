import os
os.environ["MKL_THREADING_LAYER"] = "GNU"
import theano.tensor as tt
from theano import function, shared
import pymc3
import numpy as np
import matplotlib.pyplot as plt


# Model emissivity
emisCoeffs = np.array([4.42e-05, 7.91e-4, 0.828, 2.79e-05])
def modelEmis(temp, den, a, b, c, d):
    return (a + b * den) * np.log10(temp / 10000.0) - np.log10(c + d * den)

# True values
temp_true, den_true = 13555.5, 245.5

# Observational data
obs_data = modelEmis(temp_true, den_true, *emisCoeffs)
obs_err = obs_data * 0.02

# # --------------------Treatment via model parametrisation fitting--------------------
# with pymc3.Model() as model:
#
#     # Physical condition priors
#     temp = pymc3.Normal('temp', mu=13000.0, sd=500.0)
#     den = pymc3.Normal('den', mu=300.0, sd=50.0)
#
#     # Compute model densities
#     emis_synth = modelEmis(temp, den, *emisCoeffs)
#     pymc3.Deterministic('emis_synth', emis_synth)
#
#     # Likelihood
#     Y = pymc3.Normal('Y_emision', mu=emis_synth, sd=obs_err, observed=obs_data)
#
#     # Launch model
#     trace = pymc3.sample(2000, tune=800, nchains=2, njobs=1, model=model)
#
# print pymc3.summary(trace).round(4)
#
# pymc3.traceplot(trace)
#
# plt.show()

# --------------------Treatment via grid interpolation--------------------

def biLinearInterpolation_Coeffs(x, y, matrix, x0_limit, y0_limit, x1_limit, y1_limit):

    x_in, y_in = x - x0_limit, y - y0_limit

    x0 = tt.floor(x_in)
    x1 = x0 + 1
    y0 = tt.floor(y_in)
    y1 = y0 + 1

    x0 = tt.cast(tt.clip(x0, 0, x1_limit), 'int64')
    x1 = tt.cast(tt.clip(x1, 0, x1_limit), 'int64')
    y0 = tt.cast(tt.clip(y0, 0, y1_limit), 'int64')
    y1 = tt.cast(tt.clip(y1, 0, y1_limit), 'int64')

    Q11 = matrix[y0, x0]
    Q21 = matrix[y1, x0]
    Q12 = matrix[y0, x1]
    Q22 = matrix[y1, x1]

    den12 = (x0-x1) / (y0-y1)
    den21 = (x0-x1) / (y1-y0)

    a0 = Q11 * x1*y1 / den12 + Q12 * x1*y0 / den21 + Q21 * x0*y1 / den21 + Q22 * x0*y0 / den12
    a1 = Q11 * y1 / den21 + Q12 * y0 / den12 + Q21 * y1 / den12 + Q22 * y0 / den21
    a2 = Q11 * x1 / den21 + Q12 * x1 / den12 + Q21 * x0 / den12 + Q22 * x0 / den21
    a3 = Q11 / den12 + Q12 / den21 + Q21 / den21 + Q22 / den12

    return a0, a1, a2, a3


# Temperature and density grid
temp_range = np.arange(9000, 25000, 1)
den_range = np.arange(1, 1000, 1)
temp_mesh, den_mesh = np.meshgrid(temp_range, den_range)

# Emissivity grid
emisGrid = modelEmis(temp_mesh, den_mesh, *emisCoeffs)
emisGrid_tt = shared(emisGrid)

with pymc3.Model() as model2:

    # Physical condition priors
    temp = pymc3.Normal('temp', mu=13000.0, sd=500.0)
    den = pymc3.Normal('den', mu=300.0, sd=50.0)

    # Compute model densities
    a0, a1, a2, a3 = biLinearInterpolation_Coeffs(temp, den, emisGrid_tt, temp_range[0], den_range[0], temp_range[-1], den_range[-1])
    emis_synth = a0 + a1*(temp - temp_range[0]) + a2*(den - den_range[0]) + a3*(temp - temp_range[0])*(den - den_range[0])
    pymc3.Deterministic('emis_synth', emis_synth)

    # Likelihood
    Y = pymc3.Normal('Y_emision', mu=emis_synth, sd=obs_err, observed=obs_data)

    # Launch model
    trace = pymc3.sample(2000, tune=800, nchains=2, njobs=1, model=model2)

print pymc3.summary(trace).round(4)

pymc3.traceplot(trace)

plt.show()