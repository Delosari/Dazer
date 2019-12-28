import os
os.environ["MKL_THREADING_LAYER"] = "GNU"
import theano.tensor as tt
from theano import function, shared
import pymc_examples
import numpy as np


# def v1():
#     x = tt.vector()
#     min_x = tt.vector()
#     max_x = tt.vector()
#     y = tt.clip(x, min_x, max_x)
#     f = function([x, min_x, max_x], outputs=y)
#     print f([2, 1, 4], [0, 2, 3], [1, 3, 5])
#
#
# def v2():
#     x = tt.vector()
#     min_max = tt.matrix()
#     y = tt.clip(x, min_max[:, 0], min_max[:, 1])
#     f = function([x, min_max], outputs=y)
#     print f([2, 1, 4], [[0, 1], [2, 3], [3, 5]])
#
# def v3():
#     x = tt.scalar()
#     min_x = tt.scalar()
#     max_x = tt.scalar()
#     y = tt.clip(x, min_x, max_x)
#     f = function([x, min_x, max_x], outputs=y)
#     print f(3, 0, 2)
#
# def main():
#     v1()
#     v2()
#     v3()
#
#
# main()


# Model emissivity
emisCoeffs = np.array([4.42e-05, 7.91e-4, 0.828, 2.79e-05])
def modelEmis(temp, den, a, b, c, d):
    return (a + b * den) * np.log10(temp / 10000.0) - np.log10(c + d * den)

# Temperature and density grid
temp_range = np.arange(9000, 25000, 1)
den_range = np.arange(1, 1000, 1)
temp_mesh, den_mesh = np.meshgrid(temp_range, den_range)

# True values
temp_true, den_true = 13555, 245

# Observational data
obs_data = modelEmis(temp_true, den_true, *emisCoeffs)
obs_err = obs_data * 0.05

# tempGridFlatten, denGridFlatten = temp_mesh.flatten(), den_mesh.flatten()
# print idx_temp, idx_den
# print temp_true, den_true
# print temp_range[idx_temp-1],temp_range[idx_temp],temp_range[idx_temp+1]
# print den_range[idx_den-1],den_range[idx_den],den_range[idx_den+1]

# --------------------Treatment via model parametrisation fitting--------------------

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
#     trace = pymc3.sample(2000, tune=500, nchains=2, njobs=1, model=model)
#
# print pymc3.summary(trace).round(4)

# --------------------Treatment via grid interpolation--------------------

def indexing_Qxx(x_range, x_coord, y_range, y_coord):
    x_i = np.argmin(np.abs(x_range - x_coord))
    y_i = np.argmin(np.abs(y_range - y_coord))
    return x_i, y_i

def indexing_Qxx_tt(x_range, x_coord, y_range, y_coord):
    x_i = tt.argmin(tt.abs_(x_range - x_coord))
    y_i = tt.argmin(tt.abs_(y_range - y_coord))
    return x_i, y_i

def biLinearInterpolation(x, y, matrix, x0_limit = 0, y0_limit = 0):

    x_in, y_in = x - x0_limit, y - y0_limit

    x0 = np.floor(x_in).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y_in).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, matrix.shape[1]-1)
    x1 = np.clip(x1, 0, matrix.shape[1]-1)
    y0 = np.clip(y0, 0, matrix.shape[0]-1)
    y1 = np.clip(y1, 0, matrix.shape[0]-1)

    Ia = matrix[y0, x0]
    Ib = matrix[y1, x0]
    Ic = matrix[y0, x1]
    Id = matrix[y1, x1]

    wa = (x1 - x_in) * (y1 - y_in)
    wb = (x_in - x0) * (y1 - y_in)
    wc = (x1 - x_in) * (y_in - y0)
    wd = (x_in - x0) * (y_in - y0)
    w0 = ((x1 - x0) * (y1 - y0))

    return wa*Ia + wb*Ib + wc*Ic + wd*Id / w0

def biLinearInterpolation_v2(x, y, matrix, x0_limit = 0, y0_limit = 0):

    x_in, y_in = x - x0_limit, y - y0_limit

    x0 = np.floor(x_in).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y_in).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, matrix.shape[1]-1)
    x1 = np.clip(x1, 0, matrix.shape[1]-1)
    y0 = np.clip(y0, 0, matrix.shape[0]-1)
    y1 = np.clip(y1, 0, matrix.shape[0]-1)

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

    return a0 + a1*x_in + a2*y_in + a3*x_in*y_in

def biLinearInterpolation_v2_tt(x, y, matrix, x0_limit, y0_limit, x1_limit, y1_limit):

    x_in, y_in = x - x0_limit, y - y0_limit

    x0 = tt.floor(x_in)
    x1 = x0 + 1
    y0 = tt.floor(y_in)
    y1 = y0 + 1

    x0 = tt.cast(tt.clip(x0, 1, x1_limit), 'int64')
    x1 = tt.cast(tt.clip(x1, 1, x1_limit), 'int64')
    y0 = tt.cast(tt.clip(y0, 1, y1_limit), 'int64')
    y1 = tt.cast(tt.clip(y1, 1, y1_limit), 'int64')

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

    return a0 + a1*x_in + a2*y_in + a3*x_in*y_in

def biLinearInterpolation_Coeffs_tt(x, y, matrix, x0_limit, y0_limit, x1_limit, y1_limit):

    x_in, y_in = x - x0_limit, y - y0_limit

    x0 = tt.floor(x_in)
    x1 = x0 + 1
    y0 = tt.floor(y_in)
    y1 = y0 + 1

    x0 = tt.cast(tt.clip(x0, 0, x1_limit), 'int32')
    x1 = tt.cast(tt.clip(x1, 0, x1_limit), 'int32')
    y0 = tt.cast(tt.clip(y0, 0, y1_limit), 'int32')
    y1 = tt.cast(tt.clip(y1, 0, y1_limit), 'int32')

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


# Theano operations
temp_vector, den_vector = tt.dvectors('temp_vector', 'den_vector')
temp_scalar, den_scalar = tt.dscalars('temp_scalar', 'den_scalar')
indexing_Qxx_ttfunction = function(inputs=[temp_vector, temp_scalar, den_vector, den_scalar],
                               outputs=indexing_Qxx_tt(temp_vector, temp_scalar, den_vector, den_scalar))

x, y = tt.dscalars('x', 'y')
x0_limit, y0_limit, x1_limit, y1_limit = tt.dscalars('x0_limit', 'y0_limit', 'x1_limit', 'y1_limit')
matrixGrid = tt.dmatrix('matrixGrid')

bilinearInterpolationttfunction = function(inputs=[x, y, matrixGrid, x0_limit, y0_limit, x1_limit, y1_limit],
                               outputs=biLinearInterpolation_v2_tt(x, y, matrixGrid, x0_limit, y0_limit, x1_limit, y1_limit))

biLinearInterpolation_CoeffsFunction_tt = function(inputs=[x, y, matrixGrid, x0_limit, y0_limit, x1_limit, y1_limit],
                               outputs=biLinearInterpolation_Coeffs_tt(x, y, matrixGrid, x0_limit, y0_limit, x1_limit, y1_limit))

# Numpy steps
idx_temp, idx_den = indexing_Qxx(temp_range, temp_true, den_range, den_true)
emisGrid = modelEmis(temp_mesh, den_mesh, *emisCoeffs)

# Theano steps
idx_temp_tt, idx_den_tt = indexing_Qxx_ttfunction(temp_range, temp_true, den_range, den_true)

# print 'Theano indeces', idx_temp, idx_den, temp_range[idx_temp], den_range[idx_den]
# print 'Numpy indeces', idx_temp_tt, idx_den_tt, temp_range[idx_temp], den_range[idx_den], '\n'
# print
# print 'Grid values', emisGrid[idx_den, idx_temp] == modelEmis(temp_range[idx_temp], den_range[idx_den], *emisCoeffs), '\n'
# print 'True interpolation',  modelEmis(temp_true, den_true, *emisCoeffs)
# print 'Interpolation', biLinearInterpolation(temp_true, den_true, emisGrid, temp_range[0], den_range[0])
# print 'Interpolation 2', biLinearInterpolation_v2(temp_true, den_true, emisGrid, temp_range[0], den_range[0])
# print 'Interpolation tt', bilinearInterpolationttfunction(temp_true, den_true, emisGrid, temp_range[0], den_range[0], temp_range[-1], den_range[-1])
# print 'Interpolation coeffs', biLinearInterpolation_CoeffsFunction_tt(temp_true, den_true, emisGrid, temp_range[0], den_range[0], temp_range[-1], den_range[-1])
#
emisGrid_tt = shared(emisGrid)

with pymc_examples.Model() as model2:

    # Physical condition priors
    temp = pymc_examples.Normal('temp', mu=13000.0, sd=500.0)
    den = pymc_examples.Normal('den', mu=300.0, sd=50.0)

    # Compute model densities
    a0, a1, a2, a3 = biLinearInterpolation_Coeffs_tt(temp, den, emisGrid_tt, temp_range[0], den_range[0], temp_range[-1], den_range[-1])
    emis_synth = a0 + a1*(temp - temp_range[0]) + a2*(den - den_range[0]) + a3*(temp - temp_range[0])*(den - den_range[0])
    pymc_examples.Deterministic('emis_synth', emis_synth)

    # Likelihood
    Y = pymc_examples.Normal('Y_emision', mu=emis_synth, sd=obs_err, observed=obs_data)

    # Launch model
    trace = pymc_examples.sample(2000, tune=500, nchains=2, njobs=1, model=model2)

print pymc_examples.summary(trace).round(4)