from numpy import sqrt, pi, exp, linspace, random, empty, array, copy
from scipy.optimize import curve_fit
from lmfit import Model
from timeit import default_timer as timer
from kapteyn import kmpfit

def gaussian(x, amp, cen, wid):
    return amp * exp(-(x-cen)**2 /wid)

def residuals_gauss(p, c):        # Function needed by fit routine
    x, y = c            # The values for x, y and weights
    amp, cen, wid = p
    return (y - gaussian(x, amp, cen, wid))   # An array with (weighted) residuals)


init_vals = [1, 0, 1]
x = linspace(-10,10)
n_points = len(x) 
gaussian_true = gaussian(x, 2.33, 0.21, 1.51)
gmod = Model(gaussian)
y_obs = gaussian_true + random.normal(0, 0.2, n_points)

print '--Single fit'
start = timer()
best_vals, covar = curve_fit(gaussian, x, y_obs, p0=init_vals)
end = timer()
print 'curve fit', best_vals, ' time ', (end - start) 

start = timer()
result = gmod.fit(y_obs, x=x, amp=1, cen=0, wid=1)
end = timer()
print 'lmfit', array(result.params.valuesdict().values()), (end - start) 

start = timer()
fitobj = kmpfit.Fitter(residuals=residuals_gauss, data= (x, y_obs))
fitobj.fit(params0  = [1, 0, 1])
p_1_i               = fitobj.params
end = timer()
print 'kmpfit', p_1_i, type(p_1_i), (end - start) 

print '--Bootstrap'
results_matrix_curvefit = empty([3,1000])
results_matrix_lmfit = empty([3,1000])
results_matrix_kmpfit = empty([3,1000])

start = timer()
for i in range(1000):    
    best_vals, covar = curve_fit(gaussian, x, gaussian_true + random.normal(0, 0.2, n_points), p0=init_vals)
    results_matrix_curvefit[:,i] = best_vals
end = timer()
print 'curve fit', results_matrix_curvefit.mean(1), ' time ', (end - start) 

start = timer()
x_i = copy(x)
y_synth = gaussian_true + random.normal(0, 0.2, n_points)
fitobj = gmod.fit(gaussian_true + random.normal(0, 0.2, n_points), x=x, amp=1, cen=0, wid=1)
 
for i in range(1000):
    x_i[:] = x
    y_synth[:] = gaussian_true + random.normal(0, 0.2, n_points)
    print fitobj.params.valuesdict().values()
    results_matrix_lmfit[:,i] = array(fitobj.params.valuesdict().values())
end = timer()
print 'lmfit', results_matrix_lmfit.mean(1), (end - start)

start = timer()
for i in range(1000):
    result = gmod.fit(gaussian_true + random.normal(0, 0.2, n_points), x=x, amp=1, cen=0, wid=1)
    results_matrix_lmfit[:,i] = array(result.params.valuesdict().values())
end = timer()
print 'lmfit', results_matrix_lmfit.mean(1), (end - start)

start = timer()
x_i = copy(x)
y_synth = gaussian_true + random.normal(0, 0.2, n_points)
fitobj = kmpfit.Fitter(residuals=residuals_gauss, data= (x_i, y_synth))
print 'Primero', fitobj.params
for i in range(1000):
    x_i[:] = x
    y_synth[:] = gaussian_true + random.normal(0, 0.2, n_points)
    fitobj.fit(params0  = [1, 0, 1])
    results_matrix_kmpfit[:,i] = fitobj.params
    print fitobj.params
end = timer()
print 'kmpfit fit', results_matrix_kmpfit.mean(1), ' time ', (end - start)





# #!/usr/bin/env python
# #------------------------------------------------------------
# # Purpose: Demonstrate that the scaled covariance errors for
# # unweighted fits are comparable to errors we find with
# # a bootstrap method.
# # Vog, 24 Nov 2011
# #------------------------------------------------------------
# 
# import numpy
# from matplotlib.pyplot import figure, show, rc
# from numpy.random import normal, randint
# from kapteyn import kmpfit
# 
# # Residual and model in 1 function. Model is straight line
# def residuals(p, data):
#     x, y, err = data
#     a, b = p
#     model = a + b*x
#     return (y-model)/err
# 
# # Artificial data
# N = 100
# a0 = 0; b0 = 1.2
# x = numpy.linspace(0.0, 2.0, N)
# y = a0 + b0*x + normal(0.0, 0.4, N)  # Mean,sigma,N
# err = numpy.ones(N)                  # All weights equal to 1
# 
# # Prepare fit routine
# fitobj = kmpfit.Fitter(residuals=residuals, data=(x, y, err))
# try:
#     fitobj.fit(params0=[1,1])
# except Exception, mes:
#     print "Something wrong with fit: ", mes
#     raise SystemExit
# 
# print "\n\n======== Results kmpfit unweighted fit ========="
# print "Params:        ", fitobj.params
# print "Errors from covariance matrix         : ", fitobj.xerror
# print "Uncertainties assuming reduced Chi^2=1: ", fitobj.stderr
# print "Chi^2 min:     ", fitobj.chi2_min
# print "Reduced Chi^2: ", fitobj.rchi2_min
# print "Iterations:    ", fitobj.niter
# print "Function ev:   ", fitobj.nfev
# print "Status:        ", fitobj.status
# print "Status Message:", fitobj.message
# 
# # Bootstrap method to find uncertainties
# A0, B0 = fitobj.params
# xr = x.copy()
# yr = y.copy()
# ery = err.copy()
# fitobj = kmpfit.Fitter(residuals=residuals, data=(xr, yr, ery))
# slopes = []
# offsets = []
# trials = 10000                # Number of synthetic data sets
# for i in range(trials):       # Start loop over pseudo sample
#     indx = randint(0, N, N)    # Do the resampling using an RNG
#     xr[:] = x[indx]
#     yr[:] = y[indx]
#     ery[:] = err[indx]
#     print len(xr) 
#     # Only do a regression if there are at least two different
#     # data points in the pseudo sample
#     ok = (xr != xr[0]).any()
#     
#     if (not ok):
#         print "All elements are the same. Invalid sample."
#         print xr, yr
#     else:
#         fitobj.fit(params0=[1,1])
#         print fitobj.params
#         offs, slope = fitobj.params
#         slopes.append(slope)
#         offsets.append(offs)
# 
# slopes = numpy.array(slopes) - B0
# offsets = numpy.array(offsets) - A0
# sigmaA0, sigmaB0 = offsets.std(), slopes.std()
# print "Bootstrap errors in A, B:", sigmaA0, sigmaB0






