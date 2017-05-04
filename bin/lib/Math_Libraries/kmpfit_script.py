#!/usr/bin/env python
#------------------------------------------------------------
# Purpose: Program to straight line parameters
#          to data with errors in both coordinates.
#          Use the (famous) data set from Pearson with
#          weights from York
# Vog, 12 Dec, 2011
#
# The data for x and y are from Pearson
# Pearson, K. 1901. On lines and planes of closest fit to systems 
# of points in space. Philosophical Magazine 2:559-572
# Copy of this article can be found at:
# stat.smmu.edu.cn/history/pearson1901.pdf
#
# Pearson's best fit through (3.82,3.70) ->
# a=5.7857  b=-0.546
# York added weights in 
# York, D. Can. J. Phys. 1968, 44, 1079-1086
# The Williamson Approach is also implemented.
# The steps are described in Ogren, Paul J., Norton, J. Russel
# Applying a simple Least-Squares Algorithm to Data
# with Uncertainties in Both Variables,
# J. of Chem. Education, Vol 69, Number 4, April 1992
# Best fit parameters for this method are:
# a=5.47991022403  b=-0.48053340745
#------------------------------------------------------------

# from collections        import OrderedDict
# 
# # from kapteyn            import kmpfit
# from matplotlib.pyplot  import figure, show, rc
# from numpy              import where, array, sqrt, linspace
# from scipy.odr          import Model, ODR, RealData
# 
# 
# def model(p, x):
#     
#     a, b = p
#     return a + b*x
# 
# def residuals(p, data):
#     # Residuals function for data with errors in both coordinates
#     
#     a, b            = p
#     x, y, ex, ey    = data
#     w               = ey*ey + ex*ex*b*b
#     wi              = sqrt(where(w==0.0, 0.0, 1.0/(w)))
#     d               = wi*(y-model(p,x))
#     return d
# 
# def residuals2(p, data):
#     # Residuals function for data with errors in y only
#     
#     a, b            = p
#     x, y, ey        = data
#     wi              = where(ey==0.0, 0.0, 1.0/ey)
#     d               = wi*(y-model(p,x))
#     return d
# 
# # Pearsons data with York's weights 
# 
# def kmpfit_effectivevariance(x, y, errx, erry, beta0):
#     
#     # Prepare fit routine
#     fitobj = kmpfit.Fitter(residuals=residuals, data=(x, y, errx, erry))
#     fitobj.fit(params0=beta0)
# #     print "\n\n======== Results kmpfit: weights for both coordinates ========="
# #     print "Fitted parameters:      ", fitobj.params
# #     print "Covariance errors:      ", fitobj.xerror
# #     print "Standard errors         ", fitobj.stderr
# #     print "Chi^2 min:              ", fitobj.chi2_min
# #     print "Reduced Chi^2:          ", fitobj.rchi2_min
# #     print "Iterations:             ", fitobj.niter
#     
#     return fitobj.params[1], fitobj.params[0], fitobj.stderr[1], fitobj.stderr[0], fitobj.xerror, fitobj.chi2_min, fitobj.rchi2_min
# 
# def scipy_ODR(x, y, errx, erry, beta0):
#     
#     # Compare result with ODR
#     linear      = Model(model)
#     mydata      = RealData(x, y, sx=errx, sy=erry)
#     myodr       = ODR(mydata, linear, beta0=beta0, maxit=5000, sstol=1e-14)
#     myoutput    = myodr.run()
#     
#     #WARNING: THE COVARIANCE VALUES ARE N, M (NOT M,N LIKE THE REST OF THE METHODS)
#     return myoutput.beta[1], myoutput.beta[0], myoutput.sd_beta[1], myoutput.sd_beta[0], sqrt(myoutput.cov_beta.diagonal()), myoutput.sum_square, myoutput.res_var

# x   = array([0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4])
# y   = array([5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5])
# wx  = array([1000.0,1000,500,800,200,80,60,20,1.8,1.0])
# wy  = array([1,1.8,4,8,20,20,70,70,100,500])
# errx = 1/sqrt(wx)  # We need the errors in the residuals functions
# erry = 1/sqrt(wy)
# N = len(x)
# 
# beta0 = [5.0, 1.0]         # Initial estimates
# 
# fitobj =   kmpfit_effectivevariance(x, y, errx, erry, beta0)
# 
# scipy_ODR(x, y, errx, erry,)
# 
# chi2 = (residuals(fitobj.params,(x,y,errx,erry))**2).sum()
# 
# 
# print "\n\nReference                    a                b"
# print "-----------------------------------------------------------"
# print "Literature results:"
# print "Pearson unweighted           5.7857           -0.546"
# print "Williamson                   5.47991022403    -0.48053340745"
# print "Reed                         5.47991022723    -0.48053340810"
# print "Lybanon                      5.47991025       -0.480533415"
# print 
# print "Practical results:"
# print "kmpfit effective variance    %13.11f    %13.11f"%(fitobj.params[0],fitobj.params[1])
# 
# # Some plotting
# rc('font', size=9)
# rc('legend', fontsize=8)
# fig = figure(1)
# d = (x.max() - x.min())/10
# X = linspace(x.min()-d, x.max()+d, 50)
# frame = fig.add_subplot(1,1,1, aspect=1, adjustable='datalim')
# frame.errorbar(x, y, xerr=errx, yerr=erry,  fmt='bo')
# frame.plot(X, model(fitobj.params,X), 'c', ls='--', lw=2, label="kmpfit effective variance")
# frame.plot(X, model((5.463,-0.477),X), 'm', label="York's values")
# frame.set_xlabel("X")
# frame.set_ylabel("Y")
# frame.set_title("Pearson's data with York's weights")
# frame.grid(True)
# leg = frame.legend(loc=1)
# show()


# 
# # Compare result with ODR
# linear = Model(model)
# mydata = RealData(x, y, sx=errx, sy=erry)
# myodr = ODR(mydata, linear, beta0=beta0, maxit=5000, sstol=1e-14)
# myoutput = myodr.run()
# print "\n\n======== Results ODR ========="
# print "Fitted parameters:      ", myoutput.beta
# print "Covariance errors:      ", numpy.sqrt(myoutput.cov_beta.diagonal())
# print "Standard errors:        ", myoutput.sd_beta
# print "Minimum chi^2:          ", myoutput.sum_square
# print "Minimum (reduced)chi^2: ", myoutput.res_var
# beta = myoutput.beta
