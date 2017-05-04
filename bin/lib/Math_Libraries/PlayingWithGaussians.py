'''
Created on Aug 6, 2015

@author: vital
'''
from bin.DZ_DataExplorer                import Plots_Manager

from kapteyn import kmpfit
from numpy import zeros, linspace, exp, ones
from scipy.optimize     import minimize, curve_fit, leastsq

from FittingTools import Fitting_Gaussians
import matplotlib.pyplot as plt


def SingleGaussian_Cont(Ind_Variables, A, mu, sigma):

    return A * exp(-(Ind_Variables[0]-mu)*(Ind_Variables[0]-mu)/(2.0*sigma*sigma)) + Ind_Variables[1]

def Residuals_Gaussian(p, data):
    x_true, y, zerolev = data[0], data[1], data[2]
    return (y - (SingleGaussian_Cont((x_true, zerolev), p[0], p[1], p[2])))


def funcG(p, x_true):
    # Model function is a gaussian
    A, mu, sigma, zerolev = p
    return( A * exp(-(x_true-mu)*(x_true-mu)/(2*sigma*sigma)) + zerolev )

def residualsG(p, data):
    # Return weighted residuals of Gauss
    x_true, y = data
    return y-funcG(p,x_true)

pv = Plots_Manager()


dz = Fitting_Gaussians()
mu, sigma   = 0, 1
A           = (1 / ((2*3.1415)**0.5) * sigma)
x_true           = linspace(-2, 2, 100)
y           = A*exp(-(x_true-mu)**2/(2*sigma**2))
zerolev     = zeros(len(x_true))
p_1, conv = leastsq(Residuals_Gaussian, [A, mu, sigma], [x_true, y, zerolev])
print 'predictions initial', p_1



L           = exp(-y / 2)
p0 = [-0.2, 0, 1, max(L)]
fitterG = kmpfit.Fitter(residuals=residualsG, data=(x_true,L))
fitterG.fit(params0=p0)
print fitterG.params
gaus_inv      = funcG(fitterG.params,x_true)





y_frac      = y / 3
p_2, conv   = leastsq(Residuals_Gaussian, [A/5, mu, sigma], [x_true, y_frac, zerolev])
print 'predictions by3', p_2
gaus_inv    = SingleGaussian_Cont([x_true, zerolev], p_2[0], p_2[1], p_2[2])


# dz.SingleGaussian_Cont(Ind_Variables, A, mu, sigma)

plt.plot(x_true, y, '-', color = 'black')
plt.plot(x_true, gaus_inv, '-', color = 'red')
# plt.plot(x_true, y_frac, '-', color = 'orange')
# plt.plot(x_true, y_simu, '-', color = 'red')
plt.plot(x_true, L, '-', color = 'blue')

plt.xlabel('x_true')
plt.ylabel('y');
plt.show()