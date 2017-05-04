'''
Created on Oct 15, 2015

@author: vital
'''

from numpy                                      import linspace, diag, zeros, mean, asarray, array, sqrt

from Astro_Libraries.Abundances_Class           import Parametrized_Emissivities
from CodeTools.PlottingManager                  import myPickle
from Math_Libraries import bces_script
import pymc as pm


def fonnesbeck_example(Oxygen_Flux, Ymagnitude, verbose = True):

    # Linear regression
    n   = pm.Normal('n', 0.0, tau=1e-5, value=0)
    m   = pm.Normal('m', 0.0, tau=1e-5, value=0)

    # No container needed
    mu  = n + m*Oxygen_Flux
    
    #Precision    
    tau = pm.Gamma('tau', 0.01, 0.01, value=0.01)    
        
    # Don't have to sepcify size variable
    resp = pm.Normal('resp', mu, tau=tau, value=Ymagnitude, observed=True) 
    
    mcmc = pm.MCMC(locals())

    mcmc.sample(30000, burn=10000)
    
    if verbose:
        print n.value
        print(n.summary())
    
    return m.value, n.value 

def aflaxman_example(x, y, x_err, y_err, verbose = True):
    
    #Priors coefficients
    m   = pm.Normal('m', 0.0, tau=1e-5, value=0)
    n   = pm.Normal('n', 0.0, tau=1e-5, value=0)
    
    x_err_mean = mean(x_err)
    
    x_pred = pm.Normal('x_prediction', mu=x, tau=(x_err_mean)**-2)
    
    @pm.deterministic(plot=False)
    def linear_equation(x_true = x_pred, m_coeff = m, n_coeff = n):
        return m_coeff * x_true + n_coeff
    
    y = pm.Normal('y', mu = linear_equation, tau = 1.0 / y_err**2, value = y, observed=True)
    
    mcmc2 = pm.MCMC(locals())

    mcmc2.sample(30000, burn=10000)

    if verbose:
        print n.value
        print(n.summary())

    return m.value, n.value 

def ChiSq_inference(Oxygen_Flux, Ymagnitude, Yerror, verbose = True):
    
    print 'Madre mia'
    
    return

def Python_linfit(x_true, y, y_err):
     
    Regression_Fit, Uncertainty_Matrix, Red_Chi_Sq, Residuals   = linfit(x_true, y, y_err, cov=True, relsigma=False, chisq=True, residuals=True)
    m_n_Matrix                                                  = [sqrt(Uncertainty_Matrix[t,t]) for t in range(2)] 
    R_Factor                                                    = Uncertainty_Matrix[0,1]/(m_n_Matrix[0]*m_n_Matrix[1])
    m, m_error                                                  = Regression_Fit[0], m_n_Matrix[0]
    n, n_error                                                  = Regression_Fit[1], m_n_Matrix[1]
     
    return m, m_error, n, n_error
 