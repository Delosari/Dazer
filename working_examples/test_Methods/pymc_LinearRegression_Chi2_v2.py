'''
Created on Aug 11, 2015

@author: vital
'''

import pymc

from Plotting_Libraries.bayesian_data import bayes_plotter
import matplotlib.pyplot as plt
import numpy as np


def Bayesian_Regression(x_data, y_data, m_0, n_0):
        
    #Priors
    m_coef          = pymc.Normal('m_coef', m_0, 0.01)
    n_coef          = pymc.Normal('n_coef', n_0, 0.01)
    sigma           = pymc.Uniform('sigma', 0.0, 5.0)
     
    #---Chi Square Moodel
    @pymc.stochastic(observed=True)
    def model(value = y_data, x_values = x_data, m = m_coef, n = n_coef, sigma = sigma):
        value_theo      = m*x_values + n
        chi_sq          = np.sum( np.square(value - value_theo) / np.square(sigma))
        log_ChiSq       = - chi_sq / 2.0
        return log_ChiSq
    
#     L = pymc.Normal('like', mu=model, tau=1/np.square(sigma))

    @pymc.deterministic()
    def log_ChiSq(y_data = y_data, x_values = x_data, m = m_coef, n = n_coef, sigma = sigma):
        value_theo      = m*x_values + n
        chi_sq          = np.sum( np.square(y_data - value_theo) / np.square(sigma))
        log_ChiSq       = - chi_sq / 2.0
        return log_ChiSq
        
    
    return locals()

dz = bayes_plotter()
Fig    = plt.figure(figsize = (16, 9))  
Axis1  = Fig.add_subplot(111)
dz.Import_FigConf(Fig, Axis1)

#Generating some data for the model_difference y = m * x + n
m_true, n_true  = 3, 2
sigma_true      = 2
x_true          = 25 * (np.random.random(50) - 0.5)
y_true          = m_true * x_true + n_true

#Adding some scatter
x_data, y_data  = x_true, y_true
# x_data, y_data  = np.random.normal(x_true, 2), np.random.normal(y_true, 2)

#Providiing estimates for priors
Np_lsf          = np.polyfit(x_data, y_data, 1)
m_0, n_0        = 10, 0

#Runnint the hammer
MCMC_dict       = Bayesian_Regression(x_data, y_data,m_0, n_0)
M               = pymc.MCMC(MCMC_dict)
M.sample(iter=10000, burn=1000)

#MCMC ouput
print '\nInitial m, n, estimations'
print m_0, n_0

print 'Bayesian estimation Model' 
print M.m_coef.value, M.n_coef.value, 'with Sigma:', M.sigma.value

print 'Las variables'
print M.variables

#Store pymc_tracers
MCMC_Traces = [M.trace('m_coef')[:], M.trace('n_coef')[:], M.trace('sigma')[:], M.trace('log_ChiSq')[:]]

#Display variables convergence graphs
# pymc.Matplot.plot(M.trace('m_coef'))
# plt.show()
# pymc.Matplot.plot(M.trace('n_coef'))
# plt.show()
# pymc.Matplot.plot(M.trace('sigma'))
# plt.show()

# #Sigma Maps plots
# dz.plot_SigmaLevels_MCMC(x_data, y_data, M.trace('n_coef')[:], M.trace('m_coef')[:], colors='red', linewidths=2)
# Axis1.set_ylabel('m coefficient')
# Axis1.set_xlabel('n coefficient')
# plt.show()

#Chi versus variable
# Fig.set_facecolor('w')
# Axis1.plot(M.trace('m_coef')[:100], M.trace('like')[:100], '-', color='blue', label = r'$\chi^{-2}$ evolution')
# Axis1.legend(loc='best')
# plt.show()

#Plotting regression
Fig.set_facecolor('w')
Axis1.plot(x_data, y_data, 'ok')
Axis1.plot(x_data, M.m_coef.value * x_data + M.n_coef.value, '-', color='red',  label = r'$\chi^{-2}$ Model')
# Axis1.plot(x_data, m_true * x_data + n_true, ':', color='black',  label = r'True data Model, $\sigma=$'+str(sigma_true))
Axis1.plot(x_data, m_true * x_data + n_true, ':', color='black',  label = r'True data Model')
Axis1.legend(loc='best')
plt.show()


#Saving output
# In [1]:
# %matplotlib inline
# from pymc import MCMC, database, Matplot
# from pymc.examples import gelman_bioassay
# In [2]:
# m = MCMC(gelman_bioassay, db='pickle', dbname='bioassay.pickle')
# m.sample(iter = 10000, burn = 5000, thin = 10)
# m.db.close()
#  [-----------------100%-----------------] 10000 of 10000 complete in 0.7 sec
# In [3]:
# db = database.pickle.load('bioassay.pickle')
# In [4]:
# alpha = db.trace('alpha')
# In [6]:
# Matplot.plot(alpha)
# Plotting alpha
