'''
Created on Aug 11, 2015

@author: vital
'''

import pymc

import matplotlib.pyplot as plt
import numpy as np


#Generating some data for the model_difference y = m * x + n
m_true, n_true  = 3, 2
sigma_true      = 2
x_true          = 25 * (np.random.random(50) - 0.5)
y_true          = m_true * x_true + n_true

#Adding some scatter
x_data, y_data  = x_true, y_true
x_data, y_data  = np.random.normal(x_true, 2), np.random.normal(y_true, 2)

#Providiing estimates for priors
Np_lsf          = np.polyfit(x_data, y_data, 1)
m_0, n_0        = Np_lsf[0], Np_lsf[1]

#Priors
m_coef          = pymc.Normal('m_coef', m_0, 0.01)
n_coef          = pymc.Normal('n_coef', n_0, 0.01)
sigma           = pymc.Uniform('sigma', 0.0, 5.0)
 
#---Chi Square Moodel
@pymc.stochastic(observed=True, trace=True)
def model(value = y_data, x_values = x_data, m = m_coef, n = n_coef, sigma = sigma):
    value_theo      = m*x_values + n
    chi_sq          = np.sum( np.square(value - value_theo) / np.square(sigma))
    log_ChiSq       = - chi_sq / 2.0
    return log_ChiSq


MCMC_dict2      = dict(m_coef=m_coef, n_coef=n_coef, sigma=sigma, model=model)
M               = pymc.MCMC(MCMC_dict2)
M.sample(iter=10000, burn=100)

#Code variables
print M.variables

#MCMC ouput
print '\nInitial m, n, estimations'
print m_0, n_0

print 'Bayesian estimation Model 2' 
print M.m_coef.value, M.n_coef.value, 'with Sigma:', M.sigma.value

#Store pymc_tracers
# MCMC_Traces = [M.trace('m_coef')[:], M.trace('n_coef')[:], M.trace('sigma')[:], M.trace('model')[:]]

# #Plotting regression
Fig    = plt.figure(figsize = (16, 9))  
Axis1  = Fig.add_subplot(111)
Fig.set_facecolor('w')
Axis1.plot(x_data, y_data, 'ok')
Axis1.plot(x_data, m_0 * x_data + n_0, '-', color='blue', label = 'least-squares fit')
Axis1.plot(x_data, M.m_coef.value * x_data + M.n_coef.value, '-', color='red',  label = r'$\chi^{-2}$ model')
Axis1.plot(x_data, m_true * x_data + n_true, '-', color='black',  label = r'true data model, $\sigma=$'+str(sigma_true))
# Axis1.xlabel('x')
# Axis1.ylabel('y')
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
