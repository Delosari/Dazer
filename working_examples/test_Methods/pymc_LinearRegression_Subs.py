import pymc

import matplotlib.pyplot as plt
import numpy as np


#Generating some data for the model y = m * x + n
m_true, n_true  = 3, 2
x_true          = 25 * (np.random.random(50) - 0.5)
y_true          = m_true * x_true + n_true

#Adding some scatter
x_data, y_data  = x_true, y_true
x_data, y_data  = np.random.normal(x_true, 2), np.random.normal(y_true, 2)

#Initial values for priors
Np_lsf          = np.polyfit(x_data, y_data, 1)
m_0, n_0        = Np_lsf[0], Np_lsf[1]

#Priors
m_coef, n_coef  = pymc.Normal('m_coef', n_0, 0.01), pymc.Normal('n_coef', m_0, 0.01)
sigma           = pymc.Uniform('sigma', 0.0, 10.0)

#Model
@pymc.deterministic()
def model(x_values=x_data, m=m_coef, n=n_coef):
    return m*x_values + n
 
#Likelihood
L               = pymc.Normal('L', mu=model, tau=1.0/sigma**2, value=y_data, observed=True)

#Running the MCMC
MCMC_dict       = dict(m_coef=m_coef, n_coef=n_coef, sigma=sigma, model=model, L=L)
M               = pymc.MCMC(MCMC_dict)
M.sample(iter=10000)
# S.sample(iter=100000)

#MCMC ouput
print '\nInitial m, n, ValuesEstimations'
print m_0, n_0
print 'Bayesian estimation' 
print M.m_coef.value, M.n_coef.value, 'with Sigma:', M.sigma.value

#Plotting regression
Fig    = plt.figure(figsize = (16, 9))  
Fig.set_facecolor('w')
plt.plot(x_data, y_data, 'ok')
plt.plot(x_data, M.m_coef.value * x_data + M.n_coef.value, '-', color='red',  label = r'Bayesian model')
plt.xlabel('x')
plt.ylabel('y');
plt.legend()
plt.show()