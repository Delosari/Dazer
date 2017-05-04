'''
Created on Aug 11, 2015

@author: vital
'''

from Plotting_Libraries.bayesian_data import bayes_plotter
import matplotlib.pyplot as plt
import numpy as np


def GetLower_ChiSquare(ChiArray, VariableArray, mean, min_value = 0.80, max_value = 1.20, nbins = 300):
    Bins                        = np.linspace(mean*min_value, mean*max_value, nbins)
    Points_Index                = np.digitize(VariableArray, Bins)
        
    Unique_Bins             = np.unique(Points_Index)
    Min_Index, Max_Index    = min(Unique_Bins), max(Unique_Bins) + 1
    
    Min_values_Chi              = np.zeros(len(Points_Index))
    Min_values_Variable         = np.zeros(len(Points_Index))
    
    for i in range(Min_Index, Max_Index):
        inBinIndex              = np.where(Points_Index == i)[0]

        if len(ChiArray[inBinIndex]) != 0:
            Min_values_Chi[i]       = np.min(ChiArray[inBinIndex])
            Min_values_Variable[i]  = Bins[i]
            print 'indice', np.where(ChiArray[inBinIndex] == Min_values_Chi[i])[0][0]

    return Min_values_Variable, Min_values_Chi

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
        L_prob          = - chi_sq / 2.0
        return L_prob

#   #---Chi Square Moodel
    @pymc.deterministic()
    def Chi(value = y_data, x_values = x_data, m = m_coef, n = n_coef, sigma = sigma):
        value_theo      = m*x_values + n
        ChiSq             = np.sum( np.square(value - value_theo) / np.square(sigma))
        return ChiSq

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
print M.Chi.value

print 'Las variables'
print M.variables

#Store pymc_tracers
MCMC_Traces = [M.trace('m_coef')[:], M.trace('n_coef')[:], M.trace('sigma')[:], M.trace('Chi')[:]]

#Display variables convergence graphs
# pymc.Matplot.plot(M.trace('m_coef'))
# plt.show()
# pymc.Matplot.plot(M.trace('n_coef'))
# plt.show()
# pymc.Matplot.plot(M.trace('sigma'))
# plt.show()
# pymc.Matplot.plot(M.trace('Chi'))
# plt.show()

# #Sigma Maps plots
# dz.plot_SigmaLevels_MCMC(x_data, y_data, M.trace('n_coef')[:], M.trace('m_coef')[:], colors='red', linewidths=2)
# Axis1.set_ylabel('m coefficient')
# Axis1.set_xlabel('n coefficient')
# plt.show()

#Chi versus variable
Fig.set_facecolor('w')

# Obs_BestChi = where((Obs_Wavelength < (Wmin + self.SpectraEdges_Limit)) | (Obs_Wavelength > (Wmax - self.SpectraEdges_Limit)))[0]
Obs_n_indx          = np.where((M.trace('n_coef')[:] < n_true*1.01) & (M.trace('n_coef')[:] > n_true*0.99))[0]
Obs_BestChi_indx    = np.where((M.trace('Chi')[:] < 4) & (M.trace('Chi')[:] > 2))[0]
Limit1              = (M.trace('Chi')[:] < 4) & (M.trace('Chi')[:] > 2)
Limit2              = (M.trace('n_coef')[:] < n_true*1.05) & (M.trace('n_coef')[:] > n_true*0.95)
Obs_n_Chi_indx      = np.where(Limit1  &  Limit2)[0]
# Obs_n_Chi_indx      = np.where( ((M.trace('Chi')[:] < 4) & (M.trace('Chi')[:] > 2))   &  ((M.trace('n_coef')[:] < n_true*1.10) & (M.trace('n_coef')[:] > n_true*0.90)))[0]

x_min, y_min        = GetLower_ChiSquare(M.trace('Chi')[:], M.trace('m_coef')[:], 3.0)


print 'Indexes'
print Obs_BestChi_indx
print M.trace('Chi')[:]
print type(M.trace('Chi')[:])

Axis1.plot(M.trace('m_coef')[:], M.trace('Chi')[:], 'o', color='blue', label = r'$\chi^{-2}$ evolution')
Axis1.plot(M.m_coef.value, M.Chi.value, 'o', color='red', label = r'Bayesian Predicion')
Axis1.plot(M.trace('m_coef')[Obs_BestChi_indx], M.trace('Chi')[Obs_BestChi_indx], 'o', color='green')
Axis1.plot(M.trace('m_coef')[Obs_n_indx], M.trace('Chi')[Obs_n_indx], 'o', color='orange')
Axis1.plot(M.trace('m_coef')[Obs_n_Chi_indx], M.trace('Chi')[Obs_n_Chi_indx], 'x', color='red')
Axis1.plot(x_min, y_min, 'o', color='purple')

Axis1.legend(loc='best')
Axis1.set_ylim(-0.5, 3.5)
Axis1.set_xlim(2.8,3.2)
plt.show()



#Plotting regression
# Fig.set_facecolor('w')
# Axis1.plot(x_data, y_data, 'ok')
# Axis1.plot(x_data, M.m_coef.value * x_data + M.n_coef.value, '-', color='red',  label = r'$\chi^{-2}$ Model')
# # Axis1.plot(x_data, m_true * x_data + n_true, ':', color='black',  label = r'True data Model, $\sigma=$'+str(sigma_true))
# Axis1.plot(x_data, m_true * x_data + n_true, ':', color='black',  label = r'True data Model')
# Axis1.legend(loc='best')
# plt.show()


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
