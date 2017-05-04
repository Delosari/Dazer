'''
Created on Oct 22, 2015

@author: vital
'''
from matplotlib                         import rc, image, rcParams
import pymc

from Plotting_Libraries.bayesian_data   import bayes_plotter


rc('legend', frameon = True)
rcParams['svg.fonttype'] = 'none'
#Import bayesian plotter class
bp                      = bayes_plotter()
  
#Database location
database_name           = 'he_Abundance__30000_5000_10_NoNuissanceParameters_EqwHbetaPrior5_Model3_VeryGood_4649sec'
database_folder         = '/home/vital/Workspace/X_ModelData/MCMC_databases/'
  
#Variables to plot
Traces_code             = ['He_abud', 'n_e',        'T_e',      'abs_H',    'abs_He', 'c_Hbeta',        'Xi',       'Tau']
Traces_labels           = [r'$y^{+}$', r'$n_{e}$',  r'$T_{e}$',  r'$a_{H}$', r'$a_{He}$', r'$c(H\beta)$', r'$\xi$',   r'$\tau$' ]
  
#Load database
bp.load_pymc_database(database_folder+database_name)
  
#Extract variable stats
statistics_dict         = bp.extract_traces_statistics(Traces_code)
  
#Plot with the data
bp.plot_posteriors_histagram(Traces_code, Traces_labels)
bp.savefig(database_folder + 'posteriors_hist9_best', extension='.eps')
   
bp.plot_tracers(Traces_code, Traces_labels)
bp.savefig(database_folder + 'traces_evolution9_best', extension='.eps')
   
bp.plot_acorrelation(Traces_code, Traces_labels)
bp.savefig(database_folder + 'traces_acorr9_best', extension='.eps')
  
bp.plot_chiSq_Behaviour(Traces_code, Traces_labels)
bp.savefig(database_folder + 'ChiSq_Parameters9_best', extension='.eps')
  
bp.close_database()
 
print 'All data treated'

# # extract and plot results
# y_min = MDL.stats()['gauss']['quantiles'][2.5]
# y_max = MDL.stats()['gauss']['quantiles'][97.5]
# y_fit = MDL.stats()['gauss']['mean']
# plt.plot(x_true,f_true,'b', marker='None', ls='-', lw=1, label='True')
# plt.errorbar(x_true,f,yerr=f_error, color='r', marker='.', ls='None', label='Observed')
# plt.plot(x_true,y_fit,'k', marker='+', ls='None', ms=5, mew=1, label='Fit')
# plt.fill_between(x_true, y_min, y_max, color='0.5', alpha=0.5)
# plt.legend()

# print 'Las variables...', M.variables
# # map(str, M.variables) THIS WAS COOL
# M.write_csv(StoringDataFolder + csv, variables=['m_coef', 'sigma', 'Chi', 'n_coef'])
# M.db.close()
 
# db = pymc.database.pickle.load(StoringDataFolder + db_name)
# db = py
# 
# # #MCMC ouput
# print '\nInitial m, n, estimations'
# print m_0, n_0
# #
# # db.write_csv(StoringDataFolder + csv, variables=['m_coef', 'sigma', 'Chi', 'n_coef'])


 
# chi_trace =  db.trace('Chi')
# 
# pymc.Matplot.plot(chi_trace)
# 
# db.close()

# print 'Bayesian estimation Model' 
# print db.m_coef.value, db.n_coef.value, 'with Sigma:', db.sigma.value
# print db.Chi.value

# #MCMC ouput
# print '\nInitial m, n, estimations'
# print m_0, n_0
# 
# print 'Bayesian estimation Model' 
# print M.m_coef.value, M.n_coef.value, 'with Sigma:', M.sigma.value
# print M.Chi.value
# 
# print 'Las variables'
# print M.variables
# 
# #Store pymc_tracers
# MCMC_Traces = [M.trace('m_coef')[:], M.trace('n_coef')[:], M.trace('sigma')[:], M.trace('Chi')[:]]

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

# #Chi versus variable
# Fig.set_facecolor('w')
# 
# # Obs_BestChi = where((Obs_Wavelength < (Wmin + self.SpectraEdges_Limit)) | (Obs_Wavelength > (Wmax - self.SpectraEdges_Limit)))[0]
# Obs_n_indx          = np.where((M.trace('n_coef')[:] < n_true*1.01) & (M.trace('n_coef')[:] > n_true*0.99))[0]
# Obs_BestChi_indx    = np.where((M.trace('Chi')[:] < 4) & (M.trace('Chi')[:] > 2))[0]
# Limit1              = (M.trace('Chi')[:] < 4) & (M.trace('Chi')[:] > 2)
# Limit2              = (M.trace('n_coef')[:] < n_true*1.05) & (M.trace('n_coef')[:] > n_true*0.95)
# Limit3              = (M.trace('sigma')[:] < sigma_true*1.05) & (M.trace('sigma')[:] > sigma_true*0.95)
# Obs_n_Chi_indx      = np.where(Limit1  &  Limit2 & Limit3)[0]
# # Obs_n_Chi_indx      = np.where( ((M.trace('Chi')[:] < 4) & (M.trace('Chi')[:] > 2))   &  ((M.trace('n_coef')[:] < n_true*1.10) & (M.trace('n_coef')[:] > n_true*0.90)))[0]
# 
# x_min, y_min        = GetLower_ChiSquare(M.trace('Chi')[:], M.trace('m_coef')[:], 3.0)
# 
# 
# print 'Indexes'
# print Obs_BestChi_indx
# print M.trace('Chi')[:]
# print type(M.trace('Chi')[:])
# 
# # Axis1.plot(M.trace('m_coef')[:], M.trace('Chi')[:], 'o', color='blue', label = r'$\chi^{-2}$ evolution')
# # Axis1.plot(M.m_coef.value, M.Chi.value, 'o', color='red', label = r'Bayesian Predicion')
# # Axis1.plot(M.trace('m_coef')[Obs_BestChi_indx], M.trace('Chi')[Obs_BestChi_indx], 'o', color='green')
# # Axis1.plot(M.trace('m_coef')[Obs_n_indx], M.trace('Chi')[Obs_n_indx], 'o', color='orange')
# Axis1.plot(x_min, y_min, 'o', color='purple')
# Axis1.plot(M.trace('m_coef')[Obs_n_Chi_indx], M.trace('Chi')[Obs_n_Chi_indx], 'o', color='red')
# 
# Axis1.legend(loc='best')
# # Axis1.set_ylim(-0.5, 3.5)
# # Axis1.set_xlim(2.8,3.2)
# plt.show()
# 
# 
# 
# #Plotting regression
# Fig.set_facecolor('w')
# Axis1.plot(x_data, y_data, 'ok')
# Axis1.plot(x_data, M.m_coef.value * x_data + M.n_coef.value, '-', color='red',  label = r'$\chi^{-2}$ Model')
# # Axis1.plot(x_data, m_true * x_data + n_true, ':', color='black',  label = r'True data Model, $\sigma=$'+str(sigma_true))
# Axis1.plot(x_data, m_true * x_data + n_true, ':', color='black',  label = r'True data Model')
# Axis1.legend(loc='best')
# plt.show()
# 
# 
# #Saving output
# # In [1]:
# # %matplotlib inline
# # from pymc import MCMC, database, Matplot
# # from pymc.examples import gelman_bioassay
# # In [2]:
# # m = MCMC(gelman_bioassay, db='pickle', dbname='bioassay.pickle')
# # m.sample(iter = 10000, burn = 5000, thin = 10)
# # m.db.close()
# #  [-----------------100%-----------------] 10000 of 10000 complete in 0.7 sec
# # In [3]:
# # db = database.pickle.load('bioassay.pickle')
# # In [4]:
# # alpha = db.trace('alpha')
# # In [6]:
# # Matplot.plot(alpha)
# # Plotting alpha