#!/usr/bin/python

import pymc
from CodeTools.PlottingManager                              import myPickle
from Astro_Libraries.Abundances_InferenceModel_Helium_v13   import HeAbundance_InferenceStructure


#Default address to store the data
Config              = '1000_100_NebularContinuum'
Databases_Folder    = '/home/vital/Desktop/db_Testing/'
Db_name             = 'he_Abundance_' + Config
Db_global_file_name = 'he_Abundance_Global_' + Config

#Generate dazer object
pv = myPickle()
bm = HeAbundance_InferenceStructure()

#Define plot frame and colors
pv.FigFormat_One(ColorConf = 'Night1')

testing_model = 'Model3'

# bm.Calculate_Synthetic_Fluxes(testing_model)
 
#Declare synthetic data
bm.Import_Synthetic_Fluxes(Model = testing_model)
         
#Check the observational data for Hydrogen and Helium lines which can be used for the analysis:
bm.Check_EmissionLinesObserved(bm.SyntheticData_dict)
         
#Load the data into vectors for the MCMC
bm.Preload_Data(Hbeta_Normalized = True, Deblend_Check= False)
          
#Declare the MCMC dictionary
if '_abs' not in testing_model:
    MCMC_dict = bm.Bayesian_HeliumAbundance_Analysis()
else:
    MCMC_dict = bm.Bayesian_HeliumAbundance_Analysis_abs()
          
#Run MCMC with MAP
MAP_Model               = pymc.MAP(MCMC_dict)
MAP_Model.fit(method    = 'fmin_powell') 
   
M = pymc.MCMC(MAP_Model.variables, db = 'pickle', dbname =  Databases_Folder + Db_name)
M.sample(iter=1000, burn=100, thin=2)
M.write_csv(Databases_Folder + Db_global_file_name, variables=['ChiSq', 'He_abud', 'T_e', 'n_e', 'Tau', 'c_Hbeta', 'Xi'])
M.db.close() 
     
print 'Analysis completed'
print 'He_abud map', MAP_Model.y_plus.value
print 'T_e map', MAP_Model.Te.value
print 'n_e map', MAP_Model.ne.value
print 'Tau map', MAP_Model.tau.value
print 'c_Hbeta map', MAP_Model.cHbeta.value
print 'Xi map', MAP_Model.xi.value