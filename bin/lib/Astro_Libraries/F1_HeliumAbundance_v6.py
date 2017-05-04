#!/usr/bin/python

import pymc

from Astro_Libraries.Abundances_InferenceModel_Helium_v6    import HeAbundance_InferenceStructure
from CodeTools.PlottingManager                              import myPickle


#Default address to store the data
Config              = '_30000_5000_10_NoNuissanceParameters_8_Model4_continuous'
Databases_Folder    = '/home/vital/Workspace/X_ModelData/MCMC_databases/' 
Db_name             = 'he_Abundance_' + Config
Db_global_file_name = 'he_Abundance_Global_' + Config

#Generate dazer object
pv = myPickle()
bm = HeAbundance_InferenceStructure()

#Define plot frame and colors
pv.FigFormat_One(ColorConf = 'Night1')

testing_model = 'Model1'

# bm.Calculate_Synthetic_Fluxes(testing_model)
 
# #Declare synthetic data
bm.Import_Synthetic_Fluxes(Model = testing_model)
     
#Check the observational data for Hydrogen and Helium lines which can be used for the analysis:
bm.Check_EmissionLinesObserved(bm.SyntheticData_dict)
     
#Load the data into vectors for the MCMC
bm.Preload_Data(Hbeta_Normalized = True, Deblend_Check= False)
     
#Declare the MCMC dictionary
MCMC_dict = bm.Bayesian_HeliumAbundance_Analysis()
     
#Run MCMC with MAP
MAP_Model           = pymc.MAP(MCMC_dict)
MAP_Model.fit(method = 'fmin_powell') 
M                   = pymc.MCMC(MCMC_dict, db = 'pickle', dbname =  Databases_Folder + Db_name)
M.sample(iter=30000, burn=5000, thin=10) #Optimun test
M.write_csv(Databases_Folder + Db_global_file_name, variables=['ChiSq', 'He_abud', 'T_e', 'n_e','abs_H', 'abs_He', 'Tau', 'c_Hbeta', 'Xi', 'Eqw_hbeta'])
M.db.close() 
  
print 'Analysis completed'
