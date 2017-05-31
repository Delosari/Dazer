import pymc
from bin.lib.Astro_Libraries.Abundances_InferenceModel_Helium_v14 import HeAbundance_InferenceStructure

#Default address to store the data
Config              = '_10000_1500_2_NoAbs_Model3_MAP2MCMC_Checkingv12'
Databases_Folder    = '/home/vital/Astrodata/Inferece_output/'
db_name             = 'he_Abundance_' + Config
db_global_filename  = 'he_Abundance_Global_' + Config

#Generate dazer object
bm = HeAbundance_InferenceStructure()

#Generate the synthetic data
bm.calculate_synthFluxes(model = 'Model3')
        
# #Check the observational data for Hydrogen and Helium lines which can be used for the analysis:
# bm.Check_EmissionLinesObserved(bm.SyntheticData_dict)
#         
# #Load the data into vectors for the MCMC
# bm.Preload_Data(Hbeta_Normalized = True, Deblend_Check= False)
#         
# #Declare the MCMC dictionary
# if '_abs' not in testing_model:
#     MCMC_dict = bm.Bayesian_HeliumAbundance_Analysis()
# else:
#     MCMC_dict = bm.Bayesian_HeliumAbundance_Analysis_abs()
#         
# #Run MCMC with MAP
# MAP_Model               = pymc.MAP(MCMC_dict)
# MAP_Model.fit(method    = 'fmin_powell') 
#  
# M  = pymc.MCMC(MAP_Model.variables, db = 'pickle', dbname =  Databases_Folder + Db_name)
# M.sample(iter=10000, burn=1500, thin=2) #Optimun test
# M.write_csv(Databases_Folder + Db_global_file_name, variables=['ChiSq', 'He_abud', 'T_e', 'n_e', 'Tau', 'c_Hbeta', 'Xi'])
# M.db.close() 
#    
# print 'Analysis completed'
# print 'He_abud map', MAP_Model.y_plus.value
# print 'T_e map', MAP_Model.Te.value
# print 'n_e map', MAP_Model.ne.value
# print 'Tau map', MAP_Model.tau.value
# print 'c_Hbeta map', MAP_Model.cHbeta.value
# print 'Xi map', MAP_Model.xi.value