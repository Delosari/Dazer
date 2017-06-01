import pymc
from bin.lib.Astro_Libraries.Abundances_InferenceModel_Helium_v14 import HeAbundance_InferenceStructure
# from bin.dazer_methods import Dazer

# dz = Dazer()

#Default address to store the data
Config              = '_10000_1500_2_NoAbs_Model3_MAP2MCMC_v0_CorrectingFit'
database_folder     = '/home/vital/Astrodata/Inferece_output/'
db_name             = 'he_Abundance_' + Config
db_global_filename  = 'he_Abundance_Global_' + Config
physical_model      = 'Model3'

#Generate dazer object
bm = HeAbundance_InferenceStructure()
  
#Generate the synthetic data
bm.calculate_synthFluxes(model = physical_model)
                   
#Declare the MCMC dictionary
if '_abs' not in physical_model:
    MCMC_dict = bm.helium_abundance_model()
else:
    MCMC_dict = bm.helium_abundance_model_abs()

#Run MCMC with MAP
print 'Starting prefit analysis'
MAP_Model = pymc.MAP(MCMC_dict)
MAP_Model.fit(method = 'fmin_powell') 
  
print 'Starting run'
M  = pymc.MCMC(MAP_Model.variables, db = 'pickle', dbname =  database_folder + db_name)
M.sample(iter=10000, burn=1500, thin=2)
# M.sample(iter=30000, burn=5000, thin=10) #Optimun test
 
print 'Running the data'
db  = pymc.database.pickle.load(database_folder + db_name)
M   = pymc.MCMC(MAP_Model, db=db)
M.write_csv(database_folder + db_global_filename, variables=['ChiSq','He_abud', 'T_e', 'n_e', 'Tau', 'c_Hbeta', 'Xi'])
M.db.close() 
 
print '--Analysis completed'
print 'He_abud map',    MAP_Model.y_plus.value
print 'T_e map',        MAP_Model.Te.value
print 'n_e map',        MAP_Model.ne.value
print 'Tau map',        MAP_Model.tau.value
print 'c_Hbeta map',    MAP_Model.cHbeta.value
print 'Xi map',         MAP_Model.xi.value

# dz.load_pymc_database(database_folder + db_name)
# dz.extract_traces_statistics(['He_abud', 'T_e', 'n_e', 'Tau', 'c_Hbeta', 'Xi'])
# 
# for coso in dz.statistics_dict:
#     print coso, dz.statistics_dict[coso]['mean']

