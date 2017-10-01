from numpy import array
from dazer_methods import Dazer
import matplotlib.pyplot as plt
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v20 import Run_MCMC

#Declare operation objects
dz  = Dazer()
bm  = Run_MCMC()

#MCMC configuration run
default_db_folder   = '/home/vital/Astrodata/Inference_output/'
db_name_ext         = 'pymc3_stellar'
iterat, burn, thin  = 7000, 1000, 1 #15000, 1500, 2
synth_model         = 'initial_tests'
params_list         = ['z_star','sigma_star','Av_star', 'chiSq_ssp']
synth_model         = 'ssp_synthesis_pymc3'
db_name_ext         = 'ssp_testsErrNew'

#Generate the synthetic data
bm.calculate_synthFluxes(synth_model)

print '-----Declaring the mod'         

# #Select right model according to data
# bm.select_inference_model(synth_model) 
  
#Variables to save
db_address  = '{}{}_it{}_burn{}_thin{}_{}'.format(default_db_folder, synth_model, iterat, burn, thin, db_name_ext) 

print '-----This is the runner'
# # #Run sampler
bm.run_pymc3(db_address, iterat, burn, thin, variables_list=params_list)
   
# #Load database
# pymc2_db, stat_db_dict = bm.load_pymc_database(db_address)
#    
# #Traces plot
# dz.traces_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_tracesPlot', save_pickle = False)
#     
# #Posteriors plot
# dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_posteriorPlot', save_pickle = False)
#     
# #Posteriors plot
# dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=1, n_rows=len(params_list))
# dz.save_manager(db_address + '_acorrPlot', save_pickle = False)
#    
# #Corner plot
# dz.corner_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_cornerPlot', save_pickle = False)

