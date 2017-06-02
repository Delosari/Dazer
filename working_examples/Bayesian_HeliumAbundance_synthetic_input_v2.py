from bin.dazer_methods import Dazer
from bin.lib.Astro_Libraries.Abundances_InferenceModel_Helium_v15 import Run_MCMC

#MCMC configuration run
default_db_folder   = '/home/vital/Astrodata/Inferece_output/'
db_name_ext         = 'v2_CorrectingFit'
iterat, burn, thin  = 10000, 1500, 2
synth_model         = 'Model3'
params_list         = ['ChiSq','y_plus', 'Te', 'ne', 'tau', 'cHbeta', 'xi']             

#Generate dazer object
dz = Dazer()
bm = Run_MCMC()
  
#Generate the synthetic data
bm.calculate_synthFluxes(synth_model)
     
#Select right model according to data
bm.select_inference_model(synth_model)     
            
#Variables to save
#db_address  = '{}HeAbund_pymc2_{}_it{}_burn{}_thin{}_{}'.format(default_db_folder, synth_model, iterat, burn, thin, db_name_ext) 
db_address  = '{}HeAbund_pymc2_it{}_burn{}_thin{}_{}'.format(default_db_folder, iterat, burn, thin, db_name_ext) 
                
# #Run sampler
# bm.run_pymc2(iterat, burn, thin, db_address, variables_list=params_list)

#Load database
pymc2_db, stat_db_dict = bm.load_pymc_database(db_address)

#Traces plot
dz.traces_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_tracesPlot', save_pickle = False)
 
#Posteriors plot
dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_posteriorPlot', save_pickle = False)
 
#Posteriors plot
dz.acorr_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_acorrPlot', save_pickle = False)

#Corner plot
dz.corner_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_cornerPlot', save_pickle = False)

print '\nData treated'

