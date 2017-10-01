from dazer_methods import Dazer
import matplotlib.pyplot as plt
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v18 import Run_MCMC
from numpy import array

#Declare operation objects
dz  = Dazer()
bm  = Run_MCMC()

#MCMC configuration run
default_db_folder   = '/home/vital/Astrodata/Inference_output/'
db_name_ext         = 'SSPsynth_v1'
iterat, burn, thin  = 10000, 1000, 1 #15000, 1500, 2
synth_model         = 'initial_tests'
params_list         = ['z_star','sigma_star','Av_star', 'chiSq_ssp']
synth_model         = 'ssp_synthesis'
db_name_ext         = 'ssp_testsErrNew'

#Generate the synthetic data
bm.calculate_synthFluxes(synth_model)
 
#Select right model according to data
bm.select_inference_model(synth_model) 
  
#Variables to save
db_address  = '{}{}_it{}_burn{}_thin{}_{}'.format(default_db_folder, synth_model, iterat, burn, thin, db_name_ext) 
                    
# #Run sampler
bm.run_pymc2(db_address, iterat, burn, thin, variables_list=params_list)
   
#Load database
pymc2_db, stat_db_dict = bm.load_pymc_database(db_address)
 
# #Ready data for fit
# output_data = bm.ssp_fit(bm.obj_data['z_star'], bm.obj_data['sigma_star'], bm.obj_data['Av_star'], bm.obj_data)
# fig, axis = plt.subplots(1, 1)
# axis.plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm_masked'], label='obs_flux_norm_masked')
# axis.plot(bm.obj_data['obs_wave_resam'], output_data['flux_sspFit_norm'], label='flux_sspFit_norm', linestyle='--')
# axis.plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_fluxEr_norm'], label='obs_fluxEr_norm', linestyle='--')
# # axis.plot(bm.synth_SED['obs_wave_resam'], bm.synth_SED['obs_flux_resam'], label='Synth spectrum')
# # axis.plot(bm.synth_SED['obs_wave_resam'], bm.synth_SED['obs_flux_True'], label='True synth spectrum')
# # axis.plot(bm.synth_SED['obs_wave_resam'], bm.synth_SED['obs_fluxEr_resam'], label='Err spectrum')
# # axis.plot(bm.synth_SED['obs_wave_resam'], output_data['flux_sspFit'] , label='Output fit', linestyle=':', color='tab:red')
# axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
# axis.legend()
# plt.show() 
  
#Traces plot
dz.traces_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_tracesPlot', save_pickle = False)
    
#Posteriors plot
dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_posteriorPlot', save_pickle = False)
    
#Posteriors plot
dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=1, n_rows=len(params_list))
dz.save_manager(db_address + '_acorrPlot', save_pickle = False)
   
#Corner plot
dz.corner_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_cornerPlot', save_pickle = False)

  
  
# #Generate the synthetic data
# bm.calculate_synthFluxes(synth_model)
#      
# #Select right model according to data
# bm.select_inference_model(synth_model)     
#             
# #Variables to save
# db_address  = '{}HeAbund_pymc2_it{}_burn{}_thin{}_{}'.format(default_db_folder, iterat, burn, thin, db_name_ext) 
#                 
# # #Run sampler
# bm.run_pymc2(db_address, iterat, burn, thin, variables_list=params_list)
#  
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
# dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=4, n_rows=3)
# dz.save_manager(db_address + '_acorrPlot', save_pickle = False)
#  
# #Corner plot
# dz.corner_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_cornerPlot', save_pickle = False)
#  
# print '\nData treated'

# from dazer_methods import Dazer
# import matplotlib.pyplot as plt
# from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v18 import Run_MCMC
# from numpy import array
# 
# #Declare operation objects
# dz  = Dazer()
# bm  = Run_MCMC()
# 
# #MCMC configuration run
# ssp_library_data = {}
# ssp_library_data['default_Starlight_file']   = '/home/vital/Starlight/Dani_Bases_Extra.txt'
# ssp_library_data['default_Starlight_folder'] = '/home/vital/Starlight/Bases/'
# ssp_library_data['default_Starlight_coeffs'] = '/home/vital/Starlight/Bases/coeffs_sync.txt' 
# default_db_folder   = '/home/vital/Astrodata/Inferece_output/'
# db_name_ext         = 'SSPsynth_v1'
# iterat, burn, thin  = 15000, 1500, 2 #15000, 1500, 2
# synth_model         = 'initial_tests'
# params_list         = ['z_star','sigma_star','Av_star', 'ChiSq']
# true_values         = {'z_star':0.0,'sigma_star':1.2,'Av_star':0.7}
# 
# #Load stellar libraries
# sspLib_starlight = bm.load_stellar_bases('starlight', ssp_library_data['default_Starlight_folder'], ssp_library_data['default_Starlight_file'], resample_int=1, resample_range = (4000, 6900), norm_interval = (5100,5150))
# 
# #Generate synthetic observation using default values
# synth_SED = bm.calculate_synthStellarSED(true_values['Av_star'], true_values['z_star'], true_values['sigma_star'], ssp_library_data['default_Starlight_coeffs'], sspLib_starlight, (4000, 6900))
# 
# #Prepare data for fit
# fit_data = bm.ready_data_sspFit(synth_SED, sspLib_starlight)
# 
# # #Ready data for fit
# output_data = bm.ssp_fit(true_values['z_star'], true_values['sigma_star'], true_values['Av_star'], fit_data)
# 
# fig, axis = plt.subplots(1, 1)
# axis.plot(synth_SED['obs_wave_resam'], synth_SED['obs_flux_resam'], label='Synth spectrum')
# axis.plot(synth_SED['obs_wave_resam'], synth_SED['obs_flux_True'], label='True synth spectrum')
# axis.plot(synth_SED['obs_wave_resam'], synth_SED['obs_fluxEr_resam'], label='Err spectrum')
# axis.plot(fit_data['obs_wave_resam'], output_data['flux_sspFit'] , label='Output fit', linestyle=':', color='tab:red')
# axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
# axis.legend()
# plt.show() 
