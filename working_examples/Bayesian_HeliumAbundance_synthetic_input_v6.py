from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v19 import Run_MCMC
import matplotlib.pyplot as plt

#MCMC configuration run
default_db_folder   = '/home/vital/Astrodata/Inferece_output/'
db_name_ext         = 'He_S_O_NebStellar_v1'
iterat, burn, thin  = 10000, 1500, 1 #15000, 1500, 2
synth_model         = 'complete_spec'
params_list         = ['y_plus','T_low','ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund', 'z_star','sigma_star','Av_star']             
                        #'ChiSq_Recomb','ChiSq_S','ChiSq_O',
#Generate dazer object
dz = Dazer()
bm = Run_MCMC()
  
#Generate the synthetic data
bm.calculate_synthFluxes(synth_model)

# fig, axis = plt.subplots(1, 1)
# axis.plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], label='obs_flux')
# axis.plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm_masked'], label='obs_flux_masked')
# # axis.plot(bm.synth_SED['obs_wave_resam'], bm.synth_SED['obs_flux_norm'] * bm.synth_SED['int_mask'], label='obs_flux')
# #axis.plot(bm.synth_SED['obs_wave_rest'], bm.synth_SED['neb_flux_norm'], label='neb_flux')
# axis.plot(bm.synth_SED['obs_wave_resam'], bm.synth_SED['neb_int_norm'], label='neb_int_norm')
# axis.plot(bm.synth_SED['obs_wave_resam'], bm.synth_SED['neb_flux_norm'], label='neb_flux_norm')
# axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
# axis.legend()
# plt.show() 

#Select right model according to data
bm.select_inference_model(synth_model)     
             
#Variables to save
#db_address  = '{}HeAbund_pymc2_{}_it{}_burn{}_thin{}_{}'.format(default_db_folder, synth_model, iterat, burn, thin, db_name_ext) 
db_address  = '{}{}_it{}_burn{}_thin{}_{}'.format(default_db_folder, synth_model, iterat, burn, thin, db_name_ext) 

# #Run sampler
bm.run_pymc2(db_address, iterat, burn, thin, variables_list=params_list)
  
#Load database
pymc2_db, stat_db_dict = bm.load_pymc_database(db_address)
  
#Traces plot
dz.traces_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_tracesPlot', save_pickle = False)
   
#Posteriors plot
dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_posteriorPlot', save_pickle = False)
   
#Posteriors plot
dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=4, n_rows=3)
dz.save_manager(db_address + '_acorrPlot', save_pickle = False)
  
#Corner plot
dz.corner_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_cornerPlot', save_pickle = False)
  
print '\nData treated'

