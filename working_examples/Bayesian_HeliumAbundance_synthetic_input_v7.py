from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v20 import Run_MCMC
from numpy import mean
# #MCMC configuration run
# default_db_folder   = '/home/vital/Astrodata/Inference_output/'
# iterat, burn, thin  = 10000, 1500, 2
# sim_model           = 'checking'
# sim_components      = '_He_S_O_neb_stellar_v1'
# sim_name            = sim_model + sim_components
# params_list         = ['y_plus', 'T_low','ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund','sigma_star','Av_star'] #'ChiSq_Recomb','ChiSq_S','ChiSq_O']      

default_db_folder   = '/home/vital/Astrodata/Inference_output/'
iterat, burn, thin  = 20000, 0, 1
sim_model           = 'NoPrefit_NoBurnThinning'
sim_components      = '_He_S_O_neb_stellar_v1'
sim_name            = sim_model + sim_components
params_list         = ['y_plus', 'T_low','ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund','sigma_star','Av_star'] #'ChiSq_Recomb','ChiSq_S','ChiSq_O']      
burning             = 5000
              
#Generate dazer object
dz = Dazer()
bm = Run_MCMC()

#Generate the synthetic data
bm.calculate_synthFluxes(sim_components)

# dz.FigConf()
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], label = 'obs_flux_norm')
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.nebular_SED['neb_flux_norm'], label = 'neb_flux_norm')
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.stellar_SED['stellar_flux_norm'], label = 'stellar_flux_norm')
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.stellar_SED['stellar_flux_norm'] + bm.nebular_SED['neb_flux_norm'], label = 'sumano ')
# dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
# dz.display_fig()

# #Select right model according to data
bm.select_inference_model(sim_components)     
#              
# #Variables to save
db_address  = '{}{}_it{}_burn{}_thin{}'.format(default_db_folder, sim_name, iterat, burn, thin) 
  
# #Run sampler
# bm.run_pymc2(db_address, iterat, burn, thin, variables_list=params_list)
    
#Load database
pymc2_db, stat_db_dict = bm.load_pymc_database_manual(db_address, burning=burning)

dz.FigConf()
dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], label = 'obs_flux_norm')
dz.data_plot(bm.obj_data['obs_wave_resam'], mean(stat_db_dict['ssp_Synthesis']['trace'],axis=0), label = 'fitted spectrum')
dz.data_plot(bm.obj_data['obs_wave_resam'], bm.nebular_SED['neb_flux_norm'], label = 'neb_flux_norm')
dz.data_plot(bm.obj_data['obs_wave_resam'], bm.stellar_SED['stellar_flux_norm'], label = 'stellar_flux_norm')
dz.FigWording(xlabel = 'Wavelength', ylabel = 'Normalized flux', title = '')
dz.display_fig()


# #Traces plot
# dz.traces_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_tracesPlot_Burn', save_pickle = False)
#      
# #Posteriors plot
# dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_posteriorPlot_Burn', save_pickle = False)
#      
# #Posteriors plot
# dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=4, n_rows=3)
# dz.save_manager(db_address + '_acorrPlot_Burn', save_pickle = False)
#     
# #Corner plot
# dz.corner_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_cornerPlot_Burn', save_pickle = False)
   
print '\nData treated'

