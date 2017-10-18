from numpy import empty
from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v36 import Run_MCMC

iterat, burn, thin  = 15000, 0, 1
sim_model           = 'ReDoing_WithAbs_theOtherWay'
sim_components      = '_He_S_O_neb_stellar'
obs_metals          = ['H', 'He1', 'S2', 'S3', 'O2', 'O3', 'N2', 'Ar3', 'Ar4']
sim_name            = sim_model + sim_components
params_list         = ['He1_abund', 'T_He', 'T_low', 'ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund', 'N2_abund', 'Ar3_abund', 'Ar4_abund', 'sigma_star', 'Av_star'] 
burning             = 5000
                        
#Generate dazer object
dz = Dazer()
bm = Run_MCMC()

#Generate the synthetic data
bm.calculate_simObservation(sim_components, obs_lines = obs_metals)

# obs_fluxes  = empty(bm.n_recombLines)
# emis_fluxes = empty(bm.n_recombLines)
# 
# dz.FigConf()
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], label = 'obs_flux_norm')
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'] - bm.stellar_SED['stellar_flux_norm'], label = 'emission')
# for i in bm.range_recombLines:
#     line_label = bm.Recomb_labels[i]
#     bm.Current_Label          = bm.Recomb_labels[i]
#     bm.Current_Ion            = bm.obj_data['recombLine_ions'][i]
#     bm.Current_TheoLoc        = bm.obj_data['recombLine_waves'][i]
#     selections                = bm.obj_data['lick_idcs_df'].loc[bm.Recomb_labels[i]][3:9].values
#     line_data_with_abs    = bm.measure_line(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], selections, None, Measuring_Method = 'lmfit', store_data = False)   
#     obs_fluxes[i]         = line_data_with_abs['flux_intg'] * (bm.obj_data['normFlux_obs'] / bm.Hbeta_Flux)
#     print line_label, bm.obs_recomb_fluxes[i]
# #     line_data_emis        = bm.measure_line(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'] - bm.stellar_SED['stellar_flux_norm'], selections, None, Measuring_Method = 'lmfit', store_data = False)
# #     emis_fluxes[i]        = line_data_emis['flux_intg'] * (bm.obj_data['normFlux_obs'] / bm.Hbeta_Flux)
#     
# #     print line_data_with_abs['flux_intg'], line_data_emis['flux_intg']
#     
#     dz.data_plot(line_data_with_abs['x_resample'], line_data_with_abs['y_resample'], label = line_label, linestyle=':')
#     
# # abs_vector =  emis_fluxes - obs_fluxes
# # 
# # print 'obs_fluxes', obs_fluxes
# # print 'emis_fluxes', emis_fluxes
# # 
# # print obs_fluxes
# # print emis_fluxes - abs_vector
# dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
# dz.display_fig()

# #Select right model according to data
bm.select_inference_model(sim_components)     
  
#Variables to save
db_address  = '{}{}_it{}_burn{}'.format(bm.paths_dict['inference_folder'], sim_name, iterat, burn) 
   
#Run sampler
bm.run_pymc2(db_address, iterat, variables_list=params_list, prefit=False)
         
#Load database
pymc2_db, stat_db_dict = bm.load_pymc_database_manual(db_address, burning, params_list)
          
#Traces plot
print '-Generating traces plot'
dz.traces_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_tracesPlot_Test', save_pickle = False)
            
#Posteriors plot
print '-Generating posteriors plot'
dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_posteriorPlot', save_pickle = False)
             
#Posteriors plot
print '-Generating acorrelation plot'
dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=4, n_rows=4)
dz.save_manager(db_address + '_acorrPlot', save_pickle = False)
           
#Corner plot
print '-Generating corner plot'
dz.corner_plot(params_list, pymc2_db, stat_db_dict, plot_true_values=True)
dz.save_manager(db_address + '_cornerPlot', save_pickle = False)
          
print '\nData treated'