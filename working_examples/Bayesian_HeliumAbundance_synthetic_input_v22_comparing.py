import numpy as np
from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v33_notroll import Run_MCMC

iterat, burn, thin  = 15000, 0, 1
sim_model           = 'ReDoing_WithoutAbs'
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
 
#Variables to save
db_address1  = '{}{}_it{}_burn{}'.format(bm.paths_dict['inference_folder'], sim_name, iterat, burn) 
db_address2  = '/home/vital/Astrodata/Inference_output/ReDoing_WithAbs_He_S_O_neb_stellar_it15000_burn0'

#Load database
pymc2_db1, stat_db_dict1 = bm.load_pymc_database_manual(db_address1, burning, params_list)
pymc2_db2, stat_db_dict2 = bm.load_pymc_database_manual(db_address2, burning, params_list)
               
print 'Without abs',np.mean(stat_db_dict1['calc_recomb_fluxes']['trace'], axis=0)
print 'With abs',   np.mean(stat_db_dict2['calc_recomb_fluxes']['trace'], axis=0)
print 'difference', np.mean(stat_db_dict2['calc_recomb_fluxes']['trace'], axis=0)/np.mean(stat_db_dict1['calc_recomb_fluxes']['trace'], axis=0) - 1

print 'Obs Without abs',np.mean(stat_db_dict1['calc_obs_emission']['trace'], axis=0)
print 'obs With abs',   np.mean(stat_db_dict2['calc_obs_emission']['trace'], axis=0)
print 'obs difference', np.mean(stat_db_dict2['calc_obs_emission']['trace'], axis=0)/np.mean(stat_db_dict1['calc_obs_emission']['trace'], axis=0) - 1


print 'With abs st',   np.std(stat_db_dict2['calc_recomb_fluxes']['trace'], axis=0)


dz.FigConf()
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm_masked'], label = 'obs_flux_norm_masked')
dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'] - np.mean(stat_db_dict1['stellar_continua_calculation']['trace'], axis=0), label = 'stellar without')
dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'] - np.mean(stat_db_dict2['stellar_continua_calculation']['trace'], axis=0), label = 'stellar with')

dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
dz.display_fig()

# #Traces plot
# print '-Generating traces plot'
# dz.traces_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_tracesPlot_Test', save_pickle = False)
#             
# #Posteriors plot
# print '-Generating posteriors plot'
# dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
# dz.save_manager(db_address + '_posteriorPlot', save_pickle = False)
#              
# #Posteriors plot
# print '-Generating acorrelation plot'
# dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=4, n_rows=4)
# dz.save_manager(db_address + '_acorrPlot', save_pickle = False)
#            
# #Corner plot
# print '-Generating corner plot'
# dz.corner_plot(params_list, pymc2_db, stat_db_dict, plot_true_values=True)
# dz.save_manager(db_address + '_cornerPlot', save_pickle = False)

print '\nData treated'