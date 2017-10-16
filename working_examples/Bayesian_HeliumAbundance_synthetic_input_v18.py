import numpy as np
from scipy.integrate import simps
from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v31 import Run_MCMC

iterat, burn, thin  = 10000, 0, 1
sim_model           = 'With_Absorptions'
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

# #Measure the line
# dz.FigConf()
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm_masked'], label = 'obs_flux_norm_masked')
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], label = 'emission_Spec')
# # dz.data_plot(bm.obj_data['obs_wave_resam'], bm.emission_Spec + bm.obj_data['obs_flux_norm'], label = 'Emission + stellar')
# dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
# dz.display_fig()

# dz.FigConf()
# emis_spectrum = bm.emission_Spec
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.emission_Spec, label = 'emission_Spec')
# obs_fluxes =  np.empty(bm.n_recombLines)
# for i in range(len(bm.obj_data['recombLine_labes'])):
#     line_label = bm.obj_data['recombLine_labes'][i]
#     bm.Current_Label    = bm.obj_data['lick_idcs_df'].loc[line_label].name
#     bm.Current_Ion      = bm.obj_data['lick_idcs_df'].loc[line_label].Ion
#     bm.Current_TheoLoc  = bm.obj_data['lick_idcs_df'].loc[line_label].lambda_theo * (1 + bm.obj_data['z_star'])
#     selections          = bm.obj_data['lick_idcs_df'].loc[line_label][3:9].values
#     line_data           = bm.measure_line(bm.obj_data['obs_wave_resam'], emis_spectrum, selections, None, Measuring_Method = 'lmfit', store_data = False)
#     obs_fluxes[i]       = line_data['flux_intg']
#     dz.data_plot(line_data.x_resample, line_data.y_resample, label = line_label, linestyle='--')
# print 'obs_fluxes', obs_fluxes * (bm.obj_data['normFlux_obs'] /bm.obj_data['Hbeta_Flux'])
# print 'recomb_fluxes', bm.recomb_fluxes
#   
# dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
# dz.display_fig()

#Select right model according to data
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