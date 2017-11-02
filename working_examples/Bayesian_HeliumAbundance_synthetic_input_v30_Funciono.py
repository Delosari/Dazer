import numpy as np
from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v42 import Run_MCMC
from matplotlib import rcParams
import matplotlib.pyplot as plt

iterat, burn, thin  = 15000, 0, 1
sim_model           = 'Repetimos_por_si_acaso'
sim_components      = '_neb_stars_Abunds'
obs_metals          = ['H1', 'He1', 'S2', 'S3', 'O2', 'O3', 'N2', 'Ar3', 'Ar4']
sim_name            = sim_model + sim_components
params_list         = ['He1_abund', 'T_He', 'T_low', 'ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund', 'N2_abund', 'Ar3_abund', 'Ar4_abund', 'sigma_star', 'Av_star'] 
burning             = 8000
                           
#Generate dazer object
dz = Dazer()
bm = Run_MCMC()

#Generate the synthetic data
bm.calculate_simObservation(sim_components, obs_lines = obs_metals)

# lines_fluxes_obs = bm.measure_lines(bm.Recomb_labels, bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], calibration_factor = (bm.obj_data['normFlux_obs'] / bm.obj_data['Hbeta_Flux']))
# lines_fluxes_emis = bm.measure_lines(bm.Recomb_labels, bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'] - bm.stellar_SED['stellar_flux_norm'], calibration_factor = (bm.obj_data['normFlux_obs'] / bm.obj_data['Hbeta_Flux']))

# print 'obs_fluxes',     lines_fluxes_obs
# print 'Emis_fluxes',    lines_fluxes_emis
# print 'Theo_fluxes',    bm.recomb_fluxes

# dz.FigConf()
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], label = 'obs_flux_norm')
# intg_flux = np.empty(len(bm.Recomb_labels))
# for i in range(len(bm.Recomb_labels)):
#     bm.Current_Label      = bm.Recomb_labels[i]
#     bm.Current_Ion        = bm.obj_data['recombLine_ions'][i]
#     bm.Current_TheoLoc    = bm.obj_data['recombLine_waves'][i]
#     selections            = bm.obj_data['lick_idcs_df'].loc[bm.Current_Label][3:9].values
#     line_data             = bm.measure_line(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], selections, None, Measuring_Method = 'lmfit', store_data = False)
#     intg_flux[i]          = line_data['flux_intg']
#     dz.data_plot(line_data.x_resample, line_data.y_resample, label = 'obs_flux_norm')
# dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
# dz.display_fig()

#Select right model according to data
# bm.select_inference_model(sim_components)

# #Variables to save
db_address = '{}{}_it{}_burn{}'.format(bm.paths_dict['inference_folder'], sim_name, iterat, burn) 
    
    # Run sampler
# bm.run_pymc2(db_address, iterat, variables_list=params_list, prefit=False)
            
# #Load database
pymc2_db, stat_db_dict = bm.load_pymc_database_manual(db_address, burning, params_list)

# dz.FigConf()
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'] * bm.obj_data['int_mask'], label = 'mask')
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], label = 'obs_flux_norm')
# dz.data_plot(bm.obj_data['obs_wave_resam'], np.mean(stat_db_dict['calc_continuum']['trace'], axis=0), label = 'calc_continuum')
# dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
# dz.display_fig()


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
#Corner plot

#Default sizes for computer
sizing_dict = {}
# sizing_dict['figure.figsize'] = (14, 8)
# sizing_dict['legend.fontsize'] = 15
sizing_dict['axes.labelsize'] = 50
sizing_dict['axes.titlesize'] = 50
sizing_dict['xtick.labelsize'] = 50
sizing_dict['ytick.labelsize'] = 14
rcParams.update(sizing_dict)

print '-Generating corner plot'
dz.corner_plot(params_list, pymc2_db, stat_db_dict, plot_true_values=True)
plt.savefig('E:/Cloud Storage/Dropbox/Astrophysics/Seminars/Stasinska conference/big_corner.png', dpi=100.0, bbox_inches='tight')

#                
# print '\nData treated'