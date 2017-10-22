import numpy as np
from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v48 import Run_MCMC
from collections import OrderedDict
  
iterat, burn, thin  = 7500, 0, 1
sim_model           = 'Ssp_object'
sim_components      = '_ssp_semi'
obs_metals          = ['H1', 'He1', 'S2', 'S3', 'O2', 'O3', 'N2', 'Ar3', 'Ar4']
sim_name            = sim_model + sim_components
params_list         = ['sigma_star', 'Av_star']#['He1_abund', 'T_He', 'T_low', 'ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund', 'N2_abund', 'Ar3_abund', 'Ar4_abund', 'sigma_star', 'Av_star'] 
burning             = 2500

#Lines to employed in the object
lines_to_use_dict   = OrderedDict()
lines_to_use_dict['H1']     = ['H1_4102A','H1_4341A','H1_6563A']
lines_to_use_dict['He1']    = ['He1_4026A','He1_4471A','He1_5876A','He1_6678A']
lines_to_use_dict['S2']     = ['S2_6716A', 'S2_6731A']
lines_to_use_dict['S3']     = ['S3_6312A', 'S3_9069A', 'S3_9531A']
lines_to_use_dict['O2']     = ['O2_3726A', 'O2_3729A']
lines_to_use_dict['O3']     = ['O3_4363A', 'O3_4959A', 'O3_5007A']
lines_to_use_dict['N2']     = ['N2_6548A', 'N2_6584A']
lines_to_use_dict['Ar3']    = ['Ar3_7136A']
lines_to_use_dict['Ar4']    = ['Ar4_4740A']
                    
#Generate dazer object
dz = Dazer()
bm = Run_MCMC()

#Load catalogue dataframe
catalogue_dict      = dz.import_catalogue()
catalogue_df        = dz.load_excel_DF('/home/vital/Dropbox/Astrophysics/Data/WHT_observations/WHT_Galaxies_properties.xlsx')

#Load object scientific data file
obj_code            = 'SHOC579'
obj_data            = catalogue_df.loc[obj_code]

populations_coeffs_file = '{}{}/{}_populations_prefit.txt'.format(bm.paths_dict['observations_folder'], obj_code, obj_code)
ssp_coefs_mean, ssp_coefs_std = np.loadtxt(populations_coeffs_file, unpack=True) 

idx_populations = ssp_coefs_mean > 0.0
bm.nBases       = len(np.where(idx_populations)[0]) 
bm.range_bases  = np.arange(bm.nBases)  
bm.mean_coeffs  = ssp_coefs_mean[idx_populations]
bm.population_limitss = np.vstack((bm.mean_coeffs*0, bm.mean_coeffs*2)).T

print bm.mean_coeffs
print bm.population_limitss

# idx_greater = trace_coefficients > 0.0
# self.population_limitss     = vstack((self.stellar_SED['bases_coeff'][where(idx_populations)]*0.95, self.stellar_SED['bases_coeff'][where(idx_populations)]*1.05)).T
# self.nBases                 = len(where(idx_populations)[0]) 
# self.range_bases            = arange(self.nBases)       
# columns_to_delete           = where(~idx_populations)
 



# #Load object spectrum
# fits_file           = obj_data.reduction_fits
# wave, flux, header  = dz.get_spectra_data(fits_file)
# 
# #Load object lines measurements
# lineslog_extension  = '_' + catalogue_dict['Datatype'] + '_linesLog_reduc.txt'
# ouput_folder        = '{}{}/'.format(catalogue_dict['Obj_Folder'], obj_code) 
# lineslog_address    = '{objfolder}{codeName}{lineslog_extension}'.format(objfolder=ouput_folder, codeName=obj_code, lineslog_extension=lineslog_extension)
# lineslog_frame      = dz.load_lineslog_frame(lineslog_address)
# 
# #Load observation into fit model
# bm.load_obs_data(lineslog_frame, obj_data, obj_code, wave, flux, valid_lines=lines_to_use_dict)

# dz.FigConf()
# dz.data_plot(wave, flux, label = 'observed flux ' + obj_code)
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_resam'], label = 'obs_flux_resam')
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_resam'] * bm.obj_data['int_mask'], label = 'mask')
# dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
# dz.display_fig()

# #Select right model according to data
# bm.select_inference_model(sim_components)
 
# Variables to save
# db_address = '{}{}_it{}_burn{}'.format(bm.paths_dict['inference_folder'], sim_name, iterat, burn) 
     
# Run sampler
# bm.run_pymc2(db_address, iterat, variables_list=params_list, prefit=False)
             
# #Load database


# pymc2_db, stat_db_dict = bm.load_pymc_database_manual(db_address, burning, params_list)
# 
# trace_coefficients_mean = np.mean(stat_db_dict['ssp_coefficients']['trace'], axis=0)
# trace_coefficients_std  = np.mean(stat_db_dict['ssp_coefficients']['trace'], axis=0)
# 
# populations_coeffs_file = '{}{}/{}_populations_prefit.txt'.format(bm.paths_dict['observations_folder'], obj_code, obj_code)
# np.savetxt(populations_coeffs_file, np.vstack((trace_coefficients_mean, trace_coefficients_std)).T)
# 
# print 'Bichos', np.vstack((trace_coefficients_mean, trace_coefficients_std)).T


# idx_greater = trace_coefficients > 0.0
# self.population_limitss     = vstack((self.stellar_SED['bases_coeff'][where(idx_populations)]*0.95, self.stellar_SED['bases_coeff'][where(idx_populations)]*1.05)).T
# self.nBases                 = len(where(idx_populations)[0]) 
# self.range_bases            = arange(self.nBases)       
# columns_to_delete           = where(~idx_populations)
# 
# 
# print 'Coefficients', trace_coefficients[idx_greater]
# print np.mean(stat_db_dict['ssp_coefficients']['trace'], axis=0)




# dz.FigConf()
# dz.data_plot(bm.obj_data['obs_wave_resam'], bm.obj_data['obs_flux_norm'], label = 'obs_flux_norm')
# dz.data_plot(bm.obj_data['obs_wave_resam'], np.mean(stat_db_dict['stellar_continua_calculation']['trace'], axis=0), label = 'stellar_continua_calculation')
# dz.data_plot(bm.obj_data['obs_wave_resam'], np.mean(stat_db_dict['calc_nebular_cont']['trace'], axis=0), label = 'calc_nebular_cont')
# dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
# dz.display_fig()

#Traces plot
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
# dz.corner_plot(params_list, pymc2_db, stat_db_dict, plot_true_values=False)
# dz.save_manager(db_address + '_cornerPlot', save_pickle = False)
                 
# print '\nData treated'