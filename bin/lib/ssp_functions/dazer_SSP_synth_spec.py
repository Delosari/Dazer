from ssp_synthesis_tools import ssp_fitter
from timeit import default_timer as timer
1.2import numpy as np

def ndprint(a, format_string ='{0:.3f}'):
    return [format_string.format(v,i) for i,v in enumerate(a)]

dz = ssp_fitter()
 
#Run data
config_file_FIT3D           = 'auto_ssp_V500_several_Hb.config'
stellar_library_FIT3D       = 'ssp_lib.fits'
STARLIGHT_stellar_file      = '/home/vital/Starlight/Dani_Bases_Extra.txt'
STARLIGHT_stellar_folder    = '/home/vital/Starlight/Bases/'
data_folder                 = '/home/vital/workspace/Fit_3D/example_files/'
synth_coeffs_file           = '/home/vital/Starlight/Bases/coeffs_sync.txt'
 
#Observational data
obs_Fit3D = dz.load_FIT3D_data(data_folder + config_file_FIT3D, data_folder) 
  
#Preload stellar bases fit from FIT3D
sspLib_dict = dz.load_stellar_bases('FIT3D', data_folder, obs_Fit3D['SSPs_lib'])
  
#Preparing data for the fit
ssp_fit_fit3d = dz.ready_data_sspFit(obs_Fit3D, sspLib_dict, resample_int = None, norm_interval = None)
 
#Preload stellar bases fit from starlight
sspLib_starlight = dz.load_stellar_bases('starlight', STARLIGHT_stellar_folder, STARLIGHT_stellar_file)
 
#Generate Synthetic spectrum
ssp_indeces, coeffs_synth = np.loadtxt(synth_coeffs_file, usecols=[0,1], unpack=True)
synth_wave, synth_flux = dz.generate_synthObs(sspLib_starlight['basesWave'], sspLib_starlight['fluxBases'], coeffs_synth, obs_Fit3D['input_Av'], 0.0, obs_Fit3D['input_sigma'], resample_range = (4000, 6900))
 
sync_dict = {}
sync_dict['obs_wave']         = synth_wave
sync_dict['obs_flux']         = synth_flux
sync_dict['obs_flux_err']     = np.zeros(len(synth_flux))
sync_dict['obs_fluxErrAdj']   = np.zeros(len(synth_flux))
sync_dict['nObsPix']          = len(synth_flux)
sync_dict['input_Av']         = obs_Fit3D['input_Av']
sync_dict['input_sigma']      = obs_Fit3D['input_sigma']
sync_dict['input_z']          = 0.0
sync_dict['data_type']        = 'starlight'
 
#Resample the starlight database
starlight_sspfit_data = dz.ready_data_sspFit(sync_dict, sspLib_starlight, resample_int=1, resample_range = (4000, 6900), norm_interval = (5100,5150))
 
#Reproducing FIT3D example
output_dict = dz.run_ssp_fit(0.0, obs_Fit3D['input_sigma'], obs_Fit3D['input_Av'], starlight_sspfit_data, fit_scheme = 'nnls')
 
#Show me the bases:
fig, axis = plt.subplots(1, 1)
axis.plot(synth_wave, synth_flux, label='Sync spectrum')
axis.plot(starlight_sspfit_data['obs_wave_resam'], starlight_sspfit_data['obs_flux_masked'] * starlight_sspfit_data['normFlux_obs'] , label='Resampling spectrum', linestyle=':')
axis.plot(starlight_sspfit_data['bases_wave_resam'], output_dict['flux_sspFit'], label='SSP synthesis product', linestyle='--')
axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
axis.legend()
plt.show()

 
 
# #Show me the bases:
# fig, axis = plt.subplots(1, 1)
# # axis.plot(synth_wave, synth_flux, label='Sync spectrum')
# axis.plot(starlight_sspfit_data['bases_wave_resam'], output_dict['flux_sspFit'], label='SSP synthesis product', linestyle='--')
# axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
# axis.legend()
# plt.show()

# #Show me the bases:
# fig, axis = plt.subplots(1, 1)
# axis.plot(obs_wave, synthFlux, label='Sync spectrum')
# axis.plot(obs_wave, synthFlux, label='Output fit spectrum')
# starlight_sspfit_data = dz.ready_data_sspFit(obs_Fit3D, sspLib_starlight, resample_int=1, resample_range = (4000, 6900), norm_interval = (5100,5150))
#
# plus_coeffs = coeffs_synth[coeffs_synth > 0]
# idx_true    = np.where(coeffs_synth > 0)[0]
# 
# for i in range(len(plus_coeffs)):
#     axis.plot(obs_wave, plus_coeffs[i] * bases_flux[idx_true[i]], label='Population {}, coeff {}'.format(idx_true[i], plus_coeffs[i]))
# axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
# axis.legend()
# plt.show()

# from ssp_synthesis_tools import ssp_fitter
# from timeit import default_timer as timer
# import matplotlib.pyplot as plt
# import numpy as np
# 
# def ndprint(a, format_string ='{0:.3f}'):
#     return [format_string.format(v,i) for i,v in enumerate(a)]
# 
dz = ssp_fitter()
 
#Run data
config_file_FIT3D           = 'auto_ssp_V500_several_Hb.config'
stellar_library_FIT3D       = 'ssp_lib.fits'
STARLIGHT_stellar_file      = '/home/vital/Starlight/Dani_Bases_Extra.txt'
STARLIGHT_stellar_folder    = '/home/vital/Starlight/Bases/'
data_folder                 = '/home/vital/workspace/Fit_3D/example_files/'
 
#Observational data
obs_Fit3D = dz.load_FIT3D_data(data_folder + config_file_FIT3D, data_folder) 
 
#Preload stellar bases fit from FIT3D
sspLib_dict = dz.load_stellar_bases('FIT3D', data_folder, obs_Fit3D['SSPs_lib'])
 
#Preparing data for the fit
ssp_fit_fit3d = dz.ready_data_sspFit(obs_Fit3D, sspLib_dict, resample_int = None, norm_interval = None)

 
# #Show me the bases:
# fig, axis = plt.subplots(1, 1)
# for i in [0, 10, 20, 40, 50, 60, 70, 80]:
#     axis.plot(ssp_fit_fit3d['bases_wave_resam'], ssp_fit_fit3d['bases_flux_norm'][i], label='SSP n{}'.format(i))
# axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
# axis.legend()
# plt.show()


#Reproducing FIT3D example
output_dict = dz.run_ssp_fit(obs_Fit3D['input_z'], obs_Fit3D['input_sigma'], obs_Fit3D['input_Av'], ssp_fit_fit3d, fit_scheme = 'nnls')
  
#Plot output data
fig, axis = plt.subplots(1, 1)
axis.plot(obs_Fit3D['obs_wave'], obs_Fit3D['obs_flux'], label='input spectrum')
axis.plot(obs_Fit3D['obs_wave'], output_dict['flux_sspFit'], label='SSP synthesis product', linestyle='--')
axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
axis.legend()
plt.show()

# #Show me the bases:
# fig, axis = plt.subplots(1, 1)
# for i in [0, 10, 20, 40, 50, 60, 70, 80]:
#     axis.plot(sspLib_dict['basesWave'], sspLib_dict['fluxBases'][i], label='SSP n{}'.format(i))
# axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
# axis.legend()
# plt.show()





# #Read observational data
# obs_data = np.loadtxt(data_folder + obs_Fit3D['input_spec'])
# obs_wave, obs_flux, obs_fluxVar = obs_data[:,1], obs_data[:,2], obs_data[:,3]
# obs_fluxErr = np.sqrt(abs(obs_fluxVar))

# # # #---Reading data from STARLIGHT configuration
# # #Preload stellar bases fit
# # config_starlight = {}
# # config_starlight['obs_wave'] = np.linspace(4200, 6500, 6500 - 4200 + 1)
# # config_starlight['norm_basesWave'] = 5200.0  F
# # sspLib_dict_PhD = dz.load_stellar_bases(stellar_library2, 'starlight_mode', config_starlight)
# 
# #Generate a synthetic spectrum
# coefficients_address = '/home/vital/workspace/Fit_3D/example_files/coeffs_2.txt'
# coeff_theo = np.loadtxt(coefficients_address)
# synth_obs = dz.generate_synthObs(obs_wave, coeff_theo, sspLib_dict['basesWave'], sspLib_dict['fluxBases'], obs_Fit3D['input_Av'], obs_Fit3D['input_z'], obs_Fit3D['input_sigma'])
# synth_err = obs_fluxErr
#   
# #Input data
# #dz.preload_fitData(obs_wave, synth_obs, obs_Fit3D, sspLib_dict, obs_fluxErr, data_folder)
#   
# #Run the fit
# scheme_fit = 'nnls'
# fitting_data_dict = dz.fit_ssp(dz.sspFit_dict['input_z'], dz.sspFit_dict['input_sigma'], dz.sspFit_dict['input_Av'], scheme_fit)
#   
# coeff_output = fitting_data_dict['weight_coeffs']
#   
# # idx_coefs = np.where(coeff_output != 0)
# # coeff_lsq = fitting_data_dict['coeffs_lsq']
# # idx_coefs_lsq = np.where(coeff_lsq > 0.01)
# print 'Initial coefficients', ndprint(coeff_theo[coeff_theo != 0])#, np.where(coeff_theo != 0)
# print '{}  coefficients'.format(scheme_fit), ndprint(coeff_output[coeff_output > 0.01])#, idx_coefs
# # print 'leastSq coefficients', ndprint(coeff_lsq[idx_coefs_lsq])#, idx_coefs_lsq
#   
# #Plot output data
# fig, axis = plt.subplots(1, 1)
# axis.plot(dz.sspFit_dict['obs_wave'], synth_obs, label='Synthetic spectrum', linestyle='--')
# # for idx in idx_coefs[0]:
# #     label = 'Component {}'.format(round(coeff_output[idx]),4)
# #     axis.plot(dz.sspFit_dict['obs_wave'], fitting_data_dict['flux_components'][:,idx], label=label)
# axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['flux_sspFit'], label='SSP synthesis product', linestyle='--')
#    
#    
# # axis.plot(dz.sspFit_dict['obs_wave'], obs_flux, label='Original spectrum')
# # axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['obs_fluxMasked'], label='Observed maskd spectrum')
# # axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['fluxMasked_sspFit'], label='SSP synthesis product mask')
# # axis.plot(dz.sspFit_dict['obs_wave'], dz.sspFit_dict['zero_mask'], label='Mask')
# axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
# axis.legend()
# plt.show()

