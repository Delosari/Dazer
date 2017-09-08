from ssp_synthesis_tools import ssp_fitter
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import numpy as np
import collections

#USEFULL_Function
def ndprint(a, format_string ='{0:.3f}'):
    return [format_string.format(v,i) for i,v in enumerate(a)]

#Create ssp functions object
dz = ssp_fitter()

#Run data
config_file_FIT3D           = 'auto_ssp_V500_several_Hb.config'
stellar_library_FIT3D       = 'ssp_lib.fits'
stellar_library_STARLIGHT   = '/home/vital/Starlight/Dani_Bases_Extra.txt'
data_folder                 = '/home/vital/workspace/Fit_3D/example_files/'

#Observational data
obs_Fit3D = dz.load_FIT3D_data(data_folder + config_file_FIT3D, data_folder) 

#Stellar library fits
ssp_lib_folder  = obs_Fit3D['data_folder']
ssp_lib_name    = obs_Fit3D['SSPs_lib']
sspLib_Fit3D    = dz.load_stellar_bases('FIT3D_example', ssp_lib_folder, ssp_lib_name)

#Stellar library fits
ssp_lib_folder  = '/home/vital/Starlight/Bases/'
ssp_lib_name    = '/home/vital/Starlight/Dani_Bases_Extra.txt'
sspLib_StarL    = dz.load_stellar_bases('starlight', ssp_lib_folder, ssp_lib_name)

#Resample and normalized observations and stellar bases 
#and prepare the masks


#Run the fit


#Plot the results

# #Read observational data
# obs_data = np.loadtxt(data_folder + obs_Fit3D['input_spec'])
# obs_wave, obs_flux, obs_fluxVar = obs_data[:,1], obs_data[:,2], obs_data[:,3]
# obs_fluxErr = np.sqrt(abs(obs_fluxVar))
# 
# #Preload stellar bases fit
# sspLib_dict = dz.load_stellar_bases(data_folder + stellar_library_FIT3D, 'FIT3D_example', obs_Fit3D)
# 
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

