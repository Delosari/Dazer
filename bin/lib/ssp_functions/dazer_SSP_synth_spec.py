from ssp_synthesis_tools import ssp_fitter
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import numpy as np

dz = ssp_fitter()

#Read parameters from command line
command_dict        = dz.load_command_params()

#Read parameters from config file
conf_file_address   = 'auto_ssp_V500_several_Hb.config'
config_dict         = dz.load_config_params(conf_file_address)

#Update the fit configuration giving preference to the values from the command line
config_dict.update(command_dict)

#Preload stellar bases
dz.load_stellar_bases(dz.ssp_folder, 'FIT3D_example', config_dict)

#We use prepare synthetic observation to compare the quality of the fit
coeff_theo                      = np.loadtxt(dz.ssp_folder + 'bases_coeff.txt')
spectrum_matrix                 = np.loadtxt(dz.ssp_folder + config_dict['input_spec'])
obs_wave, obs_flux, obs_fluxVar = spectrum_matrix[:,1], spectrum_matrix[:,2], spectrum_matrix[:,3]
obs_fluxErr                     = np.sqrt(abs(obs_fluxVar))

print 'Input coefficients', coeff_theo[coeff_theo != 0], np.where(coeff_theo != 0)

synth_obs           = dz.generate_synthObs(obs_wave, coeff_theo, dz.sspLib_dict['basesWave'], dz.sspLib_dict['fluxBases'], 
                     config_dict['input_Av'], 
                     config_dict['input_z'], 
                     config_dict['input_sigma'])

#Import input data: spectrum, masks, emision line loc, stellar bases...
dz.preload_fitData(obs_wave, obs_flux, config_dict, obs_fluxErr)
fitting_data_dict = dz.fit_ssp(dz.sspFit_dict['input_z'], dz.sspFit_dict['input_sigma'], dz.sspFit_dict['input_Av'])
coeff_output = fitting_data_dict['weight_coeffs']
idx_coefs = np.where(coeff_output != 0)
print 'Output New coefficients', coeff_output[coeff_output != 0], idx_coefs

#Plot output data
fig, axis = plt.subplots(1, 1)

axis.plot(dz.sspFit_dict['obs_wave'], synth_obs, label='Synthetic spectrum')

#Spectra components
for idx in idx_coefs[0]:
    label = 'Component {}'.format(round(coeff_output[idx]),4)
    axis.plot(dz.sspFit_dict['obs_wave'], fitting_data_dict['flux_components'][:,idx], label=label)
axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['flux_sspFit'], label='SSP synthesis product')
 
axis.plot(dz.sspFit_dict['obs_wave'], obs_flux, label='Original spectrum')
# axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['obs_fluxMasked'], label='Observed maskd spectrum')
# axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['fluxMasked_sspFit'], label='SSP synthesis product mask')
# axis.plot(dz.sspFit_dict['obs_wave'], dz.sspFit_dict['zero_mask'], label='Mask')
axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
axis.legend()
plt.show()

# np.savetxt(dz.ssp_folder + 'bases_coeff.txt', fitting_data_dict['weight_coeffs'])

