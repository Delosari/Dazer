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

#Import input data: spectrum, masks, emision line loc, stellar bases...
dz.load_input_data(config_dict)

#Funcom
coeff_theory = np.loadtxt(dz.ssp_folder + 'bases_coeff.txt')
print 'Initial coefficients', coeff_theory[coeff_theory != 0]

#Calculating synthetic spectrum

synth_obs = np.sum(coeff_theory * dz.sspFit_dict['fluxBases'].T, axis=1)

# #Perform SSP synthesis
# start = timer()   
# fitting_data_dict = dz.fit_ssp(dz.sspFit_dict['input_z'], dz.sspFit_dict['input_sigma'], dz.sspFit_dict['input_Av'])
# end = timer()
# print 'ssp', ' time ', (end - start)
# coeff_output = fitting_data_dict['weight_coeffs']
# print 'Output coefficients', coeff_output[coeff_output != 0]



#Plot output data
fig, axis = plt.subplots(1, 1)
axis.plot(dz.sspFit_dict['obs_wave'], dz.sspFit_dict['obs_flux'], label='Original spectrum')
axis.plot(dz.sspFit_dict['basesWave_zCor'], synth_obs, label='Synthetic spectrum')

# axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['flux_sspFit'], label='SSP synthesis product')
# axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['obs_fluxMasked'], label='Observed maskd spectrum')
# axis.plot(fitting_data_dict['obs_wave'], fitting_data_dict['fluxMasked_sspFit'], label='SSP synthesis product mask')
# axis.plot(dz.sspFit_dict['obs_wave'], dz.sspFit_dict['zero_mask'], label='Mask')
axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
axis.legend()
plt.show()


# np.savetxt(dz.ssp_folder + 'bases_coeff.txt', fitting_data_dict['weight_coeffs'])

