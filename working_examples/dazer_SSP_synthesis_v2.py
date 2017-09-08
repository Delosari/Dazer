import numpy as np
from dazer_methods import Dazer
from lib.ssp_functions.ssp_synthesis_tools import ssp_fitter
from timeit import default_timer as timer

dzp = Dazer()
dz = ssp_fitter()

#Data folder location
data_folder         = '/home/vital/workspace/Fit_3D/example_files/'
defaut_conf         = 'auto_ssp_V500_several_Hb.config'

#Read parameters from command line
command_dict        = dz.load_command_params(data_folder)

#Read parameters from config file
conf_file_address   = command_dict['config_file_address'] if 'config_file_address' in command_dict else data_folder + defaut_conf 
config_dict         = dz.load_config_params(conf_file_address)

#Update the fit configuration giving preference to the values from the command line
config_dict.update(command_dict)

#Import input data: spectrum, masks, emision line loc, stellar bases...
dz.load_input_data(config_dict)

#Perform SSP synthesis
start = timer()   
fit_products = dz.fit_ssp(config_dict['input_z'], config_dict['input_sigma'], config_dict['input_Av'])
end = timer()
print 'ssp', ' time ', (end - start)

np.savetxt(data_folder + 'coeffs_orig.txt', fit_products['weight_coeffs'].T)

#Plot the results
dzp.FigConf()
dzp.data_plot(dz.sspFit_dict['obs_wave'], dz.sspFit_dict['obs_flux'],   label='obs_flux')
dzp.data_plot(dz.sspFit_dict['obs_wave'], dz.sspFit_dict['zero_mask'],  label='my mask')
dzp.data_plot(dz.sspFit_dict['obs_wave'], fit_products['flux_sspFit'],  label='my fit')

dzp.FigWording('Wave', 'Flux', 'Input spectra')
dzp.display_fig()

