from bin.dazer_methods import Dazer
from bin.lib.ssp_functions.ssp_synthesis_tools import ssp_fitter

dzp = Dazer()
dz = ssp_fitter()

#Data folder location
data_folder         = '/home/vital/workspace/Fit_3D/example_files/'
defaut_conf         = 'auto_ssp_V500_several_Hb.config'

#Read parameters from command line
command_fit_dict    = dz.load_command_params()

#Read parameters from config file
conf_file_address   = command_fit_dict['config_file_address'] if 'config_file_address' in command_fit_dict != None else data_folder + defaut_conf 
config_fit_dict     = dz.load_config_params(conf_file_address)

#Update the fit configuration giving preference to the values from the command line
config_fit_dict.update(command_fit_dict)

#Import input data: spectrum, masks, emision line loc, stellar bases...
dz.load_input_data(config_fit_dict)

obs_fit_spectrum = dz.fit_ssp()

#Plot input spectra and regions
dzp.FigConf()
dzp.data_plot(dz.fit_conf['obs_wave'], dz.fit_conf['obs_flux'], label='obs_flux')
dzp.data_plot(dz.fit_conf['obs_wave'], dz.fit_conf['obs_flux_mask'], label='obs_flux_mask')
# dzp.data_plot(dz.fit_conf['obs_wave'], dz.fit_conf['Av_mask'], label='Av mask')
# dzp.data_plot(dz.fit_conf['obs_wave'], dz.fit_conf['obs_flux_err'], label='obs_flux_err')
# dz.data_plot(wave_unc, flux_unc_input, label='flux_unc_input')
# dz.data_plot(wave_unc, e_flux_unc, label='e_flux_unc')
# dz.data_plot(wave_unc, color_unc, label='color_unc')
# dz.data_plot(wave_unc, masked, label='masked')
# dz.data_plot(wave_unc, masked2, label='masked2')
# dz.data_plot(wave_unc, masked_Av, label='masked_Av')
# dz.data_plot(wave_unc, flux_masked, label='flux_masked')
# dz.data_plot(wave_unc, flux_masked2, label='flux_masked2')
# dz.data_plot(wave_unc, e_flux_unc_kin, label='e_flux_unc_kin')
dzp.data_plot(dz.fit_conf['obs_wave'], obs_fit_spectrum, label='fit spectrum, si')
dzp.FigWording('Wave', 'Flux', 'Input spectra')
dzp.display_fig()






#Perform unlinear parameters fit

#Perform linear fit

#Plot the data