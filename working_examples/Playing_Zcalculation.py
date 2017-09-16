from dazer_methods import Dazer
import matplotlib.pyplot as plt
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v18 import Run_MCMC
from dask.dataframe.tseries.tests.test_resample import resample
from scipy.interpolate.interpolate import interp1d

#Declare operation objects
dz  = Dazer()
bm  = Run_MCMC()

#MCMC configuration run
ssp_library_data = {}
ssp_library_data['default_Starlight_file']     = '/home/vital/Starlight/Dani_Bases_Extra.txt'
ssp_library_data['default_Starlight_folder']   = '/home/vital/Starlight/Bases/'
ssp_library_data['default_Starlight_coeffs']   = '/home/vital/Starlight/Bases/coeffs_sync.txt' 
default_db_folder   = '/home/vital/Astrodata/Inferece_output/'
db_name_ext         = 'SSPsynth_v1'
iterat, burn, thin  = 15000, 1500, 2 #15000, 1500, 2
synth_model         = 'initial_tests'
params_list         = ['z_star','sigma_star','Av_star', 'ChiSq']
true_values         = {'z_star':0.01,'sigma_star':1.2,'Av_star':0.7}

#Generate synthetic observation using default values
synth_SED = bm.calculate_synthStellarSED(true_values['Av_star'], true_values['z_star'], true_values['sigma_star'], ssp_library_data)
synth_SED2 = bm.calculate_synthStellarSED(true_values['Av_star'], 0.0, true_values['sigma_star'], ssp_library_data, resample_range=(3550, 6990))
z_guess = 0.015
wave_z_guess = synth_SED2['obs_wave'] * (1 + z_guess)
inter_spec = (interp1d(wave_z_guess, synth_SED2['obs_flux_true'], bounds_error=True)(synth_SED['obs_wave'])).T
 

#Show me the bases:
fig, axis = plt.subplots(1, 1)
axis.plot(synth_SED['obs_wave'], synth_SED['obs_flux_true'], label='Obs spectrum')
axis.plot(synth_SED2['obs_wave'], synth_SED2['obs_flux_true'], label='Bases spectrum')
axis.plot(synth_SED['obs_wave'], inter_spec, label='Interpolated_spec', linestyle=':')

axis.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'SSP synthesis fitting with 102 bases'})
axis.legend()
plt.show() 



