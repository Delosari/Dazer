from errno                          import ENOENT
from os                             import remove, getcwd
from sys                            import argv, exit
from collections                    import OrderedDict
from pandas                         import DataFrame, read_csv
from astropy.io.fits                import getdata
from scipy.signal.signaltools       import convolve2d
from scipy.interpolate.interpolate  import interp1d
from scipy.optimize                 import lsq_linear, nnls
from numpy import array, power, searchsorted, loadtxt, empty, sqrt, abs, sum as np_sum, square, isnan, ones, copy, median, arange, exp, zeros, transpose, mean, diag, linalg, dot, multiply, outer, asfortranarray, inf, random, int64
import pyneb as pn
import ntpath
from timeit import default_timer as timer
import matplotlib.pyplot as plt

#Function to delete files
def silent_remove(filename_list):
    for filename in filename_list:   
        try:
            remove(filename)
        except OSError as e:        # this would be "except OSError, e:" before Python 2.6
            if e.errno != ENOENT:   # errno.ENOENT = no such file or directory
                raise               # re-raise exception if a different error occurred

def example_data(data_folder):
    
    arguments_dict                      = OrderedDict()
    arguments_dict['script']            = 'auto_ssp_elines_rnd.py'                          #0
    arguments_dict['input_spec']        = 'NGC5947.spec_5.txt'                              #1
    arguments_dict['SSPs_lib']          = 'ssp_lib.fits,'+'ssp_lib.fits'                    #2
    arguments_dict['output_file']       = 'auto_ssp.NGC5947.cen.only.out'                   #3
    arguments_dict['mask_file']         = 'mask_elines.txt'                                 #4
    arguments_dict['conf_file']         = 'auto_ssp_V500_several_Hb.config'                 #5
    arguments_dict['plot_tag']          = 1                                                 #6
    arguments_dict['min']               = -1                                                #7
    arguments_dict['max']               = 40                                                #8
    arguments_dict['wmin']              = '3850'                                            #9
    arguments_dict['wmax']              = '6800'                                            #10
    arguments_dict['z_elines_mask']     = 'emission_lines.txt'                              #11
    arguments_dict['input_z']           = 0.02                                              #12
    arguments_dict['delta_z']           = 0.001                                             #13
    arguments_dict['min_z']             = 0.015                                             #14
    arguments_dict['max_z']             = 0.025                                             #15
    arguments_dict['input_sigma']       = 2.0                                               #16
    arguments_dict['delta_sigma']       = 0.5                                               #17
    arguments_dict['min_sigma']         = 1                                                 #18
    arguments_dict['max_sigma']         = 9                                                 #19
    arguments_dict['input_Av']          = 0.5                                               #20
    arguments_dict['delta_Av']          = 0.1                                               #21
    arguments_dict['min_Av']            = 0.0                                               #22
    arguments_dict['max_Av']            = 1.6                                               #23
     
    return arguments_dict

def check_missing_flux_values(flux):
    
    #Evaluate the nan array
    nan_idcs    = isnan(flux)
    nan_count   = np_sum(nan_idcs)
    
    #Directly save if not nan
    if nan_count > 0:
        print '--WARNING: missing flux entries'
    
    return

def CCM89_Bal07(Rv, wave):
    
    x   = 1e4 / wave            #Assuming wavelength is in Amstrongs
    ax  = zeros(len(wave))
    bx  = zeros(len(wave))
    
    idcs        = x > 1.1
    y           = (x[idcs] - 1.82)
    
    ax[idcs]    = 1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
    bx[idcs]    = 1.*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
    ax[~idcs]   = 0.574*x[~idcs]**1.61
    bx[~idcs]   = -0.527*x[~idcs]**1.61
            
    Xx  = ax+bx/Rv              #WARNING better to check this definition
    
    return Xx

def path_leaf(path):
    head, tail  = ntpath.split(path)
    file_name   = tail or ntpath.basename(head)
    folder_name = ntpath.dirname(path) + '/'
    return folder_name, file_name

class ssp_synthesis_importer():
    
    def __init__(self):
        
        #------------Configuration of Fit3D 
        self.sspSyn_commands_params = [
                                       'script',            #0 python script name
                                       'input_spec',        #1 input galactic spectrum name
                                       'SSPs_lib',          #2 fits-table to use with python
                                       'output_file',       #3 Reference name for output files
                                       'mask_file',         #4 File with the spectrum region masks
                                       'conf_file',         #5 Configuration file for the masks
                                       'plot_tag',          #6 tag to launch the plotting 
                                       'min',               #7 Min flux for ploting    
                                       'max',               #8 Max flux for ploting
                                       'wmin',              #9 Minimum wavelength for plotting
                                       'wmax',              #10 Maximum wavelength for plotting
                                       'z_elines_mask',     #11 Emission lines file
                                       'input_z',           #12 Input redshift
                                       'delta_z',           #13 Increments for redshift
                                       'min_z',             #14 Minimum redshift
                                       'max_z',             #15 Maximum redshift
                                       'input_sigma',       #16 Input velocity dispersion
                                       'delta_sigma',       #17 Increments for velocity dispersion
                                       'min_sigma',         #18 Minimum velocity dispersion
                                       'max_sigma',         #19 Maximum velocity dispersion                          
                                       'input_Av',          #20 Input reddening
                                       'delta_Av',          #21 Increments for reddening
                                       'min_Av',            #22 Minimum reddening
                                       'max_Av',            #23 Maximum reddening                                      
                                       ] 
        
        #The first 4 lines in the configuration file describe the input
        self.sspSyn_config_params   = [['input_z','delta_z','min_z','max_z','DV','RV','DS','RS','MIN_W','MAX_W'],                       #12-16
                                       ['input_sigma','delta_sigma','min_sigma','max_sigma'],                                           #17-20                          
                                       ['input_Av','delta_Av','min_Av','max_Av'],                                                       #21-24                                             
                                       ['N_Systems'],                                                                                   #Number of SSP bases
                                       ['START_W','END_W','MASK_FILE','CONFIG_FILE','NPOLY','MASK_FILE_POLY', 'N_MIN_E', 'N_MAX_E'],    #Bases config
                                       ['MIN_DELTA_CHISQ', 'MAX_NITER', 'CUT_MEDIAN_FLUX'],
                                       ['start_w_peak', 'end_w_peak'],
                                       ['wavelength_to_norm', 'width_AA', 'new_back_templates.fits']] 
        
        #Bases float indeces
        self.idcs_floats = array([0,1,4,6,7])
        
        #Emision lines mask columns headers
        self.eline_mask_header = ['start_wave', 'end_wave', 'mask_file', 'mask_config_file', 'n_poly', 'mask_file_poly', 'n_min_e', 'n_max_e']
        
        #Number of montercarlo iterations
        self.n_mc = 30
        
        #Initial value for the chiSq_min
        self.chiSq_min = 1e12
       
        return
        
    def load_FIT3D_data(self, conf_file, data_folder = None):
        
        #Check if we are executing from the folder file
        data_folder = getcwd() + '/' if data_folder == None else data_folder
        
        #Read parameters from command line
        command_dict    = self.load_FIT3D_command_params(data_folder = data_folder)
        config_dict     = self.load_FIT3D_config_file(conf_file)
        
        #Update the fit configuration giving preference to the values from the command line
        config_dict.update(command_dict)
        
        #Load observational data and masks
        config_dict = self.load_FIT3D_observational_fits(data_folder, config_dict)
        
        #Prepare output files
        output_root = config_dict['output_file'][:config_dict['output_file'].rfind('.')]
        config_dict['single_output_file']   = '{rootname}_{file_code}.{ext}'.format(rootname=output_root, file_code='single', ext='txt')
        config_dict['coeffs_output_file']   = '{rootname}_{file_code}.{ext}'.format(rootname=output_root, file_code='coeffs', ext='txt')
        config_dict['spectrum_output_file'] = '{rootname}_{file_code}.{ext}'.format(rootname=output_root, file_code='spec', ext='txt')
        config_dict['em_lines_output_file'] = '{rootname}_{file_code}.{ext}'.format(rootname=output_root, file_code='elines', ext='txt')
        
        #Delete these output files if they had been generated from a previos run #USEFULL_Function
        silent_remove([config_dict['output_file'], config_dict['single_output_file'], config_dict['coeffs_output_file'], config_dict['spectrum_output_file'], config_dict['em_lines_output_file']])
        
        #Store folder with the data and configuration folder
        config_dict['data_folder']  = data_folder
        config_dict['conf_file']    = conf_file
        config_dict['data_type']    = 'FIT3D'
        
        return config_dict

    def load_FIT3D_command_params(self,  data_folder):
                
        #Empty dictionary to store the data from the commands from the command line
        command_dict = OrderedDict()
        
        #Extract line command arguments
        self.args_list = argv
                
        #Check if the minimum parameters have been introduced (WARNING: Need to convert these to the right units)
        if len(self.args_list) > 7:
            command_dict = OrderedDict(zip(self.sspSyn_commands_params[:len(self.args_list)], self.args_list))            
        else:
            print '--Error: The input command must include all these arguments:'
            print ', '.join(self.sspSyn_commands_params[:7])
            
            #Leave the program
            #exit(0)
            
            #Currently run test example if not enought data is provided
            print '---Using example data'
            command_dict = example_data(data_folder = data_folder)
            
        return command_dict
    
    def load_FIT3D_config_file(self, config_file_address):
        
        #Empty dictionary to store the data from the config file
        fit_conf_dict = {}
        
        #Read the configuration text file
        with open(config_file_address) as conf_file:
            
            conf_lines = conf_file.readlines()
            
            #Read redshift, sigma and Av params rows
            for i in range(3):
                param_values = array(conf_lines[i].split(), dtype=float)
                fit_conf_dict.update(zip(self.sspSyn_config_params[i], param_values))
            
            #Read masks rows: 'START_W_n','END_W_n','MASK_FILE_n' ...
            nLineMasks = int(conf_lines[3])
            fit_conf_dict['nLineMasks'] = int(conf_lines[3])
            
            for i in range(4, 4 + fit_conf_dict['nLineMasks']):
                
                bases_key       = 'base_{}'.format(i - 4)
                param_values    = array(conf_lines[i].split())

                #Convert to float numerical entries
                param_values[0] = float(param_values[0])
                param_values[1] = float(param_values[1])
                param_values[4] = float(param_values[4])
                param_values[6] = float(param_values[6])
                param_values[7] = float(param_values[7])
 
                fit_conf_dict[bases_key] = param_values
                                  
            #Add ChiSq row (converting to float)
            param_values = array(conf_lines[4 + nLineMasks].split(), dtype=float)
            fit_conf_dict.update(zip(self.sspSyn_config_params[5], param_values))
            
            #Add peak wavelength row (converting to float)
            param_values = array(conf_lines[5 + nLineMasks].split(), dtype=float)
            fit_conf_dict.update(zip(self.sspSyn_config_params[6], param_values))            
            
            #Normalizing row (if available) (converting to float)
            if len(conf_lines) == 7 + nLineMasks:
                param_values = array(conf_lines[6 + nLineMasks].split(), dtype=float)
                fit_conf_dict.update(zip(self.sspSyn_config_params[7], param_values))      
            else:
                fit_conf_dict['wave_norm']      = None
                fit_conf_dict['w_wave_norm']    = None
                fit_conf_dict['new_back_file']  = None
                        
        return fit_conf_dict    

    def load_FIT3D_mask(self, config_dict, obs_flux_resam):
        
        obs_wave = config_dict['obs_wave']
        
        #--------------Generating spectrum mask
        #Load spectrum masks
        mask_xmin, mask_xmax    = loadtxt(config_dict['data_folder'] + config_dict['mask_file'], unpack = True)
        
        #Load emission lines reference to generate artificial mask
        emLine_wave             = loadtxt(config_dict['data_folder'] + config_dict['z_elines_mask'], usecols=([0]), unpack=True)
        emLine_mask_xmin        = emLine_wave * (1 + config_dict['input_z']) - 4.0 * config_dict['input_sigma']
        emLine_mask_xmax        = emLine_wave * (1 + config_dict['input_z']) + 4.0 * config_dict['input_sigma']
        
        #Firt check non zero entries
        idx_mask_zero = (obs_flux_resam != 0)
        
        #Pixels within the spectrum mask
        idx_spec_mask = ones(len(obs_wave), dtype=bool)
        for i in range(len(mask_xmin)):
            idx_cur_spec_mask   = (obs_wave > mask_xmin[i]) & (obs_wave < mask_xmax[i])
            idx_spec_mask       = idx_spec_mask & ~idx_cur_spec_mask

        #Pixels within the emline mask
        idx_emline_mask = ones(len(obs_wave), dtype=bool)
        for i in range(len(emLine_wave)):        
            idx_cur_emline_mask = (obs_wave > emLine_mask_xmin[i]) & (obs_wave < emLine_mask_xmax[i])
            idx_emline_mask     = idx_emline_mask & ~idx_cur_emline_mask

        #Recover wavelength limits for the masks
        wmin_str, wmax_str = config_dict['wmin'].split(','), config_dict['wmax'].split(',')
        wmin = float(wmin_str[0]) if len(wmin_str) == 2 else float(config_dict['wmin'])
        wmax = float(wmax_str[0]) if len(wmax_str) == 2 else float(config_dict['wmax'])        
        idx_mask_wmin, idx_mask_wmax  = (obs_wave > wmin), (obs_wave < wmax)

        #Combined individual indeces into a global mask        
        print idx_mask_zero.shape
        print idx_spec_mask.shape
        print idx_emline_mask.shape
        print idx_mask_wmax.shape
        total_masks = idx_mask_zero & idx_spec_mask & idx_emline_mask & idx_mask_wmin & idx_mask_wmax
     
        return total_masks

    def load_FIT3D_observational_fits(self, data_folder, config_dict):
        
        #--------------Read observational data
        obs_data    = loadtxt(data_folder + config_dict['input_spec'])
        obs_wave    = obs_data[:,1]
        obs_flux    = obs_data[:,2]
        obs_fluxVar = obs_data[:,3]
        
        #Issues with spectra: nan entries
        check_missing_flux_values(obs_flux)
        
        #Get the error from the library fits
        if obs_fluxVar is not None:
            obs_flux_err = sqrt(abs(obs_fluxVar))
        #Else calculate it from the spectrum
        else:
            obs_flux_err = sqrt(abs(obs_flux)/10)
        
        #Remove big error entries
        median_err                  = median(obs_flux_err)
        idx_big_err                 = (obs_flux_err > 1.5 * median_err)
        obs_fluxErrAdj              = copy(obs_flux_err)
        obs_fluxErrAdj[idx_big_err] = 1.5 * median_err        
                        
        #--------------Store data
        config_dict['obs_wave']         = obs_wave
        config_dict['obs_flux']         = obs_flux 
        config_dict['obs_flux_err']     = obs_flux_err
        config_dict['obs_fluxErrAdj']   = obs_fluxErrAdj
        config_dict['nObsPix']          = len(obs_flux)
               
        return config_dict

    def import_Fit3D_ssplibrary(self, ssp_file_address):
       
        #Dictionary to store the data
        ssp_lib_dict = {}
                
        fluxBases, hdrBases                 = getdata(ssp_file_address, 0, header=True)
        fluxBases                           = asfortranarray(fluxBases)
        nBases, nPixelsBases                = fluxBases.shape
        crpix, cdelt, crval                 = hdrBases['CRPIX1'], hdrBases['CDELT1'], hdrBases['CRVAL1']
        pixArray                            = arange(0, nPixelsBases) #WARNING should this arrange start at one?
        basesWavelength                     = (crval + cdelt * (pixArray + 1 - crpix))
                    
        #Extract age and metallicity from the bases names
        Z_vector, age_vector = empty(nBases), empty(nBases)
        for i in range(nBases):
                  
            header_code = 'NAME{}'.format(i)
             
            #Read metallicity and age from and headers list
            base_keyname    = hdrBases[header_code]
            age_str         = base_keyname[9:base_keyname.find('_z')]
            metal_str       = base_keyname[base_keyname.find('_z')+2:base_keyname.rfind('.')]
             
            age_factor      = 1000.0 if 'Myr' in age_str else 1
            age_vector[i]   = float(age_str[:-3]) / age_factor
            Z_vector[i]     = float('0.'+ metal_str)

        #Staore library data in a dictionary        
        ssp_lib_dict['crpix_bases']     = crpix
        ssp_lib_dict['cdelt_bases']     = cdelt
        ssp_lib_dict['crval_bases']     = crval  
        ssp_lib_dict['basesWave']       = basesWavelength
        ssp_lib_dict['nBases']          = nBases        
        ssp_lib_dict['nPixBases_max']   = nPixelsBases
        ssp_lib_dict['fluxBases']       = fluxBases        
        ssp_lib_dict['hdrBases']        = hdrBases
        ssp_lib_dict['ageBases']        = age_vector    
        ssp_lib_dict['zBases']          = Z_vector    
        #ssp_lib_dict['bases_one_array'] = ones(nBases)
        
        return ssp_lib_dict
        
    def import_STARLIGHT_ssplibrary(self, bases_folder, libraries_file_list):
        
        print '\n--Importing STARLIGHT library'
        print '---Bases file: {}'.format(libraries_file_list)
        print '---Bases folder: {}'.format(bases_folder)
        
        #Dictionary to store the data
        ssp_lib_dict = {}
                
        columns_names = ['file_name', 'age_yr', 'z_star', 'bases_nickname', 'f_star', 'YAV_flag', 'alpha/Fe']
        bases_df = read_csv(libraries_file_list, delim_whitespace = True, names=columns_names, skiprows=1)
        
        #Initial pass to check the biggest size
        nBases = len(bases_df.index)
        max_nPixelsBases = 0
        
        #Empty contaiores to store the data
        waveBases_orig          = []
        fluxBases_orig          = []
        Z_vector, age_vector    = empty(nBases), empty(nBases)
                
        for i in range(nBases):
        
            bases_file = bases_folder + bases_df.iloc[i]['file_name']
            wave_base_i, flux_base_i = loadtxt(bases_file, unpack=True)
                            
            #Original wavelength range and fluxes from the bases. They may have different wavelength range
            waveBases_orig.append(wave_base_i)
            fluxBases_orig.append(flux_base_i)
            
            #Interpolate the bases to observed wavelength resolution (1 angstrom per pixel is the current rule)
            age_vector[i]           = bases_df.iloc[i]['age_yr']
            Z_vector[i]             = bases_df.iloc[i]['z_star']        

        ssp_lib_dict['basesWave']           = waveBases_orig
        ssp_lib_dict['nBases']              = nBases        
        ssp_lib_dict['nPixBases_max']       = max_nPixelsBases
        ssp_lib_dict['fluxBases']           = fluxBases_orig        
        ssp_lib_dict['ageBases']            = age_vector    
        ssp_lib_dict['zBases']              = Z_vector    
        #ssp_lib_dict['bases_one_array']     = ones(nBases)
       
        print '--Library imported'
        
        return ssp_lib_dict
        
class ssp_fitter(ssp_synthesis_importer):
    
    def __init__(self):
        
        #Load tools for importing ssp libraries
        ssp_synthesis_importer.__init__(self)
        
        self.ssp_conf_dict = OrderedDict()
         
    def load_stellar_bases(self, ssp_lib_type, data_folder = None, data_file = None, resample_int = None, resample_range = None, norm_interval = (5100, 5150)):
        
        #Store stellar base type
        sspLib_dict = {}
        sspLib_dict['data_type'] = ssp_lib_type
        
        #-----------Import the base type       
        if ssp_lib_type == 'FIT3D':
            
            #Check if more files are being introduced
            if ',' in data_file:
                ssp_lib1, ssp_lib2 = data_file.split(',') #Corrently we are only using the first one (the big)
            else:
                ssp_lib1 = data_file
            
            sspLib_dict = self.import_Fit3D_ssplibrary(data_folder + ssp_lib1)
            
        elif ssp_lib_type == 'starlight':
            
            sspLib_dict = self.import_STARLIGHT_ssplibrary(data_folder, data_file)
        
        #Store stellar base type
        sspLib_dict['data_type'] = ssp_lib_type
        
        #-----------Resampling the bases
        if resample_int != None:
            bases_wave_resam = arange(int(resample_range[0]), int(resample_range[-1]), resample_int, dtype=float)
            
            #Loop through the bases. (It is assumed the bases may have different wavelength ranges)
            bases_flux_resam = empty((sspLib_dict['nBases'], len(bases_wave_resam)))
            for i in range(sspLib_dict['nBases']):
                bases_flux_resam[i,:] = interp1d(sspLib_dict['basesWave'][i], sspLib_dict['fluxBases'][i], bounds_error=True)(bases_wave_resam)
            
            sspLib_dict['basesWave_resam'] = bases_wave_resam
            sspLib_dict['basesFlux_resam'] = bases_flux_resam
        
        else:
            sspLib_dict['basesWave_resam'] = sspLib_dict['basesWave']
            sspLib_dict['basesFlux_resam'] = sspLib_dict['fluxBases']                
        
            
        #-----------Normalize the bases
        if norm_interval != None:
            normFlux_bases = empty(sspLib_dict['nBases'])
            bases_flux_norm = empty((sspLib_dict['nBases'], len(sspLib_dict['basesWave_resam'])))           
            
            #Loop through the bases (we are using the original spectra to normalize)
            normFlux_bases  = empty(sspLib_dict['nBases'])
            for i in range(sspLib_dict['nBases']):
                idx_Wavenorm_min, idx_Wavenorm_max = searchsorted(sspLib_dict['basesWave'][i], norm_interval)
                normFlux_bases[i] = median(sspLib_dict['fluxBases'][i][idx_Wavenorm_min:idx_Wavenorm_max])
                bases_flux_norm[i,:] = sspLib_dict['basesFlux_resam'][i] / normFlux_bases[i]
 
            sspLib_dict['normFlux_bases'] = normFlux_bases
            sspLib_dict['bases_flux_norm'] = bases_flux_norm 
                                                
        else:
            sspLib_dict['normFlux_bases'] = 1.0
            sspLib_dict['bases_flux_norm'] = sspLib_dict['basesFlux_resam']             
        
        
        #Display the bases limits:
        if resample_range != None:
            print '---Bases wavelength range:', sspLib_dict['basesWave_resam'][0], sspLib_dict['basesWave_resam'][-1]

             
        return sspLib_dict

    def calculate_synthStellarSED(self, Av_star, z_star, sigma_star, coeff_file, sspLib_starlight, resample_range=(4000, 6900), resample_int = 1, norm_factor = 100):

        #Dictionary to store the data   
        synth_dict = {}
        
        #Load population weights
        bases_idx, bases_coeff = loadtxt(coeff_file, usecols=[0,1], unpack=True)
        
        #Synth spectrum in bases library units
        synth_dict['obs_wave']          = sspLib_starlight['basesWave_resam']
        synth_dict['obs_flux']          = np_sum(bases_coeff.T * sspLib_starlight['basesFlux_resam'].T, axis=1) * 0.1          
                
        #Synth spectrum including physical parameters
        obs_wave_rest                   = arange(int(resample_range[0]), int(resample_range[-1]), resample_int, dtype=float)
        obs_wave_resam                  = obs_wave_rest * (1.0 + z_star)
        ssp_grid                        = self.physical_SED_model(sspLib_starlight['basesWave_resam'], obs_wave_resam, sspLib_starlight['bases_flux_norm'], Av_star, z_star, sigma_star, Rv_coeff = 3.4)
        obs_flux_norm                   = np_sum(bases_coeff.T * ssp_grid, axis=1)
        obs_flux_resam                  = obs_flux_norm * norm_factor
        nBases                          = ssp_grid.shape[1]
        nBases_one_vector               = ones(nBases) 
         
        #Generate synth spectrum
        sigma_err                       = 0.01 * median(obs_flux_resam)
        obs_fluxEr_resam                = random.normal(0.0, sigma_err, len(obs_flux_norm))        
        #obs_fluxEr_norm                 = obs_fluxEr_resam / norm_factor
        obs_fluxEr_norm                 = 0.02
                              
        #Synthetic masks
        masks_default = OrderedDict()
        
        masks_default['He1_4026A']      = (4019, 4033)
        masks_default['He1_4471A']      = (4463, 4480)
        masks_default['He1_5876A']      = (5867, 5885)
        masks_default['He1_6678A']      = (6667, 6687)
        
        masks_default['HI_delta']       = (4090,4114)
        masks_default['HI_gamma']       = (4329,4353)
        masks_default['HI_beta']        = (4840,4880)
        masks_default['HI_alpha']       = (6551,6575)
        
        #Pixels within the spectrum mask
        redshift = 1 + z_star
        boolean_mask = ones(len(obs_flux_resam), dtype=bool)
        for line in masks_default:
            wmin, wmax = masks_default[line]
            idx_cur_spec_mask   = (obs_wave_resam > wmin * redshift) & (obs_wave_resam < wmax * redshift)
            boolean_mask        = boolean_mask & ~idx_cur_spec_mask
                
        #Store the data vector
        synth_dict['normFlux_stellar']      = norm_factor
        synth_dict['stellar_flux_True']     = obs_flux_resam
        synth_dict['stellar_flux_norm']     = obs_flux_norm + obs_fluxEr_norm
        synth_dict['stellar_wave_resam']    = obs_wave_resam
        synth_dict['stellar_wave_rest']     = obs_wave_rest
        synth_dict['stellar_flux_resam']    = obs_flux_resam + obs_fluxEr_resam
        synth_dict['stellar_fluxEr_resam']  = obs_fluxEr_resam
        synth_dict['stellar_fluxEr_norm']   = obs_fluxEr_norm
        synth_dict['stellar_one_array']   = ones(nBases) 
        synth_dict['nObsPix']           = len(obs_flux_resam)       
        synth_dict['boolean_mask']      = boolean_mask
        synth_dict['int_mask']          = synth_dict['boolean_mask'] * 1
                           
        return synth_dict

    def ready_data_sspFit(self, obs_dict, ssp_lib_dict):
                    
        #-----------Load object mask
#         boolean_mask = self.load_observation_mask(obs_dict, obs_flux_resam)
#         int_mask = boolean_mask * 1
        int_mask = ones(len(obs_dict['obs_wave_resam']))
                                    
        #-----------Save the data        
        fit_data = {}
        
        #Quick naming
        fit_data['obs_wave_resam']          = obs_dict['obs_wave_resam'] 
        fit_data['normFlux_obs']            = obs_dict['normFlux_obs'] 
        fit_data['obs_flux_norm_masked']    = obs_dict['obs_flux_norm'] * int_mask
        fit_data['basesWave_resam']         = ssp_lib_dict['basesWave_resam'] 
        fit_data['bases_flux_norm']         = ssp_lib_dict['bases_flux_norm']
        fit_data['int_mask']                = int_mask
        
        return fit_data

    def physical_SED_model(self, bases_wave_rest, obs_wave, bases_flux, Av_star, z_star, sigma_star, Rv_coeff = 3.4):
                        
        #Calculate wavelength at object z
        wave_z = bases_wave_rest * (1 + z_star) 
        
        #Calculate stellar broadening kernel
        r_sigma         = sigma_star/(wave_z[1] - wave_z[0])

        box             = int64(3 * r_sigma) if int64(3 * r_sigma) < 3 else 3
        kernel_len      = 2 * box + 1
        kernel          = zeros((1, kernel_len)) 
        kernel_range    = arange(0, 2 * box + 1)

        #Generating the kernel with sigma (the norm factor is the sum of the gaussian)
        kernel[0,:]     = exp(-0.5 * ((square(kernel_range-box)/r_sigma)))
        norm            = np_sum(kernel[0,:])        
        kernel          = kernel / norm
         
        #Convove bases with respect to kernel for dispersion velocity calculation
        bases_grid_convolve = convolve2d(bases_flux, kernel, mode='same', boundary='symm')  
        
        #Interpolate bases to wavelength range
#         print 'Z', z_star
#         print wave_z[0], wave_z[-1]
#         print obs_wave[0], obs_wave[-1]
#         print
        bases_grid_interp       = (interp1d(wave_z, bases_grid_convolve, axis=1, bounds_error=True)(obs_wave)).T       
        
        #Generate final flux model including reddening
        Av_vector               = Av_star * ones(bases_grid_interp.shape[1])
        obs_wave_resam_rest     = obs_wave/(1 + z_star)
        Xx_redd                 = CCM89_Bal07(Rv_coeff, obs_wave_resam_rest) 
        dust_attenuation        = power(10, -0.4 * outer(Xx_redd, Av_vector))
        bases_grid_redd         = bases_grid_interp * dust_attenuation

        return bases_grid_redd

    def generate_synthObs(self, bases_wave, bases_flux, basesCoeff, Av_star, z_star, sigma_star, resample_range = None, resample_int = 1):
        
        '''basesWave: Bases wavelength must be at rest'''
        nbases              = basesCoeff.shape[0]        
        bases_wave_resam    = arange(int(resample_range[0]), int(resample_range[-1]), resample_int, dtype=float)
        npix_resample       = len(bases_wave_resam)
        
        #Resampling the range
        bases_flux_resam = empty((nbases, npix_resample))
        for i in range(nbases):
#             print bases_wave[i][0], bases_wave[i][-1]
#             print bases_wave_resam[0], bases_wave_resam[-1]
            bases_flux_resam[i,:] = interp1d(bases_wave[i], bases_flux[i], bounds_error=True)(bases_wave_resam)            
        
        #Display physical parameters
        synth_wave                  = bases_wave_resam * (1 + z_star)
        
        Av_vector                   = Av_star * ones(nbases)
        Xx_redd                     = CCM89_Bal07(3.4, bases_wave_resam)   
        r_sigma                     = sigma_star/(synth_wave[1] - synth_wave[0])

        #Defining empty kernel
        box                         = int(3 * r_sigma) if int(3 * r_sigma) < 3 else 3
        kernel_len                  = 2 * box + 1
        kernel                      = zeros((1, kernel_len)) 
        kernel_range                = arange(0, 2 * box + 1)
        
        #Generating the kernel with sigma (the norm factor is the sum of the gaussian)
        kernel[0,:]                 = exp(-0.5 * ((square(kernel_range-box)/r_sigma)))
        norm                        = np_sum(kernel[0,:])        
        kernel                      = kernel / norm

        #Convove bases with respect to kernel for dispersion velocity calculation
        bases_grid_convolve         = convolve2d(bases_flux_resam, kernel, mode='same', boundary='symm')  
        
        #Interpolate bases to wavelength range
        interBases_matrix           = (interp1d(bases_wave_resam, bases_grid_convolve, axis=1, bounds_error=True)(bases_wave_resam)).T       

        #Generate final flux model including dust        
        dust_attenuation            = power(10, -0.4 * outer(Xx_redd, Av_vector))
        bases_grid_model            = interBases_matrix * dust_attenuation
                    
        #Generate combined flux
        synth_flux                  = np_sum(basesCoeff.T * bases_grid_model, axis=1)
                                      
        return synth_wave, synth_flux

    def ssp_fitting(self, ssp_grid_masked, obs_flux_masked):
                
        optimize_result     = nnls(ssp_grid_masked, obs_flux_masked)
        
        coeffs_bases_norm   = optimize_result[0]
        
        return coeffs_bases_norm
    
    def ssp_fit(self, input_z, input_sigma, input_Av, obs_wave, obs_flux_masked, fit_data, fit_scheme = 'nnls'):
                
        #Quick naming
        #obs_wave        = fit_data['obs_wave_resam']
        #obs_flux_masked = fit_data['obs_flux_norm_masked']
        rest_wave       = fit_data['basesWave_resam']           #This one could be part of the self.
        bases_flux      = fit_data['bases_flux_norm']           #This one could be part of the self.
        int_mask        = fit_data['int_mask']                  
        obsFlux_mean    = fit_data['normFlux_obs']              #This one could be part of the self.
        
        #Apply physical data to stellar grid
        ssp_grid = self.physical_SED_model(rest_wave, obs_wave, bases_flux, input_Av, input_z, input_sigma, 3.4)
        ssp_grid_masked = (int_mask * ssp_grid.T).T

        #Non-negative least square fit
        if fit_scheme == 'nnls':
            #start = timer()   
            optimize_result = nnls(ssp_grid_masked, obs_flux_masked)
            #end = timer()
            #print 'nnls', ' time ', (end - start), '\n'   
            
            coeffs_bases_norm   = optimize_result[0]
            coeffs_bases        = coeffs_bases_norm * obsFlux_mean

        #---Leasts square fit
        elif fit_scheme == 'lsq':
            start = timer()   
            optimize_result = lsq_linear(ssp_grid_masked, obs_flux_masked, bounds=(0, inf))
            end = timer()
            print 'lsq', ' time ', (end - start)
            
            coeffs_bases = optimize_result.x * obsFlux_mean
          
        #---Linear fitting without restrictions
        else:
            
            start = timer()
               
            #First guess
            coeffs_bases = self.linfit1d(obs_flux_masked, obsFlux_mean, ssp_grid_masked, inv_pdl_error_i)
    
            #Count positive and negative coefficients
            idx_plus_0  = coeffs_bases[:] > 0 
            plus_coeff  = idx_plus_0.sum()
            neg_coeff   = (~idx_plus_0).sum()
            
            #Start loops
            counter = 0
            if plus_coeff > 0:
                
                while neg_coeff > 0:
                    counter += 1
                                        
                    bases_model_n = zeros([nObsPix, plus_coeff])
                                    
                    idx_plus_0                          = (coeffs_bases[:] > 0)
                    bases_model_n[:,0:idx_plus_0.sum()] = bases_grid_model_masked[:,idx_plus_0] #These are replace in order
                    coeffs_bases[~idx_plus_0]               = 0 
                    
                    #Repeat fit
                    coeffs_n = self.linfit1d(obsFlux_normMasked, obsFlux_mean, bases_model_n, inv_pdl_error_i)
                    
                    idx_plus_n  = coeffs_n[:] > 0
                    idx_min_n   = ~idx_plus_n
                    plus_coeff  = idx_plus_n.sum()
                    neg_coeff   = (idx_min_n).sum()
                    
                    #Replacing negaive by zero
                    coeffs_n[idx_min_n] = 0
                    coeffs_bases[idx_plus_0] = coeffs_n
    
                    if plus_coeff == 0:
                        neg_coeff = 0     
            else:
                plus_coeff = nBases
            
            end = timer()
            print 'FIT3D', ' time ', (end - start)
              
        #Save data to export
        fit_products                        = {}
        flux_sspFit_norm                    = np_sum(coeffs_bases_norm.T * ssp_grid_masked, axis=1)
        flux_sspFit                         = np_sum((coeffs_bases_norm * obsFlux_mean).T * ssp_grid_masked, axis=1)
        fit_products['flux_components']     = coeffs_bases.T * ssp_grid_masked
        fit_products['weight_coeffs']       = coeffs_bases
        fit_products['flux_sspFit_norm']    = flux_sspFit_norm
        fit_products['flux_sspFit']         = flux_sspFit

        return fit_products 
    
    def load_observation_mask(self, obs_dict, obs_flux_resam):

        #Declare type of masks
        data_type = obs_dict['data_type']
        
        if data_type == 'FIT3D':      

            boolean_mask = self.load_FIT3D_mask(obs_dict, obs_flux_resam)
            
        elif data_type == 'starlight':
              
            boolean_mask = ones(len(obs_flux_resam)) #self.starlight_mask(obs_dict, obs_flux_resam)
                             
        return boolean_mask
  
    def linfit1d(self, obsFlux_norm, obsFlux_mean, basesFlux, weight):
         
        nx, ny = basesFlux.shape
         
        #Case where the number of pixels is smaller than the number of bases
        if nx < ny:
            basesFlux = transpose(basesFlux)
            nx = ny
         
        A = basesFlux
        B = obsFlux_norm
                 
        #Weight definition #WARNING: Do we need to use the diag?
        if weight.shape[0] == nx:
            weight   = diag(weight)
            A       = dot(weight, A)
            B       = dot(weight, transpose(B))
        else:
            B       = transpose(B)
     
        coeffs_0    = dot(linalg.inv(dot(A.T, A)), dot(A.T, B)) * obsFlux_mean
         
        return coeffs_0         
        
        
        
        
        
        
        
        
        
        
    