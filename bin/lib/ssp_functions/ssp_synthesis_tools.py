from errno                          import ENOENT
from os                             import remove
from sys                            import argv, exit
from collections                    import OrderedDict
from pandas                         import DataFrame
from astropy.io.fits                import getdata
from scipy.signal.signaltools       import convolve2d
from scipy.interpolate.interpolate  import interp1d
from numpy import array, core, loadtxt, empty, sqrt, abs, sum as np_sum, isnan, ones, copy, median, arange, exp, zeros, transpose, mean, diag, linalg, dot
import pyneb as pn

pn.RedCorr.law

#Function to delete files
def silentremove(filename_list):
    for filename in filename_list:   
        try:
            remove(filename)
        except OSError as e:        # this would be "except OSError, e:" before Python 2.6
            if e.errno != ENOENT:   # errno.ENOENT = no such file or directory
                raise               # re-raise exception if a different error occurred

def example_data(folder_data):
    
    arguments_dict = OrderedDict()
    arguments_dict['script']            = folder_data + 'auto_ssp_elines_rnd.py'            #0
    arguments_dict['input_spec']        = folder_data + 'NGC5947.spec_5.txt'                #1
    arguments_dict['SSPs_lib']          = folder_data +'ssp_lib.fits,'+folder_data +'ssp_lib.fits'                     #2
    arguments_dict['output_file']       = folder_data + 'auto_ssp.NGC5947.cen.only.out'     #3
    arguments_dict['mask_file']         = folder_data + 'mask_elines.txt'                   #4
    arguments_dict['conf_file']         = folder_data + 'auto_ssp_V500_several_Hb.config'   #5
    arguments_dict['plot_tag']          = 1                                                 #6
    arguments_dict['min']               = -1                                                #7
    arguments_dict['max']               = 40                                                #8
    arguments_dict['wmin']              = '3850'                                              #9
    arguments_dict['wmax']              = '6800'                                              #10
    arguments_dict['z_elines_mask']     = folder_data + 'emission_lines.txt'                #11
    arguments_dict['input_z']           = 0.02                                              #12
    arguments_dict['delta_z']           = 0.001                                             #13
    arguments_dict['min_z']             = 0.015                                             #14
    arguments_dict['max_z']             = 0.025                                             #15
    arguments_dict['input_sigma']       = 2.0                                               #16 #This one used to be 2
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
    
    x   = 1e4 / wave #Assuming wavelength is in Amstrongs
    ax  = zeros(len(wave))
    bx  = zeros(len(wave))
    
    idcs    = x > 1.1
    y       = (x[idcs] - 1.82)
    
    ax[idcs]    = 1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
    bx[idcs]    = 1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
    ax[~idcs]   = 0.574*x[~idcs]**1.61
    bx[~idcs]   = -0.527*x[~idcs]**1.61
            
    Xx  = ax+bx/Rv #WARNING better to rewrite this law
    
    return Xx

class ssp_fitter():
    
    def __init__(self):
        
        self.ssp_conf_dict = OrderedDict()
        
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
        self.chiSq_min = 1e12
        
    def load_command_params(self):
        
        #Empty dictionary to store the data from the config file
        command_dict = OrderedDict()
        
        #Extract line command arguments
        self.args_list = argv
                
        #Check if the minimum parameters have been introduced (WARNING: Need to convert these to the right units)
        if len(self.args_list) > 7:
            command_dict = OrderedDict(zip(self.sspSyn_commands_params[:len(self.args_list)], self.args_list))            
        else:
            print '--Error: The input command must include all these arguments:'
            print ', '.join(self.sspSyn_commands_params[:7])
            # exit(0)
            
            print '-Using sample data'
            command_dict = example_data(folder_data = '/home/vital/workspace/Fit_3D/example_files/')
            
        return command_dict
        
    def load_config_params(self, config_file_address):
        
        #Empty dictionary to store the data from the config file
        fit_conf_dict = OrderedDict()
        
        with open(config_file_address) as conf_file:
            
            conf_lines = conf_file.readlines()
            
            #Add redshift, sigma and Av params rows
            for i in range(3):
                param_values = array(conf_lines[i].split(), dtype=float)
                fit_conf_dict.update(zip(self.sspSyn_config_params[i], param_values))
            
            #Add bases rows: 'START_W_n','END_W_n','MASK_FILE_n' ...
            self.n_line_masks = int(conf_lines[3])
            self.emlines_df   = DataFrame(columns=self.eline_mask_header)
            for i in range(4, 4 + self.n_line_masks):
                
                bases_key       = 'base_{}'.format(i - 4)
                param_values    = array(conf_lines[i].split())

                #Convert to float numerical entries
                param_values[0] = float(param_values[0])
                param_values[1] = float(param_values[1])
                param_values[4] = float(param_values[4])
                param_values[6] = float(param_values[6])
                param_values[7] = float(param_values[7])
 
                fit_conf_dict[bases_key]    = param_values
                
                #Add row to dataframe
                self.emlines_df.loc[i - 4] = param_values
                  
            #Add ChiSq row (converting to float)
            param_values = array(conf_lines[4 + self.n_line_masks].split(), dtype=float)
            fit_conf_dict.update(zip(self.sspSyn_config_params[5], param_values))
            
            #Add peak wavelength row (converting to float)
            param_values = array(conf_lines[5 + self.n_line_masks].split(), dtype=float)
            fit_conf_dict.update(zip(self.sspSyn_config_params[6], param_values))            
            
            #Normalizing row (if available) (converting to float)
            if len(conf_lines) == 7 + self.n_line_masks:
                param_values = array(conf_lines[6 + self.n_line_masks].split(), dtype=float)
                fit_conf_dict.update(zip(self.sspSyn_config_params[7], param_values))      
            else:
                fit_conf_dict['wave_norm'] = None
                fit_conf_dict['w_wave_norm'] = None
                fit_conf_dict['new_back_file'] = None
                        
        #Update the configuration dictionary giving priority a those from the command line
        return fit_conf_dict
    
    def load_input_data(self, conf_dict, data_folder = ''):
        
        #Dictionary with all the simulation configuration parameters
        self.fit_conf = conf_dict

        #WARNING: DEFT PARAMETER
        self.deft       = 1 if self.args_list == 9 else 2
        self.abs_min    = self.fit_conf['CUT_MEDIAN_FLUX']
        self.iscale     = 0 #'Pixel resolution'
        self.ymin       = 0 #Min flux
        self.ymax       = 0 #Max flux
        
        #WARNING: WAVELENGTHS LIMITS
        wmin_str, wmax_str = self.fit_conf['wmin'].split(','), self.fit_conf['wmax'].split(',')
        self.wmin   = float(wmin_str[0]) if len(wmin_str) == 2 else float(self.fit_conf['wmin'])
        self.wmin2  = float(wmin_str[1]) if len(wmin_str) == 2 else self.wmin
        self.wmax   = float(wmax_str[0]) if len(wmax_str) == 2 else float(self.fit_conf['wmax'])
        self.wmax2  = float(wmax_str[1]) if len(wmax_str) == 2 else self.wmax
       
        #--Load the stellar databases.
        if ',' in self.fit_conf['SSPs_lib']:
            #Two databases model, the second list provides de kinematic data            
            bases_list = self.fit_conf['SSPs_lib'].split(',')
            self.ssp_lib_lum_address = data_folder + bases_list[0]
            self.ssp_lib_kin_address = data_folder + bases_list[1]   
        else:
            #One bases library provides all the data
            self.ssp_lib_lum_address = data_folder + self.fit_conf['SSPs_bases']
            self.ssp_lib_kin_address = self.ssp_lib_lum_address            
    
        #--Output file
        output_root                = self.fit_conf['output_file'][:self.fit_conf['output_file'].rfind('.')]
        self.fit_conf['single']    = '{rootname}_{file_code}.{ext}'.format(rootname=output_root, file_code='single', ext='txt')
        self.fit_conf['coeffs']    = '{rootname}_{file_code}.{ext}'.format(rootname=output_root, file_code='coeffs', ext='txt')
        self.fit_conf['spectrum']  = '{rootname}_{file_code}.{ext}'.format(rootname=output_root, file_code='spec', ext='txt')
        self.fit_conf['em_lines']  = '{rootname}_{file_code}.{ext}'.format(rootname=output_root, file_code='elines', ext='txt')
        silentremove([self.fit_conf['output_file'], self.fit_conf['single'], self.fit_conf['spectrum'], self.fit_conf['em_lines']])
    
        #--Setting SFH template         
        #WARNING WE SHOULD JUST MAKE A ROUTINE TO OPEN ALL FITS INDIVIDUALLY
        sfh_base_flux, sfh_hdr          = getdata(self.ssp_lib_lum_address, 0, header=True)
        nBases_sfh, nPixels_sfh         = sfh_base_flux.shape
        coeffs_basis_sfh                = empty([nBases_sfh, 3])
        crpix_sfh, cdelt_shf, crval_sfh = sfh_hdr['CRPIX1'], sfh_hdr['CDELT1'], sfh_hdr['CRVAL1']
                
        #--Setting kinematics template
        kin_base_flux, kin_hdr          = getdata(self.ssp_lib_kin_address, 0, header=True)
        nBases_kin, nPixels_kin         = kin_base_flux.shape
        coeffs_basis_kin                = empty([nBases_sfh, 3])
        crpix_kin, cdelt_kin, crval_kin = sfh_hdr['CRPIX1'], sfh_hdr['CDELT1'], sfh_hdr['CRVAL1']
        
#         print 'Que estamos importando'
#         print '1', nBases_sfh, nPixels_sfh, self.ssp_lib_lum_address
#         print '2', nBases_kin, nPixels_kin, self.ssp_lib_kin_address
        
        #Load spectrum masks
        mask_xmin, mask_xmax = loadtxt(self.fit_conf['mask_file'], unpack = True)
        
        #Load emission lines reference to generate artificial mask
        emLine_wave         = loadtxt(self.fit_conf['z_elines_mask'], usecols=(0), unpack=True)
        emLine_mask_xmin    = emLine_wave * (1 + self.fit_conf['input_z']) - 4.0 * self.fit_conf['input_sigma']
        emLine_mask_xmax    = emLine_wave * (1 + self.fit_conf['input_z']) + 4.0 * self.fit_conf['input_sigma']
        
        #Read the data
        spectrum_matrix     = loadtxt(self.fit_conf['input_spec'])
        obs_wave, obs_flux  = spectrum_matrix[:,1], spectrum_matrix[:,2]
        
        #Get obs flux variance data if available
        if spectrum_matrix.shape[1] == 4:
            obs_flux_var = spectrum_matrix[:,3]
            obs_flux_err = sqrt(abs(obs_flux_var))
        else:
            obs_flux_err = sqrt(abs(obs_flux)/10)
        
        #Flux error vector
        median_err                      = median(obs_flux_err)
        idx_big_err                     = (obs_flux_err > 1.5 * median_err)
        obs_flux_err_adj                = copy(obs_flux_err)
        obs_flux_err_adj[idx_big_err]   = 1.5 * median_err
                
        #Issues with spectra: nan entries
        check_missing_flux_values(obs_flux)
        
        #--Generating spectrum mask
        
        #Firt check non zero entries
        idx_mask_zero = (obs_flux != 0)
        
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

        #Pixels outside wavelength edges
        idx_mask_wmin = (obs_wave > self.wmin)
        idx_mask_wmax = (obs_wave < self.wmax)
        
        #Combine flux into one array        
        total_masks = idx_mask_zero & idx_spec_mask & idx_emline_mask & idx_mask_wmin & idx_mask_wmax

        #The masked flux is set to zero
        obs_flux_mask = copy(obs_flux)
        obs_flux_mask[~total_masks] = 0 
        
        #Save data into dictionary
        self.fit_conf['obs_wave']           = obs_wave
        self.fit_conf['obs_flux']           = obs_flux
        self.fit_conf['obs_flux_mask']      = obs_flux_mask
        self.fit_conf['obs_flux_err']       = obs_flux_err
        self.fit_conf['obs_flux_err_adj']   = obs_flux_err_adj
        self.fit_conf['total_mask']         = total_masks
        self.fit_conf['crpix_kin']          = crpix_kin
        self.fit_conf['cdelt_kin']          = cdelt_kin
        self.fit_conf['crval_kin']          = crval_kin  
        
        pix_array                           = arange(0, nPixels_kin)
                
        self.fit_conf['kinWave_zCor']       = (crval_kin + cdelt_kin * (pix_array + 1 - crpix_kin)) * (1 + self.fit_conf['input_z'])

        self.fit_conf['nBases_lum']         = nBases_sfh        
        self.fit_conf['nPixels_lum']        = nPixels_sfh
         
        self.fit_conf['nBases_kin']         = nBases_kin        
        self.fit_conf['nPixels_kin']        = nPixels_kin        
               
        self.fit_conf['kin_base_flux']      = kin_base_flux        
        self.fit_conf['kin_hdr']            = kin_hdr    
        self.fit_conf['zero_mask']          = total_masks * 1.0
                    
        return

    def fit_ssp(self):
        
        #obs_wave, obs_flux, obs_flux_err, obs_flux_err_adj, input_z, input_sigma, input_Av, nBases_kin, nPixels_kin, kin_base_flux, kin_hdr, Av_mask
        obs_wave            = self.fit_conf['obs_wave']
        obs_flux            = self.fit_conf['obs_flux'] 
        obs_flux_mask       = self.fit_conf['obs_flux'] * self.fit_conf['zero_mask']
        obs_flux_err_adj    = self.fit_conf['obs_flux_err_adj']
        input_z             = self.fit_conf['input_z']
        input_sigma         = self.fit_conf['input_sigma']
        input_Av            = self.fit_conf['input_Av']
        
        #WE should remove this kin rubbish and use a standard
        nBases_kin          = self.fit_conf['nBases_kin']
        nPixels_kin         = self.fit_conf['nPixels_kin']
        kin_base_flux       = self.fit_conf['kin_base_flux']
        kin_hdr             = self.fit_conf['kin_hdr']
        zero_mask           = self.fit_conf['zero_mask']

        #Extra constants
        n_pixFlux           = len(obs_flux)
        Av_vector           = input_Av * ones(nBases_kin) #WARNING: Now this is going with the basis...
        kin_model           = empty((n_pixFlux, nBases_kin))
        kin_model_no_mask   = empty((n_pixFlux, nBases_kin))
                
        #Redshift corrected wavelength for the bases
        wave_corr           = self.fit_conf['kinWave_zCor']
        dpix_c_val          = None #WARNING: Do we need this one?
          
        #Dispersion velocity definition
        r_sigma             = input_sigma/(wave_corr[1] - wave_corr[0])
        
        #Loop through the 'small' kinematic populations
        age_vector          = empty(nBases_kin)
        Z_vector            = empty(nBases_kin) 

        #--Non-linear parameters fitting
        for i in range(nBases_kin):
            
            header_code = 'NAME{}'.format(i)
            
            #Read metallicity and age from and headers list
            base_keyname    = kin_hdr[header_code]
            age_str         = base_keyname[9:base_keyname.find('_z')]
            metal_str       = base_keyname[base_keyname.find('_z')+2:base_keyname.rfind('.')]
            
            age_factor      = 1000.0 if 'Myr' in age_str else 1
            age_vector[i]   = float(age_str[:-3]) / age_factor
            Z_vector[i]     = float('0.'+ metal_str)
            
            #Defining empty kernel
            box     = int(3 * r_sigma) if int(3 * r_sigma) < 3 else 3
            kernel  = zeros([1,2*box+1])
        
            #Generating the kernel with sigma (the norm factor is the sum of the gaussian)
            norm = 0.0
            for j in range(0, 2*box + 1):
                gaus        = exp(-0.5 * (((j-box)/r_sigma)**2))    
                kernel[0,j] = gaus
                norm        = norm+gaus
            kernel = kernel / norm
              
            #Generate kinetic model
            kin_base_flux_conv      = convolve2d(kin_base_flux, kernel, mode='same')
            kin_base_flux_conv_i    = kin_base_flux_conv[i,:]
            
            #print '-----kin_base_flux_conv_i', kin_base_flux_conv_i.shape
            
            base_kin_flux_resamp_i  = interp1d(wave_corr, kin_base_flux_conv_i, bounds_error=False, fill_value=0.)(obs_wave)
                
            
            #print '-----base_kin_flux_resamp_i', base_kin_flux_resamp_i.shape
                        
            #Adjust the bases to the mask            
            kin_model[:,i]          = base_kin_flux_resamp_i * zero_mask
            kin_model_no_mask[:,i]  = base_kin_flux_resamp_i
            kin_model_err_i         = 0.01 * abs(obs_flux_err_adj) * zero_mask #Error in the bases modelled after the object spectrum
            
#             print 'My kin_model_no_mask', np_sum(kin_model_no_mask[:,i])
#             print 'My kin_model_err_i', np_sum(kin_model_err_i)

            #Check if there are nan entries in the bases
            check_missing_flux_values(kin_model[:,i])
        
        #--Adjusting the dust attenuation        
        kin_model_dust              = empty((n_pixFlux, nBases_kin))        
        kin_model_final             = empty((n_pixFlux, nBases_kin))
        kin_model_final_nomask      = empty((n_pixFlux, nBases_kin))
        pdl_error_i                 = ones(len(obs_flux_err_adj))
        idx_zero                    = (obs_flux_err_adj == 0)
        pdl_error_i[~idx_zero]      = 1.0/(abs(obs_flux_err_adj[~idx_zero])**2)
        
        for j in range(0, nBases_kin):
            wave_res                    = obs_wave/(1 + input_z)
            dust_rat                    = CCM89_Bal07(3.1, wave_res)
            atten_dust                  = 10 ** (-0.4 * Av_vector[j] * dust_rat)  
            kin_model_dust[:, j]        = atten_dust
            kin_model_final[:,j]        = kin_model[:, j] * atten_dust
            kin_model_final_nomask[:,j] = kin_model_no_mask[:,j] * atten_dust
            #pdl_error_i                 = 1.0/(abs(obs_flux_err_adj)**2) #WARNING: This one is just rewritten     

        print 'My kin_model_final', np_sum(kin_model_final)
  
        #--Linear parameters fitting 

#         print '-----The redshift', input_z
#         print '-----The input_sigma', input_sigma
#         print '-----The input_Av', input_Av, Av_vector.shape
#         print '-----The crval2', self.fit_conf['crval_kin']
#         print '-----The cdelt2', self.fit_conf['cdelt_kin']
#         print '-----The crpix2', self.fit_conf['crpix_kin']        
#         print '-----nf2/nBases_kin', nBases_kin
#         print '-----n_c2/nPixels_kin', nPixels_kin
#         print '-----pdl_flux_c_ini2/kin_base_flux', kin_base_flux.shape
# #         print '-----kin_hdr/hdr2', kin_hdr
#         print '-----obs_wave/wave_unc', obs_wave
#         print '-----masked_Av/Av_vector', Av_vector.shape, Av_vector
#         print '-----zero_mask', zero_mask.shape, zero_mask
#         print '-----obs_flux_err_adj/e_flux_unc', obs_flux_err_adj.shape, obs_flux_err_adj
#         print '-----obs_flux/flux_unc', obs_flux.shape, obs_flux
#         print '-----self.n_mc', self.n_mc
#         print '-----self.chiSq_min', self.chiSq_min
#         print '-----wave_corr', wave_corr
#         print '-----r_sigma', r_sigma
#         print '-----wave_corr', wave_corr
#         print '-----r_sigma', r_sigma
#         print '-----kin_model', kin_model.shape, np_sum(kin_model)
#         print '-----kin_model_final', kin_model_final.shape, np_sum(kin_model_final)

    
        #Initial stellar synthesis without restrictions
        
#         print 'Inputs linfit1d'
#         print np_sum(obs_flux)
#         print np_sum(obs_flux_mask)
#         print np_sum(kin_model_final)
#         print np_sum(1/pdl_error_i)
        
        y_fit_0, coeffs_0   = self.linfit1d(obs_flux_mask, kin_model_final, 1/pdl_error_i)
        y_fit_0             = y_fit_0[:,0]
 
#         print '-----y_model_now', y_fit_0, np_sum(y_fit_0), y_fit_0.shape
#         print '-----coeffs', coeffs_0 
        
        #Count positive and negative coefficients
        plus_coeff  = (coeffs_0[:,0] > 0).sum()
        neg_coeff   = (coeffs_0[:,0] < 0).sum()

        print '----nf_new:', plus_coeff
        print '----nf_neg:', neg_coeff
        print '----n_unc:', n_pixFlux          
        print '----coeffs_0:', coeffs_0.shape          
        print '----nBases_kin:', nBases_kin          
             
        #Start loops
        counter = 0
        if plus_coeff > 0:
            
            while neg_coeff > 0:
                                                
                bases_model_n = zeros([n_pixFlux, plus_coeff])
                                
                idx_plus_0                          = (coeffs_0[:,0] > 0)
                bases_model_n[:,0:idx_plus_0.sum()] = kin_model_final[:,idx_plus_0] #These are replace in order
                coeffs_0[~idx_plus_0,0]             = 0 
                
                #Repeat fit
                y_fit_n, coeffs_n = self.linfit1d(obs_flux_mask, bases_model_n, 1/pdl_error_i)
                
                idx_plus_n  = coeffs_n[:,0] > 0
                idx_min_n   = ~idx_plus_n
                plus_coeff  = idx_plus_n.sum()
                neg_coeff   = (idx_min_n).sum()
                
                #Replacing negaive by zero
                coeffs_n[idx_min_n] = 0
                coeffs_0[idx_plus_0] = coeffs_n

                if plus_coeff == 0:
                    neg_coeff = 0     
        else:
            plus_coeff = nBases_kin

        #Calculate physical parameters from fit
        obs_fit_spectrum = np_sum(coeffs_0.T * kin_model_final_nomask, axis=1)

        print '-----obs_fit_spectrum', np_sum(obs_fit_spectrum)
       
        return obs_fit_spectrum
                       
    def linfit1d(self, flux, bases_kin_dust_model, weigh=[1]):
        
        nx, ny      = bases_kin_dust_model.shape
        flux        = array([flux])
        
        if nx < ny:
            bases_kin_dust_model = transpose(bases_kin_dust_model)
            nx = ny
        
        pdl_flux_m  = flux/mean(flux)
        A           = bases_kin_dust_model
        B           = pdl_flux_m
        
        #Weight definition
        if len(weigh) == nx: 
            weigh   = diag(weigh)
            A       = dot(weigh, A)
            B       = dot(weigh, transpose(B))
        else:
            B       = transpose(B)
    
        coeffs_0    = dot(linalg.inv(dot(A.T, A)), dot(A.T, B)) * mean(flux)
        pdl_model_0 = dot(A,coeffs_0)
    
        return pdl_model_0, coeffs_0
      
        
        
        
        
        
        
        
        
        
        
        
        
    