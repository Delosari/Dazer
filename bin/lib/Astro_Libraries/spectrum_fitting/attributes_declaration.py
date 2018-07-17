def gaussian_filter1d_ppxf(spec, sig):
    """
    Convolve a spectrum by a Gaussian with different sigma for every pixel.
    If all sigma are the same this routine produces the same output as
    scipy.ndimage.gaussian_filter1d, except for the border treatment.
    Here the first/last p pixels are filled with zeros.
    When creating a template library for SDSS data, this implementation
    is 60x faster than a naive for loop over pixels.

    :param spec: vector with the spectrum to convolve
    :param sig: vector of sigma values (in pixels) for every pixel
    :return: spec convolved with a Gaussian with dispersion sig

    """
    sig = sig.clip(0.01)  # forces zero sigmas to have 0.01 pixels
    p = int(np.ceil(np.max(3*sig)))
    m = 2*p + 1  # kernel size
    x2 = np.linspace(-p, p, m)**2

    n = spec.size
    a = np.zeros((m, n))
    # fig, ax = plt.subplots(1, 1, figsize=(16, 10))

    for j in range(m):   # Loop over the small size of the kernel
        #print j, n-m+j+1
        indices = n-m+j+1
        a[j,:] = spec
        a[j, p:-p] = spec[j:n-m+j+1]
        # ax.plot(waveData, a[j,:], label=j)

    # ax.update({'xlabel': 'Wavelength (nm)', 'ylabel': 'Flux (normalised)'})
    # ax.legend()
    # plt.show()

    gau = np.exp(-x2[:, None]/(2*sig**2))
    gau /= np.sum(gau, 0)[None, :]  # Normalize kernel

    conv_spectrum = np.sum(a*gau, 0)

    return conv_spectrum

def gaussian_filter1d_vital(sigma, wave, flux):

    # Calculate stellar broadening kernel
    r_sigma         = sigma / (wave[1] - wave[0])

    # Kernel matrix
    box             = 2 #np.int64(3 * r_sigma) if np.int64(3 * r_sigma) < 3 else 3
    kernel_len      = 2 * box + 1
    kernel          = np.zeros((1, kernel_len))
    kernel_range    = np.arange(0, 2 * box + 1)

    # Generating gaussian kernel with sigma (the norm factor is the sum of the gaussian)
    kernel[0, :]    = np.exp(-0.5 * ((np.square(kernel_range - box) / r_sigma)))
    norm            = np.sum(kernel[0, :])
    kernel          = kernel / norm

    #Perform convolution
    flux_convolved  = convolve2d(flux, kernel, mode='same', boundary='symm')

    p = int(np.ceil(np.max(3*sigma)))
    m = 2*p + 1  # kernel size
    x2 = np.linspace(-p, p, m)**2
    n = flux.size
    a = np.zeros((m, n))
    for j in range(m):   # Loop over the small size of the kernel
        a[j, p:-p] = flux[j:n-m+j+1]
    gau = np.exp(-x2[:, None]/(2*sigma**2))
    gau /= np.sum(gau, 0)[None, :]  # Normalize kernel
    myKernel = np.exp(-0.5 * ((np.square(kernel_range - box) / sigma**2)))
    myKernelN = myKernel / np.sum(myKernel)
    flux_convolved = convolve2d(np.array([flux]), myKernelN, mode='same', boundary='symm')

    print
    print 'box', box
    print 'p', p
    print
    print 'kernel (initial)', np.zeros((1, kernel_len)).shape
    print 'a (initial)', np.zeros((m, n)).shape
    print
    print 'kernel_len', kernel_len
    print 'm', m
    print
    print 'kernel_range', kernel_range
    print 'x2', x2
    print 'np.square(kernel_range - box)', np.square(kernel_range - box)
    print
    print 'x2[:, None]',x2[:, None]
    print 'np.square(kernel_range - box)', np.square(kernel_range - box)
    print
    print 'np.exp(-x2[:, None]/(2*sig**2))', np.exp(-x2[:, None]/(2*sigma**2))
    print 'np.exp(-0.5 * ((np.square(kernel_range - box) / sigma)))', myKernel
    print
    print 'myKernelN', myKernelN
    print 'gau_N', gau
    print
    print 'kernel', kernel
    print 'gau', gau
    print

    return flux_convolved

# class Attributes_SpectrumFitter():
#
#     def __init__(self):
#
#         # Temperature and density grid declaration
#         self.tem_grid_range = arange(temp_grid[0], temp_grid[1], temp_grid[2])
#         self.den_grid_range = arange(den_grid[0], den_grid[1], den_grid[2])
#         self.den_grid_range[0] = 1 #Change the minimum value for the grid to 1
#
#         # Reddening parameters
#         self.Rv_model = R_v
#         self.reddedning_curve_model = reddenig_curve
#
#         # Declare high ionization temperature ions
#         self.high_temp_ions = high_temp_ions
#
#         # Lower accepted limit for ssps
#         self.lowlimit_sspContribution = lowlimit_sspContribution
#
#         # Dictionary to store parameters
#         self.conf_dic = {'Te_neb': 10000.0, 'z_neb': 0.0, 'cHbeta_neb': 0.1, 'He1_neb': 0.085, 'He2_neb': 0.0}
#
#
#         # Dictionary to store the data (This is in gen_synth_obs)
#         self.obj_data = {}
#         self.obj_data['obj_properties_file']   = obj_properties_file
#         self.obj_data['obj_lines_file']        = obj_lines_file
#         self.obj_data['obj_ssp_coeffs_file']   = obj_ssp_coeffs_file
#         self.obj_data['obs_mask_address']      = obj_mask_file
#         self.obj_data['output_folder']         = output_folder
#         self.obj_data['flux_hbeta']            = obj_prop_df.loc['flux_hbeta'][0]
#
#         # Dictionary with synthetic abundances
#         self.abund_dict = dict(zip(abund_keys, abund_values.T))
#
#         # Import the physical data from the log
#         for param in obj_prop_df.index:
#             self.obj_data[param] = obj_prop_df.loc[param][0]
#
#         # Reddening parameters
#         self.obj_data['lineFlambda'] = self.gasExtincParams(self.obj_data['lineWaves'], self.Rv_model, self.reddedning_curve_model,Hbeta_wave = 4861.331)
#
#
#         #Calculate T_high assuming we get T_low
#         self.obj_data['T_high'] = TOIII_TSIII_relation(self.obj_data['T_low'])
#
#         # Compute lines flux
#         self.obj_data['lineFluxes'] = self.calcEmFluxes(self.obj_data['T_low'], self.obj_data['T_high'],
#                                                         obj_prop_df.loc['n_e'][0],
#                                                         obj_prop_df.loc['cHbeta'][0], obj_prop_df.loc['tau'][0],
#                                                         self.abund_dict,
#                                                         self.obj_data['lineLabes'],
#                                                         self.obj_data['lineIons'],
#                                                         self.obj_data['lineFlambda'])
#
#
#         # Use general error if this is provided
#         self.obj_data['lineErr'] = self.obj_data['lineFluxes'] * error_lines
#
#         # Save input conditions:
#         self.obj_data['wavelengh_limits'] = wavelengh_limits
#         self.obj_data['resample_inc'] = resample_inc
#         self.obj_data['norm_interval'] = norm_interval
#         self.obj_data['z_obj'] = obj_prop_df.loc['z_obj'][0]
#
#         # Rest and observed wavelength
#         obj_wave_rest = np.arange(wavelengh_limits[0], wavelengh_limits[-1], resample_inc, dtype=float)
#         obj_wave_obs = obj_wave_rest * (1.0 + self.obj_data['z_obj'])
#         self.obj_data['obs_wave_rest'] = obj_wave_rest
#         self.obj_data['obs_wave'] = obj_wave_obs
#
#         # Get Halpha flux to calibrate
#         idx_Halpha = (self.obj_data['lineLabels'] == 'H1r_6563A')
#         self.obj_data['flux_halpha'] = self.obj_data['recomb_fluxes'][idx_Halpha] * obj_prop_df.loc['flux_hbeta'][0]
#
#         # Reddening parameters for the nebular continuum
#         self.obj_data['nebFlambda'] = self.gasExtincParams(obj_wave_rest, self.Rv_model, self.reddedning_curve_model)
#
#         # Calculate the nebular continua
#         self.obj_data['obs_flux'] = self.nebFluxCont(obj_wave_rest,
#                                                      obj_prop_df.loc['cHbeta'][0], self.obj_data['nebFlambda'],
#                                                      obj_prop_df.loc['T_low'][0],
#                                                      obj_prop_df.loc['He1_abund'][0], obj_prop_df.loc['He2_abund'][0],
#                                                      self.obj_data['flux_halpha'])
#
#
#         # Save input conditions:
#         self.obj_data['wavelengh_limits'] = wavelengh_limits
#         self.obj_data['resample_inc'] = resample_inc
#         self.obj_data['norm_interval'] = norm_interval
#         self.obj_data['z_star'] = obj_prop_df.loc['z_star'][0]
#         self.obj_data['Av_star'] = obj_prop_df.loc['Av_star'][0]
#         self.obj_data['sigma_star'] = obj_prop_df.loc['sigma_star'][0]
#         self.obj_data['flux_hbeta'] = obj_prop_df.loc['flux_hbeta'][0]
#         self.obj_data['eqw_hbeta'] = obj_prop_df.loc['eqw_hbeta'][0]
#
#
#         self.sspPrefit_Idcs = np.where(idx_populations)[0]
#         self.sspPrefit_Coeffs = bases_coeff[idx_populations]
#         self.sspPrefit_Limits = np.vstack((self.sspPrefit_Coeffs * 0.8, self.sspPrefit_Coeffs * 1.2)).T
#
#         #Object data
#         self.input_continuum        = object_continuum
#         self.input_continuum_er     = obj_flux_err
#         self.input_wave             = obj_wave
#         self.int_mask               = obj_mask
#
#         #Bases parameters
#         neglected_populations       = np.where(~idx_populations)
#         self.onBasesWave            = self.ssp_lib['wave_resam']
#         self.onBasesFlux            = np.delete(self.ssp_lib['flux_resam'], neglected_populations, axis=0)
#         self.onBasesFluxNorm        = np.delete(self.ssp_lib['flux_norm'], neglected_populations, axis=0)
#         self.onBasesFluxNormCoeffs  = np.delete(self.ssp_lib['normFlux_coeff'], neglected_populations, axis=0)
#         self.range_bases            = np.arange(self.onBasesFlux.shape[0])
#
#         # Limit for bases
#         z_max_ssp                   = (self.ssp_lib['wave_resam'][0] / obj_wave[0]) - 1.0
#         z_min_ssp                   = (self.ssp_lib['wave_resam'][-1] / obj_wave[-1]) - 1.0
#         self.zMin_SspLimit          = z_min_ssp
#         self.zMax_SspLimit          = round(z_max_ssp - 0.001, 3)
#         self.z_object               = self.obj_data['z_star']
#
#         # Store the data vector
#         self.obj_data['obs_wave_rest'] = obj_wave_rest
#         self.obj_data['obs_wave'] = obj_wave_obs
#         self.obj_data['stellar_flux'] = obj_flux
#         self.obj_data['stellar_flux_err'] = obj_flux + stellar_err
#         self.obj_data['sigma_continuum'] = 0.02  # TODO put this in the object data
#
#
#
# class Attributes_ImportModelData():
#
#     def __init__(self):
#
#         self.paths_dict = {}
#
#         self.lines_df = read_excel(self.paths_dict['lines_data_file'], sheetname=0, header=0, index_col=0)
#
#         # Save the data in the dictionary (Excluding Hbeta)
#         idx_lines = (obj_lines_df.index.isin(input_lines)) & (obj_lines_df.index != 'H1_4861A')
#         self.obj_data['lineIons']   = obj_lines_df.loc[idx_lines].ion.values
#         self.obj_data['lineLabes']  = obj_lines_df.loc[idx_lines].index.values
#         self.obj_data['lineWaves']  = obj_lines_df.loc[idx_lines].obs_wavelength.values
#         self.obj_data['linePynebCode'] = obj_lines_df.loc[idx_lines].pynebCode.values
#         self.obj_data['lineFluxes'] = obj_lines_df.loc[idx_lines].obs_flux.values
#         self.obj_data['lineErr']    = obj_lines_df.loc[idx_lines].obs_fluxErr.values
#
#         # Generate the dictionary with pyneb ions
#         print('-- Loading atoms data with PyNeb')
#         self.ionDict = pn.getAtomDict(atom_list=obj_lines_df.ion.values)
#
#         print('-- Atomic sources Loaded ')
#         # for atom in self.ionDict.keys():
#         #     textPrint = '--- {}: {}'.format(atom, self.ionDict[atom].printSources())
#         #     print(textPrint)
#
#         # Establish index of lines which below to high and low ionization zones
#         self.idx_highU = np.in1d(self.obj_data['lineIons'], self.high_temp_ions)
#         self.idx_lowU = ~self.idx_highU
#
#         # Attributes to increase calculation speed
#         self.n_recombLines, self.n_colExcLines = np.sum(self.obj_data['lineType'] == 'rec'), np.sum(self.obj_data['lineType'] == 'col')
#         self.range_recombLines, self.range_colExcLines = np.arange(self.n_recombLines), np.arange(self.n_colExcLines)
#
#         # Empty Ionic dictionary for the MCMC
#         self.abund_dict = {ion: 0 for ion in self.ionDict.keys()}
#
#         # Save input conditions:
#         self.obj_data['wavelengh_limits'] = wavelengh_limits
#         self.obj_data['resample_inc'] = resample_inc
#         self.obj_data['norm_interval'] = norm_interval
#         self.obj_data['z_star'] = obj_prop_df.loc['z_star'][0]
#
#
# class Attributes_EmissionComponents():
#
#     def __init__(self):
#
#         #Hbeta configuration
#         self.Hbeta_label = 'H1_4861A'
#         self.Hbeta_wave = 4862.683
#         self.Hbeta_pynebCode = '4_2'
#
#         #Import Optical depth function
#         posHelium_Lines = ['He1_3889A','He1_4026A','He1_4387A', 'He1_4471A', 'He1_4686A','He1_4714A','He1_4922A','He1_5876A','He1_6678A','He1_7065A','He1_7281A','He1_10830A']
#         self.Coef_ftau_dict = self.import_optical_depth_coeff_table(posHelium_Lines)
#
#         self.Hbeta_xX = rc_Gas.X(self.Hbeta_wave)
#
#         self.obj_data['lineXx'] = rc_Gas.X(self.obj_data['lineWaves'])
#         self.obj_data['lineFlambda'] = self.obj_data['lineXx'] / self.Hbeta_xX - 1.0
#
#         # Emissivity grid
#         self.emis_grid = self.manage_emissivity_grids(self.obj_data['linePynebCode'], self.obj_data['lineIons'], forceGridsReset)
#
#         # Emissivity grid for collisional excited lines
#         self.Hbeta_emis_grid = self.manage_emissivity_grids(np.array([4861]), np.array(['H1r']), forceGridsReset)
#
#         # Normalizing emis_grid grids by Hbeta
#         self.emis_grid = self.recomb_emis_grid / self.Hbeta_emis_grid
#
#         # logarithmic scale log scales
#         self.log_emis_grid = np.log10(self.emis_grid)
#
#         # Dictionary to store the emissivity surface coeffients
#         self.emisCoeffs = {}
#
# class EmissionEquations():
#
#     def __init__(self):