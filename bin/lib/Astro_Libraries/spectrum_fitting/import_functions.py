import os
import sys
import numpy as np
import ConfigParser
from errno import ENOENT
from numpy import loadtxt
from pandas import read_excel, read_csv
from collections import OrderedDict
from scipy.interpolate import interp1d
from distutils.util import strtobool
from astropy.io import fits as astropyfits

# Function to create folders
def make_folder(folder_path):

    #TODO This one is only valid for 2.7
    #TODO add this one to your collection
    try:
        os.makedirs(folder_path)
    except OSError:
        if not os.path.isdir(folder_path):
            raise

    return


# Function to delete files
def silent_remove(filename_list):
    for filename in filename_list:
        try:
            os.remove(filename)
        except OSError as e:  # this would be "except OSError, e:" before Python 2.6
            if e.errno != ENOENT:  # errno.ENOENT = no such file or directory
                raise  # re-raise exception if a different error occurred


# Sample data for FIT3D compilation
def example_data(data_folder):
    arguments_dict = OrderedDict()
    arguments_dict['script'] = 'auto_ssp_elines_rnd.py'  # 0
    arguments_dict['input_spec'] = 'NGC5947.spec_5.txt'  # 1
    arguments_dict['SSPs_lib'] = 'ssp_lib.fits,' + 'ssp_lib.fits'  # 2
    arguments_dict['output_file'] = 'auto_ssp.NGC5947.cen.only.out'  # 3
    arguments_dict['mask_file'] = 'mask_elines.txt'  # 4
    arguments_dict['conf_file'] = 'auto_ssp_V500_several_Hb.config'  # 5
    arguments_dict['plot_tag'] = 1  # 6
    arguments_dict['min'] = -1  # 7
    arguments_dict['max'] = 40  # 8
    arguments_dict['wmin'] = '3850'  # 9
    arguments_dict['wmax'] = '6800'  # 10
    arguments_dict['z_elines_mask'] = 'emission_lines.txt'  # 11
    arguments_dict['input_z'] = 0.02  # 12
    arguments_dict['delta_z'] = 0.001  # 13
    arguments_dict['min_z'] = 0.015  # 14
    arguments_dict['max_z'] = 0.025  # 15
    arguments_dict['input_sigma'] = 2.0  # 16
    arguments_dict['delta_sigma'] = 0.5  # 17
    arguments_dict['min_sigma'] = 1  # 18
    arguments_dict['max_sigma'] = 9  # 19
    arguments_dict['input_Av'] = 0.5  # 20
    arguments_dict['delta_Av'] = 0.1  # 21
    arguments_dict['min_Av'] = 0.0  # 22
    arguments_dict['max_Av'] = 1.6  # 23

    return arguments_dict


# Function to check for nan entries
def check_missing_flux_values(flux):
    # Evaluate the nan array
    nan_idcs = np.isnan(flux)
    nan_count = np.sum(nan_idcs)

    # Directly save if not nan
    if nan_count > 0:
        print '--WARNING: missing flux entries'

    return


# Function to import configuration data
def parseObjData(file_address, sectionName, objData):

    # json.dump(obj_data, open('/home/vital/PycharmProjects/thesis_pipeline/spectrum_fitting/synth_objProperties_json.txt', 'w'),
    # indent=4, sort_keys=True)

    parser = ConfigParser.SafeConfigParser()
    parser.optionxform = str

    if os.path.isfile(file_address):
        parser.read(file_address)

    if not parser.has_section(sectionName):
        parser.add_section(sectionName)

    for key in objData.keys():
        value = objData[key]
        if value is not None:
            if isinstance(value, list) or isinstance(value, np.ndarray):
                value = ','.join(str(x) for x in value)
            else:
                value = str(value)
        else:
            value = ''

        parser.set(sectionName, key, value)

    with open(file_address, 'w') as f:
        parser.write(f)


    return


class SspSynthesisImporter:

    def __init__(self):

        # ------------Configuration of Fit3D
        self.sspSyn_commands_params = [
            'script',  # 0 python script name
            'input_spec',  # 1 input galactic spectrum name
            'SSPs_lib',  # 2 fits-table to use with python
            'output_file',  # 3 Reference name for output files
            'mask_file',  # 4 File with the spectrum region masks
            'conf_file',  # 5 Configuration file for the masks
            'plot_tag',  # 6 tag to launch the plotting
            'min',  # 7 Min flux for ploting
            'max',  # 8 Max flux for ploting
            'wmin',  # 9 Minimum wavelength for plotting
            'wmax',  # 10 Maximum wavelength for plotting
            'z_elines_mask',  # 11 Emission lines file
            'input_z',  # 12 Input redshift
            'delta_z',  # 13 Increments for redshift
            'min_z',  # 14 Minimum redshift
            'max_z',  # 15 Maximum redshift
            'input_sigma',  # 16 Input velocity dispersion
            'delta_sigma',  # 17 Increments for velocity dispersion
            'min_sigma',  # 18 Minimum velocity dispersion
            'max_sigma',  # 19 Maximum velocity dispersion
            'input_Av',  # 20 Input reddening
            'delta_Av',  # 21 Increments for reddening
            'min_Av',  # 22 Minimum reddening
            'max_Av',  # 23 Maximum reddening
        ]

        # The first 4 lines in the configuration file describe the input
        self.sspSyn_config_params = [['input_z', 'delta_z', 'min_z', 'max_z', 'DV', 'RV', 'DS', 'RS', 'MIN_W', 'MAX_W'],
                                     # 12-16
                                     ['input_sigma', 'delta_sigma', 'min_sigma', 'max_sigma'],
                                     # 17-20
                                     ['input_Av', 'delta_Av', 'min_Av', 'max_Av'],
                                     # 21-24
                                     ['N_Systems'],  # Number of SSP bases
                                     ['START_W', 'END_W', 'MASK_FILE', 'CONFIG_FILE', 'NPOLY', 'MASK_FILE_POLY',
                                      'N_MIN_E', 'N_MAX_E'],  # Bases config
                                     ['MIN_DELTA_CHISQ', 'MAX_NITER', 'CUT_MEDIAN_FLUX'],
                                     ['start_w_peak', 'end_w_peak'],
                                     ['wavelength_to_norm', 'width_AA', 'new_back_templates.fits']]

        # Bases float indeces
        self.idcs_floats = np.array([0, 1, 4, 6, 7])

        # Emision lines mask columns headers
        self.eline_mask_header = ['start_wave', 'end_wave', 'mask_file', 'mask_config_file', 'n_poly', 'mask_file_poly',
                                  'n_min_e', 'n_max_e']

        # Number of montercarlo iterations
        self.n_mc = 30

        # Initial value for the chiSq_min
        self.chiSq_min = 1e12

        return

    def load_FIT3D_data(self, conf_file, data_folder=None):

        # Check if we are executing from the folder file
        data_folder = os.getcwd() + '/' if data_folder is None else data_folder

        # Read parameters from command line
        command_dict = self.load_FIT3D_command_params(data_folder=data_folder)
        config_dict = self.load_FIT3D_config_file(conf_file)

        # Update the fit configuration giving preference to the values from the command line
        config_dict.update(command_dict)

        # Load observational data and masks
        config_dict = self.load_FIT3D_observational_fits(data_folder, config_dict)

        # Prepare output files
        output_root = config_dict['output_file'][:config_dict['output_file'].rfind('.')]
        config_dict['single_output_file'] = '{rootname}_{file_code}.{ext}'.format(rootname=output_root,
                                                                                  file_code='single', ext='txt')
        config_dict['coeffs_output_file'] = '{rootname}_{file_code}.{ext}'.format(rootname=output_root,
                                                                                  file_code='coeffs', ext='txt')
        config_dict['spectrum_output_file'] = '{rootname}_{file_code}.{ext}'.format(rootname=output_root,
                                                                                    file_code='spec', ext='txt')
        config_dict['em_lines_output_file'] = '{rootname}_{file_code}.{ext}'.format(rootname=output_root,
                                                                                    file_code='elines', ext='txt')

        # Delete these output files if they had been generated from a previos run #USEFULL_Function
        silent_remove([config_dict['output_file'], config_dict['single_output_file'], config_dict['coeffs_output_file'],
                       config_dict['spectrum_output_file'], config_dict['em_lines_output_file']])

        # Store folder with the data and configuration folder
        config_dict['data_folder'] = data_folder
        config_dict['conf_file'] = conf_file
        config_dict['data_type'] = 'FIT3D'

        return config_dict

    def load_FIT3D_command_params(self, data_folder):

        # Empty dictionary to store the data from the commands from the command line
        command_dict = OrderedDict()

        # Extract line command arguments
        self.args_list = sys.argv

        # Check if the minimum parameters have been introduced (WARNING: Need to convert these to the right units)
        if len(self.args_list) > 7:
            command_dict = OrderedDict(zip(self.sspSyn_commands_params[:len(self.args_list)], self.args_list))
        else:
            print '--Error: The input command must include all these arguments:'
            print ', '.join(self.sspSyn_commands_params[:7])

            # Currently run test example if not enought data is provided
            print '---Using example data'
            command_dict = example_data(data_folder=data_folder)

        return command_dict

    def load_FIT3D_config_file(self, config_file_address):

        # Empty dictionary to store the data from the config file
        fit_conf_dict = {}

        # Read the configuration text file
        with open(config_file_address) as conf_file:

            conf_lines = conf_file.readlines()

            # Read redshift, sigma and Av params rows
            for i in range(3):
                param_values = np.array(conf_lines[i].split(), dtype=float)
                fit_conf_dict.update(zip(self.sspSyn_config_params[i], param_values))

            # Read masks rows: 'START_W_n','END_W_n','MASK_FILE_n' ...
            nLineMasks = int(conf_lines[3])
            fit_conf_dict['nLineMasks'] = int(conf_lines[3])

            for i in range(4, 4 + fit_conf_dict['nLineMasks']):
                bases_key = 'base_{}'.format(i - 4)
                param_values = np.array(conf_lines[i].split())

                # Convert to float numerical entries
                param_values[0] = float(param_values[0])
                param_values[1] = float(param_values[1])
                param_values[4] = float(param_values[4])
                param_values[6] = float(param_values[6])
                param_values[7] = float(param_values[7])

                fit_conf_dict[bases_key] = param_values

            # Add ChiSq row (converting to float)
            param_values = np.array(conf_lines[4 + nLineMasks].split(), dtype=float)
            fit_conf_dict.update(zip(self.sspSyn_config_params[5], param_values))

            # Add peak wavelength row (converting to float)
            param_values = np.array(conf_lines[5 + nLineMasks].split(), dtype=float)
            fit_conf_dict.update(zip(self.sspSyn_config_params[6], param_values))

            # Normalizing row (if available) (converting to float)
            if len(conf_lines) == 7 + nLineMasks:
                param_values = np.array(conf_lines[6 + nLineMasks].split(), dtype=float)
                fit_conf_dict.update(zip(self.sspSyn_config_params[7], param_values))
            else:
                fit_conf_dict['wave_norm'] = None
                fit_conf_dict['w_wave_norm'] = None
                fit_conf_dict['new_back_file'] = None

        return fit_conf_dict

    def load_FIT3D_mask(self, config_dict, obs_flux_resam):

        obs_wave = config_dict['obs_wave']

        # --------------Generating spectrum mask
        # Load spectrum masks
        mask_xmin, mask_xmax = loadtxt(config_dict['data_folder'] + config_dict['mask_file'], unpack=True)

        # Load emission lines reference to generate artificial mask
        emLine_wave = loadtxt(config_dict['data_folder'] + config_dict['z_elines_mask'], usecols=([0]), unpack=True)
        emLine_mask_xmin = emLine_wave * (1 + config_dict['input_z']) - 4.0 * config_dict['input_sigma']
        emLine_mask_xmax = emLine_wave * (1 + config_dict['input_z']) + 4.0 * config_dict['input_sigma']

        # Firt check non zero entries
        idx_mask_zero = (obs_flux_resam != 0)

        # Pixels within the spectrum mask
        idx_spec_mask = np.ones(len(obs_wave), dtype=bool)
        for i in range(len(mask_xmin)):
            idx_cur_spec_mask = (obs_wave > mask_xmin[i]) & (obs_wave < mask_xmax[i])
            idx_spec_mask = idx_spec_mask & ~idx_cur_spec_mask

        # Pixels within the emline mask
        idx_emline_mask = np.ones(len(obs_wave), dtype=bool)
        for i in range(len(emLine_wave)):
            idx_cur_emline_mask = (obs_wave > emLine_mask_xmin[i]) & (obs_wave < emLine_mask_xmax[i])
            idx_emline_mask = idx_emline_mask & ~idx_cur_emline_mask

        # Recover wavelength limits for the masks
        wmin_str, wmax_str = config_dict['wmin'].split(','), config_dict['wmax'].split(',')
        wmin = float(wmin_str[0]) if len(wmin_str) == 2 else float(config_dict['wmin'])
        wmax = float(wmax_str[0]) if len(wmax_str) == 2 else float(config_dict['wmax'])
        idx_mask_wmin, idx_mask_wmax = (obs_wave > wmin), (obs_wave < wmax)

        # Combined individual indeces into a global mask
        print idx_mask_zero.shape
        print idx_spec_mask.shape
        print idx_emline_mask.shape
        print idx_mask_wmax.shape
        total_masks = idx_mask_zero & idx_spec_mask & idx_emline_mask & idx_mask_wmin & idx_mask_wmax

        return total_masks

    def load_FIT3D_observational_fits(self, data_folder, config_dict):

        # --------------Read observational data
        obs_data = loadtxt(data_folder + config_dict['input_spec'])
        obs_wave = obs_data[:, 1]
        obs_flux = obs_data[:, 2]
        obs_fluxVar = obs_data[:, 3]

        # Issues with spectra: nan entries
        check_missing_flux_values(obs_flux)

        # Get the error from the library fits
        if obs_fluxVar is not None:
            obs_flux_err = np.sqrt(abs(obs_fluxVar))
        # Else calculate it from the spectrum
        else:
            obs_flux_err = np.sqrt(abs(obs_flux) / 10)

        # Remove big error entries
        median_err = np.median(obs_flux_err)
        idx_big_err = (obs_flux_err > 1.5 * median_err)
        obs_fluxErrAdj = np.copy(obs_flux_err)
        obs_fluxErrAdj[idx_big_err] = 1.5 * median_err

        # --------------Store data
        config_dict['obs_wave'] = obs_wave
        config_dict['obs_flux'] = obs_flux
        config_dict['obs_flux_err'] = obs_flux_err
        config_dict['obs_fluxErrAdj'] = obs_fluxErrAdj
        config_dict['nObsPix'] = len(obs_flux)

        return config_dict

    def import_Fit3D_ssplibrary(self, ssp_file_address):

        # Dictionary to store the data
        ssp_lib_dict = {}

        fluxBases, hdrBases = astropyfits.getdata(ssp_file_address, 0, header=True)
        fluxBases = np.asfortranarray(fluxBases)
        nBases, nPixelsBases = fluxBases.shape
        crpix, cdelt, crval = hdrBases['CRPIX1'], hdrBases['CDELT1'], hdrBases['CRVAL1']
        pixArray = np.arange(0, nPixelsBases)  # WARNING should this arrange start at one?
        basesWavelength = (crval + cdelt * (pixArray + 1 - crpix))

        # Extract age and metallicity from the bases names
        Z_vector, age_vector = np.empty(nBases), np.empty(nBases)
        for i in range(nBases):
            header_code = 'NAME{}'.format(i)

            # Read metallicity and age from and headers list
            base_keyname = hdrBases[header_code]
            age_str = base_keyname[9:base_keyname.find('_z')]
            metal_str = base_keyname[base_keyname.find('_z') + 2:base_keyname.rfind('.')]

            age_factor = 1000.0 if 'Myr' in age_str else 1
            age_vector[i] = float(age_str[:-3]) / age_factor
            Z_vector[i] = float('0.' + metal_str)

        # Staore library data in a dictionary
        ssp_lib_dict['crpix_bases'] = crpix
        ssp_lib_dict['cdelt_bases'] = cdelt
        ssp_lib_dict['crval_bases'] = crval
        ssp_lib_dict['basesWave'] = basesWavelength
        ssp_lib_dict['nBases'] = nBases
        ssp_lib_dict['nPixBases_max'] = nPixelsBases
        ssp_lib_dict['fluxBases'] = fluxBases
        ssp_lib_dict['hdrBases'] = hdrBases
        ssp_lib_dict['ageBases'] = age_vector
        ssp_lib_dict['zBases'] = Z_vector
        # ssp_lib_dict['bases_one_array'] = ones(nBases)

        return ssp_lib_dict

    def import_STARLIGHT_ssplibrary(self, bases_folder, libraries_file_list):

        print '\n--Importing STARLIGHT library'
        print '---Bases file: {}'.format(libraries_file_list)
        print '---Bases folder: {}'.format(bases_folder)

        # Dictionary to store the data
        ssp_lib_dict = {}

        columns_names = ['file_name', 'age_yr', 'z_star', 'bases_nickname', 'f_star', 'YAV_flag', 'alpha/Fe']
        bases_df = read_csv(libraries_file_list, delim_whitespace=True, names=columns_names, skiprows=1)

        # Initial pass to check the biggest size
        nBases = len(bases_df.index)
        max_nPixelsBases = 0

        # Empty contaiores to store the data
        waveBases_orig = []
        fluxBases_orig = []
        Z_vector, age_vector = np.empty(nBases), np.empty(nBases)

        for i in range(nBases):
            bases_file = bases_folder + bases_df.iloc[i]['file_name']
            wave_base_i, flux_base_i = loadtxt(bases_file, unpack=True)

            # Original wavelength range and fluxes from the bases. They may have different wavelength range
            waveBases_orig.append(
                wave_base_i)  # This is not pretty but not other option if bases do not have same length...
            fluxBases_orig.append(
                flux_base_i)  # This is not pretty but not other option if bases do not have same length...

            # Interpolate the bases to observed wavelength resolution (1 angstrom per pixel is the current rule)
            age_vector[i] = bases_df.iloc[i]['age_yr']
            Z_vector[i] = bases_df.iloc[i]['z_star']

        ssp_lib_dict['basesWave'] = waveBases_orig  # This is not pretty but not other option if bases do not have same length...
        ssp_lib_dict['nBases'] = nBases
        ssp_lib_dict['nPixBases_max'] = max_nPixelsBases
        ssp_lib_dict['fluxBases'] = fluxBases_orig  # This is not pretty but not other option if bases do not have same length...
        ssp_lib_dict['ageBases'] = age_vector
        ssp_lib_dict['zBases'] = Z_vector
        # ssp_lib_dict['bases_one_array']     = ones(nBases)

        print '--Library imported'

        return ssp_lib_dict


class ImportModelData(SspSynthesisImporter):
    '''
    This class loads the binaries required for the dazer spectra synthesizer module
    '''

    def __init__(self):

        # Class with tools to import starlight bases
        SspSynthesisImporter.__init__(self)

        self.paths_dict = {}

        # Paths for linux:
        if os.name == 'posix':
            # TODO This data should be read from a text file
            self.paths_dict['inference_folder']     = '/home/vital/Astrodata/Inference_output/'
            self.paths_dict['nebular_data_folder']  = '/home/vital/Dropbox/Astrophysics/Lore/NebularContinuum/'
            self.paths_dict['Hydrogen_CollCoeff']   = '/home/vital/PycharmProjects/dazer/bin/lib/Astro_Libraries/literature_data/Neutral_Hydrogen_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_CollCoeff']     = '/home/vital/PycharmProjects/dazer/bin/lib/Astro_Libraries/literature_data/Neutral_Helium_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_OpticalDepth']  = '/home/vital/PycharmProjects/dazer/bin/lib/Astro_Libraries/literature_data/Benjamin1999_OpticalDepthFunctionCoefficients.txt'
            self.paths_dict['lines_data_file']      = '/home/vital/PycharmProjects/dazer/bin/lib/Astro_Libraries/spectrum_fitting/lines_data.xlsx'
            self.paths_dict['dazer_path']           = '/home/vital/PycharmProjects/dazer/'
            self.paths_dict['stellar_data_folder']  = '/home/vital/Starlight/'
            self.paths_dict['observations_folder']  = '/home/vital/Dropbox/Astrophysics/Data/WHT_observations/objects/'
            self.paths_dict['emisGridsPath']        = '/home/vital/PycharmProjects/thesis_pipeline/spectrum_fitting/testing_output/input_data/'

        # Paths for windows:
        elif os.name == 'nt':
            # TODO Check guides for windows and linux
            self.paths_dict['inference_folder']     = 'D:/Inference_data/'
            self.paths_dict['nebular_data_folder']  = 'E:/Cloud Storage/Dropbox/Astrophysics/Lore/NebularContinuum/'
            self.paths_dict['Hydrogen_CollCoeff']   = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Neutral_Hydrogen_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_CollCoeff']     = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Neutral_Helium_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_OpticalDepth']  = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Benjamin1999_OpticalDepthFunctionCoefficients.txt'
            self.paths_dict['stellar_data_folder']  = 'E:/Cloud Storage/Dropbox/Astrophysics/Tools/Starlight/'
            self.paths_dict['lines_data_file']      = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/spectrum_fitting/lines_data.xlsx'
            self.paths_dict['dazer_path']           = 'C:/Users/lativ/git/dazer/'
            self.paths_dict['observations_folder']  = 'E:/Cloud Storage/Dropbox/Astrophysics/Data/WHT_observations/objects/'
            self.paths_dict['emisGridsPath']        = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/spectrum_fitting/'

        # Import database with lines labels information
        self.linesDb = read_excel(self.paths_dict['lines_data_file'], sheetname=0, header=0, index_col=0)

        # Additional code to adapt to old lines log labelling
        self.linesDb['pyneb_code'] = self.linesDb['pyneb_code'].astype(str)
        idx_numeric_pynebCode = ~self.linesDb['pyneb_code'].str.contains('A+')
        self.linesDb.loc[idx_numeric_pynebCode, 'pyneb_code'] = self.linesDb.loc[
            idx_numeric_pynebCode, 'pyneb_code'].astype(int)
        self.linesDb['ion'].apply(str)

    def import_optical_depth_coeff_table(self):

        Data_dict = OrderedDict()

        opticalDepthCoeffs_df = read_csv(self.paths_dict['Helium_OpticalDepth'], delim_whitespace=True, header=0)

        opticalDepthCoeffs = {}
        for column in opticalDepthCoeffs_df.columns:
            opticalDepthCoeffs[column] = opticalDepthCoeffs_df[column].values

        return opticalDepthCoeffs

    def load_obsData(self, obsFile=None, objName=None):

        # TODO this should go into the master configuration
        list_parameters     = ['input_lines', 'Av_prefit','sigma_star_prefit', 'coeffsPop_prefit', 'coeffsPopErr_prefit', 'wavelengh_limits', 'norm_interval'] #also all 'param_prior'
        boolean_parameters  = ['Normalized_by_Hbeta']
        string_parameters   = ['address_lines_log', 'address_spectrum', 'address_obs_mask', 'obsFile', 'objName']

        # ----Load the obj data
        if obsFile is not None:
            cfg = ConfigParser.SafeConfigParser()
            cfg.optionxform = str
            cfg.read(obsFile)

            # If not section is provided we assume the file only has one and it gives us the properties of the observation
            if objName is None:
                objName = cfg.options(cfg.sections()[0])

            # Dictionary with the observation data
            obj_data = dict(cfg.items(objName))
            obj_data['obsFile'] = obsFile
            obj_data['objName'] = objName

            results_section = objName + '_results'
            #Recover data from previous fits
            if cfg.has_section(results_section):
                prefit_data = dict(cfg.items(results_section))
                obj_data.update(prefit_data)

        else:
            # Dictionary with the observation data # TODO This does not work so well
            obj_data = locals()

        # Convert to the right format # TODO Add security warnings for wrong data
        for key in obj_data.keys():

            # if key not in list_parameters + boolean_parameters + string_parameters:

            # Empty variable
            if obj_data[key] == '':
                obj_data[key] = None

            # Arrays (The last boolean overrides the parameters
            elif (key in list_parameters) or ('_prior' in key) or (',' in obj_data[key]):
                if key in ['input_lines']:
                    if obj_data[key] == 'all':
                        obj_data[key] = 'all'
                    else:
                        obj_data[key] = np.array(map(str, obj_data[key].split(',')))
                else:
                    obj_data[key] = np.array(map(float, obj_data[key].split(',')))

            # Boolean
            elif key in boolean_parameters:
                obj_data[key] = strtobool(obj_data[key]) == 1

            # Remaining are either strings (rest floats)
            elif key not in string_parameters:
                obj_data[key] = float(obj_data[key])

            # #Unrecognize object function
            # else:
            #     print 'WARNING: Parameter {} in {} not recognize. Exiting code'.format(key, obsFile)
            #     exit()

        # ----Load the obj spectrum, #TODO read this one using pandas and that way you can chek if there  is a third column for the error
        obj_data['obs_wavelength'], obj_data['obs_flux'] = loadtxt(obj_data['address_spectrum'], usecols=(0, 1), unpack=True)

        # ----Load obj lines log # TODO update code to use address_lines_log
        obj_data['obj_lines_file'] = obj_data['address_lines_log']

        return obj_data

    def load_ssp_library(self, ssp_lib_type, data_folder=None, data_file=None, wavelengh_limits=None, resample_inc=None, norm_interval=None):

        # TODO In here we need to add a test sample library

        # Store stellar base type
        sspLib_dict = {'data_type': ssp_lib_type}

        # Import the base type
        if ssp_lib_type == 'FIT3D':

            # Check if more files are being introduced
            if ',' in data_file:
                ssp_lib1, ssp_lib2 = data_file.split(',')  # Corrently we are only using the first one (the big)
            else:
                ssp_lib1 = data_file

            sspLib_dict = self.import_Fit3D_ssplibrary(data_folder + ssp_lib1)

        elif ssp_lib_type == 'starlight':

            sspLib_dict = self.import_STARLIGHT_ssplibrary(data_folder, data_file)

        # Store stellar base type
        sspLib_dict['data_type'] = ssp_lib_type

        # Trim, resample and normalized the ssp library if required
        if wavelengh_limits or resample_inc or norm_interval:
            self.treat_input_spectrum(sspLib_dict, sspLib_dict['basesWave'], sspLib_dict['fluxBases'], wavelengh_limits,
                                      resample_inc, norm_interval)

        return sspLib_dict

    def treat_input_spectrum(self, output_dict, spec_wave, spec_flux, wavelengh_limits=None, resample_inc=None, norm_interval=None):

        # TODO we should remove the nBases requirement by some style which can just read the number of dimensions

        # Store input values
        output_dict['wavelengh_limits'] = wavelengh_limits
        output_dict['resample_inc'] = resample_inc
        output_dict['norm_interval'] = norm_interval

        # Special case using 0, -1 indexing
        if wavelengh_limits is not None:
            if (wavelengh_limits[0] != 0) and (wavelengh_limits[0] != -1):
                inputWaveLimits = wavelengh_limits
            else:
                inputWaveLimits = wavelengh_limits
                if wavelengh_limits[0] == 0:
                    inputWaveLimits[0] = int(np.ceil(spec_wave[0]) + 1)
                if wavelengh_limits[-1] == -1:
                    inputWaveLimits[-1] = int(np.floor(spec_wave[-1]) - 1)

        # Resampling the spectra
        if resample_inc is not None:
            wave_resam = np.arange(inputWaveLimits[0], inputWaveLimits[-1], resample_inc, dtype=float)

            # Loop throught the fluxes (In the case of the bases it is assumed they may have different wavelength ranges)
            if isinstance(spec_flux, list):

                flux_resam = np.empty((output_dict['nBases'], len(wave_resam)))
                for i in range(output_dict['nBases']):
                    flux_resam[i, :] = interp1d(spec_wave[i], spec_flux[i], bounds_error=True)(wave_resam)

            # In case only one dimension
            elif spec_flux.ndim == 1:
                flux_resam = interp1d(spec_wave, spec_flux, bounds_error=True)(wave_resam)

            output_dict['wave_resam'] = wave_resam
            output_dict['flux_resam'] = flux_resam

        else:
            output_dict['wave_resam'] = spec_wave
            output_dict['flux_resam'] = spec_flux

        # Normalizing the spectra
        if norm_interval is not None:

            # Loop throught the fluxes (In the case of the bases it is assumed they may have different wavelength ranges)
            if isinstance(spec_flux, list):

                normFlux_coeff = np.empty(output_dict['nBases'])
                flux_norm = np.empty((output_dict['nBases'], len(wave_resam)))
                for i in range(output_dict['nBases']):
                    idx_Wavenorm_min, idx_Wavenorm_max = np.searchsorted(spec_wave[i], norm_interval)
                    normFlux_coeff[i] = np.mean(spec_flux[i][idx_Wavenorm_min:idx_Wavenorm_max])
                    flux_norm[i] = output_dict['flux_resam'][i] / normFlux_coeff[i]

            elif spec_flux.ndim == 1:
                idx_Wavenorm_min, idx_Wavenorm_max = np.searchsorted(spec_wave, norm_interval)
                normFlux_coeff = np.mean(spec_flux[idx_Wavenorm_min:idx_Wavenorm_max])
                flux_norm = output_dict['flux_resam'] / normFlux_coeff

            output_dict['flux_norm'] = flux_norm
            output_dict['normFlux_coeff'] = normFlux_coeff

        else:

            output_dict['flux_norm'] = output_dict['flux_resam']
            output_dict['normFlux_coeff'] = 1.0

        return

    def generate_object_mask(self, linesDf, wavelength, linelabels):

        # TODO This will not work for a redshifted lines log
        idcs_lineMasks = linesDf.index.isin(linelabels)
        idcs_spectrumMasks = ~linesDf.index.isin(linelabels)

        # Matrix mask for integring the emission lines
        n_lineMasks = idcs_lineMasks.sum()
        self.boolean_matrix = np.zeros((n_lineMasks, wavelength.size), dtype=bool)

        # Array with line wavelength resolution which we fill with default value (This is because there are lines beyong the continuum range)
        self.lineRes = np.ones(n_lineMasks) * (wavelength[1] - wavelength[0])

        # Total mask for valid regions in the spectrum
        n_objMasks = idcs_spectrumMasks.sum()
        self.int_mask = np.ones(wavelength.size, dtype=bool)
        self.object_mask = np.ones(wavelength.size, dtype=bool)

        # Loop through the emission lines
        wmin, wmax = linesDf['w3'].loc[idcs_lineMasks].values, linesDf['w4'].loc[idcs_lineMasks].values
        idxMin, idxMax = np.searchsorted(wavelength, [wmin, wmax])
        for i in range(n_lineMasks):
            if not np.isnan(wmin[i]) and not np.isnan(wmax[i]) and (wmax[i] < wavelength[-1]): # We need this for lines beyong continuum range #TODO propose better
                w2, w3 = wavelength[idxMin[i]], wavelength[idxMax[i]]
                idx_currentMask = (wavelength >= w2) & (wavelength <= w3)
                self.boolean_matrix[i, :] = idx_currentMask
                self.int_mask = self.int_mask & ~idx_currentMask
                self.lineRes[i] = wavelength[idxMax[i]] - wavelength[idxMax[i] - 1]

        # Loop through the object masks
        wmin, wmax = linesDf['w3'].loc[idcs_spectrumMasks].values, linesDf['w4'].loc[idcs_spectrumMasks].values
        idxMin, idxMax = np.searchsorted(wavelength, [wmin, wmax])
        for i in range(n_objMasks):
            if not np.isnan(wmin[i]) and not np.isnan(wmax[i]) and (wmax[i] < wavelength[-1]):
                w2, w3 = wavelength[idxMin[i]], wavelength[idxMax[i]]
                idx_currentMask = (wavelength >= w2) & (wavelength <= w3)
                self.int_mask = self.int_mask & ~idx_currentMask
                self.object_mask = self.object_mask & ~idx_currentMask

        return
