from os                                     import path, name, makedirs, environ
environ["MKL_THREADING_LAYER"] = "GNU"
import theano
import theano.tensor as tt
import ConfigParser
import pymc as pymc2
import pymc3
from pymc3 import math as pm3_math
from sys                                    import argv, exit
from pyneb                                  import atomicData, RecAtom, Atom
from collections                            import OrderedDict
from scipy.interpolate import interp1d
from lib.Astro_Libraries.Nebular_Continuum  import NebularContinuumCalculator
from lib.Astro_Libraries.Reddening_Corrections import ReddeningLaws
from lib.ssp_functions.ssp_synthesis_tools  import ssp_fitter, ssp_synthesis_importer
from DZ_LineMesurer                         import LineMesurer_v2
from numpy import array, savetxt, loadtxt, genfromtxt, copy, isnan, arange, insert, concatenate, mean, std, power, exp, \
    square, empty, percentile, random, median, ones, isnan, sum as np_sum, transpose, vstack, delete, where, \
    searchsorted, in1d, save as np_save, load as np_load, unique, log10
from numpy.core import defchararray
from timeit import default_timer as timer
from pandas import read_excel, read_csv
from bisect import bisect_left
from plot_tools import MCMC_printer
from distutils.util import strtobool
import matplotlib.pyplot as plt
import cPickle as pickle
from lib.Math_Libraries.functions_for_theano import BilinearInterpTheano

def make_folder(folder_path):

    #TODO This one is only valid for 2.7
    #TODO add this one to your collection
    try:
        makedirs(folder_path)
    except OSError:
        if not path.isdir(folder_path):
            raise

    return

def import_table_data(address, columns):

    Imported_Array = genfromtxt(address, dtype=float, usecols = columns, skip_header = 2).T
    Datarray = Imported_Array[:,~isnan(Imported_Array).all(0)]

    return Datarray

def correct_df_row_name(df, lines_to_cor):

    list_index = df.index.tolist()

    for i in range(len(lines_to_cor)):
        if lines_to_cor[i][0] in list_index:
            idx_i = list_index.index(lines_to_cor[i][0])
            list_index[idx_i] = lines_to_cor[i][1]

    df.index = list_index

    return df

def bilinear_interpolator_axis(x, y, x_range, y_range, data_grid):
    # TODO Check this emissivity interpolation with pyneb

    i = bisect_left(x_range, x) - 1
    j = bisect_left(y_range, y) - 1

    x1, x2 = x_range[i:i + 2]
    y1, y2 = y_range[j:j + 2]

    z11, z12 = data_grid[j][i:i + 2]
    z21, z22 = data_grid[j + 1][i:i + 2]

    inter_value = (z11 * (x2 - x) * (y2 - y) +
                   z21 * (x - x1) * (y2 - y) +
                   z12 * (x2 - x) * (y - y1) +
                   z22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1))


    return inter_value

def generate_emissivity_grid(pynebcode_list, ions_list, ion_dict, temp_interval, den_interval):

    grid_emis = empty((len(temp_interval), len(den_interval), len(pynebcode_list)))

    for i in range(len(pynebcode_list)):
        ion = ions_list[i]
        print('--- Generating grid: {}, {}'.format(ion, pynebcode_list[i]))
        grid_emis[:, :, i] =ion_dict[ion].getEmissivity(temp_interval, den_interval, wave=pynebcode_list[i],product=True)

    return grid_emis

def TOIII_TSIII_relation(TSIII):
    # TODO we should make a class with these physical relations
    return (1.0807 * TSIII / 10000.0 - 0.0846) * 10000.0

def parse_synth_objData(file_address, objName, objData):

    #json.dump(obj_data, open('/home/vital/PycharmProjects/thesis_pipeline/spectrum_fitting/synth_objProperties_json.txt', 'w'),
    # indent=4, sort_keys=True)

    parser = ConfigParser.SafeConfigParser()
    parser.optionxform = str
    parser.add_section(objName)

    for key in objData.keys():
        value = objData[key]
        if value is not None:
            if isinstance(value, list):
                value = ','.join(str(x) for x in value)
            else:
                value = str(value)
        else:
            value = ''

        parser.set(objName, key, value)

    with open(file_address, 'w') as f:
        parser.write(f)

    # parser2 = ConfigParser.SafeConfigParser()
    # parser2.read('/home/vital/PycharmProjects/thesis_pipeline/spectrum_fitting/synth_objProperties_config_parser.txt')
    #
    # print 'Aqui'
    # print parser2.options(parser2.sections()[0])
    # print 'hola', parser2.get('Object-Info', 'flux_hbeta') == ''

    return

class Import_model_data(ssp_synthesis_importer):

    '''
    This class loads the binaries required for the dazer spectra synthesizer module
    '''

    def __init__(self):

        #Class with tools to import starlight bases
        ssp_synthesis_importer.__init__(self)

        self.paths_dict = {}

        # Paths for linux:
        if name == 'posix':
            # TODO This data should be read from a text file
            self.paths_dict['inference_folder']     = '/home/vital/Astrodata/Inference_output/'
            self.paths_dict['nebular_data_folder']  = '/home/vital/Dropbox/Astrophysics/Lore/NebularContinuum/'
            self.paths_dict['Hydrogen_CollCoeff']   = '/home/vital/PycharmProjects/dazer/bin/lib/Astro_Libraries/Neutral_Hydrogen_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_CollCoeff']     = '/home/vital/PycharmProjects/dazer/bin/lib/Astro_Libraries/Neutral_Helium_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_OpticalDepth']  = '/home/vital/PycharmProjects/dazer/bin/lib/Astro_Libraries/Helium_OpticalDepthFunction_Coefficients.txt'
            self.paths_dict['lines_data_file']      = '/home/vital/PycharmProjects/dazer/bin/lib/Astro_Libraries/spectrum_fitting/lines_data.xlsx'
            self.paths_dict['dazer_path']           = '/home/vital/PycharmProjects/dazer/'
            self.paths_dict['stellar_data_folder']  = '/home/vital/Starlight/'
            self.paths_dict['observations_folder']  = '/home/vital/Dropbox/Astrophysics/Data/WHT_observations/objects/'
            self.paths_dict['emisGridsPath']        = '/home/vital/PycharmProjects/thesis_pipeline/spectrum_fitting/testing_output/input_data/'

        # Paths for windows:
        elif name == 'nt':
            # TODO Check guides for windows and linux
            self.paths_dict['inference_folder']     = 'D:/Inference_data/'
            self.paths_dict['nebular_data_folder']  = 'E:/Cloud Storage/Dropbox/Astrophysics/Lore/NebularContinuum/'
            self.paths_dict['Hydrogen_CollCoeff']   = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Neutral_Hydrogen_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_CollCoeff']     = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Neutral_Helium_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_OpticalDepth']  = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Helium_OpticalDepthFunction_Coefficients.txt'
            self.paths_dict['stellar_data_folder']  = 'E:/Cloud Storage/Dropbox/Astrophysics/Tools/Starlight/'
            self.paths_dict['lines_data_file']      = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/spectrum_fitting/lines_data.xlsx'
            self.paths_dict['dazer_path']           = 'C:/Users/lativ/git/dazer/'
            self.paths_dict['observations_folder']  = 'E:/Cloud Storage/Dropbox/Astrophysics/Data/WHT_observations/objects/'
            self.paths_dict['emisGridsPath']        = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/spectrum_fitting/'

        # Import database with lines labels information
        self.lines_df = read_excel(self.paths_dict['lines_data_file'], sheetname=0, header=0, index_col=0)

        # Additional code to adapt to old lines log labelling
        self.lines_df['pyneb_code'] = self.lines_df['pyneb_code'].astype(str)
        idx_numeric_pynebCode = ~self.lines_df['pyneb_code'].str.contains('A+')
        self.lines_df.loc[idx_numeric_pynebCode, 'pyneb_code'] = self.lines_df.loc[idx_numeric_pynebCode, 'pyneb_code'].astype(int)
        self.lines_df['ion'].apply(str)

    def import_optical_depth_coeff_table(self, helium_lines):

        Data_dict = OrderedDict()

        for i in range(len(helium_lines)):

            Line = helium_lines[i]

            # TODO addapt the table to make sure it can read both objects H13889A_label
            Data_dict[Line] = loadtxt(self.paths_dict['Helium_OpticalDepth'], dtype=float, skiprows=2, usecols=(i,))


            # if Line != self.H13889A_label:
            #     Data_dict[Line] = loadtxt(self.paths_dict['Helium_OpticalDepth'], dtype = float, skiprows = 2, usecols = (i,))
            #
            # else:
            #     Data_dict[self.He3889_label] = loadtxt(self.paths_dict['Helium_OpticalDepth'], dtype = float, skiprows = 2, usecols = (i,))

        return Data_dict

    def load_obsData(self, obsFile = None, objName = None, obs_wavelength = None, obs_flux = None):

        # TODO this should go into the master configuration
        # QUESTION Ask about this config parser best strategies
        list_parameters = ['obs_wavelength', 'obs_flux', 'wavelengh_limits', 'norm_interval']
        boolean_parameters = ['Normalized_by_Hbeta']
        string_parameters = ['address_lines_log', 'address_spectrum', 'address_obs_mask', 'obsFile', 'objName']
        print '-Reading observation {} {}'.format(obsFile if obsFile is not None else '', ' for object ' + objName if objName is not None else '')

        #----Load the obj data
        if obsFile is not None:
            cfg = ConfigParser.SafeConfigParser()
            cfg.optionxform = str
            cfg.read(obsFile)

            #If not section is provided we assume the file only has one and it gives us the properties of the observation
            if objName is None:
                objName = cfg.options(cfg.sections()[0])

            # Dictionary with the observation data
            obj_data =  dict(cfg.items(objName))
            obj_data['obsFile'] = obsFile
            obj_data['objName'] = objName

        else:
            # Dictionary with the observation data
            # TODO This does not work so well
            obj_data = locals()

        # Convert to the right format
        for key in obj_data.keys():

            # TODO Add security warnings for wrong data
            #if key not in list_parameters + boolean_parameters + string_parameters:

            #Empty variable
            if obj_data[key] == '':
                obj_data[key] = None

            #Arrays
            elif key in list_parameters:
                obj_data[key] = array(map(float, obj_data[key].split(',')))

            #Boolean
            elif key in boolean_parameters:
                obj_data[key] = strtobool(obj_data[key]) == 1

            # Remaining are either strings (rest floats)
            elif key not in string_parameters:
                obj_data[key] = float(obj_data[key])

            # #Unrecognize object function
            # else:
            #     print 'WARNING: Parameter {} in {} not recognize. Exiting code'.format(key, obsFile)
            #     exit()

        #----Load the obj spectrum
        if obj_data['obs_flux'] is None:
            if obj_data['address_spectrum'] is not None:
                obj_data['obs_wavelength'], obj_data['obs_flux'] = loadtxt(obj_data['address_spectrum'], usecols=(0,1), unpack=True)
        if obj_data['obs_flux'] is None:
            print 'WARNING: Input spectrum not found. Exiting code'
            exit()

        #----Load obj lines log
        #obj_lines_df = read_csv(obj_data['address_lines_log'], delim_whitespace=True, header=0, index_col=0)
        #obj_lines_df.index.name = 'line_label'
        #TODO update code to use address_lines_log
        obj_data['obj_lines_file'] = obj_data['address_lines_log']

        return obj_data

    def import_emission_line_data(self, obj_lines_df, input_lines):

        # If wavelengths are provided for the observation we use them, else we use the theoretical values
        if all(obj_lines_df.obs_wavelength.isnull().values):
            idx_obs_labels = self.lines_df.index.isin(obj_lines_df.index)
            obj_lines_df['obs_wavelength'] = self.lines_df.loc[idx_obs_labels].wavelength

        # Sort the dataframe by wavelength in case it isn't
        obj_lines_df.sort_values(by=['obs_wavelength'], ascending=True, inplace=True)

        # Get the references for the lines treatment
        idx_obj_lines               = self.lines_df.index.isin(obj_lines_df.index)
        obj_lines_df['emis_type']   = self.lines_df.loc[idx_obj_lines].emis_type.astype(str)
        obj_lines_df['pynebCode']   = self.lines_df.loc[idx_obj_lines].pyneb_code
        obj_lines_df['ion']         = self.lines_df.loc[idx_obj_lines].ion.astype(str)

        # Save the data in the dictionary according to recomb or colExt format
        # TODO Code to change if we remove recomb/colExt dividion
        # TODO add fluxes when we store synthetic data to a file

        #Import lines data, except Hbeta #TODO we might need to import Hbeta line individually.
        idx_recomb = (obj_lines_df.emis_type == 'rec') & (obj_lines_df.index.isin(input_lines)) & (obj_lines_df.index != 'H1_4861A')
        self.obj_data['recombLine_ions']        = obj_lines_df.loc[idx_recomb].ion.values
        self.obj_data['recombLine_labes']       = obj_lines_df.loc[idx_recomb].index.values
        self.obj_data['recombLine_waves']       = obj_lines_df.loc[idx_recomb].obs_wavelength.values
        self.obj_data['recombLine_pynebCode']   = obj_lines_df.loc[idx_recomb].pynebCode.values
        self.obj_data['recomb_fluxes']          = obj_lines_df.loc[idx_recomb].obs_flux.values
        self.obj_data['recomb_err']             = obj_lines_df.loc[idx_recomb].obs_fluxErr.values

        idx_col = (obj_lines_df.emis_type == 'col') & (obj_lines_df.index.isin(input_lines))
        self.obj_data['colLine_ions']           = obj_lines_df.loc[idx_col].ion.values
        self.obj_data['colLine_labes']          = obj_lines_df.loc[idx_col].index.values
        self.obj_data['colLine_waves']          = obj_lines_df.loc[idx_col].obs_wavelength.values
        self.obj_data['colLine_pynebCode']      = obj_lines_df.loc[idx_col].pynebCode.values
        self.obj_data['metal_fluxes']           = obj_lines_df.loc[idx_col].obs_flux.values
        self.obj_data['metal_err']              = obj_lines_df.loc[idx_col].obs_fluxErr.values

        # Generate the dictionary with pyneb ions
        self.ionDict = {}
        for raw_idx in obj_lines_df.index:
            ion = obj_lines_df.loc[raw_idx].ion
            if ion not in self.ionDict:
                element, state = ion[:-1], ion[-1]
                if obj_lines_df.loc[raw_idx].emis_type == 'rec':
                    self.ionDict[ion] = RecAtom(element, state)
                elif obj_lines_df.loc[raw_idx].emis_type == 'col':
                    self.ionDict[ion] = Atom(element, state)

        # Establish index of lines which below to high and low ionization zones
        # TODO This should also include recomb (unnecesary if we  we remove recomb/colExt dividion)
        self.idx_highU = in1d(self.obj_data['colLine_ions'], self.high_temp_ions)
        self.idx_lowU = ~self.idx_highU

        # Parameters to improve calculation speed:
        self.n_recombLines, self.n_colExcLines = len(self.obj_data['recombLine_labes']), len(self.obj_data['colLine_labes'])
        self.range_recombLines, self.range_colExcLines = arange(self.n_recombLines), arange(self.n_colExcLines)

        # Empty Ionic dictionary for the MCMC
        self.abund_dict = {ion: 0 for ion in self.ionDict.keys()}

        return

    def import_emissivity_grids(self, forceGridsReset = False):

        # TODO Improve the reseting mechanic
        # QUESTION Ask how to
        # Emissivity grid for recombination lines
        self.recomb_emis_grid = self.manage_emissivity_grids(self.obj_data['recombLine_pynebCode'], self.obj_data['recombLine_ions'], forceGridsReset)

        # Emissivity grid for collisional excited lines
        self.metals_emis_grid = self.manage_emissivity_grids(self.obj_data['colLine_pynebCode'], self.obj_data['colLine_ions'], forceGridsReset)

        # Emissivity grid for collisional excited lines
        self.Hbeta_emis_grid = self.manage_emissivity_grids([4861], ['H1'], forceGridsReset)

        #Normalizing the grids by Hbeta
        self.recomb_emis_grid = self.recomb_emis_grid / self.Hbeta_emis_grid
        self.metals_emis_grid = self.metals_emis_grid / self.Hbeta_emis_grid

        self.log_recomb_emis_grid = log10(self.recomb_emis_grid)
        self.log_metals_emis_grid = log10(self.metals_emis_grid)

        return

    def manage_emissivity_grids(self, pynebcode_list, ions_list, forceGridReset = False):

        # Empty grid
        emisGrid = empty((len(self.tem_grid_range), len(self.den_grid_range), len(ions_list)))

        for i in range(len(ions_list)):

            # Line emissivity references
            line_label = '{}_{}'.format(ions_list[i], pynebcode_list[i])
            grid_address = self.paths_dict['emisGridsPath'] + line_label + '.npy'

            # Check if grids are available
            if path.isfile(grid_address) and forceGridReset is False:
                emis_grid_i = np_load(grid_address)

            # Otherwise generate it (and save it)
            else:
                print('--- Generating grid: {}'.format(line_label))
                emis_grid_i = self.ionDict[ions_list[i]].getEmissivity(self.tem_grid_range, self.den_grid_range,
                                                                       wave=pynebcode_list[i], product=True)
                np_save(grid_address, emis_grid_i)

            # Save along the number of points
            emisGrid[:, :, i] = emis_grid_i

        return emisGrid

    def load_ssp_library(self, ssp_lib_type, data_folder=None, data_file=None, wavelengh_limits=None, resample_inc=None, norm_interval=None):

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

        #Trim, resample and normalized the ssp library if required
        if wavelengh_limits or resample_inc or norm_interval:
            self.treat_input_spectrum(sspLib_dict, sspLib_dict['basesWave'], sspLib_dict['fluxBases'], wavelengh_limits, resample_inc, norm_interval)

        return sspLib_dict

    def treat_input_spectrum(self, output_dict, spec_wave, spec_flux, wavelengh_limits = None, resample_inc = None, norm_interval = None):

        #TODO we should remove the nBases requirement by some style which can just read the number of dimensions

        #Store input values
        output_dict['wavelengh_limits']    = wavelengh_limits
        output_dict['resample_inc']        = resample_inc
        output_dict['norm_interval']       = norm_interval

        #Resampling the spectra
        if resample_inc is not None:
            wave_resam = arange(wavelengh_limits[0], wavelengh_limits[-1], resample_inc, dtype=float)

            #Loop throught the fluxes (In the case of the bases it is assumed they may have different wavelength ranges)
            if isinstance(spec_flux, list):

                flux_resam = empty((output_dict['nBases'], len(wave_resam)))
                for i in range(output_dict['nBases']):
                    flux_resam[i, :] = interp1d(spec_wave[i], spec_flux[i], bounds_error=True)(wave_resam)

            #In case only one dimension
            elif spec_flux.ndim == 1:
                flux_resam = interp1d(spec_wave, spec_flux, bounds_error=True)(wave_resam)

            output_dict['wave_resam'] = wave_resam
            output_dict['flux_resam'] = flux_resam

        else:
            output_dict['wave_resam'] = spec_wave
            output_dict['flux_resam'] = spec_flux

        #Normalizing the spectra
        if norm_interval is not None:

            #Loop throught the fluxes (In the case of the bases it is assumed they may have different wavelength ranges)
            if isinstance(spec_flux, list):

                normFlux_coeff  = empty(output_dict['nBases'])
                flux_norm       = empty((output_dict['nBases'], len(wave_resam)))
                for i in range(output_dict['nBases']):
                    idx_Wavenorm_min, idx_Wavenorm_max = searchsorted(spec_wave[i], norm_interval)
                    normFlux_coeff[i]   = mean(spec_flux[i][idx_Wavenorm_min:idx_Wavenorm_max])
                    flux_norm[i]        = output_dict['flux_resam'][i] / normFlux_coeff[i]

            elif spec_flux.ndim == 1:
                idx_Wavenorm_min, idx_Wavenorm_max = searchsorted(spec_wave, norm_interval)
                normFlux_coeff          = mean(spec_flux[idx_Wavenorm_min:idx_Wavenorm_max])
                flux_norm               = output_dict['flux_resam'] / normFlux_coeff

            output_dict['flux_norm'] = flux_norm
            output_dict['normFlux_coeff'] = normFlux_coeff

        else:

            output_dict['flux_norm'] = output_dict['flux_resam']
            output_dict['normFlux_coeff'] = 1.0

        return

    def generate_object_mask(self):

        # TODO We need to adapt this from the lick indeces system to a mask_label wmin wmax (or viceversa)
        # TODO If not mask is included generate one from emission lines
        # TODO if not mask is applied mask = all True
        obj_mask_df = read_csv(self.obj_data['obs_mask_address'], delim_whitespace=True, header=0, index_col=0)

        boolean_mask = ones(len(self.obj_data['wave_resam']), dtype=bool)

        for mask_label in obj_mask_df.index:
            wmin, wmax = obj_mask_df.loc[mask_label, 'wmin'], obj_mask_df.loc[mask_label, 'wmax']
            idx_cur_spec_mask = (self.obj_data['wave_resam'] > wmin) & (self.obj_data['wave_resam'] < wmax)
            boolean_mask = boolean_mask & ~idx_cur_spec_mask

        self.obj_data['boolean_mask'] = boolean_mask
        self.obj_data['int_mask']     = boolean_mask * 1

        self.obj_data['flux_norm_mask'] = self.obj_data['flux_norm'] * self.obj_data['int_mask']

        return

class ContinuaComponents(ssp_fitter, NebularContinuumCalculator):

    def __init__(self):

        #Load tools for continua calculation
        ssp_fitter.__init__(self)
        NebularContinuumCalculator.__init__(self)

    def nebular_Cont(self, wave_obs, z, cHbeta, Te, He1_abund, He2_abund, Halpha_Flux):

        wave_obs_rest   = wave_obs / (1.0 + z)

        neb_gCont       = self.calculate_neb_gCont(wave_obs_rest, Te, He1_abund, He2_abund)

        neb_int_norm    = self.gCont_calibration(wave_obs_rest, Te, Halpha_Flux, neb_gCont)

        neb_xX          = self.reddening_Xx(wave_obs_rest, 'G03_average', 3.4)
        flambda_neb     = neb_xX/self.Hbeta_xX - 1.0

        neb_flux_norm   = neb_int_norm * power(10, -1 * flambda_neb * cHbeta)

        return neb_flux_norm

    def calculate_nebular_SED(self, wave_obs, z, cHbeta, Te, He1_abund, He2_abund, Halpha_Flux):

        wave_obs_rest   = wave_obs / (1.0 + z)

        neb_gCont       = self.calculate_neb_gCont(wave_obs_rest, Te, He1_abund, He2_abund)

        neb_int_norm    = self.gCont_calibration(wave_obs_rest, Te, Halpha_Flux, neb_gCont)

        neb_xX          = self.reddening_Xx(wave_obs_rest, 'G03_average', 3.4) #TODO Change by global gas reddening law
        flambda_neb     = neb_xX/self.Hbeta_xX - 1.0

        neb_flux_norm   = neb_int_norm * power(10, -1 * flambda_neb * cHbeta)

        nebular_SED = {'neb_gCont': neb_gCont, 'neb_int_norm': neb_int_norm, 'neb_xX': neb_xX,
                       'flambda_neb': flambda_neb, 'neb_flux_norm': neb_flux_norm}

        return nebular_SED

class EmissionComponents(LineMesurer_v2):

    def __init__(self):

        #Load tools to measure lines
        LineMesurer_v2.__init__(self,  self.paths_dict['dazer_path'] + 'format/', 'DZT_LineLog_Headers.dz')

        #Hbeta configuration
        self.Hbeta_label            = 'H1_4861A'
        self.Hbeta_wave             = 4862.683
        self.Hbeta_pynebCode        = '4_2'

        #Import Optical depth function
        posHelium_Lines             = ['He1_3889A',  'He1_4026A',  'He1_4387A',    'He1_4471A',    'He1_4686A',    'He1_4714A',    'He1_4922A',   'He1_5876A',    'He1_6678A',   'He1_7065A',    'He1_7281A',    'He1_10830A']
        self.Coef_ftau_dict         = self.import_optical_depth_coeff_table(posHelium_Lines)

    def OpticalDepth_He(self, tau, T_4, ne, He_label):

        f_tau   = 1 + (tau/2) * (self.Coef_ftau_dict[He_label][0] + (self.Coef_ftau_dict[He_label][1] + self.Coef_ftau_dict[He_label][2]*ne + self.Coef_ftau_dict[He_label][3]*ne*ne) * T_4)

        return f_tau

    def calculate_recomb_fluxes(self, Thigh, ne, cHbeta, tau, lines_abund_dict, lines_labels, lines_ions, lines_flambda):

        #Emissivity for the two temperature layers
        t4 = Thigh / 10000.0

        #Interpolate the emissivities
        emis_module = bilinear_interpolator_axis(ne, Thigh, self.den_grid_range, self.tem_grid_range, self.recomb_emis_grid)

        #Loop throught the lines to generate the f_tau correction
        lines_ftau_vector = empty(self.n_recombLines)
        lines_abund_vector = empty(self.n_recombLines)

        for i in self.range_recombLines:
            if lines_ions[i] == 'He1':
                lines_ftau_vector[i]    = self.OpticalDepth_He(tau = tau, T_4 = t4, ne = ne, He_label = lines_labels[i])
            else:
                lines_ftau_vector[i]    = 1.0
            lines_abund_vector[i] = lines_abund_dict[lines_ions[i]]

        #Reddening factor for all the lines
        f_module = power(10, -1 * lines_flambda * cHbeta)

        #Metals flux
        recomb_fluxes = lines_abund_vector * emis_module * lines_ftau_vector * f_module

        return recomb_fluxes

    def calculate_colExcit_flux(self, Tlow, Thigh, ne, cHbeta, lines_abund_dict, lines_waves, lines_ions, lines_flambda):

        # Array to store the emissivities
        lines_emis_vector = empty(self.n_colExcLines)

        #Calculate low ionization emissivities
        lines_emis_vector[self.idx_lowU] = bilinear_interpolator_axis(ne, Tlow, self.den_grid_range,
                                                                      self.tem_grid_range,
                                                                      self.metals_emis_grid[:, :, self.idx_lowU])

        #Calculate high ionization emissivities
        lines_emis_vector[self.idx_highU] = bilinear_interpolator_axis(ne, Thigh, self.den_grid_range,
                                                                       self.tem_grid_range,
                                                                      self.metals_emis_grid[:, :, self.idx_highU])

        # Loop through the lines to assign the ionic abundances
        abund_vector = empty(self.n_colExcLines)
        for i in self.range_colExcLines:
            ion = lines_ions[i]
            abund_vector[i] = lines_abund_dict[ion]
            #print ion, abund_vector[i], lines_emis_vector[i], abund_vector[i] * lines_emis_vector[i]

        # Reddening factor for all the lines
        f_module = power(10, -1 * lines_flambda * cHbeta)

        # Metals flux
        colExcit_fluxes = abund_vector * lines_emis_vector * f_module

        return colExcit_fluxes

    def reddening_parameters(self):

        # TODO a class should be created to include reddening features for the function
        # Hbeta parameters #TODO This could go into the configuration file (change the normalizing wavelength?)
        self.Hbeta_xX   = self.reddening_Xx(array([self.Hbeta_wave]), self.reddedning_curve_model, self.Rv_model)[0]

        recomblines_xX  = self.reddening_Xx(self.obj_data['recombLine_waves'], self.reddedning_curve_model, self.Rv_model)
        colLines_xX     = self.reddening_Xx(self.obj_data['colLine_waves'], self.reddedning_curve_model, self.Rv_model)

        self.obj_data['recombLine_flambda'] = recomblines_xX / self.Hbeta_xX - 1.0
        self.obj_data['colLine_flambda'] = colLines_xX / self.Hbeta_xX - 1.0

        return

class SpectrumFitter(Import_model_data, ContinuaComponents, EmissionComponents, ReddeningLaws, MCMC_printer):

    def __init__(self):

        # QUESTION How to collapse all in class

        #Import tools to load data
        Import_model_data.__init__(self)
        ContinuaComponents.__init__(self)
        EmissionComponents.__init__(self)

        #Import sub classes
        ReddeningLaws.__init__(self)

        #For generating graphs
        MCMC_printer.__init__(self)

    def calc_emis_spectrum(self, wavelength_range, lines_waves, lines_fluxes, lines_mu, lines_sigma, redshift):

        #TODO This should go into the cinematics class

        lines_mu = lines_mu * (1 + redshift)

        A_lines = lines_fluxes / (lines_sigma * self.sqrt2pi)

        wave_matrix = wavelength_range * ones((lines_waves.shape[0], wavelength_range.shape[0]))

        emis_spec_matrix = A_lines[:, None] * exp(-(wave_matrix - lines_mu[:, None]) * (wave_matrix - lines_mu[:, None]) / (2 * lines_sigma * lines_sigma))

        combined_spectrum = emis_spec_matrix.sum(axis=0)

        return combined_spectrum

    def gen_synth_obs(self, spectra_components = None, input_ions = None,
                      wavelengh_limits = None, resample_inc = None, norm_interval = None,
                      obj_properties_file = None, obj_lines_file = None, obj_mask_file = None, output_folder = None,
                      ssp_lib_type = None, data_folder = None, data_file = None ,obj_ssp_coeffs_file = None,
                      error_stellarContinuum = None, error_recomb_lines = None, error_collexc_lines = None):

        # TODO we will need to change to recomb scheme: He1r hence, ion_dict = pn.getAtomDict(atom_list = atoms_list)
        # Dictionary to store the data
        self.obj_data = {}

        # Store input files:
        self.obj_data['obj_properties_file']   = obj_properties_file
        self.obj_data['obj_lines_file']        = obj_lines_file
        self.obj_data['obj_ssp_coeffs_file']   = obj_ssp_coeffs_file
        self.obj_data['obs_mask_address']      = obj_mask_file
        self.obj_data['output_folder']         = output_folder

        # Read simulation data from file
        obj_prop_df = read_csv(obj_properties_file, delim_whitespace=True, header=0, index_col=0)

        #--------Lines calculation
        if 'emission' in spectra_components:

            #Store input data
            self.obj_data['flux_hbeta'] = obj_prop_df.loc['flux_hbeta'][0]

            #Read lines file
            obj_lines_df = read_csv(obj_lines_file, delim_whitespace=True, header=0, index_col=0)

            # Prepare data from emission line file (trick to put all the lines)
            self.import_emission_line_data(obj_lines_df, obj_lines_df.index.values)

            # Create or load emissivity grids
            self.import_emissivity_grids()

            # Reddening parameters
            self.reddening_parameters()

            #Dictionary with synthetic abundances
            abund_df        = obj_prop_df[obj_prop_df.index.str.contains('_abund')]
            abund_keys      = abund_df.index.str.rstrip('_abund').astype(str).values
            abund_values    = abund_df.variable_magnitude.values
            self.abund_dict = dict(zip(abund_keys, abund_values.T))

            # Import the physical data from the log
            for param in obj_prop_df.index:
                self.obj_data[param] = obj_prop_df.loc[param][0]

            #Calculate T_high assuming we get T_low
            self.obj_data['T_high'] = TOIII_TSIII_relation(self.obj_data['T_low'])


            # Prepare the lines data
            self.obj_data['recomb_fluxes'] = self.calculate_recomb_fluxes(self.obj_data['T_high'], obj_prop_df.loc['n_e'][0],
                                                            obj_prop_df.loc['cHbeta'][0], obj_prop_df.loc['tau'][0],
                                                            self.abund_dict,
                                                            self.obj_data['recombLine_labes'],
                                                            self.obj_data['recombLine_ions'],
                                                            self.obj_data['recombLine_flambda'])

            self.obj_data['metal_fluxes'] = self.calculate_colExcit_flux(self.obj_data['T_low'], self.obj_data['T_high'],
                                                                obj_prop_df.loc['n_e'][0], obj_prop_df.loc['cHbeta'][0],
                                                                self.abund_dict,
                                                                self.obj_data['colLine_waves'],
                                                                self.obj_data['colLine_ions'],
                                                                self.obj_data['colLine_flambda'])

            #Use general error if this is provided
            if error_recomb_lines is not None:
                self.obj_data['recomb_err']= self.obj_data['recomb_fluxes'] * error_recomb_lines
            if error_collexc_lines is not None:
                self.obj_data['metal_err'] = self.obj_data['metal_fluxes'] * error_collexc_lines

        #--------Nebular calculation
        if 'nebular' in spectra_components:

            #Save input conditions:
            self.obj_data['wavelengh_limits']  = wavelengh_limits
            self.obj_data['resample_inc']      = resample_inc
            self.obj_data['norm_interval']     = norm_interval
            self.obj_data['z_star']            = obj_prop_df.loc['z_star'][0]

            # Rest stellar and observed wavelength
            obj_wave_rest = arange(wavelengh_limits[0], wavelengh_limits[-1], resample_inc, dtype=float)
            obj_wave_obs = obj_wave_rest * (1.0 + self.obj_data['z_star']) # TODO we should change this to z_obj

            #Get Halpha flux to calibrate
            idx_Halpha = (self.obj_data['recombLine_labes'] == 'H1_6563A')
            self.obj_data['flux_halpha'] = self.obj_data['recomb_fluxes'][idx_Halpha] * obj_prop_df.loc['flux_hbeta'][0] #/ self.stellar_SED['normFlux_stellar']

            #Calculate the nebular continua
            self.obj_data['neb_SED'] = self.calculate_nebular_SED(obj_wave_rest, # TODO We must change this to the function used on the final code
                                                     obj_prop_df.loc['z_star'][0],
                                                     obj_prop_df.loc['cHbeta'][0],
                                                     obj_prop_df.loc['T_low'][0],
                                                     obj_prop_df.loc['He1_abund'][0],
                                                     obj_prop_df.loc['He2_abund'][0],
                                                    self.obj_data['flux_halpha'])

            self.obj_data['obs_wave_rest'] = obj_wave_rest
            self.obj_data['obs_wave']      = obj_wave_obs
            self.obj_data['obs_flux']      = self.obj_data['neb_SED']['neb_int_norm']

        #--------Stellar calculation
        if 'stellar' in spectra_components:

            # Save input conditions:
            self.obj_data['wavelengh_limits']  = wavelengh_limits
            self.obj_data['resample_inc']      = resample_inc
            self.obj_data['norm_interval']     = norm_interval
            self.obj_data['z_star']            = obj_prop_df.loc['z_star'][0]
            self.obj_data['Av_star']           = obj_prop_df.loc['Av_star'][0]
            self.obj_data['sigma_star']        = obj_prop_df.loc['sigma_star'][0]
            self.obj_data['flux_hbeta']        = obj_prop_df.loc['flux_hbeta'][0]
            self.obj_data['eqw_hbeta']         = obj_prop_df.loc['eqw_hbeta'][0]

            # Import the stellar library
            self.ssp_lib = self.load_ssp_library(ssp_lib_type, data_folder, data_file, wavelengh_limits, resample_inc, norm_interval)

            # Trim bases flux to include only the stellar spectra contributing in our object
            zeros_array = ones(10)
            self.prepare_ssp_data(self.obj_data['obj_ssp_coeffs_file'], obj_wave = ones(10), obj_flux = ones(10),
                                  obj_flux_err = ones(10), obj_mask = ones(10))

            # Rest stellar and observed wavelengths
            obj_wave_rest = arange(wavelengh_limits[0], wavelengh_limits[-1], resample_inc, dtype=float)
            obj_wave_obs = obj_wave_rest * (1.0 + self.obj_data['z_star'])

            # SSP library flux at modelled Av, z_star and sigma_star
            ssp_grid_model_norm = self.physical_SED_model(self.ssp_lib['wave_resam'], obj_wave_obs, self.onBasesFluxNorm,
                                    self.obj_data['Av_star'], self.obj_data['z_star'], self.obj_data['sigma_star'], Rv_coeff=self.Rv_model)

            #Normalized object flux
            obj_flux_norm = ssp_grid_model_norm.dot(self.sspPrefit_Coeffs)

            #Instead of denormalizing the grid we use the Hbeta flux and equivalent width
            #To get a better shape of the nebular and stellar continua
            cont_hbeta  = self.obj_data['flux_hbeta'] / self.obj_data['eqw_hbeta']
            obj_flux = obj_flux_norm * cont_hbeta

            # Generate synth error
            sigma_err   = error_stellarContinuum * median(obj_flux) # TODO We should play with this error
            stellar_err = random.normal(0.0, sigma_err, len(obj_flux_norm))

            # Store the data vector
            self.obj_data['obs_wave_rest']         = obj_wave_rest
            self.obj_data['obs_wave']              = obj_wave_obs
            self.obj_data['stellar_flux']          = obj_flux
            self.obj_data['stellar_flux_err']      = obj_flux + stellar_err
            self.obj_data['sigma_continuum']       = 0.02   #TODO put this in the object data

            #Save synthetic continua according to the observations
            if 'obs_flux' not in self.obj_data:
                self.obj_data['obs_flux'] = self.obj_data['stellar_flux_err']
            else:
                self.obj_data['obs_flux'] = self.obj_data['obs_flux'] + self.obj_data['stellar_flux_err']

        #--------Store data to synthetic observation files
        idx_recomb  = (obj_lines_df.emis_type == 'rec') & (obj_lines_df.index != 'H1_4861A')
        idx_col     = (obj_lines_df.emis_type == 'col')

        obj_lines_df.loc[idx_recomb, 'obs_flux']    = self.obj_data['recomb_fluxes']
        obj_lines_df.loc[idx_recomb, 'obs_fluxErr'] = self.obj_data['recomb_err']

        obj_lines_df.loc[idx_col, 'obs_flux']       = self.obj_data['metal_fluxes']
        obj_lines_df.loc[idx_col, 'obs_fluxErr']    = self.obj_data['metal_err']

        obj_lines_df.loc['H1_4861A', 'obs_flux']    = 1.0
        obj_lines_df.loc['H1_4861A', 'obs_fluxErr'] = 1.0 * error_recomb_lines

        #Store synthetic lines log
        #TODO Add this function to your collection
        synth_lines_log = obj_lines_file.replace('.txt', '_new.txt')
        with open(synth_lines_log, 'w') as f:
            f.write(obj_lines_df.ix[:,:'blended_label'].to_string(float_format = lambda x: "{:15.8f}".format(x), index = True, index_names=False))

        #Store the spectrum as a text file
        synth_spectrum_address = '{}/{}'.format(path.split(obj_lines_file)[0], 'synth_spectrum.txt')
        savetxt(synth_spectrum_address, transpose([self.obj_data['obs_wave'], self.obj_data['obs_flux']]))

        #Store synthetic object data log
        #TODO Add here the comments of the function
        obj_dict = OrderedDict()
        obj_dict['address_lines_log']   = synth_lines_log
        obj_dict['address_spectrum']    = synth_spectrum_address
        obj_dict['address_obs_mask']    = self.obj_data['obs_mask_address']
        obj_dict['obs_wavelength']      = None
        obj_dict['obs_flux']            = None
        obj_dict['Normalized_by_Hbeta'] = True
        obj_dict['flux_hbeta']          = self.obj_data['flux_hbeta']
        obj_dict['flux_halpha']         = ''
        obj_dict['eqw_hbeta']           = self.obj_data['eqw_hbeta']
        obj_dict['T_low']               = self.obj_data['T_low']
        obj_dict['T_high']              = self.obj_data['T_high']
        obj_dict['cHbeta']              = None
        obj_dict['sigma_gas']           = ''
        obj_dict['z_star']              = 0.0
        obj_dict['continuum_sigma']     = 0.02
        obj_dict['resample_inc']        = 1
        obj_dict['wavelengh_limits']    = [4000, 6900]
        obj_dict['norm_interval']       = [5100, 5150]

        #Save the data into an "ini" configuration file
        conf_address = '{}/{}'.format(path.split(obj_lines_file)[0], 'synth_obj.txt')
        parse_synth_objData(conf_address, 'synth_obj', obj_dict)

        #return self.obj_data.copy()
        return

    def ready_simulation(self, output_folder, obs_data, ssp_data, fitting_components, overwrite_grids=False,
                         input_lines=None, wavelengh_limits = None, resample_inc=None, norm_interval=None):

        # Dictionary to store the data
        self.obj_data = obs_data.copy()

        self.fitting_components = fitting_components #TODO I could change this one locaation

        if output_folder is None: # TODO need to establish a global path to save the data
            self.output_folder = self.obj_data['output_folder']

        # Prepare emission data
        if 'emission' in fitting_components:

            obj_lines_df = read_csv(obs_data['obj_lines_file'], delim_whitespace=True, header=0, index_col=0)

            #Prepare data from emission line file #TODO this could be cleaner
            self.import_emission_line_data(obj_lines_df, input_lines = input_lines)

            #Create or load emissivity grids
            self.import_emissivity_grids()

            # Reddening parameters
            self.reddening_parameters()

        # Prepare stellar data
        if 'stellar' in fitting_components:

            # Create dictionary with the stellar library configuration which this object will use
            self.ssp_lib = ssp_data.copy()

            # Trim, resample and normalize according to the input object
            self.treat_input_spectrum(self.obj_data, self.obj_data['obs_wave'], self.obj_data['obs_flux'], wavelengh_limits, resample_inc, norm_interval)

            # Generate object masks:
            self.generate_object_mask()

        #Special operations for synthetic data
        if 'T_low' in self.obj_data:
            self.obj_data['T_He'] = TOIII_TSIII_relation(self.obj_data['T_low'])

        return

    def prepare_gas_data(self, recomb_fluxes, recomb_err, colExc_fluxes, colExc_err, TSIII, TSIII_err, nSII, nSII_err, cHbeta, cHbeta_err):

        # QUESTION pycharm give me all the arguments in a function
        # QUESTION faster documentation

        self.abund_iter_dict    = {}

        # Emission line fluxes and errors
        self.recomb_fluxes      = recomb_fluxes
        self.recomb_err         = recomb_err
        self.colExc_fluxes      = colExc_fluxes
        self.colExc_err         = colExc_err

        self.log_colExc_fluxes  = log10(colExc_fluxes)
        self.log_colExc_Err     = log10(colExc_err)**-2


        # Physical parameters
        self.TSIII              = TSIII
        self.TSIII_err          = TSIII_err
        self.nSII               = nSII
        self.nSII_err           = nSII_err
        self.cHbeta             = cHbeta
        self.cHbeta_err         = cHbeta_err

        #Lines data
        if 'nebular' in self.fitting_components:
            self.f_HalphaNorm = self.obj_data['flux_halpha'] / self.obj_data['normFlux_coeff']

        #Determine abundaces present
        abund_list = []
        if 'colLine_ions' in self.obj_data:
            abund_list += list(unique(self.obj_data['colLine_ions']))
        if 'recombLine_ions' in self.obj_data:
            recomb_exc_list = list(unique(self.obj_data['recombLine_ions']))
            if 'H1' in recomb_exc_list:
                recomb_exc_list.remove('H1')
            abund_list += recomb_exc_list

        abund_list = defchararray.add(array(abund_list), '_abund')

        self.abund_list = list(abund_list)

        # Vector temperature regions keys
        self.tempRegionRecomb = array(['T_high']*self.n_recombLines)
        self.tempRegionCol    = array(['T_low ']*self.n_colExcLines)

        #Assign values
        self.tempRegionCol[self.idx_highU] = 'T_high'

        return

    def prepare_ssp_data(self, input_data, obj_wave, obj_flux, obj_flux_err, obj_mask):

        # QUESTION Show me keys in dictionary

        #Case of a prefit run
        if input_data is None:

            #Default parameters for the nebular continuum
            Te_neb, z_neb, cHbeta_neb = self.conf_dic['Te_neb'], self.conf_dic['z_neb'], self.conf_dic['cHbeta_neb']
            He1_neb, He2_neb = self.conf_dic['He1_neb'], self.conf_dic['He2_neb']

            #Generate synthetic nebular emission to remove from object
            self.obj_data['synth_neb_flux'] = self.nebular_Cont(obj_wave, z_neb, cHbeta_neb, Te_neb,
                                            He1_neb, He2_neb, self.obj_data['flux_halpha'] / self.obj_data['normFlux_coeff'])

            object_continuum = (obj_flux - self.obj_data['synth_neb_flux']) * obj_mask

            idx_populations = ones(self.ssp_lib['flux_resam'].shape[0], dtype=bool)

        else:

            bases_idx, bases_coeff = loadtxt(input_data, usecols=[0, 1], unpack=True)

            idx_populations = bases_coeff >= self.lowlimit_sspContribution # TODO we should check this limit thing currentl it must be 0.001

            object_continuum = obj_flux * obj_mask

            self.sspPrefit_Idcs = where(idx_populations)[0]
            self.sspPrefit_Coeffs = bases_coeff[idx_populations]
            self.sspPrefit_Limits = vstack((self.sspPrefit_Coeffs * 0.8, self.sspPrefit_Coeffs * 1.2)).T

        #Object data
        self.input_continuum        = object_continuum
        self.input_continuum_er     = obj_flux_err
        self.input_wave             = obj_wave
        self.int_mask               = obj_mask

        #Bases parameters
        neglected_populations       = where(~idx_populations)
        self.onBasesWave            = self.ssp_lib['wave_resam']
        self.onBasesFlux            = delete(self.ssp_lib['flux_resam'], neglected_populations, axis=0)
        self.onBasesFluxNorm        = delete(self.ssp_lib['flux_norm'], neglected_populations, axis=0)
        self.onBasesFluxNormCoeffs  = delete(self.ssp_lib['normFlux_coeff'], neglected_populations, axis=0)
        self.range_bases            = arange(self.onBasesFlux.shape[0])

        # Limit for bases
        z_max_ssp                   = (self.ssp_lib['wave_resam'][0] / obj_wave[0]) - 1.0
        z_min_ssp                   = (self.ssp_lib['wave_resam'][-1] / obj_wave[-1]) - 1.0
        self.zMin_SspLimit          = z_min_ssp
        self.zMax_SspLimit          = round(z_max_ssp - 0.001, 3)
        self.z_object               = self.obj_data['z_star']

        return

    def save_prefit_data(self, obj_wave):

        # Read input data # TODO this should not be necessary, we should get it in rum_pymc2
        stellarPrefit_db, stat_db_dict      = self.load_pymc_database_manual(self.prefit_db, burning=1000,
                                             params_list=['Av_star', 'sigma_star', 'ssp_coefficients'])

        # Default parameters for the nebular continuum #TODO we need to change this by the actual prefit data the redshift here
        Te_neb, z_neb, cHbeta_neb           = self.conf_dic['Te_neb'], self.conf_dic['z_neb'], self.conf_dic['cHbeta_neb']
        He1_neb, He2_neb                    = self.conf_dic['He1_neb'], self.conf_dic['He2_neb']

        # Generate synthetic nebular emission to remove from object
        self.ssp_lib['synth_neb_flux']      = self.nebular_Cont(obj_wave, z_neb, cHbeta_neb, Te_neb, He1_neb, He2_neb,
                                                    self.obj_data['flux_halpha'] / self.obj_data['normFlux_coeff'])

        # Save population coefficients      # TODO Default files names should go into the configuration
        sspCoeffs_file                      = self.input_folder + 'prefit_populations.txt'
        coeffs_vector                       = stat_db_dict['ssp_coefficients']['trace'].mean(axis=0)
        pops_vector                         = arange(len(coeffs_vector), dtype=int)
        savetxt(sspCoeffs_file, transpose(array([pops_vector, coeffs_vector])), fmt="%4i %10.8f")

        #Store data
        self.ssp_lib['db_sspPrefit']        = stellarPrefit_db
        self.ssp_lib['dbDict_sspPrefit']    = stat_db_dict

        # Mean parameter values #TODO we need to add the redshift here
        self.sspPrefit_Coeffs               = coeffs_vector
        self.ssp_lib['Av_sspPrefit']        = stat_db_dict['Av_star']['mean']
        self.ssp_lib['sigma_sspPrefit']     = stat_db_dict['sigma_star']['mean']

    def load_prefit_data(self, obj_wave):

        # TODO are an extension/code/folder for each run
        #-----Observation versus prefit comparison # TODO do we need to load this
        stellarPrefit_db, stat_db_dict      = self.load_pymc_database_manual(self.prefit_db, burning=1000, #TODO Add global burning
                                              params_list=['Av_star', 'sigma_star', 'ssp_coefficients'])

        # Default parameters for the nebular continuum
        Te_neb, z_neb, cHbeta_neb           = self.conf_dic['Te_neb'], self.conf_dic['z_neb'], self.conf_dic['cHbeta_neb']
        He1_neb, He2_neb                    = self.conf_dic['He1_neb'], self.conf_dic['He2_neb']

        # Generate synthetic nebular emission to remove from object
        self.ssp_lib['synth_neb_flux']      = self.nebular_Cont(obj_wave, z_neb, cHbeta_neb, Te_neb, He1_neb, He2_neb,
                                                    self.obj_data['flux_halpha'] / self.obj_data['normFlux_coeff'])

        #Store data
        self.ssp_lib['db_sspPrefit']        = stellarPrefit_db
        self.ssp_lib['dbDict_sspPrefit']    = stat_db_dict

        # Mean parameter values #TODO we need to add the redshift here
        self.sspPrefit_Coeffs               = stat_db_dict['ssp_coefficients']['trace'].mean(axis=0)
        self.ssp_lib['Av_sspPrefit']        = stat_db_dict['Av_star']['mean']
        self.ssp_lib['sigma_sspPrefit']     = stat_db_dict['sigma_star']['mean']

        return

class SpecSynthesizer(SpectrumFitter, BilinearInterpTheano):

    def __init__(self, atomic_ref=None, nebular_data = None, ssp_type=None, ssp_folder=None, ssp_conf_file=None,
                 temp_grid=None, den_grid=None, high_temp_ions=None, R_v=3.4, reddenig_curve='G03_average', lowlimit_sspContribution=0.0):

        #TODO in here we give the choice: Configuration from text file but it can also be imported from text file create a default conf to upload from there if missing

        # Import parent classes
        SpectrumFitter.__init__(self)
        BilinearInterpTheano.__init__(self)

        # Load nebular constants data
        if nebular_data is None:
            self.load_neb_constants(self.paths_dict['nebular_data_folder'])

        # Define atomic data:
        if atomic_ref == 'default':
            atomicData.setDataFile('s_iii_coll_HRS12.dat') #TODO actually this should be uplodaded from a file

        # Temperature and density grid declaration
        self.tem_grid_range = arange(temp_grid[0], temp_grid[1], temp_grid[2])
        self.den_grid_range = arange(den_grid[0], den_grid[1], den_grid[2])
        self.den_grid_range[0] = 1 #Change the minimum value for the grid to 1

        # Reddening parameters
        self.Rv_model = R_v
        self.reddedning_curve_model = reddenig_curve

        # Declare high ionization temperature ions
        self.high_temp_ions = high_temp_ions

        # Lower accepted limit for ssps
        self.lowlimit_sspContribution = lowlimit_sspContribution

        # Dictionary to store parameters
        self.conf_dic = {'Te_neb': 10000.0, 'z_neb': 0.0, 'cHbeta_neb': 0.1, 'He1_neb': 0.085, 'He2_neb': 0.0}

    def fit_observation(self, model_name, model_type, obs_data, ssp_data, iterations=15000, burn=7000, thin=1, output_folder = None,
                        fitting_components=None, input_lines=None, params_list = None, prefit_SSP = True, prefit_model = False,
                        wavelengh_limits = None, resample_inc = None, norm_interval=None):

        # Declare model type
        self.model_type = model_type

        # Folders to store inputs and outputs
        self.input_folder = output_folder + 'input_data/'
        self.output_folder = output_folder + 'output_data/'

        # Create them if not available
        make_folder(self.input_folder)
        make_folder(self.output_folder)

        # Prepare data for the MCMC fitting
        self.ready_simulation(output_folder, obs_data, ssp_data, fitting_components, input_lines=input_lines,
                              wavelengh_limits=wavelengh_limits, resample_inc=resample_inc, norm_interval=norm_interval)

        # Prefit stellar continua
        if 'stellar' in fitting_components:

            # Run prefit output #TODO these prefit parameters should go into the master conf
            self.prefit_db = self.input_folder + 'SSP_prefit_database'

            # Perform a new SPP synthesis otherwise use available data
            if prefit_SSP:

                # Ready data for a prefit on our stellar continua
                self.prepare_ssp_data(None, obj_wave=self.obj_data['wave_resam'], obj_flux=self.obj_data['flux_norm'],
                                      obj_flux_err=self.obj_data['sigma_continuum'], obj_mask=self.obj_data['int_mask'])

                # Select model
                self.select_inference_model('stelar_prefit')

                # Run stellar continua prefit and store the results
                self.run_pymc(self.prefit_db, iterations = 5000, variables_list = ['Av_star', 'sigma_star'], prefit = True)
                self.save_prefit_data(self.obj_data['wave_resam'])

                # Plot input simulation data
                self.plot_prefit_data()

            else:
                # Load prefit data
                self.load_prefit_data(self.obj_data['wave_resam'])

            # Ready stellar
            self.prepare_ssp_data(self.input_folder + 'prefit_populations.txt',
                                  obj_wave=self.obj_data['wave_resam'], obj_flux=self.obj_data['flux_norm'],
                                  obj_flux_err=self.obj_data['sigma_continuum'], obj_mask=self.obj_data['int_mask'])

        # Ready gas data
        self.prepare_gas_data(self.obj_data['recomb_fluxes'], self.obj_data['recomb_err'],
                            self.obj_data['metal_fluxes'], self.obj_data['metal_err'],
                            self.obj_data['TSIII'], self.obj_data['TSIII_error'],
                            self.obj_data['nSII'], self.obj_data['nSII_error'],
                            self.obj_data['cHbeta_prior'], self.obj_data['cHbeta_error_prior'])

        # Select model
        self.select_inference_model(model_type, params_list + self.abund_list)

        # Run model
        output_address = self.output_folder + model_name
        self.run_pymc(output_address, iterations=iterations, variables_list=params_list, prefit=prefit_model, model_type = model_type)

        # Load the results
        output_db, output_db_dict = self.load_pymc_database_manual(output_address, burning=burn, params_list=params_list)

        # Plot output data
        self.plot_ouput_data(output_address, output_db, output_db_dict, params_list)

        return

    def select_inference_model(self, model = 'complete', params_list = None):

        #Parameter list with the parameter space
        self.params_list = params_list

        # TODO Esta no es elegante, adjuntar a la configuracion
        # Check if we must analyse recombination and/or collisional excited lines
        possible_Recomb = ['He1_abund', 'He2_abund']
        recombination_check = set(self.params_list) & set(possible_Recomb)
        if recombination_check:
            self.fitting_components.append('recomb_lines')

        # Check if we must analyse recombination and/or collisional excited lines
        possible_colExcit = ['S2_abund','S3_abund','O2_abund','O3_abund', 'N2_abund', 'Ar3_abund', 'Ar4_abund']
        collisional_check = set(self.params_list) & set(possible_colExcit)

        if collisional_check:
            self.fitting_components.append('colExcited_lines')

        # Declare the simulation type
        if model == 'complete':
            self.inf_dict = self.complete_model()

        elif model == 'detect':
            self.inf_dict = self.detect_model()

        elif model == 'stelar_prefit':
            print 'Limites', self.zMin_SspLimit, self.zMax_SspLimit
            self.inf_dict = self.stellarContinua_model()

    def stellarContinua_model(self):

        #Stellar parameters
        Av_star     = pymc2.Uniform('Av_star', 0.0, 5.00)
        sigma_star  = pymc2.Uniform('sigma_star', 0.0, 5.00)
        #z_star      = pymc2.Uniform('z_star', self.z_min_ssp_limit, self.z_max_ssp_limit)

        @pymc2.deterministic #Shift, multiply and convolve by a factor given by the model parameters
        def ssp_coefficients(z_star=0.0, Av_star=Av_star, sigma_star=sigma_star, input_flux=self.input_continuum):

            ssp_grid_i = self.physical_SED_model(self.onBasesWave, self.input_wave, self.onBasesFluxNorm, Av_star, z_star, sigma_star, self.Rv_model)

            self.ssp_grid_i_masked = (self.int_mask * ssp_grid_i.T).T

            ssp_coeffs_norm = self.ssp_fitting(self.ssp_grid_i_masked, input_flux)

            return ssp_coeffs_norm

        @pymc2.deterministic # Theoretical normalized flux
        def stellar_continua_calculation(ssp_coeffs=ssp_coefficients):

            flux_sspFit_norm = np_sum(ssp_coeffs.T * self.ssp_grid_i_masked, axis=1)

            return flux_sspFit_norm

        @pymc2.stochastic(observed=True) # Likelihood
        def likelihood_ssp(value = self.input_continuum, stellarTheoFlux=stellar_continua_calculation, sigmaContinuum=self.input_continuum_er):
            chi_F = sum(square(stellarTheoFlux - value) / square(sigmaContinuum))
            return - chi_F / 2

        return locals()

    def complete_model(self):

        # TODO Priors data should go into configuration file

        # Gas parameters
        ne          = pymc2.TruncatedNormal('ne', self.obj_data['nSII'], self.obj_data['nSII_error'] ** -2, a=50.0, b=1000.0)
        cHbeta      = pymc2.TruncatedNormal('cHbeta', 0.15, 0.05 ** -2, a=0.0, b=3.0)
        T_low       = pymc2.TruncatedNormal('T_low', self.obj_data['TSIII'], self.obj_data['TSIII_error'] ** -2, a=7000.0, b=20000.0)

        # Metals abundances
        S2_abund    = pymc2.Uniform('S2_abund', 0.000001, 0.001)
        S3_abund    = pymc2.Uniform('S3_abund', 0.000001, 0.001)
        O2_abund    = pymc2.Uniform('O2_abund', 0.000001, 0.001)
        O3_abund    = pymc2.Uniform('O3_abund', 0.000001, 0.001)
        N2_abund    = pymc2.Uniform('N2_abund', 0.000001, 0.001)
        Ar3_abund   = pymc2.Uniform('Ar3_abund', 0.000001, 0.001)
        Ar4_abund   = pymc2.Uniform('Ar4_abund', 0.000001, 0.001)

        # Helium parameters
        He1_abund   = pymc2.Uniform('He1_abund', 0.050, 0.15)
        tau         = pymc2.TruncatedNormal('tau', 0.75, 0.5 ** -2, a=0.0, b=7.0)
        cHbeta      = pymc2.TruncatedNormal('cHbeta', 0.15, 0.05 ** -2, a=0.0, b=3.0)
        T_He        = pymc2.TruncatedNormal('T_He', self.obj_data['TSIII'], self.obj_data['TSIII_error'] ** -2, a=7000.0, b=20000.0, value=14500.0)

        #Stellar parameters
        Av_star     = pymc2.Uniform('Av_star', 0.0, 5.00)
        sigma_star  = pymc2.Uniform('sigma_star', 0.0, 5.00)
        # z_star    = pymc2.Uniform('z_star', self.z_min_ssp_limit, self.z_max_ssp_limit)
        ssp_coefs   = [pymc2.Uniform('ssp_coefs_%i' % i, self.sspPrefit_Limits[i][0], self.sspPrefit_Limits[i][1]) for i in self.range_bases]

        @pymc2.deterministic()
        def calc_Thigh(Te=T_low):
            return (1.0807 * Te / 10000.0 - 0.0846) * 10000.0

        @pymc2.deterministic()
        def calc_abund_dict(He1_abund=He1_abund, S2_abund=S2_abund, S3_abund=S3_abund, O2_abund=O2_abund, O3_abund=O3_abund, N2_abund=N2_abund, Ar3_abund=Ar3_abund, Ar4_abund=Ar4_abund):

            self.abund_iter_dict['H1']  = He1_abund
            self.abund_iter_dict['He1'] = He1_abund
            self.abund_iter_dict['S2']  = S2_abund
            self.abund_iter_dict['S3']  = S3_abund
            self.abund_iter_dict['O2']  = O2_abund
            self.abund_iter_dict['O3']  = O3_abund
            self.abund_iter_dict['N2']  = N2_abund
            self.abund_iter_dict['Ar3'] = Ar3_abund
            self.abund_iter_dict['Ar4'] = Ar4_abund

            return self.abund_iter_dict

        @pymc2.deterministic
        def calc_colExcit_fluxes(abund_dict=calc_abund_dict, T_low=T_low, T_High=calc_Thigh, ne=ne, cHbeta=cHbeta):

            colExcit_fluxes = self.calculate_colExcit_flux(T_low,
                                                           T_High,
                                                           ne,
                                                           cHbeta,
                                                           abund_dict,
                                                           self.obj_data['colLine_waves'],
                                                           self.obj_data['colLine_ions'],
                                                           self.obj_data['colLine_flambda'])

            return colExcit_fluxes

        @pymc2.deterministic
        def calc_nebular_cont(z_star=self.z_object, cHbeta=self.cHbeta, Te=self.TSIII, He1_abund=He1_abund, He2_abund=0.0, Halpha_Flux = self.f_HalphaNorm):

            neb_flux_norm = self.nebular_Cont(self.input_wave,
                                              z_star,
                                              cHbeta,
                                              Te,
                                              He1_abund,
                                              He2_abund,
                                              Halpha_Flux)

            return neb_flux_norm

        @pymc2.deterministic
        def calc_continuum(z_star=self.z_object, Av_star=Av_star, sigma_star=sigma_star, ssp_coefs=ssp_coefs, nebular_flux=calc_nebular_cont):

            ssp_grid_i = self.physical_SED_model(self.onBasesWave,
                                                 self.input_wave,
                                                 self.onBasesFluxNorm,
                                                 Av_star,
                                                 z_star,
                                                 sigma_star,
                                                 self.Rv_model)

            fit_continuum = ssp_grid_i.dot(ssp_coefs) + nebular_flux

            return fit_continuum

        @pymc2.deterministic
        def calc_recomb_fluxes(abund_dict=calc_abund_dict, T_He=T_He, ne=ne, cHbeta=cHbeta, tau=tau):

            recomb_fluxes = self.calculate_recomb_fluxes(T_He,
                                                         ne,
                                                         cHbeta,
                                                         tau,
                                                         abund_dict,
                                                         self.obj_data['recombLine_labes'],
                                                         self.obj_data['recombLine_ions'],
                                                         self.obj_data['recombLine_flambda'])

            return recomb_fluxes

        #QUESTION Issues with more than one likelihood
        @pymc2.stochastic(observed=True)  # Likelihood
        def likelihood_ssp(value=self.input_continuum, fit_continuum=calc_continuum, sigmaContinuum=self.input_continuum_er):
            calc_continuum_masked = fit_continuum * self.obj_data['int_mask']
            chi_F = sum(square(calc_continuum_masked - value) / square(sigmaContinuum))
            return - chi_F / 2

        @pymc2.stochastic(observed=True)  # Likelihood
        def likelihood_recomb(value=self.recomb_fluxes, H_He_TheoFlux=calc_recomb_fluxes, sigmaLines=self.recomb_err):
            chi_F = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2

        @pymc2.stochastic(observed=True)  # Likelihood
        def likelihood_colExcited(value=self.colExc_fluxes, theo_metal_fluzes=calc_colExcit_fluxes, sigmaLines=self.colExc_fluxes):
            chi_F = sum(square(theo_metal_fluzes - value) / square(sigmaLines))
            return - chi_F / 2

        return locals()

    def detect_model(self):

        # TODO Priors data should go into configuration file

        # Gas parameters
        if ('emission' in self.fitting_components) or ('nebular' in self.fitting_components):
            n_e     = pymc2.TruncatedNormal('n_e', self.obj_data['nSII'], self.obj_data['nSII_error'] ** -2, a=1.0, b=1000.0)
            cHbeta  = pymc2.TruncatedNormal('cHbeta', 0.15, 0.05 ** -2, a=0.0, b=3.0)
            T_low   = pymc2.TruncatedNormal('T_low', self.obj_data['TSIII'], self.obj_data['TSIII_error'] ** -2, a=7000.0, b=20000.0)

        # Metals abundances # QUESTION prior from dictionary
        S2_abund = pymc2.Uniform('S2_abund', 0.000001, 0.001) if 'S2_abund' in self.params_list else 0.0
        S3_abund = pymc2.Uniform('S3_abund', 0.000001, 0.001) if 'S3_abund' in self.params_list else 0.0
        O2_abund = pymc2.Uniform('O2_abund', 0.000001, 0.001) if 'O2_abund' in self.params_list else 0.0
        O3_abund = pymc2.Uniform('O3_abund', 0.000001, 0.001) if 'O3_abund' in self.params_list else 0.0
        N2_abund = pymc2.Uniform('N2_abund', 0.000001, 0.001) if 'N2_abund' in self.params_list else 0.0
        Ar3_abund = pymc2.Uniform('Ar3_abund', 0.000001, 0.001) if 'Ar3_abund' in self.params_list else 0.0
        Ar4_abund = pymc2.Uniform('Ar4_abund', 0.000001, 0.001) if 'Ar4_abund' in self.params_list else 0.0

        # Helium parameters
        He1_abund   = pymc2.Uniform('He1_abund', 0.050, 0.15) if 'He1_abund' in self.params_list else 0.0
        tau         = pymc2.TruncatedNormal('tau', 0.75, 0.5 ** -2, a=0.0, b=7.0)  if 'tau' in self.params_list else 0.0
        #T_He       = pymc2.TruncatedNormal('T_He', self.obj_data['TSIII'], self.obj_data['TSIII_error'] ** -2, a=7000.0, b=20000.0, value=self.obj_data['TSIII'])

        #Stellar parameters
        if 'stellar' in self.fitting_components:
            Av_star     = pymc2.Uniform('Av_star', 0.0, 5.00)
            sigma_star  = pymc2.Uniform('sigma_star', 0.0, 5.00)
            # z_star    = pymc2.Uniform('z_star', self.z_min_ssp_limit, self.z_max_ssp_limit)
            ssp_coefs   = [pymc2.Uniform('ssp_coefs_%i' % i, self.sspPrefit_Limits[i][0], self.sspPrefit_Limits[i][1]) for i in self.range_bases]

        @pymc2.deterministic()
        def calc_Thigh(Te=T_low):
            return TOIII_TSIII_relation(Te)

        @pymc2.deterministic()
        def calc_abund_dict(He1_abund=He1_abund, S2_abund=S2_abund, S3_abund=S3_abund, O2_abund=O2_abund, O3_abund=O3_abund, N2_abund=N2_abund, Ar3_abund=Ar3_abund, Ar4_abund=Ar4_abund):

            self.abund_iter_dict['H1']  = 1.0
            self.abund_iter_dict['He1'] = He1_abund
            self.abund_iter_dict['S2']  = S2_abund
            self.abund_iter_dict['S3']  = S3_abund
            self.abund_iter_dict['O2']  = O2_abund
            self.abund_iter_dict['O3']  = O3_abund
            self.abund_iter_dict['N2']  = N2_abund
            self.abund_iter_dict['Ar3'] = Ar3_abund
            self.abund_iter_dict['Ar4'] = Ar4_abund

            return self.abund_iter_dict

        if 'nebular' in self.fitting_components:
            @pymc2.deterministic
            def calc_nebular_cont(z_star=self.z_object, cHbeta=self.cHbeta, Te=self.TSIII, He1_abund=He1_abund, He2_abund=0.0, Halpha_Flux = self.f_HalphaNorm):

                neb_flux_norm = self.nebular_Cont(self.input_wave,
                                                  z_star,
                                                  cHbeta,
                                                  Te,
                                                  He1_abund,
                                                  He2_abund,
                                                  Halpha_Flux)

                return neb_flux_norm

        if 'stellar' in self.fitting_components:
            @pymc2.deterministic
            def calc_continuum(z_star=self.z_object, Av_star=Av_star, sigma_star=sigma_star, ssp_coefs=ssp_coefs, nebular_flux=calc_nebular_cont):

                ssp_grid_i = self.physical_SED_model(self.onBasesWave,
                                                     self.input_wave,
                                                     self.onBasesFluxNorm,
                                                     Av_star,
                                                     z_star,
                                                     sigma_star,
                                                     self.Rv_model)

                fit_continuum = ssp_grid_i.dot(ssp_coefs) + nebular_flux

                return fit_continuum


        if 'recomb_lines' in self.fitting_components:

            @pymc2.deterministic
            def calc_recomb_fluxes(abund_dict=calc_abund_dict, T_He=calc_Thigh, ne=n_e, cHbeta=cHbeta, tau=tau):

                recomb_fluxes = self.calculate_recomb_fluxes(T_He,
                                                             ne,
                                                             cHbeta,
                                                             tau,
                                                             abund_dict,
                                                             self.obj_data['recombLine_labes'],
                                                             self.obj_data['recombLine_ions'],
                                                             self.obj_data['recombLine_flambda'])

                return recomb_fluxes


        if 'colExcited_lines' in self.fitting_components:

            @pymc2.deterministic
            def calc_colExcit_fluxes(abund_dict=calc_abund_dict, T_low=T_low, T_High=calc_Thigh, ne=n_e, cHbeta=cHbeta):

                colExcit_fluxes = self.calculate_colExcit_flux(T_low,
                                                               T_High,
                                                               ne,
                                                               cHbeta,
                                                               abund_dict,
                                                               self.obj_data['colLine_waves'],
                                                               self.obj_data['colLine_ions'],
                                                               self.obj_data['colLine_flambda'])

                return colExcit_fluxes


        # QUESTION Issues with more than one likelihood
        if 'stellar' in self.fitting_components:
            @pymc2.stochastic(observed=True)  # Likelihood
            def likelihood_ssp(value=self.input_continuum, fit_continuum=calc_continuum, sigmaContinuum=self.input_continuum_er):
                calc_continuum_masked = fit_continuum * self.obj_data['int_mask']
                chi_F = sum(square(calc_continuum_masked - value) / square(sigmaContinuum))
                return - chi_F / 2

        if 'emission' in self.fitting_components:

            if 'recomb_lines' in self.fitting_components:

                @pymc2.stochastic(observed=True)  # Likelihood
                def likelihood_recomb(value=self.recomb_fluxes, H_He_TheoFlux=calc_recomb_fluxes, sigmaLines=self.recomb_err):
                    chi_F = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
                    return - chi_F / 2


            if 'colExcited_lines' in self.fitting_components:
                @pymc2.stochastic(observed=True)  # Likelihood
                def likelihood_colExcited(value=self.colExc_fluxes, theo_metal_fluzes=calc_colExcit_fluxes, sigmaLines=self.colExc_fluxes):
                    chi_F = sum(square(theo_metal_fluzes - value) / square(sigmaLines))
                    return - chi_F / 2

        return locals()

    def pymc3_model(self):

        with pymc3.Model() as model:

            # Physical conditions priors
            T_low   = pymc3.Normal('T_low', mu=self.obj_data['TSIII'], tau=self.obj_data['TSIII_error']**-2, shape=self.n_colExcLines)
            T_high  = (1.0807 * T_low / 10000.0 - 0.0846) * 10000.0

            n_e     = pymc3.Normal('n_e', mu=self.obj_data['nSII'], tau=self.obj_data['nSII_error']**-2, shape=self.n_colExcLines)
            cHbeta  = pymc3.Uniform('cHbeta', lower=0, upper=2)

            # Abundances priors
            abund_dict = {'S2': pymc3.Uniform('S2_abund', lower=-5, upper=-2),
                          'S3': pymc3.Uniform('S3_abund', lower=-5, upper=-2),
                          'O2': pymc3.Uniform('O2_abund', lower=-5, upper=-2),
                          'O3': pymc3.Uniform('O3_abund', lower=-5, upper=-2)}

            # Compute higher region temperature
            T_high = (1.0807 * T_low / 10000.0 - 0.0846) * 10000.0

            # Temperature dict
            Temp_dict = {'T_low ':T_low, 'T_high':T_high}

            # Loop through the lines to assign the ionic abundances
            log_abund_vector = empty(self.n_colExcLines)
            for i in self.range_colExcLines:
                line_label = self.obj_data['colLine_labes'][i]
                ion = self.obj_data['colLine_ions'][i]
                T_ion = Temp_dict[self.tempRegionCol[i]]
                abund = abund_dict[ion]

                #Interpolate the emissivity
                emis_ion = bilinear_interpolator_axis(n_e[i].tag.test_value, T_ion[i].tag.test_value,
                                                                            self.den_grid_range,
                                                                            self.tem_grid_range,
                                                                            self.log_metals_emis_grid[:, :, i])
                # Compute synthetic fluxes
                colExcit_fluxes = abund + emis_ion - cHbeta * self.obj_data['colLine_flambda'][i]

                # Likelihood
                chiSq = pymc3.Normal(line_label + '_chiSq', mu=colExcit_fluxes, tau=self.log_colExc_Err[i], observed = self.log_colExc_fluxes[i])

            for RV in model.basic_RVs:
                print(RV.name, RV.logp(model.test_point))

            # Launch model
            trace = pymc3.sample(10000, tune=2000)

        return trace, model

    def pymc3_model2(self):

        with pymc3.Model() as model:

            # Physical conditions priors
            T_low   = pymc3.Normal('T_e', mu=self.obj_data['TSIII'], tau=self.obj_data['TSIII_error']**-2)
            T_high  = (1.0807 * T_low / 10000.0 - 0.0846) * 10000.0

            n_e     = pymc3.Uniform('d_e', lower=1, upper=500)#pymc3.Normal('n_e', mu=self.obj_data['nSII'], tau=self.obj_data['nSII_error']**-2)
            cHbeta  = pymc3.Uniform('cHbeta', lower=0, upper=2)

            # Abundances priors
            abund_dict = {'S2': pymc3.Uniform('S2_abund', lower=0, upper=15),
                          'S3': pymc3.Uniform('S3_abund', lower=0, upper=15),
                          'O2': pymc3.Uniform('O2_abund', lower=0, upper=15),
                          'O3': pymc3.Uniform('O3_abund', lower=0, upper=15)}

            # Compute higher region temperature
            T_high = (1.0807 * T_low / 10000.0 - 0.0846) * 10000.0

            # Temperature dict
            Temp_dict = {'T_low ':T_low, 'T_high':T_high}

            ne_abs_dif  = pm3_math.abs_(self.den_grid_range - n_e)
            ne_idx0_tt  = tt.argmin(ne_abs_dif)

            # Get indeces from den and temp value for the emissivity grid
            ne_idx0     = self.indexing_bilinear_pymc3(self.den_grid_range, n_e)
            T_low_idx0  = self.indexing_bilinear_pymc3(self.tem_grid_range, T_low)
            T_high_idx0 = self.indexing_bilinear_pymc3(self.tem_grid_range, T_high)

            # Temperature index dict
            temp_dict       = {'T_low ':T_low, 'T_high':T_high}
            temp_idx_dict   = {'T_low ':T_low_idx0, 'T_high':T_high_idx0}

            # Loop through the lines to assign the ionic abundances
            log_abund_vector = empty(self.n_colExcLines)
            for i in self.range_colExcLines:
                line_label  = self.obj_data['colLine_labes'][i]
                ion         = self.obj_data['colLine_ions'][i]
                T_ion       = Temp_dict[self.tempRegionCol[i]]
                T_ion_idx0  = Temp_dict[self.tempRegionCol[i]]
                abund       = abund_dict[ion]

                #Interpolate the emissivity
                emis_ion = bilinear_interpolator_axis(n_e, T_ion, ne_idx0, T_ion_idx0,
                                                        self.den_grid_range,self.tem_grid_range, self.log_metals_emis_grid[:, :, i])

                # Compute synthetic fluxes
                colExcit_fluxes = abund + emis_ion - cHbeta * self.obj_data['colLine_flambda'][i] - 12

                # Likelihood
                chiSq = pymc3.Normal(line_label + '_chiSq', mu=colExcit_fluxes, tau=self.log_colExc_Err[i], observed = self.log_colExc_fluxes[i])

            for RV in model.basic_RVs:
                print(RV.name, RV.logp(model.test_point))

            # Launch model
            trace = pymc3.sample(10000, tune=2000)

        return trace, model

    def chemical_model(self):

        with pymc3.Model() as model:

            a_cord = pymc3.Uniform('d_e', lower=1, upper=500)
            a_idx0 = pymc3.math.argmin(pymc3.math.abs_(a_array-a_cord))

            # Physical conditions priors
            T_e   = pymc3.Normal('T_e', mu=self.obj_data['TSIII'], tau=self.obj_data['TSIII_error']**-2)
            d_e     = pymc3.Uniform('d_e', lower=1, upper=500)
            cHbeta  = pymc3.Uniform('cHbeta', lower=0, upper=2)

            # Abundances priors
            abund_dict = {'S2': pymc3.Uniform('S2_abund', lower=-5, upper=-2), 'S3': pymc3.Uniform('S3_abund', lower=-5, upper=-2),
                          'O2': pymc3.Uniform('O2_abund', lower=-5, upper=-2), 'O3': pymc3.Uniform('O3_abund', lower=-5, upper=-2)}

            # Temperature dict
            Temp_dict = {'T_low ':T_e.tag.test_value}

            # Loop through the lines to assign the ionic abundances
            for i in self.range_line_peaks:
                line_label = self.obj_data['line_labes'][i]
                ion = self.obj_data['line_ion'][i]
                abund = abund_dict[ion]
                T_ion = T_e.tag.test_value
                de_ion = d_e.tag.test_value

                # Interpolate the emissivity
                emis_Ion = bilinear_interpolator_axis(de_ion, T_ion, self.den_grid_range, self.tem_grid_range, self.emis_grid[:, :, i])

                # Compute synthetic fluxes
                synthFluxes = abund + emis_Ion - cHbeta * self.obj_data['colLine_flambda'][i]

                # Likelihood
                chiSq = pymc3.Normal(line_label + '_like', mu=synthFluxes, tau=self.logErrFlux[i], observed=self.logFluxes[i])

            # Launch model
            trace = pymc3.sample(10000, tune=2000)

        return trace

    def detect_model_porsi(self):

        # TODO Priors data should go into configuration file

        # Gas parameters
        if ('emission' in self.fitting_components) or ('nebular' in self.fitting_components):
            n_e          = pymc2.TruncatedNormal('n_e', self.obj_data['nSII'], self.obj_data['nSII_error'] ** -2, a=1.0, b=1000.0)
            cHbeta      = pymc2.TruncatedNormal('cHbeta', 0.15, 0.05 ** -2, a=0.0, b=3.0)
            T_low       = pymc2.TruncatedNormal('T_low', self.obj_data['TSIII'], self.obj_data['TSIII_error'] ** -2, a=7000.0, b=20000.0)

        # Metals abundances
        # S2_abund    = pymc2.Uniform('S2_abund', 0.000001, 0.001)
        # S3_abund    = pymc2.Uniform('S3_abund', 0.000001, 0.001)
        # O2_abund    = pymc2.Uniform('O2_abund', 0.000001, 0.001)
        # O3_abund    = pymc2.Uniform('O3_abund', 0.000001, 0.001)
        # N2_abund    = pymc2.Uniform('N2_abund', 0.000001, 0.001)
        # Ar3_abund   = pymc2.Uniform('Ar3_abund', 0.000001, 0.001)
        # Ar4_abund   = pymc2.Uniform('Ar4_abund', 0.000001, 0.001)

        # QUESTION prior from dictionary
        S2_abund    = 0
        S3_abund    = 0
        O2_abund    = 0
        O3_abund    = 0
        N2_abund    = 0
        Ar3_abund   = 0
        Ar4_abund   = 0

        # Helium parameters
        He1_abund   = pymc2.Uniform('He1_abund', 0.050, 0.15)
        tau         = pymc2.TruncatedNormal('tau', 0.75, 0.5 ** -2, a=0.0, b=7.0)
        cHbeta      = pymc2.TruncatedNormal('cHbeta', 0.15, 0.05 ** -2, a=0.0, b=3.0)
        #T_He        = pymc2.TruncatedNormal('T_He', self.obj_data['TSIII'], self.obj_data['TSIII_error'] ** -2, a=7000.0, b=20000.0, value=self.obj_data['TSIII'])

        #Stellar parameters
        if 'stellar' in self.fitting_components:
            Av_star     = pymc2.Uniform('Av_star', 0.0, 5.00)
            sigma_star  = pymc2.Uniform('sigma_star', 0.0, 5.00)
            # z_star    = pymc2.Uniform('z_star', self.z_min_ssp_limit, self.z_max_ssp_limit)
            ssp_coefs   = [pymc2.Uniform('ssp_coefs_%i' % i, self.sspPrefit_Limits[i][0], self.sspPrefit_Limits[i][1]) for i in self.range_bases]

        @pymc2.deterministic()
        def calc_Thigh(Te=T_low):
            return TOIII_TSIII_relation(Te)

        @pymc2.deterministic()
        def calc_abund_dict(He1_abund=He1_abund, S2_abund=S2_abund, S3_abund=S3_abund, O2_abund=O2_abund, O3_abund=O3_abund, N2_abund=N2_abund, Ar3_abund=Ar3_abund, Ar4_abund=Ar4_abund):

            self.abund_iter_dict['H1']  = 1.0
            self.abund_iter_dict['He1'] = He1_abund
            self.abund_iter_dict['S2']  = S2_abund
            self.abund_iter_dict['S3']  = S3_abund
            self.abund_iter_dict['O2']  = O2_abund
            self.abund_iter_dict['O3']  = O3_abund
            self.abund_iter_dict['N2']  = N2_abund
            self.abund_iter_dict['Ar3'] = Ar3_abund
            self.abund_iter_dict['Ar4'] = Ar4_abund

            return self.abund_iter_dict

        # @pymc2.deterministic
        # def calc_colExcit_fluxes(abund_dict=calc_abund_dict, T_low=T_low, T_High=calc_Thigh, ne=n_e, cHbeta=cHbeta):
        #
        #     colExcit_fluxes = self.calculate_colExcit_flux(T_low,
        #                                                    T_High,
        #                                                    ne,
        #                                                    cHbeta,
        #                                                    abund_dict,
        #                                                    self.obj_data['colLine_waves'],
        #                                                    self.obj_data['colLine_ions'],
        #                                                    self.obj_data['colLine_flambda'])
        #
        #     return colExcit_fluxes

        if 'nebular' in self.fitting_components:

            @pymc2.deterministic
            def calc_nebular_cont(z_star=self.z_object, cHbeta=self.cHbeta, Te=self.TSIII, He1_abund=He1_abund, He2_abund=0.0, Halpha_Flux = self.f_HalphaNorm):

                neb_flux_norm = self.nebular_Cont(self.input_wave,
                                                  z_star,
                                                  cHbeta,
                                                  Te,
                                                  He1_abund,
                                                  He2_abund,
                                                  Halpha_Flux)

                return neb_flux_norm

        if 'stellar' in self.fitting_components:

            @pymc2.deterministic
            def calc_continuum(z_star=self.z_object, Av_star=Av_star, sigma_star=sigma_star, ssp_coefs=ssp_coefs, nebular_flux=calc_nebular_cont):

                ssp_grid_i = self.physical_SED_model(self.onBasesWave,
                                                     self.input_wave,
                                                     self.onBasesFluxNorm,
                                                     Av_star,
                                                     z_star,
                                                     sigma_star,
                                                     self.Rv_model)

                fit_continuum = ssp_grid_i.dot(ssp_coefs) + nebular_flux

                return fit_continuum

        @pymc2.deterministic
        def calc_recomb_fluxes(abund_dict=calc_abund_dict, T_He=calc_Thigh, ne=n_e, cHbeta=cHbeta, tau=tau):

            recomb_fluxes = self.calculate_recomb_fluxes(T_He,
                                                         ne,
                                                         cHbeta,
                                                         tau,
                                                         abund_dict,
                                                         self.obj_data['recombLine_labes'],
                                                         self.obj_data['recombLine_ions'],
                                                         self.obj_data['recombLine_flambda'])

            return recomb_fluxes

        # QUESTION Issues with more than one likelihood

        if 'stellar' in self.fitting_components:
            @pymc2.stochastic(observed=True)  # Likelihood
            def likelihood_ssp(value=self.input_continuum, fit_continuum=calc_continuum, sigmaContinuum=self.input_continuum_er):
                calc_continuum_masked = fit_continuum * self.obj_data['int_mask']
                chi_F = sum(square(calc_continuum_masked - value) / square(sigmaContinuum))
                return - chi_F / 2

        if 'emission' in self.fitting_components:
            @pymc2.stochastic(observed=True)  # Likelihood
            def likelihood_recomb(value=self.recomb_fluxes, H_He_TheoFlux=calc_recomb_fluxes, sigmaLines=self.recomb_err):
                chi_F = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
                return - chi_F / 2

        # @pymc2.stochastic(observed=True)  # Likelihood
        # def likelihood_colExcited(value=self.colExc_fluxes, theo_metal_fluzes=calc_colExcit_fluxes, sigmaLines=self.colExc_fluxes):
        #     chi_F = sum(square(theo_metal_fluzes - value) / square(sigmaLines))
        #     return - chi_F / 2

        return locals()

    def run_pymc(self, db_address, iterations = 10000, variables_list = None, prefit = True, model_type = 'pymc2'):

        if 'pymc3' not in model_type:

            # Define MCMC model
            self.MAP_Model = pymc2.MAP(self.inf_dict)

            # Prefit:
            if prefit is not False:
                fit_method = prefit if prefit is str else 'fmin_powell'
                print '\n--Starting {} prefit'.format(fit_method)
                start = timer()
                self.MAP_Model.fit(method = fit_method)
                end = timer()
                print 'prefit interval', (end - start)

            # Print prefit data
            print 'Initial conditions:'
            self.display_run_data(self.MAP_Model, variables_list)

            # Launch sample
            print '\nInitiating fit:'
            self.pymc2_M = pymc2.MCMC(self.MAP_Model.variables, db = 'pickle', dbname =  db_address)
            self.pymc2_M.sample(iter=iterations)

            # Save the output csv mean data
            if variables_list != None:
                print '--Saving results in csv'
                self.csv_address = db_address + '_Parameters'
                self.pymc2_M.write_csv(self.csv_address, variables=variables_list)

            #Print again the output prediction for the entire trace
            self.display_run_data(self.MAP_Model, variables_list)

            #Close the database
            self.pymc2_M.db.close()

        else:
            # Address to store the pickle
            model_address = db_address + '.pkl'

            # Launch sample
            trace, pymc3_model = self.pymc3_model2()

            # # Load he data
            # with open(model_address, 'rb') as buff:
            #     pickleData = pickle.load(buff)
            #
            # pymc3_model, trace = pickleData['model'], pickleData['trace']

            #Save the data
            print pymc3.summary(trace)
            print model_address
            pymc3.traceplot(trace)

            plt.show()

            # Save the data
            with open(model_address, 'wb') as trace_pickle:
                pickle.dump({'model': pymc3_model, 'trace': trace}, trace_pickle)

    def load_pymc_database_manual(self, db_address, burning = 0, params_list = None):

        #Load the pymc output textfile database
        pymc_database = pymc2.database.pickle.load(db_address)

        #Create a dictionaries with the traces and statistics
        traces_dic = {}
        stats_dic = OrderedDict()
        stats_dic['true_values'] = empty(len(params_list))

        #This variable contains all the traces from the MCMC (stochastic and deterministic)
        traces_list = pymc_database.trace_names[0]

        #Get statistics from the run
        for i in range(len(traces_list)):

            trace               = traces_list[i]
            stats_dic[trace]    = OrderedDict()
            trace_array         = pymc_database.trace(trace)[burning:]
            traces_dic[trace]   = trace_array

            if 'dict' not in trace:
                stats_dic[trace]['mean']                    = mean(trace_array)
                stats_dic[trace]['median']                  = median(trace_array)
                stats_dic[trace]['standard deviation']      = std(trace_array)
                stats_dic[trace]['n']                       = trace_array.shape[0]
                stats_dic[trace]['16th_p']                  = percentile(trace_array, 16)
                stats_dic[trace]['84th_p']                  = percentile(trace_array, 84)
                stats_dic[trace]['95% HPD interval']        = (stats_dic[trace]['16th_p'], stats_dic[trace]['84th_p'])
                stats_dic[trace]['trace']                   = trace_array

                if trace in params_list:

                    if trace in self.obj_data: #TODO we need a better structure fo this
                        stats_dic[trace]['true_value'] = self.obj_data[trace]

        #Generate a MCMC object to recover all the data from the run
        dbMCMC = pymc2.MCMC(traces_dic, pymc_database)

        return dbMCMC, stats_dic

    def display_run_data(self, database, variables_list):

        for param in variables_list:
            param_entry = getattr(database, param, None)
            if param_entry is not None:
                try:
                    print '-{} {}'.format(param, param_entry.value)
                except:
                    print 'I could not do it ', param