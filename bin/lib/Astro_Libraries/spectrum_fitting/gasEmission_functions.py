import numpy as np
import pyneb as pn
from os import path
from pandas import read_excel
from bisect import bisect_left
from scipy.optimize import curve_fit
from DZ_LineMesurer import LineMesurer_v2
from numpy import ones, exp, zeros, log10
from tensor_tools import EmissivitySurfaceFitter_tensorOps, EmissionEquations_tensorOps, EmissionTensors, storeValueInTensor
from inspect import getargspec
from collections import OrderedDict
from bin.lib.Astro_Libraries.spectrum_fitting.import_functions import parseDataFile
from bin.lib.Astro_Libraries.spectrum_fitting.abundance_tools import ChemicalModel

def TOIII_TSIII_relation(TSIII):
    # TODO we should make a class with these physical relations
    return (1.0807 * TSIII / 10000.0 - 0.0846) * 10000.0


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


class EmissivitySurfaceFitter(EmissivitySurfaceFitter_tensorOps):

    def __init__(self):

        # Class with the tensor operations of this class
        EmissivitySurfaceFitter_tensorOps.__init__(self)

        # TODO we should read this from the xlsx file
        self.ionEmisEq_fit = {'S2_6716A'  : self.emisEquation_TeDe,
                            'S2_6731A'  : self.emisEquation_TeDe,
                            'S3_6312A'  : self.emisEquation_Te,
                            'S3_9069A'  : self.emisEquation_Te,
                            'S3_9531A'  : self.emisEquation_Te,
                            'Ar4_4740A' : self.emisEquation_Te,
                            'Ar3_7136A' : self.emisEquation_Te,
                            'Ar3_7751A' : self.emisEquation_Te,
                            'O3_4363A'  : self.emisEquation_Te,
                            'O3_4959A'  : self.emisEquation_Te,
                            'O3_5007A'  : self.emisEquation_Te,
                            'O2_7319A'  : self.emisEquation_TeDe,
                            'O2_7330A'  : self.emisEquation_TeDe,
                            'O2_7319A_b': self.emisEquation_TeDe,
                            'N2_6548A'  : self.emisEquation_Te,
                            'N2_6584A'  : self.emisEquation_Te,
                            'H1_4102A'  : self.emisEquation_HI,
                            'H1_4341A'  : self.emisEquation_HI,
                            'H1_6563A'  : self.emisEquation_HI,
                            'He1_3889A' : self.emisEquation_HeI_fit,
                            'He1_4026A' : self.emisEquation_HeI_fit,
                            'He1_4471A' : self.emisEquation_HeI_fit,
                            'He1_5876A' : self.emisEquation_HeI_fit,
                            'He1_6678A' : self.emisEquation_HeI_fit,
                            'He1_7065A' : self.emisEquation_HeI_fit,
                            'He2_4686A' : self.emisEquation_HeII_fit}

        self.EmisRatioEq_fit = {'RSII'    : self.emisEquation_TeDe,
                            'RSIII'     : self.emisEquation_Te,
                            'ROIII'     : self.emisEquation_Te}

        # Second dictionary for the line fluxes
        self.ionEmisEq = dict(self.ionEmisEq_fit)

        # for line in ['He1_3889A', 'He1_4026A', 'He1_4471A', 'He1_5876A', 'He1_5876A', 'He1_6678A', 'He1_7065A']:
        #     self.ionEmisEq[line] = self.emisEquation_HeI
        # self.ionEmisEq['He2_4686A'] = self.emisEquation_HeII

        # Initial coeffient values to help with the fitting
        self.epm2017_emisCoeffs = {'He1_3889A': np.array([0.173, 0.00054, 0.904, 1e-5]),
                            'He1_4026A': np.array([-0.09, 0.0000063, 4.297, 1e-5]),
                            'He1_4471A': np.array([-0.1463, 0.0005, 2.0301, 1.5e-5]),
                            'He1_5876A': np.array([-0.226, 0.0011, 0.745, -5.1e-5]),
                            'He1_6678A': np.array([-0.2355, 0.0016, 2.612, 0.000146]),
                            'He1_7065A': np.array([0.368, 0.0017, 4.329, 0.0024])}

        return

    def emisEquation_Te(self, xy_space, a, b, c):
        temp_range, den_range = xy_space
        return a + b / (temp_range/10000.0) + c * np.log10(temp_range/10000)

    def emisEquation_TeDe(self, xy_space, a, b, c, d, e):
        temp_range, den_range = xy_space
        return a + b / (temp_range/10000.0) + c * np.log10(temp_range/10000) + np.log10(1 + e * den_range)

    def emisEquation_HI(self, xy_space, a, b, c):
        temp_range, den_range = xy_space
        return a + b * np.log10(temp_range) + c * np.log10(temp_range) * np.log10(temp_range)

    def emisEquation_HeI_fit(self, xy_space, a, b, c, d):
        temp_range, den_range = xy_space
        return (a + b * den_range) * np.log10(temp_range / 10000.0) - np.log10(c + d * den_range)

    def emisEquation_HeII_fit(self, xy_space, a, b):
        temp_range, den_range = xy_space
        return a + b * np.log(temp_range/10000)

    def fitEmis(self, func_emis, xy_space, line_emis, p0 = None):
        p1, p1_cov = curve_fit(func_emis, xy_space, line_emis, p0)
        return p1, p1_cov

class EmissionEquations(EmissionEquations_tensorOps, EmissionTensors):

    def __init__(self):

        # Class with the tensor operations of this class
        EmissionEquations_tensorOps.__init__(self)
        EmissionTensors.__init__(self)

        self.fluxEq = {}

    def assignFluxEq2Label(self, labelsList, ionsList, recombLabels=['H1r', 'He1r', 'He2r']):

        eqLabelArray = np.copy(ionsList)

        for i in range(eqLabelArray.size):
            if eqLabelArray[i] not in recombLabels:
                # TODO integrate a dictionary to add corrections
                if labelsList[i] != 'O2_7319A_b':
                    eqLabelArray[i] = 'metals'
                else:
                    eqLabelArray[i] = 'O2_7319A_b'

        return eqLabelArray

    def H1_lines(self, emis_ratio, cHbeta, flambda, abund=None, ftau=None, continuum=None):
        return emis_ratio - flambda * cHbeta

    def He1_lines(self, emis_ratio, cHbeta, flambda, abund, ftau=None, continuum=None):
        return log10(abund) + emis_ratio + log10(ftau) - flambda * cHbeta

    def He2_lines(self, emis_ratio, cHbeta, flambda, abund, ftau=None, continuum=None):
        return log10(abund) + emis_ratio - flambda * cHbeta

    def He1_linesFlux(self, emis_ratio, cHbeta, flambda, abund, ftau=None, continuum=None):
        return abund * emis_ratio * ftau * np.power(10, -1 * flambda * cHbeta)

    def He2_linesFlux(self, emis_ratio, cHbeta, flambda, abund, ftau=None, continuum=None):
        return abund * emis_ratio * np.power(10, -1 * flambda * cHbeta)

    def metal_lines(self, emis_ratio, cHbeta, flambda, abund, ftau=None, continuum=None):
        return abund + emis_ratio - flambda * cHbeta - 12

class EmissionComponents(EmissivitySurfaceFitter, EmissionEquations, ChemicalModel, LineMesurer_v2):

    def __init__(self, tempGrid = None, denGrid = None):

        # TODO this should go to the text file
        self.blendedLinesDict = dict(O2_7319A_b = ['O2_7319A', 'O2_7330A'], O2_3721A_b = ['O2_3726A', 'O2_3729A'])

        # Tools to fit a the emissivity to a plane
        EmissivitySurfaceFitter.__init__(self)

        # Functions to compute fluxes
        EmissionEquations.__init__(self)

        # Functions to compute elemental abundances
        ChemicalModel.__init__(self)

        # Load tools to measure lines # TODO there should be a better way to declare these folders
        LineMesurer_v2.__init__(self, '/home/vital/PycharmProjects/dazer/format/', 'DZT_LineLog_Headers.dz')
        # folderHeaders = path.join( 'C:\\Users\\Vital\\PycharmProjects\\dazer\\format\\', '')
        # LineMesurer_v2.__init__(self, folderHeaders, 'DZT_LineLog_Headers.dz')

        # Load default lines and their properties
        linesDataFile = path.join(self.externalDataFolder, self.config['linesData_file'])
        self.linesDb = read_excel(linesDataFile, sheet_name=0, header=0, index_col=0)

        # Defining temperature and density grids
        if (tempGrid is not None) and (denGrid is not None):
            self.tempRange = np.linspace(tempGrid[0], tempGrid[1], tempGrid[2])
            self.denRange = np.linspace(denGrid[0], denGrid[1], denGrid[2])
            X, Y = np.meshgrid(self.tempRange, self.denRange)
            self.tempGridFlatten, self.denGridFlatten = X.flatten(), Y.flatten()

    def calc_emis_spectrum(self, wavelength_range, lines_waves, lines_fluxes, lines_mu, lines_sigma, redshift):

        #TODO This should go into the cinematics class

        lines_mu = lines_mu * (1 + redshift)

        A_lines = lines_fluxes / (lines_sigma * self.sqrt2pi)

        wave_matrix = wavelength_range * ones((lines_waves.shape[0], wavelength_range.shape[0]))

        emis_spec_matrix = A_lines[:, None] * exp(-(wave_matrix - lines_mu[:, None]) * (wave_matrix - lines_mu[:, None]) / (2 * lines_sigma * lines_sigma))

        combined_spectrum = emis_spec_matrix.sum(axis=0)

        return combined_spectrum

    def ftau_func(self, tau, temp, den, a, b, c, d):
        return 1 + tau/2.0 * (a + (b + c * den + d * den * den) * temp/10000.0)

    def import_atomic_data(self, atomic_data = None, ftau_coefficients = None, tempGrid = None, denGrid = None):

        # Load nebular continuum constants
        self.loadNebCons(self.externalDataFolder)

        # # Import default database with line labels if not provided # TODO not the right place for loading this file
        # if atomic_data is not None:
        #     linesDataFile = atomic_data if atomic_data is not None else path.join(self.externalDataFolder, self.config['linesData_file'])
        #     self.linesDb = read_excel(linesDataFile, sheet_name=0, header=0, index_col=0)

        # Additional code to adapt to old lines log labelling
        self.linesDb['pyneb_code'] = self.linesDb['pyneb_code'].astype(str)
        idx_numeric_pynebCode = ~self.linesDb['pyneb_code'].str.contains('A+')
        self.linesDb.loc[idx_numeric_pynebCode, 'pyneb_code'] = self.linesDb.loc[idx_numeric_pynebCode, 'pyneb_code']#.astype(int)
        self.linesDb['ion'].apply(str)

        print('\n-- Loading atoms data with PyNeb')

        # Load user atomic data references
        # TODO is the atomic data really changing?
        for reference in self.config['atomic_data_list']:
            pn.atomicData.setDataFile(reference)

        # Generate the dictionary with pyneb ions
        self.ionDict = pn.getAtomDict(atom_list=np.unique(self.linesDb.ion.values.astype(str)))

        # Preparing Hbeta references
        self.Hbeta_label = 'H1_4861A'
        self.Hbeta_wave = 4862.683
        self.Hbeta_pynebCode = '4_2'

        # Import Optical depth function
        ftauFile = ftau_coefficients if ftau_coefficients is not None else path.join(self.externalDataFolder, self.config['ftau_file'])
        self.ftau_coeffs = self.import_optical_depth_coeff_table(ftauFile)

        # Declare diagnostic ratios
        self.diagnosDict = {'ROIII': {'num' :np.array([4363]), 'den' :np.array([4959, 5007])},
                            'RSIII': {'num' :np.array([6312]), 'den' :np.array([9069, 9531])},
                            'RSII' : {'num' :np.array([6716]), 'den' :np.array([6731])}}

        # Defining temperature and density grids(These are already define at init)
        if (tempGrid is not None) and (denGrid is not None):
            self.tempRange = np.linspace(tempGrid[0], tempGrid[1], tempGrid[2])
            self.denRange = np.linspace(denGrid[0], denGrid[1], denGrid[2])
            X, Y = np.meshgrid(self.tempRange, self.denRange)
            self.tempGridFlatten, self.denGridFlatten = X.flatten(), Y.flatten()

        return

    def import_emission_line_data(self, obj_lines_df, input_lines = 'all'):

        # If wavelengths are provided for the observation we use them, else we use the theoretical values
        if all(obj_lines_df.obs_wavelength.isnull().values):
            idx_obs_labels = self.linesDb.index.isin(obj_lines_df.index)
            obj_lines_df['obs_wavelength'] = self.linesDb.loc[idx_obs_labels].wavelength

        # Sort the dataframe by wavelength in case it isn't
        obj_lines_df.sort_values(by=['obs_wavelength'], ascending=True, inplace=True)

        # Get the references for the lines treatment
        idx_obj_lines = self.linesDb.index.isin(obj_lines_df.index)
        obj_lines_df['ion'] = self.linesDb.loc[idx_obj_lines].ion.astype(str)
        obj_lines_df['emis_type'] = self.linesDb.loc[idx_obj_lines].emis_type.astype(str)
        obj_lines_df['pynebCode'] = self.linesDb.loc[idx_obj_lines].pyneb_code
        obj_lines_df['emis_type'] = self.linesDb.loc[idx_obj_lines].emis_type

        # Get lines for the analysis (Excluding Hbeta) #  TODO add check of possible lines w.r.t spectrum
        if input_lines is 'all':
            idx_lines = (obj_lines_df.index != 'H1_4861A')
        else:
            idx_lines = (obj_lines_df.index.isin(input_lines)) & (obj_lines_df.index != 'H1_4861A')

        # Get blended lines

        # Declare lines data
        self.obj_data['lineLabels'] = obj_lines_df.loc[idx_lines].index.values
        self.obj_data['lineIons'] = obj_lines_df.loc[idx_lines].ion.values
        self.obj_data['lineWaves'] = obj_lines_df.loc[idx_lines].obs_wavelength.values
        self.obj_data['lineFluxes'] = obj_lines_df.loc[idx_lines].obs_flux.values
        self.obj_data['lineErr'] = obj_lines_df.loc[idx_lines].obs_fluxErr.values
        self.obj_data['lineType'] = obj_lines_df.loc[idx_lines].emis_type.values
        self.obj_data['linePynebCode'] = obj_lines_df.loc[idx_lines].pynebCode.values

        # Assign flux equation to each line
        self.eqLabelArray = self.assignFluxEq2Label(self.obj_data['lineLabels'], self.obj_data['lineIons'])

        # Read emissivity coefficients if they are available
        # TODO are check for missing lines
        self.emisCoeffs = {}
        for line in self.obj_data['lineLabels']:
            if line + '_coeffs' in self.config:
                self.emisCoeffs[line] = self.config[line + '_coeffs']

        return

    def gasSamplerVariables(self, lineIons, high_temp_ions, linesFlux=None, linesErr=None, lineLabels=None, lineFlambda=None, normalized_by_Hbeta=True, linesMinimumError=None):

        # TODO Add normalization check here
        # Emission line fluxes and errors normalization
        if linesFlux is not None:
            if normalized_by_Hbeta:
                self.obsLineFluxes = linesFlux
                self.obsLineFluxErr = linesErr
            else:
                self.obsLineFluxes = linesFlux/self.obj_data['flux_hbeta']
                self.obsLineFluxErr = linesErr/self.obj_data['flux_hbeta']

            # Setting minimin error to 0.01 of the line flux
            self.fitLineFluxErr = np.copy(self.obsLineFluxErr)
            if linesMinimumError is not None:
                err_fraction = linesErr/linesFlux
                idcs_smallErr = err_fraction < linesMinimumError
                self.fitLineFluxErr[idcs_smallErr] = linesMinimumError * self.obsLineFluxes[idcs_smallErr]

        # Logic to distinguish the lines
        if lineIons is not None:

            # Type of lines
            self.H1_lineIdcs = (lineIons == 'H1r')
            self.He1_lineIdcs = (lineIons == 'He1r')
            self.He2_lineIdcs = (lineIons == 'He2r')

            #Lines with special flux equations
            self.corrFluxLines = (lineLabels == 'O2_7319A_b')

            # Establish index of lines which below to high and low ionization zones
            self.idx_highU = np.in1d(lineIons, high_temp_ions)

            # Unique atoms in the chemical analysis
            uniqueAtoms = np.unique(lineIons)
            self.obsAtoms = uniqueAtoms[uniqueAtoms!='H1r']

            # Attributes to increase calculation speed
            self.rangeLines = np.arange(lineIons.size)
            self.rangeObsAtoms = np.arange(self.obsAtoms.size)

            # Check gas components
            self.H1rCheck = self.checkIonObservance('H1r', self.obsAtoms)
            self.He1rCheck, self.He2rCheck = self.checkIonObservance('He1r', self.obsAtoms), self.checkIonObservance('He2r', self.obsAtoms)
            self.Ar3Check, self.Ar4Check = self.checkIonObservance('Ar3', self.obsAtoms), self.checkIonObservance('Ar4', self.obsAtoms)
            self.O2Check, self.O3Check = self.checkIonObservance('O2', self.obsAtoms), self.checkIonObservance('O3', self.obsAtoms)
            self.N2Check = self.checkIonObservance('N2', self.obsAtoms)
            self.S2Check, self.S3Check = self.checkIonObservance('S2', self.obsAtoms), self.checkIonObservance('S3', self.obsAtoms)

        # Lines data
        self.lineLabels = lineLabels
        self.lineIons = lineIons
        self.lineFlambda = lineFlambda

        # # Define diagnostic ratios
        # self.obsDiagRatios = np.empty(self.diagnosRatios.size)
        # for i in range(self.diagnosRatios.size):
        #
        #     ratio = self.diagnosRatios[i]
        #
        #     # Get the numerator and denominator wave indeces
        #     num_waves = self.diagnosDict[ratio]['num']
        #     den_waves = self.diagnosDict[ratio]['den']
        #
        #     # New indices
        #     num_idcs = np.argwhere(np.in1d(self.obj_data['linePynebCode'], num_waves))
        #     den_idcs = np.argwhere(np.in1d(self.obj_data['linePynebCode'], den_waves))
        #
        #     # New diagnostics
        #     self.obsDiagRatios[i] = np.log10(self.obsLineFluxes[num_idcs].sum() / self.obsLineFluxes[den_idcs].sum())

        return

    def fitEmissivityPlane(self, ions_list, labels_list, configuration_file = None):

        # Dictionary to store the emissivity surface coeffients
        self.emisCoeffs = {}
        for i in range(labels_list.size):

            lineLabel = labels_list[i]

            # Get equation type to fit the emissivity
            line_func  = self.ionEmisEq_fit[lineLabel]
            n_args = len(getargspec(line_func).args) - 2 #TODO Not working in python 2.7 https://stackoverflow.com/questions/847936/how-can-i-find-the-number-of-arguments-of-a-python-function

            # Compute emissivity functions coefficients
            emis_grid_i = self.emis_dict[lineLabel]
            p0 = self.epm2017_emisCoeffs[lineLabel] if lineLabel in self.epm2017_emisCoeffs else np.zeros(n_args)
            p1, cov1 = self.fitEmis(line_func, (self.tempGridFlatten, self.denGridFlatten), emis_grid_i, p0=p0)
            self.emisCoeffs[lineLabel] = p1

        # Save the coefficients to the configuration file
        coeffs_dict = OrderedDict()
        for label_i in range(labels_list.size):

            lineLabel = self.obj_data['lineLabels'][label_i]

            coeffs_dict[lineLabel] = self.emisCoeffs[lineLabel]

        parseDataFile(configuration_file, 'line_emissivities', coeffs_dict, 'lists', '_coeffs')

        return

    def fitEmissivityDiagnosPlane(self, diagnRatios, emissivityRatiosGrid):

        # Temperature and density meshgrids #
        # TODO this XX YY should go to the self.
        X, Y = np.meshgrid(self.tem_grid_range, self.den_grid_range)
        XX, YY = X.flatten(), Y.flatten()

        # Dictionary to store the emissivity surface coeffients
        emisRatioCoeffs = {}
        for i in range(diagnRatios.size):

            ratioLabel = diagnRatios[i]

            # Get equation type to fit the emissivity
            diagnoFunct = self.EmisRatioEq_fit[ratioLabel]

            # Compute emissivity functions coefficients
            p1, cov1 = self.fitEmis(diagnoFunct, (XX, YY), emissivityRatiosGrid[:, i], p0=np.array([0,0,0]))
            emisRatioCoeffs[ratioLabel] = p1

            print ratioLabel, p1

        return emisRatioCoeffs

    def computeEmissivityGrids(self, pynebcode_list, ions_list, grids_folder, loadgrids_check = False):

        # Empty grid array
        emisGrid = np.empty((self.tempGridFlatten.size, ions_list.size))
        Hbeta_emis_grid = self.ionDict['H1r'].getEmissivity(self.tempGridFlatten, self.denGridFlatten, wave=4861, product=False)

        for i in range(len(ions_list)):

            # Line emissivity references
            line_label = '{}_{}'.format(ions_list[i], pynebcode_list[i])
            grid_address = path.join(grids_folder, line_label + '.' + 'npy')

            print line_label, pynebcode_list[i]

            # Check if grids are available # TODO check if we are reloading the data
            if path.isfile(grid_address) and loadgrids_check:
                emis_grid_i = np.load(grid_address)

            # Otherwise generate it (and save it)
            else:
                emis_grid_i = self.ionDict[ions_list[i]].getEmissivity(self.tempGridFlatten, self.denGridFlatten, wave=pynebcode_list[i], product=False)
                np.save(grid_address, emis_grid_i)

            # Save along the number of points
            emisGrid[:, i] = np.log10(emis_grid_i/Hbeta_emis_grid)

        return emisGrid

    def computeEmissivityDict(self, pynebcode_list, ions_list, labels_list, grids_folder, loadgrids_check = False):

        # Empty grid array
        emisDict = OrderedDict()
        Hbeta_emis_grid = self.ionDict['H1r'].getEmissivity(self.tempGridFlatten, self.denGridFlatten, wave=4861, product=False)

        for i in range(len(labels_list)):

            # Line emissivity references
            line_label = labels_list[i]
            grid_address = path.join(grids_folder, line_label + '.' + 'npy')

            # Check if grids are available # TODO check if we are reloading the data
            if path.isfile(grid_address) and loadgrids_check:
                emis_grid_i = np.load(grid_address)

            # Otherwise generate it (and save it)
            else:

                # Check if it is a blended line:
                if '_b' not in line_label:
                    # TODO I should change wave by label
                    emis_grid_i = self.ionDict[ions_list[i]].getEmissivity(self.tempGridFlatten, self.denGridFlatten,
                                                                           wave=float(pynebcode_list[i]), product=False)
                else:
                    # TODO upgrade this mechanic for different ions blending
                    for component in self.blendedLinesDict[line_label]:
                        component_wave = float(self.linesDb.loc[component].pyneb_code)
                        emis_grid_i += self.ionDict[ions_list[i]].getEmissivity(self.tempGridFlatten, self.denGridFlatten,
                                                                           wave=component_wave, product=False)

                np.save(grid_address, emis_grid_i)

            # Save along the number of points
            emisDict[line_label] = np.log10(emis_grid_i/Hbeta_emis_grid)

        return emisDict

    def computeDiagnosGrids(self, pynebcodeWaves, diagnosDict, emisRatioGrid):

        diagnosGridDict = {}
        for diagnos in diagnosDict:

            # Get the numerator and denominator wave indeces
            num_waves = diagnosDict[diagnos]['num']
            den_waves = diagnosDict[diagnos]['den']

            # Check if all the ratio lines are available
            num_set = np.in1d(pynebcodeWaves, num_waves)
            den_set = np.in1d(pynebcodeWaves, den_waves)

            #Sum up elements in numerator and denominator and divide
            if (num_waves.size == num_set.sum()) and (den_waves.size == den_set.sum()):

                num_idcs = np.argwhere(num_set)
                den_idcs = np.argwhere(den_set)

                num_array = np.power(10, emisRatioGrid[:, num_idcs])
                den_array = np.power(10, emisRatioGrid[:, den_idcs])

                diagnosGridDict[diagnos] = np.log10(num_array.sum(axis=1) / den_array.sum(axis=1))

        # Protocol in case no line ratio is available
        if len(diagnosGridDict) > 0:

            # Move data from dict to an array
            diagnosList = np.array(diagnosGridDict.keys())
            diagnosGrid = np.empty((emisRatioGrid.shape[0], diagnosList.size))

            for i in range(diagnosList.size):
                diagnosGrid[:, i] = np.squeeze(diagnosGridDict[diagnosList[i]])

        else:
            diagnosList, diagnosGrid = None, None

        return diagnosList, diagnosGrid

    def calcEmFluxes(self, Tlow, Thigh, ne, cHbeta, tau, abund_dict, fluxEqDict, fluxContainer, tensorCheck=False):

        for i in self.rangeLines:

            # Line data
            line_label = self.lineLabels[i]
            line_ion = self.lineIons[i]
            line_flambda = self.lineFlambda[i]
            line_fluxEqLabel = self.eqLabelArray[i]

            # Parameters to compute the emissivity
            line_coeffs = self.emisCoeffs[line_label]
            emis_func = self.ionEmisEq[line_label]

            # Appropriate data for the ion
            Te_calc = Thigh if self.idx_highU[i] else Tlow

            # Line Emissivity
            line_emis = emis_func((Te_calc, ne), *line_coeffs)

            # Atom abundance
            line_abund = 0.0 if self.H1_lineIdcs[i] else abund_dict[line_ion]

            # ftau correction for HeI lines
            line_ftau = self.ftau_func(tau, Te_calc, ne, *self.ftau_coeffs[line_label]) if self.He1_lineIdcs[i] else None

            # Line flux with special correction:
            if self.corrFluxLines[i] :
                fluxEq_i = fluxEqDict[line_fluxEqLabel](line_emis, cHbeta, line_flambda, line_abund, abund_dict['O3'], Thigh)

            # Classical line flux
            else:
                fluxEq_i = fluxEqDict[line_fluxEqLabel](line_emis, cHbeta, line_flambda, line_abund, line_ftau, continuum=0.0)

            # Store in container
            if tensorCheck:
                fluxContainer = storeValueInTensor(i, fluxEq_i, fluxContainer)
            else:
                fluxContainer[i] = fluxEq_i

        return fluxContainer

    def calcEmFluxes_pyneb(self, Tlow, Thigh, ne, cHbeta, tau, abund_dict):

        fluxContainer = np.zeros(self.lineLabels.size)

        for i in self.rangeLines:

            # Line data
            line_label = self.lineLabels[i]
            line_ion = self.lineIons[i]
            line_flambda = self.lineFlambda[i]
            line_PynebCode = self.linePnCode[i]

            # Appropriate data for the ion
            Te_calc = Thigh if self.idx_highU[i] else Tlow

            # Hbeta emissivity

            # Line relative emissivity
            Hbeta_emis = self.ionDict['H1r'].getEmissivity(tem=Te_calc, den=ne, wave=4861)
            line_emis = self.ionDict[line_ion].getEmissivity(tem=Te_calc, den=ne, wave=line_PynebCode)
            emisRatio = line_emis/Hbeta_emis

            # Atom abundance (it assumes metal abundance is in 12 + log format)
            line_abund = 1.0 if self.H1_lineIdcs[i] else abund_dict[line_ion]
            if line_ion not in ['H1r', 'He1r', 'He2r']:
                line_abund = np.power(10, line_abund - 12)

            # ftau correction for HeI lines
            line_ftau = self.ftau_func(tau, Te_calc, ne, *self.ftau_coeffs[line_label]) if self.He1_lineIdcs[i] else 1.0

            # Line synthetic flux
            fluxContainer[i] = line_abund * emisRatio * line_ftau * np.power(10, -cHbeta * line_flambda)

        return fluxContainer

    def calculate_recomb_fluxes(self, Thigh, ne, cHbeta, tau, lines_abund_dict, lines_labels, lines_ions, lines_flambda):

        #Emissivity for the two temperature layers
        t4 = Thigh / 10000.0

        #Interpolate the emissivities
        emis_module = bilinear_interpolator_axis(ne, Thigh, self.den_grid_range, self.tem_grid_range, self.recomb_emis_grid)

        #Loop throught the lines to generate the f_tau correction
        lines_ftau_vector = np.empty(self.n_recombLines)
        lines_abund_vector = np.empty(self.n_recombLines)

        for i in self.range_recombLines:
            if lines_ions[i] == 'He1':
                lines_ftau_vector[i]    = self.OpticalDepth_He(tau = tau, T_4 = t4, ne = ne, He_label = lines_labels[i])
            else:
                lines_ftau_vector[i]    = 1.0
            lines_abund_vector[i] = lines_abund_dict[lines_ions[i]]

        #Reddening factor for all the lines
        f_module = np.power(10, -1 * lines_flambda * cHbeta)

        #Metals flux
        recomb_fluxes = lines_abund_vector * emis_module * lines_ftau_vector * f_module

        return recomb_fluxes

    def calculate_colExcit_flux(self, Tlow, Thigh, ne, cHbeta, lines_abund_dict, lines_waves, lines_ions, lines_flambda):

        # Array to store the emissivities
        lines_emis_vector = np.empty(self.n_colExcLines)

        #Calculate low ionization emissivities
        lines_emis_vector[self.idx_lowU] = bilinear_interpolator_axis(ne, Tlow, self.den_grid_range,
                                                                      self.tem_grid_range,
                                                                      self.metals_emis_grid[:, :, self.idx_lowU])

        #Calculate high ionization emissivities
        lines_emis_vector[self.idx_highU] = bilinear_interpolator_axis(ne, Thigh, self.den_grid_range,
                                                                       self.tem_grid_range,
                                                                      self.metals_emis_grid[:, :, self.idx_highU])

        # Loop through the lines to assign the ionic abundances
        abund_vector = np.empty(self.n_colExcLines)
        for i in self.range_colExcLines:
            ion = lines_ions[i]
            abund_vector[i] = lines_abund_dict[ion]
            #print ion, abund_vector[i], lines_emis_vector[i], abund_vector[i] * lines_emis_vector[i]

        # Reddening factor for all the lines
        f_module = np.power(10, -1 * lines_flambda * cHbeta)

        # Metals flux
        colExcit_fluxes = abund_vector * lines_emis_vector * f_module

        return colExcit_fluxes