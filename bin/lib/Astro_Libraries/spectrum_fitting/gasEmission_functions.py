import numpy as np
import pyneb as pn
from os import path
from bisect import bisect_left
from scipy.optimize import curve_fit
from DZ_LineMesurer import LineMesurer_v2
from numpy import ones, exp, zeros


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


class EmissivitySurfaceFitter():

    def __init__(self):

        # TODO we should read this from the xlsx file
        self.ionEmisEq = {'S2_6716A'  : self.emisEquation_TeDe,
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
                            'N2_6548A'  : self.emisEquation_Te,
                            'N2_6584A'  : self.emisEquation_Te,
                            'H1_4102A'  : self.emisEquation_HI,
                            'H1_4341A'  : self.emisEquation_HI,
                            'H1_6563A'  : self.emisEquation_HI,
                            'He1_4026A' : self.emisEquation_HeI,
                            'He1_4471A' : self.emisEquation_HeI,
                            'He1_5876A' : self.emisEquation_HeI,
                            'He1_6678A' : self.emisEquation_HeI,
                            'He1_7065A' : self.emisEquation_HeI,
                            'He2_4686A' : self.emisEquation_HeII}

        return

    def emisEquation_Te(self, xy_space, a, b, c):
        temp_range, den_range = xy_space
        return a + b / temp_range + c * np.log10(temp_range)

    def emisEquation_TeDe(self, xy_space, a, b, c, d, e):
        temp_range, den_range = xy_space
        return a + b / temp_range + c * np.log10(temp_range) + np.log10(1 + e * den_range)

    def emisEquation_HI(self, xy_space, a, b, c):
        temp_range, den_range = xy_space
        return a + b * np.log10(temp_range) + c * np.log10(temp_range) * np.log10(temp_range)

    def emisEquation_HeI(self, xy_space, a, b, c, d):
        temp_range, den_range = xy_space
        return (c + b * den_range) * np.log10(temp_range) - np.log10(a + b * den_range)

    def emisEquation_HeII(self, xy_space, a, b):
        temp_range, den_range = xy_space
        return a * np.power(temp_range, b)

    def fitEmis(self, func_emis, xy_space, line_emis, p0 = None):
        p1, p1_cov = curve_fit(func_emis, xy_space, line_emis, p0)
        return p1, p1_cov

class EmissionEquations():

    def __init__(self):

        self.fluxEq = {}

    def assignFluxEq2Label(self, labelsList):

        for i in range(len(labelsList)):

            lineLabel = labelsList[i]

            if 'H1r' in lineLabel:
                eqFlux = self.H1_lines

            elif 'He1r' in lineLabel:
                eqFlux = self.He1_lines

            elif 'He2r' in lineLabel:
                eqFlux = self.He2_lines

            else:
                eqFlux = self.metal_lines

            self.fluxEq[lineLabel] = eqFlux

    def H1_lines(self, emis_ratio, cHbeta, flambda, abund=None, ftau=None, continuum=None):
        return emis_ratio - flambda * cHbeta

    def He1_lines(self, emis_ratio, cHbeta, flambda, abund, ftau=None, continuum=None):
        return abund + emis_ratio + ftau - flambda * cHbeta

    def He2_lines(self,emis_ratio, cHbeta, flambda, abund, ftau=None, continuum=None):
        return abund + emis_ratio - flambda * cHbeta

    def metal_lines(self, abund, emis_ratio, flambda, cHbeta, ftau=None, continuum=None):
        return abund + emis_ratio - flambda * cHbeta - 12

class EmissionComponents(EmissivitySurfaceFitter, EmissionEquations, LineMesurer_v2):

    def __init__(self, dazer_path):

        # Tools to fit a the emissivity to a plane
        EmissivitySurfaceFitter.__init__(self)

        # Functions to compute fluxes
        EmissionEquations.__init__(self)

        # Load tools to measure lines
        LineMesurer_v2.__init__(self,  dazer_path, 'DZT_LineLog_Headers.dz')

        #Hbeta configuration
        self.Hbeta_label = 'H1_4861A'
        self.Hbeta_wave = 4862.683
        self.Hbeta_pynebCode = '4_2'

        #Import Optical depth function
        self.ftau_coeffs = self.import_optical_depth_coeff_table()

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

    def import_emission_line_data(self, obj_lines_df, input_lines = 'all'):

        # If wavelengths are provided for the observation we use them, else we use the theoretical values
        if all(obj_lines_df.obs_wavelength.isnull().values):
            idx_obs_labels = self.lines_df.index.isin(obj_lines_df.index)
            obj_lines_df['obs_wavelength'] = self.lines_df.loc[idx_obs_labels].wavelength

        # Sort the dataframe by wavelength in case it isn't
        obj_lines_df.sort_values(by=['obs_wavelength'], ascending=True, inplace=True)

        # Get the references for the lines treatment
        idx_obj_lines = self.lines_df.index.isin(obj_lines_df.index)
        obj_lines_df['emis_type'] = self.lines_df.loc[idx_obj_lines].emis_type.astype(str)
        obj_lines_df['pynebCode'] = self.lines_df.loc[idx_obj_lines].pyneb_code
        obj_lines_df['ion'] = self.lines_df.loc[idx_obj_lines].ion.astype(str)
        obj_lines_df['emis_type'] = self.lines_df.loc[idx_obj_lines].emis_type

        # Save the data in the dictionary (Excluding Hbeta)
        print input_lines
        if input_lines is not 'all':
            idx_lines = (obj_lines_df.index.isin(input_lines)) & (obj_lines_df.index != 'H1_4861A')
        else:
            idx_lines = (obj_lines_df.index != 'H1_4861A')

        self.obj_data['lineLabels'] = obj_lines_df.loc[idx_lines].index.values
        self.obj_data['lineIons'] = obj_lines_df.loc[idx_lines].ion.values
        self.obj_data['lineWaves'] = obj_lines_df.loc[idx_lines].obs_wavelength.values
        self.obj_data['linePynebCode'] = obj_lines_df.loc[idx_lines].pynebCode.values
        self.obj_data['lineFluxes'] = obj_lines_df.loc[idx_lines].obs_flux.values
        self.obj_data['lineErr'] = obj_lines_df.loc[idx_lines].obs_fluxErr.values
        self.obj_data['lineType'] = obj_lines_df.loc[idx_lines].emis_type.values

        # Logic to distinguish the lines
        self.H1_lines = (obj_lines_df.loc[idx_lines].ion == 'H1r').values
        self.He1_lines = (obj_lines_df.loc[idx_lines].ion == 'He1r').values
        self.He2_lines = (obj_lines_df.loc[idx_lines].ion == 'He2r').values

        # Assign flux equation to each line
        self.assignFluxEq2Label(self.obj_data['lineLabels'])

        # Generate the dictionary with pyneb ions
        print('-- Loading atoms data with PyNeb') #TODO this must go to the lowest order to the location where we get the configuration data
        self.ionDict = pn.getAtomDict(atom_list=obj_lines_df.ion.values)

        print('-- Atomic sources Loaded ')
        # for atom in self.ionDict.keys():
        #     textPrint = '--- {}: {}'.format(atom, self.ionDict[atom].printSources())
        #     print(textPrint)

        # Establish index of lines which below to high and low ionization zones
        self.idx_highU = np.in1d(self.obj_data['lineIons'], self.high_temp_ions)
        self.idx_lowU = ~self.idx_highU

        # Attributes to increase calculation speed
        self.n_recombLines, self.n_colExcLines = np.sum(self.obj_data['lineType'] == 'rec'), np.sum(self.obj_data['lineType'] == 'col')
        self.range_recombLines, self.range_colExcLines = np.arange(self.n_recombLines), np.arange(self.n_colExcLines)
        self.n_lines, self.range_lines = self.obj_data['lineLabels'].size, np.arange(self.obj_data['lineLabels'].size)

        # Empty Ionic dictionary for the MCMC
        self.abund_dict = {ion: 0 for ion in self.ionDict.keys()}

        return

    def import_emissivity_grids(self, forceGridsReset=True):

        # Emissivity grid
        self.emis_grid = self.manage_emissivity_grids(self.obj_data['linePynebCode'], self.obj_data['lineIons'], forceGridsReset)

        return

    def fit_emissivity_surface(self):

        # Temperature and density meshgrids
        X, Y = np.meshgrid(self.tem_grid_range, self.den_grid_range)
        XX, YY = X.flatten(), Y.flatten()

        # Dictionary to store the emissivity surface coeffients
        self.emisCoeffs = {}
        for i in range(self.obj_data['lineLabels'].size):

            lineLabel = self.obj_data['lineLabels'][i]

            # Get equation type to fit the emissivity
            line_func  = self.ionEmisEq[lineLabel]

            # Compute emissivity functions coefficients
            p1, cov1 = self.fitEmis(line_func, (XX, YY), self.emis_grid[:, i])
            self.emisCoeffs[lineLabel] = p1

        return

    def manage_emissivity_grids(self, pynebcode_list, ions_list, forceGridReset=False):

        # Temperature and density meshgrids
        X, Y = np.meshgrid(self.tem_grid_range, self.den_grid_range)
        XX, YY = X.flatten(), Y.flatten()

        # Empty grid array
        emisGrid = np.empty((XX.size, ions_list.size))
        Hbeta_emis_grid = self.ionDict['H1r'].getEmissivity(XX, YY, wave=4861, product=False)

        for i in range(len(ions_list)):

            # Line emissivity references
            line_label = '{}_{}'.format(ions_list[i], pynebcode_list[i])
            grid_address = self.paths_dict['emisGridsPath'] + line_label + '.npy'

            # Check if grids are available
            if path.isfile(grid_address) and forceGridReset is False:
                emis_grid_i = np.load(grid_address)

            # Otherwise generate it (and save it)
            else:
                print('--- Generating grid: {}'.format(line_label))
                emis_grid_i = self.ionDict[ions_list[i]].getEmissivity(XX, YY, wave=pynebcode_list[i], product=False)
                np.save(grid_address, emis_grid_i)

            # Save along the number of points
            emisGrid[:, i] = np.log10(emis_grid_i/Hbeta_emis_grid)

        return emisGrid

    def calcEmFluxes(self, Tlow, Thigh, ne, cHbeta, tau, abund_dict, lineLabels, lineIons, lineFlambda):

        linesFlux = zeros(lineLabels.size)

        for i in self.range_lines:

            # Line data
            line_label      = lineLabels[i]
            line_ion        = lineIons[i]
            line_flambda    = lineFlambda[i]

            # Parameters to compute the emissivity
            line_coeffs     = self.emisCoeffs[line_label]
            emis_func       = self.ionEmisEq[line_label]

            # Appropiate data for the ion
            Te_calc = Thigh if self.idx_highU[i] else Tlow

            # Line Emissivitiy
            line_emis = emis_func((Te_calc, ne), *line_coeffs)

            # Atom abundance
            line_abund = 1.0 if self.H1_lines[i] else abund_dict[line_ion]

            # ftau correction for HeI lines
            line_ftau = self.ftau_func(tau, Te_calc, ne, *self.ftau_coeffs[line_label]) if self.He1_lines[i] else None

            # Line synthetic flux
            fluxEq_i = self.fluxEq[line_label](line_emis, cHbeta, line_flambda, line_abund, line_ftau, continuum=None)

            # Store in container
            linesFlux[i] = fluxEq_i

        return linesFlux

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