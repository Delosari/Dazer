import numpy as np
from os import path
from pandas import read_csv
from lib.Astro_Libraries.spectrum_fitting.plot_tools import MCMC_printer
from collections import OrderedDict
from lib.Astro_Libraries.spectrum_fitting.import_functions import ImportModelData, parseObjData, make_folder
from lib.Astro_Libraries.spectrum_fitting.starContinuum_functions import SspFitter, CCM89_Bal07
from lib.Astro_Libraries.spectrum_fitting.gasContinuum_functions import NebularContinuaCalculator
from lib.Astro_Libraries.spectrum_fitting.gasEmission_functions import TOIII_TSIII_relation, EmissionComponents
from lib.Astro_Libraries.spectrum_fitting.extinction_tools import ReddeningLaws


class ModelIngredients(ImportModelData, SspFitter, NebularContinuaCalculator, EmissionComponents, ReddeningLaws, MCMC_printer):

    def __init__(self):

        # Load the default configuration
        ImportModelData.__init__(self, path.dirname(path.realpath(__file__)))

        # Load tools for spectra calculation
        SspFitter.__init__(self)
        EmissionComponents.__init__(self, self.config['temp_grid'], self.config['den_grid'])
        NebularContinuaCalculator.__init__(self)

        # Import extinction classes
        ReddeningLaws.__init__(self, self.config['R_v'], self.config['reddenig_curve'])

        # For generating graphs
        MCMC_printer.__init__(self)

    def gen_synth_obs(self, obs_name, output_folder,
                      obj_properties_file=None, obj_lines_file = None,
                      wavelengh_limits = None, resample_inc = None, norm_interval = None,
                      ssp_lib_type = None, ssp_folder = None, ssp_file = None,
                      obj_ssp_coeffs_file = None, error_stellarContinuum = None, error_lines = None,
                      atomic_data=None, ftau_coeffs=None):

        # Dictionary to store the data
        self.obj_data = {}

        # Store input files:
        self.obj_data['obj_properties_file'] = obj_properties_file
        self.obj_data['obj_lines_file'] = obj_lines_file
        self.obj_data['obj_ssp_coeffs_file'] = obj_ssp_coeffs_file
        self.obj_data['output_folder'] = output_folder

        # Read simulation data from file
        obj_prop_df = read_csv(obj_properties_file, delim_whitespace=True, header=0, index_col=0)

        # Read lines file
        obj_lines_df = read_csv(obj_lines_file, delim_whitespace=True, header=0, index_col=0)

        # Import the stellar library and use it to generate the observation continuum
        self.ssp_lib = self.load_ssp_library(ssp_lib_type, ssp_folder, ssp_file, wavelengh_limits, resample_inc, norm_interval)

        # Declare wavelength for the object
        z_obj, obj_WaveRest = obj_prop_df.loc['z_obj'][0],  self.ssp_lib['wave_resam']
        obj_WaveObs = obj_WaveRest * (1.0 + z_obj)

        # Generate masks for the object from the lines log
        self.generate_object_mask(obj_lines_df, obj_WaveRest, obj_lines_df.index.values)

        # ---------------------------------------- Emission lines data -------------------------------------------------
        # --------------------------------------------------------------------------------------------------------------

        # Store input data
        self.obj_data['flux_hbeta'] = obj_prop_df.loc['flux_hbeta'][0]

        # Load atomic data
        self.import_atomic_data(atomic_data, ftau_coeffs, self.config['temp_grid'], self.config['den_grid'])

        # Prepare data from emission line file (trick to put all the lines)
        self.import_emission_line_data(obj_lines_df, obj_lines_df.index.values)

        # Reddening parameters
        self.obj_data['lineFlambda'] = self.gasExtincParams(self.obj_data['lineWaves'], self.config['R_v'], self.config['reddenig_curve'])

        # # Create or load emissivity grids
        # self.emis_dict = self.computeEmissivityDict(self.obj_data['linePynebCode'], self.obj_data['lineIons'],
        #                                             self.obj_data['lineLabels'], output_folder)
        #
        # # Fit emissivity grids to a surface
        # self.fitEmissivityPlane(self.obj_data['lineIons'], self.obj_data['lineLabels'], self.configFolder)
        #
        # # Plot fits of emissivity grids
        # self.plot_emisFits(self.obj_data['lineLabels'], self.emisCoeffs, self.emis_dict, output_folder)
        #
        # # Create or emissivity diagnostic ratios
        # self.diagnosRatios, self.diagnosGrid = self.computeDiagnosGrids(self.obj_data['linePynebCode'], self.diagnosDict, self.emis_grid)
        #
        # # Fit emissivity ratios a surface
        # self.emisRatioCoeffs = self.fitEmissivityDiagnosPlane(self.diagnosRatios, self.diagnosGrid)
        #
        # Plot fits of emissivity ratios
        # self.plot_emisRatioFits(self.diagnosRatios, self.emisRatioCoeffs, self.diagnosGrid, self.paths_dict['emisGridsPath'])


        # Reddening parameters
        self.obj_data['lineFlambda'] = self.gasExtincParams(self.obj_data['lineWaves'], self.Rv_model,
                                                            self.reddedning_curve_model)

        # Variables to make the iterations simpler
        self.obj_data['T_low_true'], self.obj_data['n_e_true'] = obj_prop_df.loc['T_low_true'][0], obj_prop_df.loc['n_e_true'][0]
        self.gasSamplerVariables(self.obj_data['lineIons'], self.config['high_temp_ions'],
                                 lineLabels=self.obj_data['lineLabels'], lineFlambda=self.obj_data['lineFlambda'])

        #Dictionary with synthetic abundances
        abund_df = obj_prop_df[obj_prop_df.index.str.contains('_abund')]
        abund_keys = abund_df.index.str.rstrip('_abund').astype(str).values
        abund_values = abund_df.variable_magnitude.values
        self.abund_dict = dict(zip(abund_keys, abund_values.T))

        #Gas physical parameters
        T_low = obj_prop_df.loc['T_low_true'][0]
        n_e = obj_prop_df.loc['n_e_true'][0]
        tau = obj_prop_df.loc['tau'][0]
        cHbeta = obj_prop_df.loc['cHbeta'][0]

        #Calculate T_high assuming we get T_low
        T_high = TOIII_TSIII_relation(T_low)

        # # Prepare gas data # TODO update the prepare_gas_data to use run in this part
        # self.lineLabels = self.obj_data['lineLabels']
        # self.lineIons = self.obj_data['lineIons']
        # self.lineFlambda = self.obj_data['lineFlambda']

        # Compute lines flux
        self.lineLabels = self.obj_data['lineLabels']
        lineFluxes = self.calcEmFluxes(T_low, T_high, n_e, cHbeta, tau, self.abund_dict, self.emFluxTensors, np.zeros(len(self.obj_data['lineLabels'])))
        self.obj_data['lineFluxes'] = lineFluxes

        # Use general error if this is provided
        self.obj_data['lineErr'] = self.obj_data['lineFluxes'] * error_lines

        # Store data to synthetic observation files
        idx_lines = (obj_lines_df.index != 'H1_4861A')
        obj_lines_df.loc[idx_lines, 'obs_flux']      = self.obj_data['lineFluxes']
        obj_lines_df.loc[idx_lines, 'obs_fluxErr']   = self.obj_data['lineErr']
        obj_lines_df.loc['H1_4861A', 'obs_flux']     = 1.0
        obj_lines_df.loc['H1_4861A', 'obs_fluxErr']  = 1.0 * error_lines

        # Assign line region
        obj_lines_df['w3'] = np.floor(obj_lines_df['obs_wavelength'].values * 0.998)
        obj_lines_df['w4'] = np.ceil(obj_lines_df['obs_wavelength'].values * 1.002)

        # Create txt lines log
        synth_lines_log = '{}{}_lineslog.txt'.format(output_folder, obs_name)
        with open(synth_lines_log, 'w') as f:
            f.write(obj_lines_df.ix[:, :'blended_label'].to_string(float_format=lambda x: "{:15.8f}".format(x), index=True, index_names=False))


        # ---------------------------------------- Continuum calculation -----------------------------------------------
        # --------------------------------------------------------------------------------------------------------------

        # Get Halpha flux to calibrate
        idx_Halpha = (self.obj_data['lineLabels'] == 'H1_6563A')
        flux_halpha = self.obj_data['lineFluxes'][idx_Halpha][0] * obj_prop_df.loc['flux_hbeta'][0]

        # Reddening parameters for the nebular continuum
        self.obj_data['nebFlambda'] = self.gasExtincParams(obj_WaveRest, self.config['R_v'], self.config['reddenig_curve'])

        # Calculate the nebular continua
        self.obj_data['nebularFlux'] = self.nebFluxCont(obj_WaveRest,
                                                    obj_prop_df.loc['cHbeta'][0], self.obj_data['nebFlambda'],
                                                    obj_prop_df.loc['T_low_true'][0],
                                                    obj_prop_df.loc['He1r_abund'][0], obj_prop_df.loc['He2r_abund'][0],
                                                    flux_halpha)

        # Save input conditions:
        Av_star = obj_prop_df.loc['Av_star'][0]
        sigma_star = obj_prop_df.loc['sigma_star'][0]
        flux_hbeta = obj_prop_df.loc['flux_hbeta'][0]
        eqw_hbeta = obj_prop_df.loc['eqw_hbeta'][0]

        #Get populations for the stellar continua
        bases_idx, bases_coeff, bases_coeff_err = np.loadtxt(self.obj_data['obj_ssp_coeffs_file'], usecols=[0, 1, 2], unpack=True)

        # SSP library flux at modelled Av, z_star and sigma_star
        ssp_grid_model_norm = self.physical_SED_model(self.ssp_lib['wave_resam'], obj_WaveObs, self.ssp_lib['flux_norm'],
                                Av_star, z_obj, sigma_star, Rv_coeff=self.config['R_v'])

        # Normalized object flux
        stellarFluxNorm = ssp_grid_model_norm.dot(bases_coeff)

        # Instead of denormalizing the grid we use the Hbeta flux and equivalent width
        # To get a better shape of the nebular and stellar continua
        flux_hbeta, eqw_hbeta = obj_prop_df.loc['flux_hbeta'][0], obj_prop_df.loc['eqw_hbeta'][0]
        cont_hbeta = flux_hbeta / eqw_hbeta
        stellarFlux = stellarFluxNorm * cont_hbeta

        # Generate synth error # TODO We should play with this error
        sigma_err = error_stellarContinuum * np.median(stellarFlux)
        self.obj_data['stellarFlux'] = stellarFlux + np.random.normal(0.0, sigma_err, stellarFlux.size)

        #Save synthetic continua according to the observations
        if 'nebularFlux' not in self.obj_data:
            self.obj_data['obsFlux'] = self.obj_data['stellarFlux']
        else:
            self.obj_data['obsFlux'] = self.obj_data['nebularFlux'] + self.obj_data['stellarFlux']

        # Store the spectrum as a text file
        synth_spectrum_address = '{}{}_spectrum.txt'.format(output_folder, obs_name)
        np.savetxt(synth_spectrum_address, np.transpose(np.array([obj_WaveObs, self.obj_data['obsFlux']])), fmt="%7.1f %10.4e")

        #Store synthetic object data log # TODO Add here the comments of the function and make it automatic
        obj_dict = OrderedDict()
        obj_dict['address_lines_log']   = synth_lines_log
        obj_dict['address_spectrum']    = synth_spectrum_address
        obj_dict['flux_hbeta']          = flux_hbeta
        obj_dict['flux_halpha']         = flux_halpha
        obj_dict['Normalized_by_Hbeta'] = True
        obj_dict['eqw_hbeta']           = eqw_hbeta
        obj_dict['sigma_gas']           = obj_prop_df.loc['sigma_gas'][0]
        obj_dict['T_low']               = T_low
        obj_dict['T_high']              = T_high
        obj_dict['n_e']                 = n_e
        obj_dict['cHbeta']              = cHbeta
        obj_dict['Av_star']             = Av_star
        obj_dict['sigma_star']          = sigma_star
        obj_dict['z_obj']               = z_obj
        obj_dict['continuum_sigma']     = error_stellarContinuum
        obj_dict['resample_inc']        = resample_inc
        obj_dict['wavelengh_limits']    = wavelengh_limits
        obj_dict['norm_interval']       = norm_interval
        obj_dict['Te_prior']            = [10000.0, 1000.0]
        obj_dict['ne_prior']            = [200, 100]
        obj_dict['cHbeta_prior']        = [0.125, 0.02]
        obj_dict['T_low_true']          = T_low
        obj_dict['T_high_true']         = T_high
        obj_dict['n_e_true']            = n_e
        obj_dict['cHbeta_true']         = cHbeta
        obj_dict['tau_true']            = tau
        obj_dict['He1r_true']           = self.abund_dict['He1r']
        obj_dict['He2r_true']           = self.abund_dict['He2r']
        obj_dict['S2_true']             = self.abund_dict['S2']
        obj_dict['S3_true']             = self.abund_dict['S3']
        obj_dict['Ar3_true']            = self.abund_dict['Ar3']
        obj_dict['Ar4_true']            = self.abund_dict['Ar4']
        obj_dict['O2_true']             = self.abund_dict['O2']
        obj_dict['O3_true']             = self.abund_dict['O3']
        obj_dict['N2_true']             = self.abund_dict['N2']

        #Save the data into an "ini" configuration file
        conf_address = '{}{}_objParams.txt'.format(output_folder, obs_name)
        parseObjData(conf_address, obs_name, obj_dict)

        return

    def ready_simulation(self, output_folder, obs_data, ssp_data, fitting_components, overwrite_grids=False,
                         input_lines=None, wavelengh_limits = None, resample_inc=None, norm_interval=None):

        # Dictionary to store the data
        self.obj_data = {}
        self.obj_data = obs_data.copy()

        # Prepare emission data
        if 'emission' in fitting_components:

            # Read log with observational features and masks
            obj_lines_df = read_csv(obs_data['obj_lines_file'], delim_whitespace=True, header=0, index_col=0)

            # Load atomic data
            self.import_atomic_data() # TODO Do we need this one for the nebular continuum?

            # Prepare data from emission line file (trick to put all the lines)
            self.import_emission_line_data(obj_lines_df, input_lines=input_lines)

            # Reddening parameters
            self.obj_data['lineFlambda'] = self.gasExtincParams(self.obj_data['lineWaves'], self.Rv_model, self.reddedning_curve_model)

            # # Create or load emissivity grids
            # self.emis_dict = self.computeEmissivityDict(self.obj_data['linePynebCode'], self.obj_data['lineIons'],
            #                                             self.obj_data['lineLabels'], output_folder)
            #
            # # Fit emissivity grids to a surface
            # self.fitEmissivityPlane(self.obj_data['lineIons'], self.obj_data['lineLabels'], self.configFolder)
            #
            # # Plot fits of emissivity grids
            # self.plot_emisFits(self.obj_data['lineLabels'], self.emisCoeffs, self.emis_dict, output_folder)
            #
            # # Create or emissivity diagnostic ratios
            # self.diagnosRatios, self.diagnosGrid = self.computeDiagnosGrids(self.obj_data['linePynebCode'], self.diagnosDict, self.emis_grid)
            #
            # # Fit emissivity ratios a surface
            # self.emisRatioCoeffs = self.fitEmissivityDiagnosPlane(self.diagnosRatios, self.diagnosGrid)
            #
            # # Plot fits of emissivity ratios
            # self.plot_emisRatioFits(self.diagnosRatios, self.emisRatioCoeffs, self.diagnosGrid, self.paths_dict['emisGridsPath'])

        # Trim, resample and normalize according to the input object
        if 'obs_wavelength' in self.obj_data: # TODO this is a bit dirty
            self.treat_input_spectrum(self.obj_data, self.obj_data['obs_wavelength'], self.obj_data['obs_flux'], wavelengh_limits, resample_inc, norm_interval)

        # Reddening parameters for the nebular continuum
        if 'nebular' in fitting_components:
            obj_wave_rest = self.obj_data['wave_resam'] / (1.0 + self.obj_data['z_obj'])
            self.obj_data['nebFlambda'] = self.gasExtincParams(obj_wave_rest, self.Rv_model,  self.reddedning_curve_model)

        # Declare stellar data
        if 'stellar' in fitting_components:
            self.ssp_lib = {}
            self.ssp_lib = ssp_data.copy()

        # Generate object masks
        if ('wave_resam' in self.obj_data) and ('lineLabels' in self.obj_data): # TODO this is a bit dirty
            self.generate_object_mask(obj_lines_df, self.obj_data['wave_resam'], self.obj_data['lineLabels'])

        return

    def prepareSimulation(self, obs_data, ssp_data=None, output_folder = None, storage_folder = None,
                        spectra_components=None, input_lines='all', normalized_by_Hbeta=True,
                        excludeReddening = False, T_high_prior = False, prefit_ssp = True,
                        wavelengh_limits = None, resample_inc = None, norm_interval=None):

        # Store components fit
        self.spectraComponents = spectra_components

        # Folders to store inputs and outputs
        self.input_folder       = output_folder + 'input_data/' # TODO sure?
        self.output_folder      = output_folder + 'output_data/'
        self.dataFolder         = self.output_folder if storage_folder is None else storage_folder # TODO sure?
        self.configFile         = obs_data['obsFile'] #TODO this one has to go into the configuration
        self.objName            = str(obs_data['objName'])
        self.prefit_db          = '{}{}_sspPrefitDB'.format(self.dataFolder, self.objName)
        self.sspCoeffsPrefit_file = '{}{}_prefitSSPpopulations.txt'.format(self.input_folder, self.objName)
        self.sspCoeffs_file     = '{}{}_SSPpopulations.txt'.format(self.input_folder, self.objName)

        # Create them if not available
        make_folder(self.input_folder)
        make_folder(self.output_folder)

        # Prepare spectrum components for fitting
        self.emissionCheck, self.stellarCheck, self.emissionCheck = False, False, False

        # Pre-analysis emission spectrum
        if 'emission' in self.spectraComponents:

            # Get emission data from input files
            self.ready_simulation(output_folder, obs_data, ssp_data, spectra_components, input_lines=input_lines,
                                  wavelengh_limits=wavelengh_limits, resample_inc=resample_inc, norm_interval=norm_interval)

            # Declare gas sampler variables
            self.gasSamplerVariables(self.obj_data['lineIons'], self.config['high_temp_ions'],
                                     self.obj_data['lineFluxes'], self.obj_data['lineErr'],
                                     self.obj_data['lineLabels'], self.obj_data['lineFlambda'],
                                     normalized_by_Hbeta, self.config['linesMinimumError'])

            # Prios definition # TODO this must go to the a special section
            self.priorsDict = {'T_low': self.obj_data['Te_prior'], 'n_e': self.obj_data['ne_prior']}

            # Confirm inputs are valid
            self.emissionCheck = True

        # Prefit stellar continua
        if 'stellar' in self.spectraComponents:

            self.stellarCheck = True

            # Perform a new SPP synthesis otherwise use available data
            if prefit_ssp:

                # Compute nebular continuum using normalise Halpha and standard conditions
                self.computeDefaultNeb(self.nebDefault['Te_neb'], self.obj_data['nebFlambda'], self.nebDefault['cHbeta_neb'],
                                       self.nebDefault['He1_neb'], self.nebDefault['He2_neb'],
                                       self.nebDefault['flux_halpha'] / self.obj_data['normFlux_coeff'], self.nebDefault['z_neb'])

                # Ready continuum data
                self.prepareContinuaData(self.ssp_lib['wave_resam'], self.ssp_lib['flux_norm'], self.ssp_lib['normFlux_coeff'],
                                         self.obj_data['wave_resam'], self.obj_data['flux_norm'], self.obj_data['continuum_sigma'],
                                         self.int_mask, nebularFlux=self.nebDefault['synth_neb_flux'])

                # Select model
                self.select_inference_model('stelar_prefit')

                # Plot input simulation data
                self.plotInputSSPsynthesis()

                # Run stellar continua prefit and store/print the results
                #self.run_pymc(self.prefit_db, iterations=8000, variables_list=['Av_star', 'sigma_star'], prefit = True)
                self.savePrefitData(self.sspCoeffsPrefit_file, self.prefit_db)

            # Compute nebular continuum using prior physical data
            self.computeDefaultNeb(self.nebDefault['Te_neb'], self.obj_data['nebFlambda'], self.nebDefault['cHbeta_neb'],
                                   self.nebDefault['He1_neb'], self.nebDefault['He2_neb'],
                                   self.obj_data['flux_halpha'] / self.obj_data['normFlux_coeff'], self.nebDefault['z_neb'])

            # Compute nebular continuum using normalise Halpha and standard conditions
            # TODO I think I need to remove nebular continuum here if I want to add it later
            self.prepareContinuaData(self.ssp_lib['wave_resam'], self.ssp_lib['flux_norm'], self.ssp_lib['normFlux_coeff'],
                                     self.obj_data['wave_resam'], self.obj_data['flux_norm'],
                                     self.obj_data['continuum_sigma'],
                                     self.int_mask,
                                     nebularFlux=None,#self.nebDefault['synth_neb_flux'],
                                     mainPopulationsFile=self.sspCoeffsPrefit_file)

        return


    def prepareContinuaData(self, basesWave, basesFlux, basesFluxCoeffs, obsWave, obsFlux, obsFluxEr, objMask, nebularFlux = None, mainPopulationsFile = None):

        # Remove nebular contribution from observed continuum
        if nebularFlux is not None:
            inputContinuum = obsFlux - nebularFlux
        else:
            inputContinuum = obsFlux

        # Trim the total SSP library to the main populations stored in a text file
        if mainPopulationsFile is not None:

            # Three column file with idcs populations, weight and mainPopulationsFile
            bases_idx, bases_coeff, bases_coeff_err = np.loadtxt(mainPopulationsFile, usecols=[0, 1, 2], unpack=True)

            # Include ssps above minimum value TODO we should check this limit thing currentl it must be 0.001
            idx_populations = bases_coeff >= self.lowlimit_sspContribution

            # Data for the analysis # TODO not sure if this one should be here
            self.stellarAv_prior = self.obj_data['Av_prefit'][0],  self.obj_data['Av_prefit'][1]
            self.stellarSigma_prior = self.obj_data['sigma_star_prefit'][0],  self.obj_data['sigma_star_prefit'][1]

        else:
            bases_idx = np.arange(basesFlux.shape[0])
            bases_coeff = np.ones(basesFlux.shape[0], dtype=bool)
            bases_coeff_err = np.zeros(basesFlux.shape[0], dtype=bool)
            idx_populations = np.ones(basesFlux.shape[0], dtype=bool)

        # Input Object data
        self.obsFluxNorm            = inputContinuum
        self.inputContinuum         = inputContinuum * objMask
        self.inputContinuumEr       = obsFluxEr # TODO this will not work if the error is a scalar.. need to rethink
        self.inputWave              = obsWave

        # Populations parameters
        self.sspPrefitIdcs          = np.where(idx_populations)[0]
        self.sspPrefitCoeffs        = bases_coeff[idx_populations]
        self.sspPrefitErr           = bases_coeff_err[idx_populations]
        self.sspPrefitLimits        = np.vstack((self.sspPrefitCoeffs * 0.8, self.sspPrefitCoeffs * 1.2)).T # Theoretical limits

        #Bases parameters
        neglected_populations       = np.where(~idx_populations)
        self.onBasesWave            = basesWave
        self.onBasesFluxNorm        = np.delete(basesFlux, neglected_populations, axis=0)
        self.onBasesFluxNormCoeffs  = np.delete(basesFluxCoeffs, neglected_populations, axis=0)
        self.nBases                 = self.onBasesFluxNorm.shape[0] #self.onBasesFlux.shape[0]
        self.range_bases            = np.arange(self.nBases)

        # Limit for bases
        self.zMin_SspLimit          = np.around((obsWave[-1] / basesWave[-1]), decimals=2 - 1)
        self.zMax_SspLimit          = np.around((obsWave[0] / basesWave[0]), decimals=2 - 1)

        # Static reddening curve for the stellar continuum
        self.Xx_stellar             = CCM89_Bal07(self.Rv_model, basesWave)

        return

    def savePrefitData(self, sspCoeffs_file, ssp_db_file):

        # Read input data
        stellarPrefit_db, stat_db_dict = self.load_pymc_database_manual(ssp_db_file, burning=5000,
                                params_list=['Av_star', 'sigma_star', 'ssp_coefficients'])

        # Get mean and uncertainty values
        Av_mean, Av_std = stat_db_dict['Av_star']['trace'].mean(axis=0), stat_db_dict['Av_star']['trace'].std(axis=0)
        sigma_mean, sigma_std = stat_db_dict['sigma_star']['trace'].mean(axis=0), stat_db_dict['sigma_star']['trace'].std(axis=0)
        coeffs_mean, coeffs_std = stat_db_dict['ssp_coefficients']['trace'].mean(axis=0), stat_db_dict['ssp_coefficients']['trace'].std(axis=0)

        # File for saving the population coefficients # TODO Default files names should go into the configuration
        pops_vector = np.arange(coeffs_mean.size, dtype=int)
        np.savetxt(sspCoeffs_file, np.transpose(np.array([pops_vector, coeffs_mean, coeffs_std])), fmt="%4i %10.8f %10.8f")

        # Add results to object config file
        sectionName = self.objName + '_results'
        objData = {'Av_prefit':[Av_mean,Av_std], 'sigma_star_prefit':[sigma_mean, sigma_std], 'coeffsPop_prefit':coeffs_mean, 'coeffsPopErr_prefit':coeffs_std}
        parseObjData(self.configFile, sectionName, objData)

        # Compute mean output spectrum
        ssp_grid_i_norm = self.physical_SED_model(self.onBasesWave, self.inputWave, self.onBasesFluxNorm,
                                                  Av_mean, 0.0, sigma_mean, self.Rv_model)

        obj_ssp_fit_flux = ssp_grid_i_norm.dot(coeffs_mean)

        # Print prefit output
        self.plotOutputSSPsynthesis(stellarPrefit_db, stat_db_dict, obj_ssp_fit_flux, coeffs_mean)

        return

    def load_prefit_data(self, obj_wave):

        # Mean parameter values #TODO we need to add the redshift here
        self.sspPrefit_Coeffs               = stat_db_dict['ssp_coefficients']['trace'].mean(axis=0)
        self.sspPrefit_err                  = stat_db_dict['ssp_coefficients']['trace'].std(axis=0)
        self.ssp_lib['Av_sspPrefit']        = stat_db_dict['Av_star']['mean']
        self.ssp_lib['sigma_sspPrefit']     = stat_db_dict['sigma_star']['mean']

        return

    def computeDefaultNeb(self, Te, nebFlambda, cHbeta, He1_abund, He2_abund, fluxHalpha_norm, z_obj):

        # Generate synthetic nebular emission to remove from object
        self.nebDefault['wave_neb'] = self.obj_data['wave_resam'] / (1 + z_obj)

        self.nebDefault['synth_neb_flux'] = self.nebFluxCont(self.nebDefault['wave_neb'], cHbeta, nebFlambda, Te, He1_abund, He2_abund, fluxHalpha_norm)

        return