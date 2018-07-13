import numpy as np
from pandas import read_csv
from plot_tools import MCMC_printer
from collections import OrderedDict
from import_functions import ImportModelData, parse_synth_objData, make_folder
from starContinuum_functions import SspFitter, CCM89_Bal07
from gasContinuum_functions import NebularContinuaCalculator
from gasEmission_functions import TOIII_TSIII_relation, EmissionComponents
from extinction_tools import ReddeningLaws

class ModelIngredients(ImportModelData, SspFitter, NebularContinuaCalculator, EmissionComponents, ReddeningLaws, MCMC_printer):

    def __init__(self):

        #Import tools to load data
        ImportModelData.__init__(self)

        #Load tools for spectra calculation
        SspFitter.__init__(self)
        NebularContinuaCalculator.__init__(self, self.paths_dict['nebular_data_folder'])
        EmissionComponents.__init__(self, self.paths_dict['dazer_path'] + 'format/')

        #Import sub classes
        ReddeningLaws.__init__(self)

        #For generating graphs
        MCMC_printer.__init__(self)

        # Default configuration to read from file # TODO actually this should be uplodaded from a file
        temp_grid = [5000, 25000, 20]
        den_grid = [0, 500, 20]
        high_temp_ions = ['He1r', 'He2r', 'O3', 'Ar4']
        R_v = 3.4
        reddenig_curve = 'G03 LMC'
        lowlimit_sspContribution =  0.001

        # Default nebular continuum configuration # TODO confirm this setup
        self.nebDefault = {'Te_neb': 10000.0, 'z_neb': 0.0, 'cHbeta_neb': 0.1, 'He1_neb': 0.1, 'He2_neb': 0.001}

        # Temperature and density grid declaration
        self.tem_grid_range = np.linspace(temp_grid[0], temp_grid[1], temp_grid[2])
        self.den_grid_range = np.linspace(den_grid[0], den_grid[1], den_grid[2])
        self.den_grid_range[0] = 1 #Change the minimum value for the grid to 1

        # Reddening parameters
        self.Rv_model = R_v
        self.reddedning_curve_model = reddenig_curve

        # Declare high ionization temperature ions
        self.high_temp_ions = high_temp_ions

        # Lower accepted limit for ssps
        self.lowlimit_sspContribution = lowlimit_sspContribution

    def gen_synth_obs(self, spectra_components = None, input_ions = None,
                      wavelengh_limits = None, resample_inc = None, norm_interval = None,
                      obj_properties_file = None, obj_lines_file = None,
                      output_folder = None, obs_name = 'synthObj', ssp_lib_type = None, data_folder = None,
                      data_file = None ,obj_ssp_coeffs_file = None,
                      error_stellarContinuum = None, error_lines = None):

        # Dictionary to store the data
        self.obj_data = {}

        # Store input files:
        self.obj_data['obj_properties_file']   = obj_properties_file
        self.obj_data['obj_lines_file']        = obj_lines_file
        self.obj_data['obj_ssp_coeffs_file']   = obj_ssp_coeffs_file
        self.obj_data['output_folder']         = output_folder

        # Read simulation data from file
        obj_prop_df = read_csv(obj_properties_file, delim_whitespace=True, header=0, index_col=0)

        # Read lines file
        obj_lines_df = read_csv(obj_lines_file, delim_whitespace=True, header=0, index_col=0)

        # Import the stellar library and use it to generate the observation continuum
        self.ssp_lib = self.load_ssp_library(ssp_lib_type, data_folder, data_file, wavelengh_limits, resample_inc, norm_interval)

        # Declare wavelength for the object
        z_obj = obj_prop_df.loc['z_obj'][0]
        obj_WaveRest = self.ssp_lib['wave_resam']
        obj_WaveObs = obj_WaveRest * (1.0 + z_obj)

        # Generate masks for the object from the lines log
        self.generate_object_mask(obj_lines_df, obj_WaveRest)

        # -------- Emission lines data
        #Store input data
        self.obj_data['flux_hbeta'] = obj_prop_df.loc['flux_hbeta'][0]

        # Prepare data from emission line file (trick to put all the lines)
        self.import_emission_line_data(obj_lines_df, obj_lines_df.index.values)

        # Variables to make the iterations simpler
        self.gasSamplerVariables(self.obj_data['lineIons'], self.high_temp_ions)

        # Create or load emissivity grids
        self.import_emissivity_grids()

        # Fit emissivity grids to a surface
        self.fit_emissivity_surface()

        # Reddening parameters
        self.obj_data['lineFlambda'] = self.gasExtincParams(self.obj_data['lineWaves'], self.Rv_model, self.reddedning_curve_model)

        # Plot fits of emissivity grids
        #self.plot_emisFits(self.obj_data['lineLabels'], self.emisCoeffs, self.emis_grid, self.paths_dict['emisGridsPath'])

        #Dictionary with synthetic abundances
        abund_df        = obj_prop_df[obj_prop_df.index.str.contains('_abund')]
        abund_keys      = abund_df.index.str.rstrip('_abund').astype(str).values
        abund_values    = abund_df.variable_magnitude.values
        self.abund_dict = dict(zip(abund_keys, abund_values.T))

        #Gas physical parameters
        T_low = obj_prop_df.loc['T_low'][0]
        n_e = obj_prop_df.loc['n_e'][0]
        tau = obj_prop_df.loc['tau'][0]
        cHbeta = obj_prop_df.loc['cHbeta'][0]

        #Calculate T_high assuming we get T_low
        T_high = TOIII_TSIII_relation(T_low)

        # Compute lines flux
        logLine_fluxes = self.calcEmFluxes(T_low, T_high, n_e, cHbeta, tau, self.abund_dict,
                                    self.obj_data['lineLabels'], self.obj_data['lineIons'], self.obj_data['lineFlambda'])

        self.obj_data['lineFluxes'] = np.power(10, logLine_fluxes)

        # Use general error if this is provided
        self.obj_data['lineErr'] = self.obj_data['lineFluxes'] * error_lines

        # Store data to synthetic observation files
        idx_lines = (obj_lines_df.index != 'H1_4861A')
        obj_lines_df.loc[idx_lines,'obs_flux']      = self.obj_data['lineFluxes']
        obj_lines_df.loc[idx_lines,'obs_fluxErr']   = self.obj_data['lineErr']
        obj_lines_df.loc['H1_4861A','obs_flux']     = 1.0
        obj_lines_df.loc['H1_4861A','obs_fluxErr']  = 1.0 * error_lines

        # Assign line region
        obj_lines_df['w3'] = np.floor(obj_lines_df['obs_wavelength'].values * 0.998)
        obj_lines_df['w4'] = np.ceil(obj_lines_df['obs_wavelength'].values * 1.002)

        # Create txt lines log
        synth_lines_log = '{}{}_lineslog.txt'.format(output_folder, obs_name)
        with open(synth_lines_log, 'w') as f:
            f.write(obj_lines_df.ix[:, :'blended_label'].to_string(float_format=lambda x: "{:15.8f}".format(x), index=True, index_names=False))


        # -------- Nebular continuum calculation

        # Get Halpha flux to calibrate
        idx_Halpha = (self.obj_data['lineLabels'] == 'H1_6563A')
        flux_halpha = self.obj_data['lineFluxes'][idx_Halpha][0] * obj_prop_df.loc['flux_hbeta'][0]

        # Reddening parameters for the nebular continuum
        self.obj_data['nebFlambda'] = self.gasExtincParams(obj_WaveRest, self.Rv_model,  self.reddedning_curve_model)

        # Calculate the nebular continua
        self.obj_data['nebularFlux'] = self.nebFluxCont(obj_WaveRest,
                                                    obj_prop_df.loc['cHbeta'][0], self.obj_data['nebFlambda'],
                                                    obj_prop_df.loc['T_low'][0],
                                                    obj_prop_df.loc['He1r_abund'][0], obj_prop_df.loc['He2r_abund'][0],
                                                    flux_halpha)


        # Save input conditions:
        Av_star = obj_prop_df.loc['Av_star'][0]
        sigma_star = obj_prop_df.loc['sigma_star'][0]
        flux_hbeta = obj_prop_df.loc['flux_hbeta'][0]
        eqw_hbeta = obj_prop_df.loc['eqw_hbeta'][0]

        #Getting all the fluxes to analyse the data
        bases_idx, bases_coeff, bases_coeff_err = np.loadtxt(self.obj_data['obj_ssp_coeffs_file'], usecols=[0, 1, 2], unpack=True)

        # SSP library flux at modelled Av, z_star and sigma_star
        ssp_grid_model_norm = self.physical_SED_model(self.ssp_lib['wave_resam'], obj_WaveObs, self.ssp_lib['flux_norm'],
                                Av_star, z_obj, sigma_star, Rv_coeff=self.Rv_model)

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
        if ('nebular' or 'stellar') in spectra_components:
            synth_spectrum_address = '{}{}_spectrum.txt'.format(output_folder, obs_name)
            np.savetxt(synth_spectrum_address, np.transpose(np.array([obj_WaveObs, self.obj_data['obsFlux']])), fmt="%7.1f %10.4e")

        # #Store synthetic object data log # TODO Add here the comments of the function
        # obj_dict = OrderedDict()
        # obj_dict['address_lines_log']   = synth_lines_log
        # obj_dict['address_spectrum']    = synth_spectrum_address
        # obj_dict['flux_hbeta']          = flux_hbeta
        # obj_dict['flux_halpha']         = flux_halpha
        # obj_dict['Normalized_by_Hbeta'] = True
        # obj_dict['eqw_hbeta']           = eqw_hbeta
        # obj_dict['sigma_gas']           = obj_prop_df.loc['sigma_gas'][0]
        # obj_dict['T_low']               = T_low
        # obj_dict['T_high']              = T_high
        # obj_dict['n_e']                 = n_e
        # obj_dict['cHbeta']              = cHbeta
        # obj_dict['Av_star']             = Av_star
        # obj_dict['sigma_star']          = sigma_star
        # obj_dict['z_obj']               = z_obj
        # obj_dict['continuum_sigma']     = error_stellarContinuum
        # obj_dict['resample_inc']        = resample_inc
        # obj_dict['wavelengh_limits']    = wavelengh_limits
        # obj_dict['norm_interval']       = norm_interval
        # obj_dict['Te_prior']            = [10000.0, 1000.0]
        # obj_dict['ne_prior']            = [80, 50]
        # obj_dict['cHbeta_prior']        = [0.125, 0.02]
        #
        # #Save the data into an "ini" configuration file
        # conf_address = '{}{}_objParams.txt'.format(output_folder, obs_name)
        # parse_synth_objData(conf_address, obs_name, obj_dict)

        return

    def ready_simulation(self, output_folder, obs_data, ssp_data, fitting_components, overwrite_grids=False,
                         input_lines=None, wavelengh_limits = None, resample_inc=None, norm_interval=None):

        # Dictionary to store the data
        self.obj_data = obs_data.copy()

        if output_folder is None: # TODO need to establish a global path to save the data
            self.output_folder = self.obj_data['output_folder']

        #Read log with observational features and masks
        obj_lines_df = read_csv(obs_data['obj_lines_file'], delim_whitespace=True, header=0, index_col=0)

        # Prepare emission data
        if 'emission' in fitting_components:

            #Prepare data from emission line file
            self.import_emission_line_data(obj_lines_df, input_lines = input_lines)

            #Create or load emissivity grids
            self.import_emissivity_grids()

            # Fit emissivity grids to a surface
            self.fit_emissivity_surface()

            # Reddening parameters
            self.obj_data['lineFlambda'] = self.gasExtincParams(self.obj_data['lineWaves'], self.Rv_model, self.reddedning_curve_model)

        # Trim, resample and normalize according to the input object
        if ('nebular' in fitting_components) or ('stellar' in fitting_components):
            self.treat_input_spectrum(self.obj_data, self.obj_data['obs_wavelength'], self.obj_data['obs_flux'], wavelengh_limits, resample_inc, norm_interval)

        # Reddening parameters for the nebular continuum
        if 'nebular' in fitting_components:
            obj_wave_rest = self.obj_data['wave_resam'] / (1.0 + self.obj_data['z_obj'])
            self.obj_data['nebFlambda'] = self.gasExtincParams(obj_wave_rest, self.Rv_model,  self.reddedning_curve_model)

        # Declare stellar data
        if 'stellar' in fitting_components:
            self.ssp_lib = ssp_data.copy()

        # Generate object masks
        if 'wave_resam' in self.obj_data: # TODO this is a bit dirty
            self.generate_object_mask(obj_lines_df, self.obj_data['wave_resam'])

        return

    def prepareSimulation(self, obs_data, ssp_data=None, output_folder = '', spectra_components=None, input_lines='all',
                        prefit_ssp = True, prefit_data = False, wavelengh_limits = None, resample_inc = None, norm_interval=None):

        # Store components fit
        self.spectra_components = spectra_components

        # Folders to store inputs and outputs
        self.input_folder = output_folder + 'input_data/'
        self.output_folder = output_folder + 'output_data/'

        # Create them if not available
        make_folder(self.input_folder)
        make_folder(self.output_folder)

        # Prepare data for the MCMC fitting
        self.ready_simulation(output_folder, obs_data, ssp_data, spectra_components, input_lines=input_lines,
                              wavelengh_limits=wavelengh_limits, resample_inc=resample_inc, norm_interval=norm_interval)

        # Prefit stellar continua
        if 'stellar' in self.spectra_components:

            self.stellarCheck = True

            # Run prefit output #TODO these prefit parameters should go into the master conf
            self.prefit_db = self.input_folder + 'SSP_prefit_database'

            # Perform a new SPP synthesis otherwise use available data
            if prefit_ssp:

                # Ready data for a prefit on our stellar continua
                self.prepare_ssp_data(None, obj_wave=self.obj_data['wave_resam'], obj_flux=self.obj_data['flux_norm'],
                                      obj_flux_err=self.obj_data['continuum_sigma'], obj_mask=self.int_mask)

                # Select model
                self.select_inference_model('stelar_prefit')

                # Run stellar continua prefit and store the results
                self.run_pymc(self.prefit_db, iterations = 5000, variables_list = ['Av_star', 'sigma_star'], prefit = True)
                self.save_prefit_data(self.obj_data['wave_resam'])

                # Plot input simulation data
                self.plot_prefit_data()

            # Load prefit data
            self.load_prefit_data(self.obj_data['wave_resam'])

            # Ready stellar
            self.prepare_ssp_data(self.input_folder + 'prefit_populations.txt',
                                  obj_wave=self.obj_data['wave_resam'], obj_flux=self.obj_data['flux_norm'],
                                  obj_flux_err=self.obj_data['continuum_sigma'], obj_mask=self.int_mask)

        else:
            self.stellarCheck = False

        # Ready gas data
        if 'emission' in self.spectra_components:
            self.emissionCheck = True
            self.prepare_gas_data()
        else:
            self.emissionCheck = False

        return

    def prepare_gas_data(self):

        self.Te_prior           = np.array(self.obj_data['Te_prior'])
        self.ne_prior           = np.array(self.obj_data['ne_prior'])
        self.cHbeta_prior       = np.array(self.obj_data['cHbeta_prior'])

        # Emission line fluxes and errors #TODO Add normalization check here
        self.obsLineFluxes      = self.obj_data['lineFluxes']
        self.obsLineFluxErr     = self.obj_data['lineErr']
        self.logObsLineFluxes   = np.log10(self.obsLineFluxes)
        self.logObsLineFluxErr  = np.abs(np.log10(self.obsLineFluxErr)) # TODO IS THIS FUCKING RIGHT?

        # Lines data
        self.lineLabels         = self.obj_data['lineLabels']
        self.lineIons           = self.obj_data['lineIons']
        self.lineFlambda        = self.obj_data['lineFlambda']

        # Variables to make the iterations simpler
        self.gasSamplerVariables(self.lineIons, self.high_temp_ions)

        return

    def prepare_ssp_data(self, input_data, obj_wave, obj_flux, obj_flux_err, obj_mask): # Remove nebularCheck, mainPopulationsFile = None

        # prepareContinuaData(basesWave, basesFlux, obsWave, obsFlux, obsFlux, objMask, nebularCheck = False, mainPopulationsFile = None)

        # if nebularCheck (remove nebular continua)

        # if mainPopulationsFile is not None (slice the bases)

        #Case of a prefit run
        if input_data is None:

            #Default parameters for the nebular continuum # TODO check this mechanic
            Te_neb, z_neb, cHbeta_neb = self.nebDefault['Te_neb'], self.nebDefault['z_neb'], self.nebDefault['cHbeta_neb']
            He1_neb, He2_neb = self.nebDefault['He1_neb'], self.nebDefault['He2_neb']

            obj_wave_rest = obj_wave / (1 + z_neb)

            #Generate synthetic nebular emission to remove from object
            self.obj_data['synth_neb_flux'] = self.nebFluxCont(obj_wave_rest, cHbeta_neb, self.obj_data['nebFlambda'],
                                                    Te_neb, He1_neb, He2_neb, self.obj_data['flux_halpha'])

            object_continuum = (obj_flux - self.obj_data['synth_neb_flux']) * obj_mask

            idx_populations = np.ones(self.ssp_lib['flux_resam'].shape[0], dtype=bool)

        # Case of a final run
        else:

            bases_idx, bases_coeff, bases_coeff_err = np.loadtxt(input_data, usecols=[0, 1, 2], unpack=True)

            idx_populations = bases_coeff >= self.lowlimit_sspContribution # TODO we should check this limit thing currentl it must be 0.001

            object_continuum = obj_flux * obj_mask

            self.sspPrefit_Idcs = np.where(idx_populations)[0]
            self.sspPrefit_Coeffs = bases_coeff[idx_populations]
            self.sspPrefit_err = bases_coeff_err[idx_populations]
            self.sspPrefit_Limits = np.vstack((self.sspPrefit_Coeffs * 0.8, self.sspPrefit_Coeffs * 1.2)).T

        #Object data
        self.input_continuum        = object_continuum
        self.input_continuum_er     = obj_flux_err
        self.input_wave             = obj_wave

        #Bases parameters
        neglected_populations       = np.where(~idx_populations)
        self.onBasesWave            = self.ssp_lib['wave_resam']
        self.onBasesFlux            = np.delete(self.ssp_lib['flux_resam'], neglected_populations, axis=0)
        self.onBasesFluxNorm        = np.delete(self.ssp_lib['flux_norm'], neglected_populations, axis=0)
        self.onBasesFluxNormCoeffs  = np.delete(self.ssp_lib['normFlux_coeff'], neglected_populations, axis=0)
        self.nBases                 = self.onBasesFlux.size #self.onBasesFlux.shape[0]
        self.range_bases            = np.arange(self.nBases)

        # Limit for bases
        self.zMin_SspLimit = np.around((obsWave[-1] / basesWave[-1]), decimals=2)
        self.zMax_SspLimit = np.around((obsWave[0] / basesWave[0]), decimals=2)

        # Static reddening curve for the stellar continuum
        self.Xx_stellar = CCM89_Bal07(self.Rv_model, self.ssp_lib['wave_resam'])

        return

    def save_prefit_data(self, obj_wave):

        # Read input data # TODO this should not be necessary, we should get it in rum_pymc2
        stellarPrefit_db, stat_db_dict      = self.load_pymc_database_manual(self.prefit_db, burning=1000,
                                             params_list=['Av_star', 'sigma_star', 'ssp_coefficients'])

        # Default parameters for the nebular continuum #TODO we need to change this by the actual prefit data the redshift here
        Te_neb, z_neb, cHbeta_neb           = self.nebDefault['Te_neb'], self.nebDefault['z_neb'], self.nebDefault['cHbeta_neb']
        He1_neb, He2_neb                    = self.nebDefault['He1_neb'], self.nebDefault['He2_neb']

        # Generate synthetic nebular emission to remove from object
        self.ssp_lib['synth_neb_flux']      = self.nebFluxCont(obj_wave, cHbeta_neb, self.obj_data['nebFlambda'],
                                              Te_neb, He1_neb, He2_neb, self.obj_data['flux_halpha'])

        # Save population coefficients      # TODO Default files names should go into the configuration
        sspCoeffs_file                      = self.input_folder + 'prefit_populations.txt'
        coeffs_mean                         = stat_db_dict['ssp_coefficients']['trace'].mean(axis=0)
        coeffs_std                          = stat_db_dict['ssp_coefficients']['trace'].std(axis=0)
        pops_vector                         = np.arange(len(coeffs_mean), dtype=int)
        np.savetxt(sspCoeffs_file, np.transpose(np.array([pops_vector,coeffs_mean, coeffs_std])), fmt="%4i %10.8f %10.8f")

        #Store data
        self.ssp_lib['db_sspPrefit']        = stellarPrefit_db
        self.ssp_lib['dbDict_sspPrefit']    = stat_db_dict

        # Mean parameter values #TODO we need to add the redshift here
        self.sspPrefit_Coeffs               = coeffs_mean
        self.ssp_lib['Av_sspPrefit']        = stat_db_dict['Av_star']['mean']
        self.ssp_lib['sigma_sspPrefit']     = stat_db_dict['sigma_star']['mean']

    def load_prefit_data(self, obj_wave):

        # TODO are an extension/code/folder for each run
        #-----Observation versus prefit comparison # TODO do we need to load this
        stellarPrefit_db, stat_db_dict      = self.load_pymc_database_manual(self.prefit_db, burning=1000, #TODO Add global burning
                                              params_list=['Av_star', 'sigma_star', 'ssp_coefficients'])

        # Default parameters for the nebular continuum
        Te_neb, z_neb, cHbeta_neb           = self.nebDefault['Te_neb'], self.nebDefault['z_neb'], self.nebDefault['cHbeta_neb']
        He1_neb, He2_neb                    = self.nebDefault['He1_neb'], self.nebDefault['He2_neb']

        # Generate synthetic nebular emission to remove from object
        self.ssp_lib['synth_neb_flux']      = self.nebFluxCont(obj_wave, cHbeta_neb, self.obj_data['nebFlambda'],
                                              Te_neb, He1_neb, He2_neb, self.obj_data['flux_halpha'])

        #Store data
        self.ssp_lib['db_sspPrefit']        = stellarPrefit_db
        self.ssp_lib['dbDict_sspPrefit']    = stat_db_dict

        # Mean parameter values #TODO we need to add the redshift here
        self.sspPrefit_Coeffs               = stat_db_dict['ssp_coefficients']['trace'].mean(axis=0)
        self.sspPrefit_err                  = stat_db_dict['ssp_coefficients']['trace'].std(axis=0)
        self.ssp_lib['Av_sspPrefit']        = stat_db_dict['Av_star']['mean']
        self.ssp_lib['sigma_sspPrefit']     = stat_db_dict['sigma_star']['mean']

        return

