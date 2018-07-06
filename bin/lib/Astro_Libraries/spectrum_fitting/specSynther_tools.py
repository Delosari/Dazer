import numpy as np
from os import path
from sys import exit
from collections import OrderedDict
from import_functions import ImportModelData
from starContinuum_functions import SspFitter
from plot_tools import MCMC_printer
from numpy import ones, exp
from pandas import read_csv
from import_functions import parse_synth_objData
from gasEmission_functions import TOIII_TSIII_relation, EmissionComponents
from gasContinuum_functions import NebularContinuaCalculator
from extinction_tools import ReddeningLaws

class SpectrumFitter(ImportModelData, SspFitter, NebularContinuaCalculator, EmissionComponents, ReddeningLaws, MCMC_printer):

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

    def gen_synth_obs(self, spectra_components = None, input_ions = None,
                      wavelengh_limits = None, resample_inc = None, norm_interval = None,
                      obj_properties_file = None, obj_lines_file = None, obj_mask_file = None,
                      output_folder = None, obs_name = 'synthObj', ssp_lib_type = None, data_folder = None,
                      data_file = None ,obj_ssp_coeffs_file = None,
                      error_stellarContinuum = None, error_lines = None):

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

            # Import the physical data from the log
            for param in obj_prop_df.index:
                self.obj_data[param] = obj_prop_df.loc[param][0]

            #Calculate T_high assuming we get T_low
            self.obj_data['T_high'] = TOIII_TSIII_relation(self.obj_data['T_low'])

            # Compute lines flux
            logLine_fluxes = self.calcEmFluxes(self.obj_data['T_low'], self.obj_data['T_high'], obj_prop_df.loc['n_e'][0],
                                                            obj_prop_df.loc['cHbeta'][0], obj_prop_df.loc['tau'][0],
                                                            self.abund_dict,
                                                            self.obj_data['lineLabels'],
                                                            self.obj_data['lineIons'],
                                                            self.obj_data['lineFlambda'])

            self.obj_data['lineFluxes'] = np.power(10, logLine_fluxes)

            # Use general error if this is provided
            self.obj_data['lineErr'] = self.obj_data['lineFluxes'] * error_lines

            # Store data to synthetic observation files
            idx_lines = (obj_lines_df.index != 'H1_4861A')
            obj_lines_df.loc[idx_lines,'obs_flux']      = self.obj_data['lineFluxes']
            obj_lines_df.loc[idx_lines,'obs_fluxErr']   = self.obj_data['lineErr']
            obj_lines_df.loc['H1_4861A','obs_flux']     = 1.0
            obj_lines_df.loc['H1_4861A','obs_fluxErr']  = 1.0 * error_lines

            # Create txt log
            synth_lines_log = '{}{}_lineslog.txt'.format(output_folder, obs_name)
            with open(synth_lines_log, 'w') as f:
                f.write(obj_lines_df.ix[:, :'blended_label'].to_string(float_format=lambda x: "{:15.8f}".format(x), index=True, index_names=False))

        # -------- Nebular calculation
        if 'nebular' in spectra_components:

            #Save input conditions:
            self.obj_data['wavelengh_limits']  = wavelengh_limits
            self.obj_data['resample_inc']      = resample_inc
            self.obj_data['norm_interval']     = norm_interval
            self.obj_data['z_obj']             = obj_prop_df.loc['z_obj'][0]

            # Rest and observed wavelength
            obj_wave_rest = np.arange(wavelengh_limits[0], wavelengh_limits[-1], resample_inc, dtype=float)
            obj_wave_obs = obj_wave_rest * (1.0 + self.obj_data['z_obj'])
            self.obj_data['obs_wave_rest'] = obj_wave_rest
            self.obj_data['obs_wave'] = obj_wave_obs

            #Get Halpha flux to calibrate
            idx_Halpha = (self.obj_data['lineLabels'] == 'H1_6563A')
            self.obj_data['flux_halpha'] = self.obj_data['lineFluxes'][idx_Halpha][0] * obj_prop_df.loc['flux_hbeta'][0]

            # Reddening parameters for the nebular continuum
            self.obj_data['nebFlambda'] = self.gasExtincParams(obj_wave_rest, self.Rv_model,  self.reddedning_curve_model)

            #Calculate the nebular continua
            self.obj_data['obs_flux'] = self.nebFluxCont(obj_wave_rest,
                                                        obj_prop_df.loc['cHbeta'][0], self.obj_data['nebFlambda'],
                                                        obj_prop_df.loc['T_low'][0],
                                                        obj_prop_df.loc['He1r_abund'][0], obj_prop_df.loc['He2r_abund'][0],
                                                        self.obj_data['flux_halpha'])

        #--------Stellar calculation
        if 'stellar' in spectra_components:

            # Save input conditions:
            self.obj_data['wavelengh_limits']   = wavelengh_limits
            self.obj_data['resample_inc']       = resample_inc
            self.obj_data['norm_interval']      = norm_interval
            self.obj_data['z_obj']              = obj_prop_df.loc['z_obj'][0]
            self.obj_data['Av_star']            = obj_prop_df.loc['Av_star'][0]
            self.obj_data['sigma_star']         = obj_prop_df.loc['sigma_star'][0]
            self.obj_data['flux_hbeta']         = obj_prop_df.loc['flux_hbeta'][0]
            self.obj_data['eqw_hbeta']          = obj_prop_df.loc['eqw_hbeta'][0]

            # Import the stellar library
            self.ssp_lib = self.load_ssp_library(ssp_lib_type, data_folder, data_file, wavelengh_limits, resample_inc, norm_interval)

            # Trim bases flux to include only the stellar spectra contributing in our object
            self.prepare_ssp_data(self.obj_data['obj_ssp_coeffs_file'], obj_wave = ones(10), obj_flux = ones(10),
                                  obj_flux_err = ones(10), obj_mask = ones(10))

            # Rest stellar and observed wavelengths
            obj_wave_rest = np.arange(wavelengh_limits[0], wavelengh_limits[-1], resample_inc, dtype=float)
            obj_wave_obs = obj_wave_rest * (1.0 + self.obj_data['z_obj'])

            # SSP library flux at modelled Av, z_star and sigma_star
            ssp_grid_model_norm = self.physical_SED_model(self.ssp_lib['wave_resam'], obj_wave_obs, self.onBasesFluxNorm,
                                    self.obj_data['Av_star'], self.obj_data['z_obj'], self.obj_data['sigma_star'], Rv_coeff=self.Rv_model)

            #Normalized object flux
            obj_flux_norm = ssp_grid_model_norm.dot(self.sspPrefit_Coeffs)

            #Instead of denormalizing the grid we use the Hbeta flux and equivalent width
            #To get a better shape of the nebular and stellar continua
            cont_hbeta  = self.obj_data['flux_hbeta'] / self.obj_data['eqw_hbeta']
            obj_flux = obj_flux_norm * cont_hbeta

            # Generate synth error
            sigma_err   = error_stellarContinuum * np.median(obj_flux) # TODO We should play with this error
            stellar_err = np.random.normal(0.0, sigma_err, len(obj_flux_norm))

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

        # Store the spectrum as a text file
        if ('nebular' or 'stellar') in spectra_components:
            synth_spectrum_address = '{}{}_spectrum.txt'.format(output_folder, obs_name)
            np.savetxt(synth_spectrum_address, np.transpose([self.obj_data['obs_wave'], self.obj_data['obs_flux']]))

        #Store synthetic object data log # TODO Add here the comments of the function
        obj_dict = OrderedDict()
        obj_dict['address_lines_log']   = synth_lines_log
        obj_dict['address_spectrum']    = synth_spectrum_address
        obj_dict['address_obs_mask']    = self.obj_data['obs_mask_address']
        obj_dict['flux_hbeta']          = self.obj_data['flux_hbeta'] if 'flux_hbeta' in self.obj_data else None
        obj_dict['flux_halpha']         = self.obj_data['flux_halpha'] if 'flux_halpha' in self.obj_data else None
        obj_dict['Normalized_by_Hbeta'] = True
        obj_dict['eqw_hbeta']           = self.obj_data['eqw_hbeta'] if 'eqw_hbeta' in self.obj_data else None
        obj_dict['T_low']               = self.obj_data['T_low'] if 'T_low' in self.obj_data else None
        obj_dict['T_high']              = self.obj_data['T_high'] if 'T_high' in self.obj_data else None
        obj_dict['n_e']                 = self.obj_data['n_e'] if 'n_e' in self.obj_data else None
        obj_dict['cHbeta']              = self.obj_data['cHbeta'] if 'cHbeta' in self.obj_data else None
        obj_dict['sigma_gas']           = self.obj_data['sigma_gas'] if 'sigma_gas' in self.obj_data else None
        obj_dict['z_obj']               = self.obj_data['z_obj'] if 'z_obj' in self.obj_data else None
        obj_dict['continuum_sigma']     = 0.02
        obj_dict['resample_inc']        = 1
        obj_dict['wavelengh_limits']    = [4000, 6900]
        obj_dict['norm_interval']       = [5100, 5150]
        obj_dict['Te_prior']            = [15000, 700]
        obj_dict['ne_prior']            = [100, 50]
        obj_dict['cHbeta_prior']        = [0.075, 0.05]

        #Save the data into an "ini" configuration file
        conf_address = '{}{}_objParams.txt'.format(output_folder, obs_name)
        parse_synth_objData(conf_address, obs_name, obj_dict)

        return

    def ready_simulation(self, output_folder, obs_data, ssp_data, fitting_components, overwrite_grids=False,
                         input_lines=None, wavelengh_limits = None, resample_inc=None, norm_interval=None):

        # Dictionary to store the data
        self.obj_data = obs_data.copy()

        self.fitting_components = fitting_components #TODO I could change this one location

        if output_folder is None: # TODO need to establish a global path to save the data
            self.output_folder = self.obj_data['output_folder']

        # Prepare emission data
        if 'emission' in fitting_components:

            obj_lines_df = read_csv(obs_data['obj_lines_file'], delim_whitespace=True, header=0, index_col=0)

            #Prepare data from emission line file
            self.import_emission_line_data(obj_lines_df, input_lines = input_lines)

            #Create or load emissivity grids
            self.import_emissivity_grids()

            # Fit emissivity grids to a surface
            self.fit_emissivity_surface()

            # Reddening parameters
            self.obj_data['lineFlambda'] = self.gasExtincParams(self.obj_data['lineWaves'], self.Rv_model, self.reddedning_curve_model)

        if ('nebular' or 'stellar') in fitting_components:

            # Trim, resample and normalize according to the input object
            self.treat_input_spectrum(self.obj_data, self.obj_data['obs_wavelength'], self.obj_data['obs_flux'], wavelengh_limits, resample_inc, norm_interval)

        if 'nebular' in fitting_components:

            # Reddening parameters for the nebular continuum
            obj_wave_rest = self.obj_data['wave_resam'] / (1.0 + self.obj_data['z_obj'])
            self.obj_data['nebFlambda'] = self.gasExtincParams(obj_wave_rest, self.Rv_model,  self.reddedning_curve_model)

        # Prepare stellar data
        if 'stellar' in fitting_components:

            # Create dictionary with the stellar library configuration which this object will use
            self.ssp_lib = ssp_data.copy()

            # Generate object masks:
            self.generate_object_mask()

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

        self.log_colExc_fluxes  = np.log10(colExc_fluxes)
        self.log_colExc_Err     = np.log10(colExc_err)**-2

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
            abund_list += list(np.unique(self.obj_data['colLine_ions']))
        if 'recombLine_ions' in self.obj_data:
            recomb_exc_list = list(np.unique(self.obj_data['recombLine_ions']))
            if 'H1' in recomb_exc_list:
                recomb_exc_list.remove('H1')
            abund_list += recomb_exc_list

        abund_list = np.core.defchararray.add(np.array(abund_list), '_abund')

        self.abund_list = list(abund_list)

        # Vector temperature regions keys
        self.tempRegionRecomb = np.array(['T_high']*self.n_recombLines)
        self.tempRegionCol    = np.array(['T_low ']*self.n_colExcLines)

        #Assign values
        self.tempRegionCol[self.idx_highU] = 'T_high'

        return

    def prepare_ssp_data(self, input_data, obj_wave, obj_flux, obj_flux_err, obj_mask):

        #Case of a prefit run
        if input_data is None:

            #Default parameters for the nebular continuum # TODO check this mechanic
            Te_neb, z_neb, cHbeta_neb = 10000.0, 0.0, 0.1
            He1_neb, He2_neb = 0.1, 0.001

            obj_wave_rest = obj_wave / (1 + z_neb)

            #Generate synthetic nebular emission to remove from object
            self.obj_data['synth_neb_flux'] = self.nebFluxCont(obj_wave_rest, cHbeta_neb, self.obj_data['nebFlambda'],
                                                    Te_neb, He1_neb, He2_neb, self.obj_data['flux_halpha'])

            object_continuum = (obj_flux - self.obj_data['synth_neb_flux']) * obj_mask

            idx_populations = ones(self.ssp_lib['flux_resam'].shape[0], dtype=bool)

        # Case of a final run
        else:

            bases_idx, bases_coeff = np.loadtxt(input_data, usecols=[0, 1], unpack=True)

            idx_populations = bases_coeff >= self.lowlimit_sspContribution # TODO we should check this limit thing currentl it must be 0.001

            object_continuum = obj_flux * obj_mask

            self.sspPrefit_Idcs = np.where(idx_populations)[0]
            self.sspPrefit_Coeffs = bases_coeff[idx_populations]
            self.sspPrefit_Limits = np.vstack((self.sspPrefit_Coeffs * 0.8, self.sspPrefit_Coeffs * 1.2)).T

        #Object data
        self.input_continuum        = object_continuum
        self.input_continuum_er     = obj_flux_err
        self.input_wave             = obj_wave
        self.int_mask               = obj_mask

        #Bases parameters
        neglected_populations       = np.where(~idx_populations)
        self.onBasesWave            = self.ssp_lib['wave_resam']
        self.onBasesFlux            = np.delete(self.ssp_lib['flux_resam'], neglected_populations, axis=0)
        self.onBasesFluxNorm        = np.delete(self.ssp_lib['flux_norm'], neglected_populations, axis=0)
        self.onBasesFluxNormCoeffs  = np.delete(self.ssp_lib['normFlux_coeff'], neglected_populations, axis=0)
        self.range_bases            = np.arange(self.onBasesFlux.shape[0])

        # Limit for bases
        z_max_ssp                   = (self.ssp_lib['wave_resam'][0] / obj_wave[0]) - 1.0
        z_min_ssp                   = (self.ssp_lib['wave_resam'][-1] / obj_wave[-1]) - 1.0
        self.z_object               = self.obj_data['z_obj']
        self.zMin_SspLimit          = z_min_ssp
        self.zMax_SspLimit          = round(z_max_ssp - 0.001, 3)

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
        pops_vector                         = np.arange(len(coeffs_vector), dtype=int)
        np.savetxt(sspCoeffs_file, np.transpose(np.array([pops_vector, coeffs_vector])), fmt="%4i %10.8f")

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

