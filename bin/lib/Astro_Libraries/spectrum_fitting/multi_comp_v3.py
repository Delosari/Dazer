from os import path, makedirs, environ
environ["MKL_THREADING_LAYER"] = "GNU"
import theano.tensor as tt
import pymc as pymc2
import pymc3
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
from pymc3 import math as pm3_math
from pyneb import atomicData, RecAtom, Atom
from collections import OrderedDict
from numpy import arange,mean, std, square, empty, percentile, median, sum as np_sum
from specSynther_tools import SpectrumFitter
from gasEmission_functions import TOIII_TSIII_relation
from import_functions import make_folder
from lib.Math_Libraries.functions_for_theano import BilinearInterpTheano


class SpecSynthesizer(SpectrumFitter, BilinearInterpTheano):

    def __init__(self, atomic_ref=None, nebular_data = None, ssp_type=None, ssp_folder=None, ssp_conf_file=None,
                 temp_grid=None, den_grid=None, high_temp_ions=None, R_v=3.4, reddenig_curve='G03_average', lowlimit_sspContribution=0.0):

        #TODO in here we give the choice: Configuration from text file but it can also be imported from text file create a default conf to upload from there if missing

        # Import parent classes
        SpectrumFitter.__init__(self)
        BilinearInterpTheano.__init__(self)

        # Define atomic data:
        if atomic_ref == 'default':
            atomicData.setDataFile('s_iii_coll_HRS12.dat') #TODO actually this should be uplodaded from a file

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

        # Dictionary to store parameters
        #self.conf_dic = {'Te_neb': 10000.0, 'z_neb': 0.0, 'cHbeta_neb': 0.1, 'He1_neb': 0.085, 'He2_neb': 0.0}

    def fit_observation(self, model_name, model_type, obs_data, ssp_data, iterations=15000, burn=7000, thin=1, output_folder = None,
                        fitting_components=None, input_lines='all', params_list = None, prefit_SSP = True, prefit_model = False,
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
                                      obj_flux_err=self.obj_data['continuum_sigma'], obj_mask=self.obj_data['int_mask'])

                # Select model
                self.select_inference_model('stelar_prefit')

                # Run stellar continua prefit and store the results
                self.run_pymc(self.prefit_db, iterations = 5000, variables_list = ['Av_star', 'sigma_star'], prefit = True)
                self.save_prefit_data(self.obj_data['wave_resam'])

                # Plot input simulation data
                self.plot_prefit_data()

                print 'Acabamos'
                exit()

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
        # possible_Recomb = ['H1r', 'He1r', 'He2r']
        # recombination_check = set(self.params_list) & set(possible_Recomb)
        #
        # if recombination_check and ('stellar' in self.fitting_components):
        #     self.fitting_components.append('recomb_lines')
        # Check if we must analyse recombination and/or collisional excited lines
        # possible_colExcit = ['S2_abund','S3_abund','O2_abund','O3_abund', 'N2_abund', 'Ar3_abund', 'Ar4_abund']
        # collisional_check = set(self.params_list) & set(possible_colExcit)
        #
        # if collisional_check:
        #     self.fitting_components.append('colExcited_lines')

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

    def run_pymc(self, db_address, iterations = 10000, variables_list = None, prefit = True, model_type = 'pymc2'):

        if 'pymc3' not in model_type:

            # Define MCMC model
            self.MAP_Model = pymc2.MAP(self.inf_dict)

            # Prefit:
            if prefit is not False:
                fit_method = prefit if prefit is str else 'fmin_powell'
                print '\n--Starting {} prefit'.format(fit_method)
                self.MAP_Model.fit(method = fit_method)

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