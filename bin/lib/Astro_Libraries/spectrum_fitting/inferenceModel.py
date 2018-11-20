from os import path, makedirs, environ
environ["MKL_THREADING_LAYER"] = "GNU"
import pymc3
import pymc as pymc2
import cPickle as pickle
import theano
import theano.tensor as tt
from collections import OrderedDict
from numpy import mean, std, square, percentile, median, sum as np_sum, array, ones, empty, any
from lib.Astro_Libraries.spectrum_fitting.specSynthesizer_tools import ModelIngredients
from gasEmission_functions import TOIII_TSIII_relation
import matplotlib.pyplot as plt
from import_functions import ImportModelData, parseObjData, make_folder

# Line to avoid the compute_test_value error
theano.config.compute_test_value = "ignore"

class SpectraSynthesizer(ModelIngredients):

    def __init__(self):

        ModelIngredients.__init__(self)

    def fitSpectra(self, model_name, hammer = 'nuts', iterations=8000, tuning=2000, output_folder=''):

        # Declare the sampler #TODO we need to differentiate the model from the sampler
        self.select_inference_model(model = hammer)

        # Run the sampler
        db_address = output_folder + model_name + '.db' # TODO Deberiamos poder quitar este .db
        self.normContants = {'He1r': 0.1, 'He2r': 0.001}  # TODO need to decide where to place this
        self.run_pymc(db_address, iterations=iterations, tuning=tuning, model_type=hammer)

        # # Load the results
        inferenceTrace, interenceParamsDict = self.load_pymc_database_manual(db_address, sampler='pymc3')

        # # # Plot output data
        self.plotOuputData(self.output_folder + model_name, inferenceTrace, interenceParamsDict, array(interenceParamsDict.keys()))

        return

    def run_pymc(self, db_address, iterations = 10000, tuning = 0, variables_list = None, prefit = True, model_type = 'pymc2'):

        #TODO this part is very dirty it is not clear where it goes
        if 'nuts' not in model_type:

            # Define MCMC model
            MAP_Model = pymc2.MAP(self.inf_dict)

            # Prefit:
            if prefit is not False:
                fit_method = prefit if prefit is str else 'fmin_powell'
                print '\n--Starting {} prefit'.format(fit_method)
                MAP_Model.fit(method = fit_method)

            # Print prefit data
            print 'Initial conditions:'
            self.display_run_data(MAP_Model, variables_list)

            # Launch sample
            print '\nInitiating fit:'
            self.pymc2_M = pymc2.MCMC(MAP_Model.variables, db = 'pickle', dbname =  db_address)
            self.pymc2_M.sample(iter=iterations)

            # Save the output csv mean data
            if variables_list != None:
                print '--Saving results in csv'
                csv_address = db_address + '_Parameters'
                self.pymc2_M.write_csv(csv_address, variables=variables_list)

            #Print again the output prediction for the entire trace
            self.display_run_data(MAP_Model, variables_list)

            #Close the database
            self.pymc2_M.db.close()

        else:
            print '-- Starting sampling'
            # Launch sample
            trace, model = self.inf_dict(iterations, tuning)

            # Save the data
            with open(db_address, 'wb') as trace_pickle:
                pickle.dump({'model': model, 'trace': trace}, trace_pickle)

    def select_inference_model(self, model = 'nuts', params_list = None):

        # Declare the simulation type
        if model == 'nuts':
            #self.inf_dict = self.nuts_model
            self.inf_dict = self.nuts_TwoTemps

        elif model == 'stelar_prefit':
            print 'Limites', self.zMin_SspLimit, self.zMax_SspLimit
            self.inf_dict = self.stellarContinua_model()

        elif model == 'complete':
            self.inf_dict = self.complete_model()

    def nuts_model(self, iterations, tuning):

        # Container to store the synthetic line fluxes
        if self.emissionCheck:
            lineFlux_tt = tt.zeros(self.lineLabels.size)
            continuum = tt.zeros(self.obj_data['wave_resam'].size)
            # idx_N2_6548A = self.lineLabels == 'N2_6548A'
            # idx_N2_6584A = self.lineLabels == 'N2_6584A'
            # self.obsLineFluxErr[idx_N2_6548A], self.obsLineFluxErr[idx_N2_6584A] = 0.1* self.obsLineFluxes[idx_N2_6548A], 0.1 * self.obsLineFluxes[idx_N2_6584A]

        # Stellar bases tensor
        if self.stellarCheck:
            Xx_tt = theano.shared(self.Xx_stellar)
            basesFlux_tt = theano.shared(self.onBasesFluxNorm)
            nebular_continuum_tt = theano.shared(self.nebDefault['synth_neb_flux'])
            err_Continuum = 0.10 * ones(self.inputContinuum.size) # TODO really need to check this
            # err_Continuum = self.obsFluxNorm * 0.05
            # err_Continuum[err_Continuum < 0.001] = err_Continuum.mean()

        with pymc3.Model() as model:

            if self.stellarCheck:

                # Stellar continuum priors
                Av_star = pymc3.Normal('Av_star', mu=self.stellarAv_prior[0], sd=self.stellarAv_prior[0] * 0.10) #pymc3.Lognormal('Av_star', mu=1, sd=0.75)
                w_i = pymc3.Normal('w_i', mu=self.sspPrefitCoeffs, sd=self.sspPrefitCoeffs*0.10, shape=self.nBases)

                # Compute stellar continuum
                stellar_continuum = w_i.dot(basesFlux_tt)

                # Apply extinction
                spectrum_reddened = stellar_continuum * tt.pow(10, -0.4 * Av_star * Xx_tt)

                # Add nebular component
                continuum = spectrum_reddened + nebular_continuum_tt #pymc3.Deterministic('continuum_Op', spectrum_reddened + nebular_continuum)

                # Apply mask
                continuum_masked = continuum * self.int_mask

                # Likelihood continuum components
                Y_continuum = pymc3.Normal('Y_continuum', mu=continuum_masked, sd=err_Continuum, observed=self.inputContinuum)

            if self.emissionCheck:

                # Gas Physical conditions priors
                T_low = pymc3.Normal('T_low', mu=self.Te_prior[0], sd=1000.0)
                cHbeta = pymc3.Lognormal('cHbeta', mu=0, sd=1) if self.NoReddening is False else self.obj_data['cHbeta_true']

                # High temperature
                T_high = TOIII_TSIII_relation(T_low)

                if self.emissionCheck:

                    # Emission lines density
                    n_e =  pymc3.Normal('n_e', mu=self.ne_prior[0], sd=self.ne_prior[1])
                    #n_e = self.normContants['n_e'] * pymc3.Lognormal('n_e', mu=0, sd=1)

                    # Helium abundance priors
                    if self.He1rCheck:
                        tau = pymc3.Lognormal('tau', mu=1, sd=0.75)

                    # Composition priors
                    abund_dict = {'H1r':1.0}
                    for j in self.rangeObsAtoms:
                        if self.obsAtoms[j] == 'He1r':
                            abund_dict[self.obsAtoms[j]] = self.normContants['He1r'] * pymc3.Lognormal(self.obsAtoms[j], mu=0, sd=1)#pymc3.Uniform(self.obsAtoms[j], lower=0, upper=1)
                        elif self.obsAtoms[j] == 'He2r':
                            abund_dict[self.obsAtoms[j]] = self.normContants['He2r'] * pymc3.Lognormal(self.obsAtoms[j], mu=0, sd=1)#pymc3.Uniform(self.obsAtoms[j], lower=0, upper=1)
                        else:
                            abund_dict[self.obsAtoms[j]] = pymc3.Normal(self.obsAtoms[j], mu=5, sd=5)

                    # Loop through the lines
                    for i in self.rangeLines:

                        # Line data
                        line_label = self.lineLabels[i]
                        line_ion = self.lineIons[i]
                        line_flambda = self.lineFlambda[i]

                        # Parameters to compute the emissivity
                        line_coeffs = self.emisCoeffs[line_label]
                        emis_func = self.ionEmisEq_tt[line_label]

                        # Appropiate data for the ion
                        Te_calc = T_high if self.idx_highU[i] else T_low

                        # Line Emissivitiy
                        line_emis = emis_func((Te_calc, n_e), *line_coeffs)

                        # Atom abundance
                        line_abund = 1.0 if self.H1_lineIdcs[i] else abund_dict[line_ion]

                        # Line continuum
                        line_continuum = tt.sum(continuum * self.boolean_matrix[i]) * self.lineRes[i]

                        # ftau correction for HeI lines
                        line_ftau = self.ftau_func(tau, Te_calc, n_e, *self.ftau_coeffs[line_label]) if self.He1_lineIdcs[i] else None

                        # Line synthetic flux
                        flux_i = self.fluxEq_tt[line_label](line_emis, cHbeta, line_flambda, line_abund, line_ftau, continuum=line_continuum)

                        # Store in container
                        lineFlux_tt = tt.inc_subtensor(lineFlux_tt[i], flux_i)

                    # Store computed fluxes
                    lineFlux_ttarray = pymc3.Deterministic('calcFluxes_Op',lineFlux_tt)

                    # Likelihood gas components
                    Y_emision = pymc3.Normal('Y_emision', mu=lineFlux_ttarray, sd=self.obsLineFluxErr, observed=self.obsLineFluxes)

            # Get energy traces in model
            for RV in model.basic_RVs:
                print(RV.name, RV.logp(model.test_point))

            # Launch model
            trace = pymc3.sample(iterations, tune=tuning, nchains=2, njobs=2)

        return trace, model

    def nuts_TwoTemps(self, iterations, tuning):

        # Container to store the synthetic line fluxes
        if self.emissionCheck:
            lineFlux_tt = tt.zeros(self.lineLabels.size)
            continuum = tt.zeros(self.obj_data['wave_resam'].size)
            # idx_N2_6548A = self.lineLabels == 'N2_6548A'
            # idx_N2_6584A = self.lineLabels == 'N2_6584A'
            # self.obsLineFluxErr[idx_N2_6548A], self.obsLineFluxErr[idx_N2_6584A] = 0.1* self.obsLineFluxes[idx_N2_6548A], 0.1 * self.obsLineFluxes[idx_N2_6584A]

        # Stellar bases tensor
        if self.stellarCheck:
            Xx_tt = theano.shared(self.Xx_stellar)
            basesFlux_tt = theano.shared(self.onBasesFluxNorm)
            nebular_continuum_tt = theano.shared(self.nebDefault['synth_neb_flux'])
            err_Continuum = 0.10 * ones(self.inputContinuum.size) # TODO really need to check this
            # err_Continuum = self.obsFluxNorm * 0.05
            # err_Continuum[err_Continuum < 0.001] = err_Continuum.mean()

        with pymc3.Model() as model:

            if self.stellarCheck:

                # Stellar continuum priors
                Av_star = pymc3.Normal('Av_star', mu=self.stellarAv_prior[0], sd=self.stellarAv_prior[0] * 0.10) #pymc3.Lognormal('Av_star', mu=1, sd=0.75)
                w_i = pymc3.Normal('w_i', mu=self.sspPrefitCoeffs, sd=self.sspPrefitCoeffs*0.10, shape=self.nBases)

                # Compute stellar continuum
                stellar_continuum = w_i.dot(basesFlux_tt)

                # Apply extinction
                spectrum_reddened = stellar_continuum * tt.pow(10, -0.4 * Av_star * Xx_tt)

                # Add nebular component
                continuum = spectrum_reddened + nebular_continuum_tt #pymc3.Deterministic('continuum_Op', spectrum_reddened + nebular_continuum)

                # Apply mask
                continuum_masked = continuum * self.int_mask

                # Likelihood continuum components
                Y_continuum = pymc3.Normal('Y_continuum', mu=continuum_masked, sd=err_Continuum, observed=self.inputContinuum)

            if self.emissionCheck:

                # Gas Physical conditions priors
                print 'Este prior es', self.Te_prior[0]
                T_low = pymc3.Normal('T_low', mu=self.Te_prior[0], sd=2000.0)
                cHbeta = pymc3.Lognormal('cHbeta', mu=0, sd=1) if self.NoReddening is False else self.obj_data['cHbeta_true']

                # # Declare a High temperature prior if ions are available, else use the empirical relation.
                # if any(self.idx_highU):
                #     T_high = pymc3.Normal('T_high', mu=10000.0, sd=1000.0)
                # else:
                #     T_high = TOIII_TSIII_relation(self.Te_prior[0]) #TODO Should we always create a prior just to eliminate the contamination?

                if self.emissionCheck:

                    # Emission lines density
                    n_e =  255.0#pymc3.Normal('n_e', mu=self.ne_prior[0], sd=self.ne_prior[1])
                    #n_e = self.normContants['n_e'] * pymc3.Lognormal('n_e', mu=0, sd=1)

                    # Helium abundance priors
                    if self.He1rCheck:
                        tau = pymc3.Lognormal('tau', mu=1, sd=0.75)

                    # Composition priors
                    abund_dict = {'H1r':1.0}
                    for j in self.rangeObsAtoms:
                        if self.obsAtoms[j] == 'He1r':
                            abund_dict[self.obsAtoms[j]] = self.normContants['He1r'] * pymc3.Lognormal(self.obsAtoms[j], mu=0, sd=1)#pymc3.Uniform(self.obsAtoms[j], lower=0, upper=1)
                        elif self.obsAtoms[j] == 'He2r':
                            abund_dict[self.obsAtoms[j]] = self.normContants['He2r'] * pymc3.Lognormal(self.obsAtoms[j], mu=0, sd=1)#pymc3.Uniform(self.obsAtoms[j], lower=0, upper=1)
                        else:
                            abund_dict[self.obsAtoms[j]] = pymc3.Normal(self.obsAtoms[j], mu=5, sd=5)

                    # Loop through the lines
                    for i in self.rangeLines:

                        # Line data
                        line_label = self.lineLabels[i]
                        line_ion = self.lineIons[i]
                        line_flambda = self.lineFlambda[i]

                        # Parameters to compute the emissivity
                        line_coeffs = self.emisCoeffs[line_label]
                        emis_func = self.ionEmisEq_tt[line_label]

                        # Appropiate data for the ion
                        #Te_calc = T_high if self.idx_highU[i] else T_low
                        Te_calc =  T_low

                        # Line Emissivitiy
                        line_emis = emis_func((Te_calc, n_e), *line_coeffs)

                        # Atom abundance
                        line_abund = 1.0 if self.H1_lineIdcs[i] else abund_dict[line_ion]

                        # Line continuum
                        line_continuum = tt.sum(continuum * self.boolean_matrix[i]) * self.lineRes[i]

                        # ftau correction for HeI lines
                        line_ftau = self.ftau_func(tau, Te_calc, n_e, *self.ftau_coeffs[line_label]) if self.He1_lineIdcs[i] else None

                        # Line synthetic flux
                        flux_i = self.fluxEq_tt[line_label](line_emis, cHbeta, line_flambda, line_abund, line_ftau, continuum=line_continuum)

                        # Store in container
                        lineFlux_tt = tt.inc_subtensor(lineFlux_tt[i], flux_i)

                    # Store computed fluxes
                    lineFlux_ttarray = pymc3.Deterministic('calcFluxes_Op',lineFlux_tt)

                    # Likelihood gas components
                    Y_emision = pymc3.Normal('Y_emision', mu=lineFlux_ttarray, sd=self.obsLineFluxErr, observed=self.obsLineFluxes)

            # Get energy traces in model
            for RV in model.basic_RVs:
                print(RV.name, RV.logp(model.test_point))

            # Launch model
            trace = pymc3.sample(iterations, tune=tuning, nchains=2, njobs=2)

        return trace, model

    def stellarContinua_model(self):

        #Stellar parameters
        Av_star     = pymc2.Uniform('Av_star', 0.0, 5.00)
        sigma_star  = pymc2.Uniform('sigma_star', 0.0, 5.00)
        #z_star      = pymc2.Uniform('z_star', self.z_min_ssp_limit, self.z_max_ssp_limit)

        # Shift, multiply and convolve by a factor given by the model parameters
        @pymc2.deterministic
        def ssp_coefficients(z_star=0.0, Av_star=Av_star, sigma_star=sigma_star, input_flux=self.inputContinuum):

            ssp_grid_i = self.physical_SED_model(self.onBasesWave, self.inputWave, self.onBasesFluxNorm, Av_star, z_star, sigma_star, self.Rv_model)

            self.ssp_grid_i_masked = (self.int_mask * ssp_grid_i.T).T

            ssp_coeffs_norm = self.ssp_fitting(self.ssp_grid_i_masked, input_flux)

            return ssp_coeffs_norm

        # Theoretical normalized flux
        @pymc2.deterministic
        def stellar_continua_calculation(ssp_coeffs=ssp_coefficients):

            flux_sspFit_norm = np_sum(ssp_coeffs.T * self.ssp_grid_i_masked, axis=1)

            return flux_sspFit_norm

        # Likelihood
        @pymc2.stochastic(observed=True)
        def likelihood_ssp(value = self.inputContinuum, stellarTheoFlux=stellar_continua_calculation, sigmaContinuum=self.inputContinuumEr):
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

    def load_pymc_database_manual(self, db_address, burning = 0, params_list = None, sampler = 'pymc2'):

        if sampler is 'pymc3':

            # Restore the trace
            with open(db_address, 'rb') as trace_restored:
                data = pickle.load(trace_restored)
            basic_model, trace = data['model'], data['trace']

            # Save the parameters you want in a dictionary of dicts
            param_dict, stats_dict = OrderedDict(), OrderedDict()
            for parameter in trace.varnames:
                if ('_log__' not in parameter) and ('interval' not in parameter):

                    trace_norm = self.normContants[parameter] if parameter in self.normContants else 1.0
                    trace_i = trace_norm * trace[parameter]
                    stats_dict[parameter] = OrderedDict()
                    param_dict[parameter] = array([trace[parameter].mean(), trace[parameter].std()])

                    stats_dict[parameter]['mean']                    = mean(trace_i)
                    stats_dict[parameter]['median']                  = median(trace_i)
                    stats_dict[parameter]['standard deviation']      = std(trace_i)
                    stats_dict[parameter]['n']                       = trace_i.size
                    stats_dict[parameter]['16th_p']                  = percentile(trace_i, 16)
                    stats_dict[parameter]['84th_p']                  = percentile(trace_i, 84)
                    stats_dict[parameter]['95% HPD interval']        = (stats_dict[parameter]['16th_p'], stats_dict[parameter]['84th_p'])
                    stats_dict[parameter]['trace']                   = trace_i

            # Save parameters into the object log
            parseObjData(self.configFile, self.objName + '_results', param_dict)

            return trace, stats_dict

        else:

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

            #Generate a pymc2 database to recover all the data from the run
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