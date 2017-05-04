from collections        import OrderedDict
from numpy              import argmax, append, exp, zeros, pi, asarray, argsort, diff, concatenate, int8, where, array, sqrt, square, ones, power, sum, mean, std, linspace, max, round, median, polyfit, vstack, random, greater, less, searchsorted, sort
from numpy              import float as npfloat, log as nplog, round as np_round, around, float32, invert
from numpy.linalg       import lstsq
from pymc               import deterministic, stochastic, Normal, Uniform, MCMC, Bernoulli, stochastic_from_dist
from scipy.integrate    import simps
from scipy.signal       import argrelextrema
from uncertainties      import ufloat
from bces_script        import bces
from linfit_script      import linfit
from lnr_script         import kelly
from lmfit              import Parameters, minimize as lmfit_minimize, fit_report
from pyneb              import Atom
from lmfit.models import GaussianModel

class Fitting_Gaussians():
    
    def __init__(self):
            
        self.Combined                           = None
        self.MontecarloCheck                    = True
        self.MC_Iterations                      = 10
                
        self.NComps                             = 0
        self.GaussianSampling                   = 101
                
        self.Fitting_dict                       = OrderedDict()
        self.Fitting_dict['Deblend check']      = False    #This logic is true when a blended group is observed
        self.Fitting_dict['Fitting method']     = None
        self.Fitting_dict['start treatment']    = False    #This logic is true to force the deblending process
        self.Fitting_dict['line group']         = None     #This integer provides the index for the group of lines describes the blended lines
        self.Fitting_dict['line number']        = None     
        self.Fitting_dict['line label']         = None     #This string describes the label for the current blended line we are mesuring: H1_6563A
        self.Fitting_dict['blended label']      = None     #This string should be used for the blended label for the lines log
        self.Fitting_dict['blended number']     = None     #This integer describes the number of blended lines. Currently it must be equal to the number in the blended list elements.
        self.Fitting_dict['blended wavelengths']= None     #This list contains the wavelengths of the blended lines
        self.Fitting_dict['blended index']      = None     #This integer provides the index of the current line in a blended group
        self.Fitting_dict['kmpfit_dict']        = None     #This dictionaries containes the parameters for kmpfit

        self.Fitting_dict['Wide component']     = False    #This keyword informs if there is a wide component on the emission line
        
        self.Fitting_dict['line type']          = None     #Spectral feature type string: 'Absorption', 'Emission'
        self.Fitting_dict['peak_waves']         = None     #This list contains the wavelenghts of all the peaks  

        self.Fitting_dict['y_scaler']           = None     #This is the magnitude used to scale the y flux (normaly the line peak or line higher peak)
        self.Fitting_dict['x_scaler']           = None     #This is the magnitude used to scale the wavelength values (normaly the line middle wavelength)
        
        self.Fitting_dict['x_resample']         = None     #x Gaussian resampling of the line required for plotting
        self.Fitting_dict['y_resample']         = None     #y Gaussian resampling of the line required for plotting
        self.Fitting_dict['y_resample_total']   = None     #y Gaussian resampling of the blended group line required for plotting
        
        self.Fitting_dict['ContinuumFlux']      = None     #Continuum intensity across the whole plotting region 
        self.Fitting_dict['m_Continuum']        = None     #Assuming a line this is the gradient (m)
        self.Fitting_dict['n_Continuum']        = None     #Assuming a line this is the y axis interception point (n) 
        self.Fitting_dict['zerolev_resample']   = None     #Continuum level resampling at the line
        self.Fitting_dict['zerolev_median']     = None     #Mean level at the line center
        self.Fitting_dict['zerolev_sigma']      = None     #Continuum dispersion assuming linear shape 
        self.Fitting_dict['ContinuumWidth']     = None     #This is the number of pixels manually selected
        
        self.Fitting_dict['x_norm']             = None
        self.Fitting_dict['y_norm']             = None
        self.Fitting_dict['zerolev_norm']       = None
        self.Fitting_dict['sig_zerolev_norm']   = None
        
        self.Fitting_dict['lmfit_params']       = None     #lmfit parameters dict
        self.Fitting_dict['lmfit_params_wide']  = None     #lmfit parameters dict
        self.Fitting_dict['MC_iteration']       = None     #Number of iterations for calculation (1 for normal 1000 for MC)
        
        self.Fitting_dict['lmfit_output']       = None
        self.Fitting_dict['FluxI_N_vector']     = None     #This vectors holds all the fluxes calculated for the case of a gaussian fit
        self.Fitting_dict['parameters_list']    = None     #This list contains all the parameters which are fitted for a line: A1, mu1, sigma1, A2, mu2, sigma2 ...
        self.Fitting_dict['FluxI']              = None
        self.Fitting_dict['Add_wideComponent']  = None    #Extra step to procced to the wide component calculation 
        self.Fitting_dict['WC_theowavelength']  = None        #Gaussian_Coefficients coefficients
        
        self.Fitting_dict['wide mask']          = None
        self.Fitting_dict['Maxima_Waves']       = None      #To store the data from the best peaks location
        self.Fitting_dict['Maxima_peaks']       = None

        self.GHcoeffs = {}
        self.GHcoeffs['c0'] = sqrt(6.0) / 4.0
        self.GHcoeffs['c1'] = -sqrt(3.0)
        self.GHcoeffs['c2'] = -sqrt(6.0)
        self.GHcoeffs['c3'] = 2.0 * sqrt(3.0) / 3.0
        self.GHcoeffs['c4'] = sqrt(6.0) / 3.0
        
        self.skeness_limit = {'fixed':(False)}
        self.kutorsis_limit = {'fixed':(False)}
        
        self.skeness_Glimit = {'fixed':(True)}
        self.kutorsis_Glimit = {'fixed':(True)}        

        N2 = Atom('N', 2)
        N2_6548A = N2.getEmissivity(tem=10000, den=100, wave=6548)
        N2_6584A = N2.getEmissivity(tem=10000, den=100, wave=6584)
        
        self.N2_Ratio = N2_6584A / N2_6548A
        
        self.sqrt2pi = sqrt(2*pi)

    def lmfit_gaussian_Residual(self, params, x, y, zerolev, Ncomps, err):
     
        return (self.gaussian_curve_SingleMixture(params.valuesdict(), x, zerolev, Ncomps) - y) / err

    def lmfit_gaussian_Residual_wide(self, params, x, y, zerolev, Ncomps, err):
     
        return (self.gaussian_curve_SingleMixture_wide(params.valuesdict(), x, zerolev, Ncomps) - y) / err
   
    def gaussian_curve_SingleMixture(self, params, x, zerolev, Ncomps):
        
        y_model = 0.0
        for i in range(Ncomps):
            index   = str(i)
            A       = params['A' + index] 
            mu      = params['mu' + index] 
            sigma   = params['sigma' + index] 
            y_model = y_model + A * exp(-(x-mu)*(x-mu)/(2*sigma*sigma))

         
        return y_model + zerolev
    
    def gaussian_curve_SingleMixture_wide(self, params, x, zerolev, Ncomps):
        
        y_model = 0.0
        for index in Ncomps:
            A       = params['A' + index] 
            mu      = params['mu' + index] 
            sigma   = params['sigma' + index] 
            y_model = y_model + A * exp(-(x-mu)*(x-mu)/(2*sigma*sigma))

         
        return y_model + zerolev

    def gaussian_curve_SingleMixture_from_dict(self, dict, x, zerolev, Ncomps):

        y_model = 0.0
        for comps in Ncomps:
            A       = dict['A' + comps].nominal_value
            mu      = dict['mu' + comps].nominal_value 
            sigma   = dict['sigma' + comps].nominal_value 
            y_model = y_model + A * exp(-(x-mu)*(x-mu)/(2*sigma*sigma))
         
        return y_model + zerolev       
        
    def Load_lmfit_parameters(self, x, y, zerolev, err_zerolev, n_comps, wide_component = False, A_limits = 0.30, mu_precission = 1, sigma_limit = 5):
        
        #Scale parameters
        ind_max = argmax(y)
        self.Fitting_dict['x_scaler'] = x[ind_max]
        self.Fitting_dict['y_scaler'] = y[ind_max]
        
        
        #Scale the range
        self.Fitting_dict['x_norm']             = x - self.Fitting_dict['x_scaler']
        self.Fitting_dict['y_norm']             = y / self.Fitting_dict['y_scaler']
        self.Fitting_dict['zerolev_norm']       = zerolev / self.Fitting_dict['y_scaler']
        self.Fitting_dict['sig_zerolev_norm']   = err_zerolev / self.Fitting_dict['y_scaler']
        
        #Get line maxima and minima
        peak_wave, peak_flux, minima_wave, minima_flux = self.get_lines_peaks(ind_max, n_comps)
        
        #Store peaks location for log
        self.Fitting_dict['peak_waves']  = peak_wave + self.Fitting_dict['x_scaler']
        self.Fitting_dict['peak_Maxima'] = peak_flux * self.Fitting_dict['y_scaler']
         
        #Lmfit dictionary        
        params = Parameters()
        for i in range(n_comps):  
            index = str(i)
            params.add('A'     + index, value = peak_flux[i] - mean(self.Fitting_dict['zerolev_norm']), min = 0.0)
            params.add('mu'    + index, value = peak_wave[i], min = peak_wave[i] - mu_precission, max = peak_wave[i] + mu_precission)
            params.add('sigma' + index, value = 1, min = 0)
            #params.add('fwhm'  + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
            params.add('FluxG' + index, expr = '{A} * {sigma} * {sqrt2pi}'.format(A = 'A'  + index, sigma = 'sigma' + index, sqrt2pi = self.sqrt2pi))
            #params.add('FluxG' + index, expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + index, fwhm = 'fwhm' + index))
            
        #For blended components we set the same wavelength:
        if n_comps > 1:
                        
            Highest_index       = argmax(self.Fitting_dict['peak_Maxima'])
            small_components    = range(n_comps)    
            del small_components[Highest_index]
            
            for indx in small_components:
                
                #We set the same sigma
                index_small = str(indx)
                expresion = 'sigma{index_big} * ( (mu{index_small} + {scaller}) / (mu{index_big} + {scaller}) )'.format(index_big = Highest_index, index_small = index_small, scaller = self.Fitting_dict['x_scaler'])
                params['sigma' + index_small].set(expr = expresion) 
            
                #We force the theoretical - biggest mu
                #expresion = '{mu_small} - mu{index_big}'.format(mu_small = self.Fitting_dict['blended wavelengths'][indx] - self.Fitting_dict['x_scaler'], index_big = Highest_index)
                #params['mu' + index_small].set(expr = expresion) 
          
        #Special condition: Wide componentine in Halpha
        Wide_params_list = []        
        if self.Fitting_dict['Add_wideComponent']:
            
            #Additional fitter
            params_W = Parameters()
            
            #TRICK TO ADD AN ADDITIONAL VALUE
            n_nindex = str(n_comps)               
            params_W.add('A'  + n_nindex,       value  = 0.2, min = 0)
            params_W.add('mu' + n_nindex,       value = 0.0)
            params_W.add('sigma' + n_nindex,    value = 6, min = 3, max = 20.0)
            #params_W.add('fwhm' + n_nindex,     expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + n_nindex))
            #params_W.add('FluxG'+ n_nindex,     expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + n_nindex, fwhm = 'fwhm' + n_nindex))
            params_W.add('FluxG' + n_nindex, expr = '{A} * {sigma} * {sqrt2pi}'.format(A = 'A'  + n_nindex, sigma = 'sigma' + n_nindex, sqrt2pi = self.sqrt2pi))
            Wide_params_list = params_W.keys()
            
            #Update for Nitrogen relation: Mode 1 adjuxt the fluxes
            params['FluxG0'].set(expr = 'FluxG2 / {N2_ratio}'.format(N2_ratio = 2.94))
            #Update for Nitrogen relation: Mode 2 adjust the amplitudes
#             params['A0'].set(expr = '(A2*sigma2) / ({N2_ratio}*sigma0) '.format(N2_ratio = 2.94)) 
            
            self.Fitting_dict['lmfit_params_wide']  = params_W

        #Store the data
        self.Fitting_dict['lmfit_params']       = params
        self.Fitting_dict['parameters_list']    = array(params.keys() + Wide_params_list) 
 
        return

    def get_lines_peaks(self, ind_max, Ncomps):
        
        target_wavelengths = None
                    
        #Get feature geometry:
        #-- Single line
        if self.Fitting_dict['Deblend check'] == False:    
            peak_flux = array([self.Fitting_dict['y_norm'][ind_max]]) 
            peak_wave = array([self.Fitting_dict['x_norm'][ind_max]])    
            minima_wave = 0.0
            minima_flux = 0.0
            
        #--Blended line
        else:
            max_index, min_index    = argrelextrema(self.Fitting_dict['y_norm'], greater)[0], argrelextrema(self.Fitting_dict['y_norm'], less)[0]
            List_blended_lines      = self.Fitting_dict['Blended list'][1][self.Fitting_dict['line group']]
            
            #With wide component #ONLY WORKS FOR THE BLENDED HALPHA SCHEME
            if self.Fitting_dict['Add_wideComponent'] == False:
                target_wavelengths = array(List_blended_lines) - self.Fitting_dict['x_scaler']
            else:
                target_wavelengths = array(List_blended_lines + [List_blended_lines[1]]) -  self.Fitting_dict['x_scaler']
            
            maxima_wavelengths = sort(self.Fitting_dict['x_norm'][max_index])
            minima_wavelengths = sort(self.Fitting_dict['x_norm'][min_index])
                        
            if len(max_index) == Ncomps:
                peak_flux, minima_flux  = self.Fitting_dict['y_norm'][max_index], self.Fitting_dict['y_norm'][min_index]
                peak_wave, minima_wave  = maxima_wavelengths, minima_wavelengths  
            else:
                closest_indeces = self.search_targets_in_array(maxima_wavelengths, target_wavelengths)
                peak_wave, peak_flux  = self.Fitting_dict['x_norm'][max_index][closest_indeces], self.Fitting_dict['y_norm'][max_index][closest_indeces]
                    
            if self.Fitting_dict['Add_wideComponent']:
                if len(peak_wave) ==  len(target_wavelengths) - 1:
                    peak_wave = append(peak_wave, [0])
                    peak_flux = append(peak_flux, [0.1])
            
            minima_wave, minima_flux =  self.Fitting_dict['x_norm'][min_index], self.Fitting_dict['y_norm'][min_index]
        
        self.Fitting_dict['Maxima_Waves'] = peak_wave + self.Fitting_dict['x_scaler']
        self.Fitting_dict['Maxima_peaks'] = peak_flux * self.Fitting_dict['y_scaler']
            
        return peak_wave, peak_flux, minima_wave, minima_flux

    def search_targets_in_array(self, known_array, test_array):
        
        #This function gives the indeces of the closest values within a sorted array
        
        index_sorted        = argsort(known_array)
        known_array_sorted  = known_array[index_sorted]
        known_array_middles = known_array_sorted[1:] - diff(known_array_sorted.astype('f'))/2
        idx1                = searchsorted(known_array_middles, test_array)
        indices             = index_sorted[idx1]
        
        return indices

    def fit_line(self, x, y, zero_lev, err_continuum, Ncomps, fitting_parameters, fitting_parameters_wide, iterations): 
        
        #Number of x points in the spectral line
        x_Grid_Length = len(x)
                
        #Generate empty containers to store the data
        self.Fitting_dict['FluxI_N_vector'] = zeros(iterations)
        for key in self.Fitting_dict['parameters_list']:
            self.Fitting_dict[key + '_norm']    = zeros(iterations)
            self.Fitting_dict[key + '_norm_er'] = zeros(iterations)
                   
        #Loop through the iterations (Only 1 if it is not a bootstrap)
        for i in range(iterations):
            
            #Fit narrow component           
            if i == 0:
                y_new = y
            else:
                noise_array = random.normal(0.0, err_continuum, x_Grid_Length).astype(float32)
                y_new       = y + noise_array
                
            fit_Output      = lmfit_minimize(self.lmfit_gaussian_Residual, fitting_parameters, args=(x, y_new, zero_lev, Ncomps, err_continuum))
            output_params   = fit_Output.params
              
            #Case with a wide component
            if self.Fitting_dict['Add_wideComponent']:
                
                #Case for halpha
                sigma_limit = fit_Output.params['sigma1'].value
                limit_0     = 6548.05 - self.Fitting_dict['x_scaler']  - sigma_limit * 1.5
                limit_1     = 6548.05 - self.Fitting_dict['x_scaler']  + sigma_limit * 1.5
                limit_2     = 0  - sigma_limit * 4
                limit_3     = 0  + sigma_limit * 4 
                limit_4     = 6583.46 - self.Fitting_dict['x_scaler']  - sigma_limit * 3
                limit_5     = 6583.46 - self.Fitting_dict['x_scaler']  + sigma_limit * 3
                
                indeces = ((x >= limit_0) & (x <= limit_1)) + ((x >= limit_2) & (x <= limit_3)) + ((x >= limit_4) & (x <= limit_5))    
                mask    = invert(indeces)
                self.Fitting_dict['wide mask'] = mask
                
                x_wide      = x[mask]
                y_wide      = y_new[mask]
                zero_wide   = zero_lev[mask]
                Ncomps_wide = ['3']
                                
                fit_Output_wide     = lmfit_minimize(self.lmfit_gaussian_Residual_wide, fitting_parameters_wide, args=(x_wide, y_wide, zero_wide, Ncomps_wide, err_continuum))

                y_wide              = self.gaussian_curve_SingleMixture_wide(fit_Output_wide.params.valuesdict(), x, zero_lev, Ncomps_wide)

                y_pure_emission     = y_new - y_wide + zero_lev 
                self.Fitting_dict['emis_limpio'] = y_pure_emission 
                
                fit_Output_emission = lmfit_minimize(self.lmfit_gaussian_Residual, fitting_parameters, args=(x, y_pure_emission, zero_lev, Ncomps, err_continuum))

                output_params       = fit_Output_emission.params + fit_Output_wide.params
                                
            #Store the integrated flux
            self.Fitting_dict['FluxI_N_vector'][i]  = simps(y_new, x) - simps(zero_lev, x)
          
            #Store the fitting parameters
            for key in self.Fitting_dict['parameters_list']:
                self.Fitting_dict[key + '_norm'][i]     = output_params[key].value
                self.Fitting_dict[key + '_norm_er'][i]  = output_params[key].stderr
                
        #Store the output fit (Only used for single line output)        
        self.Fitting_dict['lmfit_output'] = fit_Output
        
        #Finally increase the number of components in case of a wide component
        if self.Fitting_dict['Add_wideComponent']:
            self.Fitting_dict['blended number'] = self.Fitting_dict['blended number'] + 1
        
        return

#     def fit_line_together(self, x, y, zero_lev, err_continuum, Ncomps, fitting_parameters, fitting_parameters_wide, iterations): 
#         
#         #Number of x points in the spectral line
#         x_Grid_Length = len(x)
#                 
#         #Generate empty containers to store the data
#         self.Fitting_dict['FluxI_N_vector'] = zeros(iterations)
#         for key in self.Fitting_dict['parameters_list']:
#             self.Fitting_dict[key + '_norm']    = zeros(iterations)
#             self.Fitting_dict[key + '_norm_er'] = zeros(iterations)
#                    
#         #Loop through the iterations (Only 1 if it is not a bootstrap)
#         for i in range(iterations):
#             
#             #Fit narrow component           
#             if i == 0:
#                 y_new = y
#             else:
#                 noise_array = random.normal(0.0, err_continuum, x_Grid_Length).astype(float32)
#                 y_new       = y + noise_array
#                 
#             fit_Output      = lmfit_minimize(self.lmfit_gaussian_Residual, fitting_parameters, args=(x, y_new, zero_lev, Ncomps, err_continuum))
#             output_params   = fit_Output.params
#                               
#             #Store the integrated flux
#             self.Fitting_dict['FluxI_N_vector'][i]  = simps(y_new, x) - simps(zero_lev, x)
#           
#             #Store the fitting parameters
#             for key in self.Fitting_dict['parameters_list']:
#                 self.Fitting_dict[key + '_norm'][i]     = output_params[key].value
#                 self.Fitting_dict[key + '_norm_er'][i]  = output_params[key].stderr
#                 
#         #Store the output fit (Only used for single line output)        
#         self.Fitting_dict['lmfit_output'] = fit_Output
#         
#         return

    def scale_lmfit_params(self, line_wave, line_flux, fit_output, Ncomps, x_scale, y_scale, fitting_method):
        
        #Scale the integrated flux (The same for all schemes)
        flux_I_N = mean(self.Fitting_dict['FluxI_N_vector'])
        
        #Scale the gaussian parameters (divided since they have different error calculation)
        #--Simple fitting
        if 'MC' not in fitting_method:
            
            #For simple fitting the integrated flux has no error, but the gaussian components do
            self.Fitting_dict['FluxI']              = ufloat(flux_I_N, 0.0) * y_scale
                        
            for i in range(Ncomps):
                index = str(i)
                self.Fitting_dict['A' + index]      = ufloat(self.Fitting_dict['A' + index + '_norm'],      self.Fitting_dict['A' + index + '_norm_er'])      * y_scale
                self.Fitting_dict['mu' + index]     = ufloat(self.Fitting_dict['mu' + index + '_norm'],     self.Fitting_dict['mu' + index + '_norm_er'])     + x_scale
                self.Fitting_dict['sigma' + index]  = ufloat(self.Fitting_dict['sigma' + index + '_norm'],  self.Fitting_dict['sigma' + index + '_norm_er'])
                self.Fitting_dict['FluxG' + index]  = ufloat(self.Fitting_dict['FluxG' + index + '_norm'],   self.Fitting_dict['FluxG' + index + '_norm_er'])   * y_scale
               
        #--Bootstrap
        else:
            
            #For MC fitting the integrated flux has error
            self.Fitting_dict['FluxI']              = ufloat(flux_I_N, std(self.Fitting_dict['FluxI_N_vector'])) * y_scale

            for i in range(Ncomps):
                index = str(i)
                self.Fitting_dict['A' + index]      = ufloat(mean(self.Fitting_dict['A' + index + '_norm']),        std(self.Fitting_dict['A' + index + '_norm'])) * y_scale
                self.Fitting_dict['mu' + index]     = ufloat(mean(self.Fitting_dict['mu' + index + '_norm']),       std(self.Fitting_dict['mu' + index + '_norm'])) + x_scale
                self.Fitting_dict['sigma' + index]  = ufloat(mean(self.Fitting_dict['sigma' + index + '_norm']),    std(self.Fitting_dict['sigma' + index + '_norm']))
                self.Fitting_dict['FluxG' + index]  = ufloat(mean(self.Fitting_dict['FluxG' + index + '_norm']),     std(self.Fitting_dict['FluxG' + index + '_norm'])) * y_scale
                self.Fitting_dict['FluxI']          = ufloat(flux_I_N, std(self.Fitting_dict['FluxI_N_vector'])) * y_scale
        
        #Calculate the gaussian curve for plotting
        self.Fitting_dict['x_resample']             = linspace(line_wave[0], line_wave[-1], 50 * Ncomps)
        self.Fitting_dict['zerolev_resample']       = self.Fitting_dict['m_Continuum'] * self.Fitting_dict['x_resample'] + self.Fitting_dict['n_Continuum']
        self.Fitting_dict['y_resample']             = [self.gaussian_curve_SingleMixture_from_dict(self.Fitting_dict, self.Fitting_dict['x_resample'], self.Fitting_dict['zerolev_resample'], Ncomps = map(str, range(Ncomps)))]

        #We store all the gaussians for plotting alongside the big component
        for j in range(Ncomps):
            Comp = [str(j)]
            self.Fitting_dict['y_resample'].append(self.gaussian_curve_SingleMixture_from_dict(self.Fitting_dict, self.Fitting_dict['x_resample'], self.Fitting_dict['zerolev_resample'], Ncomps = Comp))

        return

class Bayesian_regressions():    
    
    def __init__(self):
        
        self.Methodology = None
        self.prob_threshold = 0.40
        
    def lr_ChiSq(self, x_array, y_array, m_0, n_0):
        
        m                       = Normal('m', m_0, 0.01)
        n                       = Normal('n', n_0, 0.01)
        sigma                   = Uniform('sigma', 0.0, 5.0)
         
        @stochastic(observed=True)
        def model(value = self.y_error, x_values = self.x_array, m = m, n = n, sigma = sigma):
            
            value_theo      = m*x_values + n
            chi_sq          = sum(square(value - value_theo) / square(sigma))
            log_ChiSq       = - chi_sq / 2.0
            return log_ChiSq
                
        return locals()
    
    def inference_outliers(self, x_array, y_array, m_0, n_0, spread_vector): 
                
        outlier_points  = Uniform('outlier_points', 0, 1.0, value=0.1)
        mean_outliers   = Uniform('mean_outliers', -100, 100, value=0)
        spread_outliers = Uniform('spread_outliers', -100, 100, value=0)
        
        @stochastic
        def slope_and_intercept(slope = m_0):
            prob_slope = nplog(1. / (1. + slope ** 2))
            return prob_slope
        
        @deterministic
        def model_(x=x_array, slope_and_intercept=slope_and_intercept):
            slope, intercept = slope_and_intercept
            fit = slope * x + intercept 
            return fit
        
        inlier = Bernoulli('inlier', p=1 - outlier_points, value=zeros(x_array.size))
        
        def log_posterior_likelihood_of_outlier(y_with_outlier, mu, spread_vector, inlier, mean_outliers, spread_outliers):
            inlier_posterior = sum(inlier * (nplog(2 * pi * spread_vector ** 2) + (y_with_outlier - mu) ** 2 / (spread_vector ** 2)))
            outlier_posterior = sum((1 - inlier) * (nplog(2 * pi * ((spread_vector ** 2) + (spread_outliers ** 2))) + (y_with_outlier - mean_outliers) ** 2 / ((spread_vector ** 2) + (spread_outliers ** 2))))
            return -0.5 * (inlier_posterior + outlier_posterior)
        
        outlier_distribution = stochastic_from_dist('outlier_distribution', logp=log_posterior_likelihood_of_outlier, dtype=npfloat, mv=True)
        
        outlier_dist = outlier_distribution('outlier_dist',  mu=model_,  spread_vector=spread_vector, mean_outliers=mean_outliers,  spread_outliers=spread_outliers, inlier=inlier, value=y_array, observed=True)
    
        return locals()
    
class Linear_Regressions(Bayesian_regressions):
    
    def __init__(self):
        
        Bayesian_regressions.__init__(self)
        
        self.x_array = None
        self.x_error = None
        
        self.y_array = None
        self.y_error = None

    def load_obs_data(self, x_values, y_values, x_errors = None, y_errors = None):
        
        #Default case we input all the values manually
        self.x_array = x_values
        self.x_error = x_errors
        
        self.y_array = y_values
        self.y_error = y_errors
         
    def perform_regression(self, Methodology):
        
        if Methodology      == 'bces':
            fit_dict        = self.bces_regression()
            
        elif Methodology    == 'Max_Likelihood':
            fit_dict        = self.max_likelihood_regression()
            
        elif Methodology    == 'linfit':
            fit_dict        = self.linfit_regression()
        
        elif Methodology    == 'scipy':
            fit_dict        = self.scipy_regression()
            
        elif Methodology    == 'kmpfit':
            fit_dict        = self.kmpfit_regression()            

        elif Methodology    == 'kelly':
            fit_dict        = self.kellyBces_regression()          

        elif 'Inference' in Methodology:
            fit_dict        = self.inference_model(Methodology)            
    
        return fit_dict

    def inference_model(self, Methodology):

        if Methodology == 'Inference - ChiSq':
            
            Inf_dict = self.inference_ChiSq()
            
        if Methodology == 'Inference - Outliers':

            Inf_dict = self.Outliers_Krough()
            
        return Inf_dict
        
    def linfit_regression(self):
        
        fit_dict    = OrderedDict()
        
        fit_dict['methodology'] = 'Linfit'
        
        Regression_Fit, Uncertainty_Matrix, fit_dict['red_ChiSq'], fit_dict['residuals'] = linfit(x_true = self.x_array, y = self.y_array, sigmay = self.y_error, relsigma = False, cov = True, chisq = True, residuals = True)
        
        m_n_Matrix                              = [sqrt(Uncertainty_Matrix[t,t]) for t in range(2)] 
        fit_dict['R_factor']                    = Uncertainty_Matrix[0,1]/(m_n_Matrix[0]*m_n_Matrix[1])
        fit_dict['m'], fit_dict['m_error']      = Regression_Fit[0], m_n_Matrix[0]
        fit_dict['n'], fit_dict['n_error']      = Regression_Fit[1], m_n_Matrix[1]
        
        return fit_dict
    
    def bces_regression(self, cov = None):
        #Rodrigo Nemmen, http://goo.gl/8S1Oo
        
        fit_dict    = OrderedDict()

        if cov == None:
            #This is the covariance between the measurements. If none provided it is assume error is independent between measurments for the 
            cov = zeros(len(self.x_array))
            
        fit_dict['methodology'] = (r'OLS(Y|X)$_{bces}$', r'OLS(X|Y)$_{bces}$', r'bisector$_{bces}$', r'Orthogonal$_{bces}$')
        
        fit_dict['m'],fit_dict['n'],fit_dict['m_error'],fit_dict['n_error'],fit_dict['cov'] = bces(self.x_array, self.x_error, self.y_array, self.y_error, cov)
        
        return fit_dict
        
    def scipy_regression(self):
        #ODR Method
        
        fit_dict    = OrderedDict()
        
        fit_dict['methodology'] = r'ODR$_{ScIPy}$'
                
        beta_0 = (0, 1)

        fit_dict['m'], fit_dict['n'], fit_dict['m_error'], fit_dict['n_error'], fit_dict['cov'], fit_dict['chiSq'],  fit_dict['red_ChiSq'] = scipy_ODR(self.x_array, self.y_array, self.y_array, self.y_error, beta_0)

        return fit_dict
    
    def kmpfit_regression(self):
        
        #Kmpfit methodology using an effective variance method
        fit_dict    = OrderedDict()
        
        fit_dict['methodology'] = r'Effective Variance$_{kmpfit}$'

        scipy_guess_dict = self.scipy_regression()
        
        beta_0 = (scipy_guess_dict['n'], scipy_guess_dict['m'])
                
        fit_dict['m'], fit_dict['n'], fit_dict['m_error'], fit_dict['n_error'], fit_dict['cov'], fit_dict['chiSq'],  fit_dict['red_ChiSq'] = kmpfit_effectivevariance(self.x_array, self.y_array, self.x_error, self.y_error, beta_0)
        
        return fit_dict

    def bayesian_regression(self, Methodology):
        
        fit_dict                = OrderedDict()
        
        fit_dict['methodology'] = r'Inference $\chi^{2}$ model'
        
        #Initial guess for the fitting:
        Np_lsf                  = polyfit(self.x_array, self.y_array, 1)
        m_0, n_0                = Np_lsf[0], Np_lsf[1]
                
        MCMC_dict               = self.lr_ChiSq(self.x_array, self.y_array, m_0, n_0)
        
        myMCMC                  = MCMC(MCMC_dict)
        
        myMCMC.sample(iter=10000, burn=1000)

        fit_dict['m'], fit_dict['n'], fit_dict['m_error'], fit_dict['n_error'] = myMCMC.stats()['m']['mean'], myMCMC.stats()['n']['mean'], myMCMC.stats()['m']['standard deviation'], myMCMC.stats()['n']['standard deviation']
        
        return fit_dict

    def kellyBces_regression(self):
        
        fit_dict                    = OrderedDict()
            
        fit_dict['methodology']   = (r'Inferences$_{bces}$')
        
        n_tuple, m_tuple, cov       = kelly(x1=self.x_array, x2=self.y_array, x1err=self.x_error, x2err=self.y_error)
        
        fit_dict['m'],fit_dict['n'],fit_dict['m_error'],fit_dict['n_error'],fit_dict['cov'] = m_tuple[0], n_tuple[0], m_tuple[1], n_tuple[1], cov
        
        return fit_dict 
        
    def Outliers_Krough(self):
        
        fit_dict                = OrderedDict()
        
        fit_dict['methodology'] = r'Outliers Krough'
        
        #Initial Guess for fitting
        Bces_guess              = self.bces_regression()
        m_0, n_0                = Bces_guess['m'][0], Bces_guess['n'][0]
                
        Spread_vector           = ones(len(self.x_array))
        
        #Model for outliers detection
        Outliers_dect_dict      = self.inference_outliers(self.x_array, self.y_array, m_0, n_0, Spread_vector)
        
        mcmc = MCMC(Outliers_dect_dict)
        mcmc.sample(100000, 20000)
        
        #Extract the data with the outliers coordinates
        probability_of_points           = mcmc.trace('inlier')[:].astype(float).mean(0)
        fit_dict['x_coords_outliers']   = self.x_array[probability_of_points < self.prob_threshold]
        fit_dict['y_coords_outliers']   = self.y_array[probability_of_points < self.prob_threshold]
                
        return fit_dict

def Python_linfit(x_true, y, y_err, errors_output = True):
     
    Regression_Fit, Uncertainty_Matrix, Red_Chi_Sq, Residuals   = linfit(x_true, y, y_err, cov=True, relsigma=False, chisq=True, residuals=True)
    m_n_Matrix                                                  = [sqrt(Uncertainty_Matrix[t,t]) for t in range(2)] 
    R_Factor                                                    = Uncertainty_Matrix[0,1]/(m_n_Matrix[0]*m_n_Matrix[1])
    m, m_error                                                  = Regression_Fit[0], m_n_Matrix[0]
    n, n_error                                                  = Regression_Fit[1], m_n_Matrix[1]
    
    if errors_output: 
        return m, m_error, n, n_error
    
    else:
        return m, n    

def NumpyRegression(x, y):
     
    Matrix_Coefficient = vstack([x, ones(len(x))]).T
    m, n = lstsq(Matrix_Coefficient, y)[0]
 
    return m, n

# from collections        import OrderedDict
# 
# from kapteyn            import kmpfit
# from numpy              import exp, zeros, asarray, argsort, diff, concatenate, int8, where, array, sqrt, square, ones, power, sum, mean, std, linspace, max, round, median, polyfit
# from numpy              import float as npfloat
# from numpy              import log as nplog
# from numpy              import pi, sum, zeros
# from numpy.linalg       import lstsq
# from pymc               import deterministic, stochastic, Normal, Uniform, MCMC, Bernoulli, stochastic_from_dist
# from scipy.integrate    import simps
# from scipy.interpolate  import interp1d
# from scipy.optimize     import minimize, curve_fit, leastsq
# from uncertainties      import unumpy
# 
# from bces_script        import bces
# from kmpfit_script      import kmpfit_effectivevariance, scipy_ODR
# from linfit_script      import linfit
# from lnr_script         import kelly
# import random           as random
# 
# 
# class Fitting_Gaussians():
#     '''This class contains the methods to measure gaussins in a scheme similar to IRAF'''
#     
#     def __init__(self):
#     
#         self.x_Norm                         = None
#         self.y_Norm                         = None
#         
#         self.A_Norm                         = None
#         self.mu_Norm                        = None
#         self.sigma_Norm                     = None
#         
#         self.zerolev_Norm                   = None
#         self.sigma_zerolev_Norm             = None
#         
#         self.Combined                       = None
#         self.MontecarloCheck                = True
#         self.Iterations                     = 1000
# 
#         self.A_Vector                       = None 
#         self.mu_Vector                      = None
#         self.sigma_Vector                   = None
#         self.y_New                          = None
#         self.Area_Norm                      = None
#         
#         self.Kmpfit_Dictionary              = []
#                 
#         self.Fitting_dict                   = OrderedDict()
#         self.Fitting_dict['Deblend check']  = False                         #This logic is true when a blended group is observed
#         self.Fitting_dict['start deblend']  = False                         #This logic is true when a blended group is observed
#         self.Fitting_dict['line group']     = None                          #This integer provides the index for the group of lines describes the blended lines
#         self.Fitting_dict['line number']    = None                          #This integer describes which line from the current blended group we are measuring. It starts at 0
#         self.Fitting_dict['line label']     = None                          #This string describes the label for the current blended line we are mesuring: H1_6563A
#         self.Fitting_dict['blended label']  = None                          #This string should be used for the blended label for the lines log
#         self.Fitting_dict['blended number'] = None                          #This integer describes the number of blended lines. Currently it must be equal to the number in the blended list elements.
#         self.Fitting_dict['p1_norm']        = {}                            #List with the normalized estimations for the gaussian fitting: [A, mu, sigma]
#         self.Fitting_dict['kmpfit_dict']    = {}                            #This dictionaries containes the parameters for kmpfit
#         self.Fitting_dict['p0']             = None
#         self.Fitting_dict['p1']             = None                          
#         self.Fitting_dict['p_1_std']        = None
#         self.Fitting_dict['p_1_Area']       = None
#         self.Fitting_dict['p_1_Area_Std']   = None
#         
#         self.NComps                         = 0
#         
#         self.GaussianSampling               = 101
#         
#         self.Max_y                  = None
#         self.Mean_x                 = None
# 
#     def Calculate_EW(self, SubWave, SubInt, Flux_Brute, No_StellarContribution, Current_Wavelength, TableAddress, TableHeaderSize):
# 
#         #Case in which we are measuring the line intensity from a zero level continuum (that means without the stellar continuum)               
#         if No_StellarContribution == True: 
#             LocalMedian = float(self.GetDataInTable("TheoWavelength", "Continuum_Median",Current_Wavelength, TableAddress,TableHeaderSize))
#         else:
#             LocalMedian = self.LocalMedian
# 
#         #Single line
#         if self.Fitting_dict['Deblend check'] == False:
#             self.EqW = Flux_Brute / LocalMedian
#             self.Classic_EmLine_Error(SubWave, SubInt, self.Flux_Brute, self.EqW)
#             self.SigmaEW_MCMC = 0.0       
# 
#         #Blended line
#         else:
#             self.Fitting_dict['Eqw']    = zeros(self.Fitting_dict['blended number'])
#             for i in range(len(self.Fitting_dict['blended number'])):
#                 self.Fitting_dict['Eqw'][i] = self.Fitting_dict['Flux_Gauss'][i] / LocalMedian
#             
#         return 
#  
# 
#        
#     def SingleGaussian_Cont(self, Ind_Variables, A, mu, sigma):
#     
#         #In order to increase somehow the speed we simplify the code by no asigning many variables
#     
#         #x_true           = Ind_Variables[0]
#         #continuum   = Ind_Variables[1]
#     
#         #y = A * exp(-(x_true-mu)*(x_true-mu)/(2.0*sigma*sigma)) + continuum
#     
#         return A * exp(-(Ind_Variables[0]-mu)*(Ind_Variables[0]-mu)/(2.0*sigma*sigma)) + Ind_Variables[1]
# 
#     def MixtureGaussian_Cont_Short(self, Ind_Variables, p, NComps):
#         y = 0.0
#         for i in range(NComps):
#             A, mu, sigma = p[i*3:(i+1)*3]
#             y += A * exp(-(Ind_Variables[0]-mu)*(Ind_Variables[0]-mu)/(2.0*sigma*sigma))
#         
#         return y + Ind_Variables[1]
#     
#     def MixtureGaussian_Cont(self, Ind_Variables, p_0):
#         
#         y       = 0.0
#         x_true       = Ind_Variables[0]
#         zerolev = Ind_Variables[1]
#                
#         N_Comps = int(len(p_0) / 3)
#         
#         for i in range(N_Comps):
#             A, mu, sigma = p_0[i*3:(i+1)*3]
#             y += A * exp(-(x_true-mu)*(x_true-mu)/(2.0*sigma*sigma))
#                 
#         return y + zerolev
#     
#     def Residuals_Gaussian(self, p, data):
#         x_true, y, zerolev, sigma_zerolev = data[0], data[1], data[2], data[3]
#         return (y - (self.SingleGaussian_Cont((x_true, zerolev), p[0], p[1], p[2]))) / sigma_zerolev
# 
#     def Residuals_Gaussian_MCMC(self, p, x_true, y, zerolev):
#         return y - (self.SingleGaussian_Cont((x_true, zerolev), p[0], p[1], p[2]))
# 
#     def Residuals_GaussianMixture(self, p, data):
#         x_true, y, zerolev, sigma_zerolev = data[0], data[1], data[2], data[3]
#         return (y - (self.MixtureGaussian_Cont_v2((x_true, zerolev), p))) / sigma_zerolev
#     
#     def MixtureGaussian_Cont_v2(self, Ind_Variables, p):
#         y = 0.0
#         for i in range( int(len(p) / 3)):
#             A, mu, sigma = p[i*3:(i+1)*3]
#             y += A * exp(-(Ind_Variables[0]-mu)*(Ind_Variables[0]-mu)/(2.0*sigma*sigma))
#         return y + Ind_Variables[1]
#     
#     def Chi2(self, p_0, x_true, y, cont_Vector, sig_Cont):
# 
#         #In order to increase somehow the speed we simplify the code by no asigning many variables
# #       Chi2 = sum(power((y - self.SingleGaussian_Cont((x_true, cont_Vector), p_0[0], p_0[1], p_0[2]))/sig_Cont, 2))
#         
#         return sum(power((y - self.SingleGaussian_Cont((x_true, cont_Vector), p_0[0], p_0[1], p_0[2]))/sig_Cont, 2))
# 
#     def Scale_Parameters(self, x_true, y, p_0, zerolev, sigma_zerolev, SingleLine = True):
#         
#         x_new = zeros(len(x_true))
#         x_list = x_true.tolist()
#         
#         #This rounding is necesary but how does this affect the precission?
#         for i in range(len(x_list)):
#             x_new[i] = round(x_list[i], 3) 
#         
#         if SingleLine:
#                 
#             p_0[1] = round(p_0[1], 3)
#             
#             self.x_Norm              = x_new             - p_0[1]
#             self.y_Norm              = y                 / p_0[0]
#             self.A_Norm              = 1.0
#             self.mu_Norm             = 0.0
#             self.sigma_Norm          = round(p_0[2], 3)
#             self.zerolev_Norm        = zerolev           / p_0[0]
#             self.sigma_zerolev_Norm  = sigma_zerolev     / p_0[0]
#                     
#         else:
#             
#             self.Max_y   = max(p_0[0])
#             self.Mean_x  = round(median(p_0[1]), 3)
#             
#             self.x_Norm             = round(x_new               - self.Mean_x, 3)
#             self.y_Norm             = round(y                   / self.Max_y,  4)
#             self.A_Norm             = round(p_0[0]              / self.Max_y,  4)
#             self.mu_Norm            = round(p_0[1]              - self.Mean_x, 3)
#             self.sigma_Norm         = round(p_0[2]                          , 3)
#             self.zerolev_Norm       = round(zerolev             / self.Max_y, 4)
#             self.sigma_zerolev_Norm = round(sigma_zerolev       / self.Max_y, 4)
#             
#         return
#     
#     def Calculate_G_Parameters(self, x_true, y, p_0, zerolev, sigma_zerolev, Method, line_type = True):
#                 
#         #Check for the line type 
#         if line_type == 'Absorption': 
#             y = 2 * zerolev - y
#         
#         #Normalize the parameters to improve the fitting quality
#         self.Scale_Parameters(x_true, y, p_0, zerolev, sigma_zerolev)
#                 
#         if Method == 'Min_Chi2': 
#             Minimize_Output                             = minimize(self.Chi2, x0 = [self.A_Norm, self.mu_Norm, self.sigma_Norm], args=(self.x_Norm, self.y_Norm, self.zerolev_Norm, self.sigma_zerolev_Norm), method = 'Nelder-Mead')
#             self.Fitting_dict['p1_norm']                = Minimize_Output['x_true']
#             self.Fitting_dict['p_1_std_Norm']           = None
#             self.Fitting_dict['p_1_Area_Norm']          = None
#             
#         elif Method == 'CurveFit':
#             self.Fitting_dict['p1_norm'], conv_curFit   = curve_fit(self.SingleGaussian_Continuum, (self.x_Norm, self.zerolev_Norm), self.y_Norm, [self.A_Norm, self.mu_Norm, self.sigma_Norm])
#             self.Fitting_dict['p_1_std_Norm']           = None
#             self.Fitting_dict['p_1_Area_Norm']          = None
# 
#         elif Method == 'leastsqr':
#             self.Fitting_dict['p1_norm'], conv          = leastsq(self.Residuals_Gaussian, [self.A_Norm, self.mu_Norm, self.sigma_Norm], args=(self.x_Norm, self.y_Norm, self.zerolev_Norm, [self.sigma_zerolev_Norm] * len(self.x_Norm)))
#             self.Fitting_dict['p_1_std_Norm']           = None
#             self.Fitting_dict['p_1_Area_Norm']          = None
# 
#         elif Method == 'kmpfit':
#             fitobj                                      = kmpfit.Fitter(residuals=self.Residuals_Gaussian, data=([self.x_Norm, self.y_Norm, self.zerolev_Norm, self.sigma_zerolev_Norm]))
#             fitobj.parinfo                              = self.Kmpfit_Dictionary
#             fitobj.fit(params0                          = [self.A_Norm, self.mu_Norm, self.sigma_Norm])
#             self.Fitting_dict['p1_norm']                = fitobj.params 
#             self.Fitting_dict['p_1_std_Norm']           = None
#             self.Fitting_dict['p_1_Area_Norm']          = None
#             
#         elif 'MCMC' in Method:      
#             self.MCMC_Scheme()
#         
#         #De normalize the fitting parameters to the physical units 
#         self.ReScale_Parameters(self.Fitting_dict['p1_norm'], p_0, self.Fitting_dict['p_1_std_Norm'], self.Fitting_dict['p_1_Area_Norm']) # ADAPT THE RESCALING FOR AN ABSORPTION LINE
#   
#         return
#     
#     def Calculate_GM_Parameters_v3(self, p_0, subWave, zerolev, Method):
#                                 
#         #Reshaping the initial guesses vector to: [A_0, A_1, A_2..., mu_0, mu_1, mu_2..., sigma_0, sigma_1, sigma_2...]
#         p_O_Reshape = zeros(self.Fitting_dict['blended number'] * 3)
#         for i in range(self.Fitting_dict['blended number']):
#             p_O_Reshape[i*3:(i+1)*3] = self.A_Norm[i], self.mu_Norm[i],  self.sigma_Norm[i]
#         
#         #Perfom fit
#         if Method == 'leastsqr':
#             self.Fitting_dict['p1_norm'], conv  = leastsq(self.Residuals_GaussianMixture, args=([self.x_Norm, self.y_Norm, self.zerolev_Norm, self.sigma_zerolev_Norm]),   x0 = p_O_Reshape)
#             self.Fitting_dict['p_1_std_Norm']   = None
#             self.Fitting_dict['p_1_Area_Norm']  = None       
#        
#             
#         if Method == 'kmpfit':
#             fitobj                                  = kmpfit.Fitter(residuals=self.Residuals_GaussianMixture, data=([self.x_Norm, self.y_Norm, self.zerolev_Norm, self.sigma_zerolev_Norm]))
#             fitobj.parinfo                          = self.Kmpfit_Dictionary
#             fitobj.fit(params0                      = p_O_Reshape)
#             self.Fitting_dict['p1_norm']            = fitobj.params 
#             self.Fitting_dict['p_1_std_Norm']       = None
#             self.Fitting_dict['p_1_Area_Norm']      = None
#             self.Fitting_dict['p_1_Area_Std_Norm']  = None 
#                      
#         elif 'MCMC' in Method:
#             self.MCMC_Scheme_Blend(Method)
#             
#         #De normalize the fitting parameters to the physical units 
#         self.ReScale_Parameters(self.Fitting_dict['p1_norm'], p_0, self.Fitting_dict['p_1_std_Norm'], self.Fitting_dict['p_1_Area_Norm'], SingleLine=False)
#         
#         return
# 
#     def ReScale_Parameters(self, p_1_Norm, p_0, p_1_stdev_Norm = None, p_1_Area_Norm = None, p_1_Area_std_Norm=None, SingleLine = True):
#         
#         #Single Line
#         if SingleLine:
#             A       = p_1_Norm[0] * p_0[0]
#             mu      = p_1_Norm[1] + p_0[1]
#             sigma   = p_1_Norm[2] 
#                           
#             p_1     = [A, mu[0], sigma]
#             
#             if (p_1_stdev_Norm == None) and (p_1_Area_Norm == None):
#                 self.Fitting_dict['p1'] = p_1
#                 return 
#     
#             else:
#                 p_1_Area                        = p_1_Area_Norm * p_0[0]
#                 p_1_stdev                       = [p_1_stdev_Norm[0] *  p_0[0], p_1_stdev_Norm[1], p_1_stdev_Norm[2]]
#                 
#                 self.Fitting_dict['p1']         = p_1
#                 self.Fitting_dict['p_1_std']    = p_1_stdev
#                 self.Fitting_dict['p_1_Area']   = p_1_Area            
#                 return
#         
#         #Blended line    
#         else:
#             p_1 = ones(len(p_1_Norm))
#             
#             for i in range(self.Fitting_dict['blended number']):
#                 A, mu, sigma                = p_1_Norm[i*3:(i+1)*3]
#                 p_1[i*3:(i+1)*3]            = A * self.Max_y, mu + self.Mean_x, sigma
# 
#             if (p_1_stdev_Norm == None) and (p_1_Area_Norm == None):
#                 self.Fitting_dict['p1']     = p_1
#                 return
# 
#             else:
#                 p_1_stdev = zeros(len(p_1_stdev_Norm))
#                 
#                 for i in range(self.Fitting_dict['blended number']):
#                     A_std, mu_std, sigma_std = p_1_stdev_Norm[i*3:(i+1)*3]
#                     p_1_stdev[i*3:(i+1)*3]   = A_std * self.Max_y, mu_std, sigma_std
#                 
#                 p_1_Area                            = p_1_Area_Norm * self.Max_y
#                 p_1_Area_Std                        = p_1_Area_std_Norm * self.Max_y
#                 
#                 self.Fitting_dict['p1']             = p_1
#                 self.Fitting_dict['p_1_std']        = p_1_stdev
#                 self.Fitting_dict['p_1_Area']       = p_1_Area
#                 self.Fitting_dict['p_1_Area_Std']   = p_1_Area
# 
#                 return
# 
#     def MCMC_Scheme(self):
#                 
#         #Number of x points in the spectral line
#         x_Grid_Length   = len(self.x_Norm)
# 
#         #Define vector to store MCMC predictions        
#         A_Vector        = zeros(self.MC_Iterations)
#         mu_Vector       = zeros(self.MC_Iterations)
#         sigma_Vector    = zeros(self.MC_Iterations)
#         AreaB_Norm      = zeros(self.MC_Iterations)
#         AreaG_Norm      = zeros(self.MC_Iterations)
#         
#         #Resample x grid for the gaussian flux
#         self.x_Norm_Resample    = linspace(self.x_Norm[0], self.x_Norm[-1], self.GaussianSampling, endpoint=True)
#         
#         #Resampling the linear continuum for the fitting #CONFLITH WITH MARTAS APPROACH?
#         Interpolation           = interp1d(self.x_Norm, self.zerolev_Norm, kind = 'slinear')
#         self.zerolev_Resample   = Interpolation(self.x_Norm_Resample)
# 
#         #Choose MC Gaussian scheme
#         if self.Fitting_dict['Fitting method']  == 'MCMC_kmpfit':
#         
#             #Loop through the bootstrap method
#             for i in range(self.MC_Iterations):
#                 
#                 #Array of simulated line
#                 y_New_Norm               = self.sigma_zerolev_Norm * random.randn(x_Grid_Length) + self.y_Norm
#                 
#                 #Run line Gaussian fitting     
#                 fitobj              = kmpfit.Fitter(residuals=self.Residuals_GaussianMixture, data=([self.x_Norm, y_New_Norm, self.zerolev_Norm, self.sigma_zerolev_Norm]))
#                 fitobj.parinfo      = self.Fitting_dict['kmpfit_dict']
#                 fitobj.fit(params0  = [self.A_Norm, self.mu_Norm, self.sigma_Norm])
#                 p_1_i               = fitobj.params 
#                 
#                 #Save parameters fitting
#                 A_Vector[i]         = p_1_i[0]
#                 mu_Vector[i]        = p_1_i[1]
#                 sigma_Vector[i]     = p_1_i[2]               
#                 
#                 #Calculate Gaussian flux
#                 y_resample          = self.SingleGaussian_Cont((self.x_Norm_Resample, self.zerolev_Resample), p_1_i[0], p_1_i[1], p_1_i[2])
#                 AreaG_Norm[i]       = self.Integrate_Gaussian(self.x_Norm_Resample, y_resample, self.zerolev_Resample)
#                 
#                 #Calculate Brute flux
#                 AreaB_Norm[i]       = self.Integrate_Gaussian(self.x_Norm, y_New_Norm, self.zerolev_Norm)
# 
#         #Case the fitting method is not recognized
#         else:
#             print 'WARNING: MCMC could not determine gaussian of emission line'
#             print '--Method', self.Methodology
#             return
#     
#         #Store the mean values from the bootstrap into the dictionary
#         self.Fitting_dict['p1_norm']            = [mean(A_Vector),     mean(mu_Vector),       mean(sigma_Vector)]      
#         self.Fitting_dict['p1_std_Norm']        = [std(A_Vector),      std(mu_Vector),        std(sigma_Vector)]
#         self.Fitting_dict['p1_Area_Norm']       = mean(AreaG_Norm)
#         self.Fitting_dict['p1_Area_Std_Norm']   = std(AreaG_Norm)
#         self.Fitting_dict['p1_AreaB_Norm']      = mean(AreaB_Norm)
#         self.Fitting_dict['p1_AreaB_Std_Norm']  = std(AreaB_Norm)
# 
#         return
#  
#     def MCMC_Scheme_Blend(self, Method):
#         
#         p_1_Norm                = []
#         p_1_std_Norm            = []
#         p_1_Area_Norm           = []
#         p_1_Area_std_Norm       = []
#                 
#         self.A_Vector           = zeros((self.Fitting_dict['blended number'], self.Iterations))
#         self.mu_Vector          = zeros((self.Fitting_dict['blended number'], self.Iterations))
#         self.sigma_Vector       = zeros((self.Fitting_dict['blended number'], self.Iterations))
#         
#         self.Area_Norm          = zeros((self.Fitting_dict['blended number'], self.Iterations))       
#         Length_y                = len(self.x_Norm)
#         self.y_New              = zeros(Length_y)
#         
#         #We put this here to improve speed
#         self.x_Norm_Resample    = linspace(self.x_Norm[0], self.x_Norm[-1], self.GaussianSampling, endpoint=True)
#         Interpolation           = interp1d(self.x_Norm, self.x_Norm, kind = 'slinear')        
#         self.zerolev_Resample   = Interpolation(self.x_Norm_Resample) 
#     
#         p_O_Reshape = zeros(self.Fitting_dict['blended number'] * 3)
#         for i in range(self.Fitting_dict['blended number']):
#             p_O_Reshape[i*3:(i+1)*3]    = self.A_Norm[i], self.mu_Norm[i],  self.sigma_Norm[i]
#     
#         if Method == 'MCMC_leastsqr':
#             for i in range(self.Iterations):
#                 for z in range(Length_y):
#                     self.y_New[z] = random.gauss(self.y_Norm[z], self.sigma_zerolev_Norm)
# 
#                 p_1_i, conv_mcmc_leastsqr   = leastsq(self.Residuals_GaussianMixture, args=([self.x_Norm, self.y_New, self.zerolev_Norm, self.sigma_zerolev_Norm]),  x0 = p_O_Reshape)
#  
#                 y_resample_c = zeros(Length_y)
#                 for j in range(self.Fitting_dict['blended number']):
#                     self.A_Vector[j][i], self.mu_Vector[j][i], self.sigma_Vector[j][i] = p_1_i[j*3:(j+1)*3]
#                     y_resample_c            = self.A_Vector[j][i] * exp(-(self.x_Norm_Resample-self.mu_Vector[j][i])*(self.x_Norm_Resample-self.mu_Vector[j][i])/(2.0*self.sigma_Vector[j][i]*self.sigma_Vector[j][i])) + self.zerolev_Resample
#                     self.Area_Norm[j][i]    = self.Integrate_Gaussian(self.x_Norm_Resample, y_resample_c, self.zerolev_Resample)
# 
#         if Method == 'MCMC_kmpfit':
#             for i in range(self.Iterations):
#                 for z in range(Length_y):
#                     self.y_New[z] = random.gauss(self.y_Norm[z], self.sigma_zerolev_Norm)
#                 
#                 fitobj          = kmpfit.Fitter(residuals=self.Residuals_GaussianMixture, data=([self.x_Norm, self.y_New, self.zerolev_Norm, self.sigma_zerolev_Norm]))
#                 fitobj.parinfo  = self.Kmpfit_Dictionary
#                 fitobj.fit(params0 = p_O_Reshape)
# 
#                 p_1_i           = fitobj.params 
#                 
#                 y_resample_c = zeros(Length_y)
#                 for j in range(self.Fitting_dict['blended number']):
#                     self.A_Vector[j][i], self.mu_Vector[j][i], self.sigma_Vector[j][i] = p_1_i[j*3:(j+1)*3]
#                     y_resample_c            = self.A_Vector[j][i] * exp(-(self.x_Norm_Resample-self.mu_Vector[j][i])*(self.x_Norm_Resample-self.mu_Vector[j][i])/(2.0*self.sigma_Vector[j][i]*self.sigma_Vector[j][i])) + self.zerolev_Resample
#                     self.Area_Norm[j][i]    = self.Integrate_Gaussian(self.x_Norm_Resample, y_resample_c, self.zerolev_Resample)
#   
#                                            
#         else:
#             print 'WARNING: MCMC could not determine gaussian of emission line'
#             print '--Method', self.Methodology
#             self.Fitting_dict['p1_norm'] = None
#             return
#         
#         for j in range(self.Fitting_dict['blended number']):
#             
#             p_1_Norm.append(mean(self.A_Vector[j]))
#             p_1_Norm.append(mean(self.mu_Vector[j]))
#             p_1_Norm.append(mean(self.sigma_Vector[j]))
# 
#             p_1_std_Norm.append(std(self.A_Vector[j]))
#             p_1_std_Norm.append(std(self.mu_Vector[j]))
#             p_1_std_Norm.append(std(self.sigma_Vector[j]))
#             
#             p_1_Area_Norm.append(mean(self.Area_Norm[j]))
#             p_1_Area_std_Norm.append(std(self.Area_Norm[j]))
#              
#         self.Fitting_dict['p1_norm']            = array(p_1_Norm)
#         self.Fitting_dict['p_1_std_Norm']       = array(p_1_std_Norm)
#         self.Fitting_dict['p_1_Area_Norm']      = array(p_1_Area_Norm)
#         self.Fitting_dict['p_1_Area_Std_Norm']  = array(p_1_Area_std_Norm)
#              
#         return
#         
#     def Resample_Gaussian(self, x_true, zerolev, sampling = None, Emission_Gaussian = True, SingleLine = True):
#         
#         #Case of a single line     
#         if SingleLine == True:
#             if sampling == None:
#                 sampling = self.GaussianSampling
#             
#             x_resample                          = linspace(x_true[0], x_true[-1], self.GaussianSampling, endpoint=True)
#             
#             #THIS interpolation is kind of stupid
#             Interpolation                       = interp1d(x_true, zerolev, kind = 'slinear')        
#             zerolev_resample                    = Interpolation(x_resample)
#             
#             if Emission_Gaussian == 'Emission':
#                 y_resample                      = self.SingleGaussian_Cont((x_resample, zerolev_resample), self.Fitting_dict['p1'][0], self.Fitting_dict['p1'][1], self.Fitting_dict['p1'][2])
#                 self.Fitting_dict['y_plotting'] = self.SingleGaussian_Cont((x_resample, zerolev_resample), self.Fitting_dict['p1'][0], self.Fitting_dict['p1'][1], self.Fitting_dict['p1'][2])
# 
#             else:
#                 y_resample                      = 2 * zerolev_resample - self.SingleGaussian_Cont((x_resample, zerolev_resample), self.Fitting_dict['p1'][0], self.Fitting_dict['p1'][1], self.Fitting_dict['p1'][2])
#                 self.Fitting_dict['y_plotting'] = self.SingleGaussian_Cont((x_resample, zerolev_resample), self.Fitting_dict['p1'][0], self.Fitting_dict['p1'][1], self.Fitting_dict['p1'][2])
#                 
#             self.Fitting_dict['x_resample']     = x_resample
#             self.Fitting_dict['y_resample']     = y_resample
#         
#         #Case of a blended line
#         else: 
#             #UPDATE THIS FOR ABSORPTIONS
#             x_resample_complete             = linspace(x_true[0], x_true[-1], self.GaussianSampling * self.NComps, endpoint=True)       
#             self.zerolev_Resample           = self.Continuum_Gradient * x_resample_complete + self.Continuum_n
#             y_resample_Complete             = self.MixtureGaussian_Cont_Short((x_resample, self.zerolev_Resample), self.p_1, self.NComps)
#         
#             self.Fitting_dict['x_resample'] = [x_resample_complete]
#             self.Fitting_dict['y_resample'] = [y_resample_Complete]
#             self.Fitting_dict['y_plotting'] = [y_resample_Complete]
#         
#             #Individual regions
#             x_resample              = linspace(x_true[0], x_true[-1], self.GaussianSampling, endpoint=True)
#             for i in range(len(self.Fitting_dict['blended number'])):
#                 A_i, mu_i, sigma_i  = self.Fitting_dict['p1'][i*3:(i+1)*3]
#                 y_resample_i        = self.SingleGaussian_Cont((x_resample, zerolev_resample), A_i, mu_i, sigma_i)
#                 
#                 self.Fitting_dict['x_resample'].append(x_resample_complete)
#                 self.Fitting_dict['y_resample'].append(y_resample_i)
#                 self.Fitting_dict['y_plotting'].append(y_resample_i)
#                       
#         return
#     
#     def Resample_Gaussian_MCMC(self, x_true, A, mu, sigma, zerolev):
#                             
#         return self.SingleGaussian_Continuum((x_true, zerolev), A, mu, sigma)   
#      
#     def Integrate_Gaussian(self, x_true, y, zerolev):
#         
#         return simps(y, x_true) - simps(zerolev,x_true)
#     
#     def Integrate_EmissionLine(self, x_true, y, zerolev):
#     
#         return simps(y, x_true) - simps(zerolev,x_true)
#     
#     def FindMaxima(self, xval, yval, MinLevel, ListLines, Deblend_Check):
#         
#         xval                = asarray(xval)
#         yval                = asarray(yval)
#         GroupLines          = asarray(ListLines)
#         
#         sort_idx            = argsort(xval)
#         yval                = yval[sort_idx]
#         gradient            = diff(yval)
#         maxima              = diff((gradient > 0).view(int8))
#         ListIndeces         = concatenate((([0],) if gradient[0] < 0 else ()) + (where(maxima == -1)[0] + 1,) + (([len(yval)-1],) if gradient[-1] > 0 else ()))
#         self.X_Maxima, self.Y_Maxima = [], []
#             
#         for index in ListIndeces:
#             if yval[index] > MinLevel:
#                 self.X_Maxima.append(xval[index])
#                 self.Y_Maxima.append(yval[index])
#                 
#         A_List              = []
#         mu_List             = []
#     
#         for i in range(len(GroupLines)):
#             TheoWave        = GroupLines[i]
#             Closest_Index   = abs(self.X_Maxima-TheoWave).argmin()
#             A_List.append(self.Y_Maxima[Closest_Index])
#             mu_List.append(self.X_Maxima[Closest_Index])
#             
#         if len(ListIndeces) == 1:
#             Deblend_Check   = True #Not sure what this is supposed to mean
#             A_List          = list(set(A_List))
#             mu_List         = list(set(mu_List))
#             self.Fitting_dict['blended number'] = 1
#  
#         A_max               = max(A_List)
#         Mean_x              = round(median(mu_List))
#                 
#         for i in range(self.Fitting_dict['blended number']):
#             if A_List[i] < 0.40 * A_max:
#                 self.Kmpfit_Dictionary.append({})
#                 self.Kmpfit_Dictionary.append({'limits':(mu_List[i]-Mean_x-0.25,mu_List[i]-Mean_x+0.25)})
#                 self.Kmpfit_Dictionary.append({'limits':(0,1.2)})
#             else:
#                 self.Kmpfit_Dictionary.append({})
#                 self.Kmpfit_Dictionary.append({})
#                 self.Kmpfit_Dictionary.append({'limits':(0,10)})
#  
#         return array(mu_List), array(A_List), Deblend_Check
#     
# class Bayesian_regressions():    
#     
#     def __init__(self):
#         
#         self.Methodology = None
#         self.prob_threshold = 0.40
#         
#     def lr_ChiSq(self, x_array, y_array, m_0, n_0):
#         
#         m                       = Normal('m', m_0, 0.01)
#         n                       = Normal('n', n_0, 0.01)
#         sigma                   = Uniform('sigma', 0.0, 5.0)
#          
#         @stochastic(observed=True)
#         def model(value = self.y_error, x_values = self.x_array, m = m, n = n, sigma = sigma):
#             
#             value_theo      = m*x_values + n
#             chi_sq          = sum(square(value - value_theo) / square(sigma))
#             log_ChiSq       = - chi_sq / 2.0
#             return log_ChiSq
#                 
#         return locals()
#     
#     def inference_outliers(self, x_array, y_array, m_0, n_0, spread_vector): 
#         
#         print m_0
#         
#         outlier_points  = Uniform('outlier_points', 0, 1.0, value=0.1)
#         mean_outliers   = Uniform('mean_outliers', -100, 100, value=0)
#         spread_outliers = Uniform('spread_outliers', -100, 100, value=0)
#         
#         @stochastic
#         def slope_and_intercept(slope = m_0):
#             prob_slope = nplog(1. / (1. + slope ** 2))
#             print 'y este', prob_slope
#             return prob_slope
#         
#         @deterministic
#         def model_(x=x_array, slope_and_intercept=slope_and_intercept):
#             slope, intercept = slope_and_intercept
#             fit = slope * x + intercept 
#             return fit
#         
#         inlier = Bernoulli('inlier', p=1 - outlier_points, value=zeros(x_array.size))
#         
#         def log_posterior_likelihood_of_outlier(y_with_outlier, mu, spread_vector, inlier, mean_outliers, spread_outliers):
#             inlier_posterior = sum(inlier * (nplog(2 * pi * spread_vector ** 2) + (y_with_outlier - mu) ** 2 / (spread_vector ** 2)))
#             outlier_posterior = sum((1 - inlier) * (nplog(2 * pi * ((spread_vector ** 2) + (spread_outliers ** 2))) + (y_with_outlier - mean_outliers) ** 2 / ((spread_vector ** 2) + (spread_outliers ** 2))))
#             return -0.5 * (inlier_posterior + outlier_posterior)
#         
#         outlier_distribution = stochastic_from_dist('outlier_distribution', logp=log_posterior_likelihood_of_outlier, dtype=npfloat, mv=True)
#         
#         outlier_dist = outlier_distribution('outlier_dist',  mu=model_,  spread_vector=spread_vector, mean_outliers=mean_outliers,  spread_outliers=spread_outliers, inlier=inlier, value=y_array, observed=True)
#     
#         return locals()
#     
# class Linear_Regressions(Bayesian_regressions):
#     
#     def __init__(self):
#         
#         Bayesian_regressions.__init__(self)
#         
#         self.x_array = None
#         self.x_error = None
#         
#         self.y_array = None
#         self.y_error = None
# 
#     def load_obs_data(self, x_values, y_values, x_errors = None, y_errors = None):
#         
#         #Default case we input all the values manually
#         self.x_array = x_values
#         self.x_error = x_errors
#         
#         self.y_array = y_values
#         self.y_error = y_errors
#          
#     def perform_regression(self, Methodology):
#         
#         if Methodology      == 'bces':
#             fit_dict        = self.bces_regression()
#             
#         elif Methodology    == 'Max_Likelihood':
#             fit_dict        = self.max_likelihood_regression()
#             
#         elif Methodology    == 'linfit':
#             fit_dict        = self.linfit_regression()
#         
#         elif Methodology    == 'scipy':
#             fit_dict        = self.scipy_regression()
#             
#         elif Methodology    == 'kmpfit':
#             fit_dict        = self.kmpfit_regression()            
# 
#         elif Methodology    == 'kelly':
#             fit_dict        = self.kellyBces_regression()          
# 
#         elif 'Inference' in Methodology:
#             fit_dict        = self.inference_model(Methodology)            
#     
#         return fit_dict
# 
#     def inference_model(self, Methodology):
# 
#         if Methodology == 'Inference - ChiSq':
#             
#             Inf_dict = self.inference_ChiSq()
#             
#         if Methodology == 'Inference - Outliers':
# 
#             Inf_dict = self.Outliers_Krough()
#             
#         return Inf_dict
#         
#     def linfit_regression(self):
#         
#         fit_dict    = OrderedDict()
#         
#         fit_dict['methodology'] = 'Linfit'
#         
#         Regression_Fit, Uncertainty_Matrix, fit_dict['red_ChiSq'], fit_dict['residuals'] = linfit(x_true = self.x_array, y = self.y_array, sigmay = self.y_error, relsigma = False, cov = True, chisq = True, residuals = True)
#         
#         m_n_Matrix                              = [sqrt(Uncertainty_Matrix[t,t]) for t in range(2)] 
#         fit_dict['R_factor']                    = Uncertainty_Matrix[0,1]/(m_n_Matrix[0]*m_n_Matrix[1])
#         fit_dict['m'], fit_dict['m_error']      = Regression_Fit[0], m_n_Matrix[0]
#         fit_dict['n'], fit_dict['n_error']      = Regression_Fit[1], m_n_Matrix[1]
#         
#         return fit_dict
#     
#     def bces_regression(self, cov = None):
#         #Rodrigo Nemmen, http://goo.gl/8S1Oo
#         
#         fit_dict    = OrderedDict()
# 
#         if cov == None:
#             #This is the covariance between the measurements. If none provided it is assume error is independent between measurments for the 
#             cov = zeros(len(self.x_array))
#             
#         fit_dict['methodology'] = (r'OLS(Y|X)$_{bces}$', r'OLS(X|Y)$_{bces}$', r'bisector$_{bces}$', r'Orthogonal$_{bces}$')
#         
#         fit_dict['m'],fit_dict['n'],fit_dict['m_error'],fit_dict['n_error'],fit_dict['cov'] = bces(self.x_array, self.x_error, self.y_array, self.y_error, cov)
#         
#         return fit_dict
#         
#     def scipy_regression(self):
#         #ODR Method
#         
#         fit_dict    = OrderedDict()
#         
#         fit_dict['methodology'] = r'ODR$_{ScIPy}$'
#                 
#         beta_0 = (0, 1)
# 
#         fit_dict['m'], fit_dict['n'], fit_dict['m_error'], fit_dict['n_error'], fit_dict['cov'], fit_dict['chiSq'],  fit_dict['red_ChiSq'] = scipy_ODR(self.x_array, self.y_array, self.y_array, self.y_error, beta_0)
# 
#         return fit_dict
#     
#     def kmpfit_regression(self):
#         
#         #Kmpfit methodology using an effective variance method
#         fit_dict    = OrderedDict()
#         
#         fit_dict['methodology'] = r'Effective Variance$_{kmpfit}$'
# 
#         scipy_guess_dict = self.scipy_regression()
#         
#         beta_0 = (scipy_guess_dict['n'], scipy_guess_dict['m'])
#         
#         print 'input values', beta_0
#         
#         fit_dict['m'], fit_dict['n'], fit_dict['m_error'], fit_dict['n_error'], fit_dict['cov'], fit_dict['chiSq'],  fit_dict['red_ChiSq'] = kmpfit_effectivevariance(self.x_array, self.y_array, self.x_error, self.y_error, beta_0)
#         
#         return fit_dict
# 
#     def bayesian_regression(self, Methodology):
#         
#         fit_dict                = OrderedDict()
#         
#         fit_dict['methodology'] = r'Inference $\chi^{2}$ model'
#         
#         #Initial guess for the fitting:
#         Np_lsf                  = polyfit(self.x_array, self.y_array, 1)
#         m_0, n_0                = Np_lsf[0], Np_lsf[1]
#         
#         print 'Bayesian guesses:', m_0, n_0
#         
#         MCMC_dict               = self.lr_ChiSq(self.x_array, self.y_array, m_0, n_0)
#         
#         myMCMC                  = MCMC(MCMC_dict)
#         
#         myMCMC.sample(iter=10000, burn=1000)
# 
#         fit_dict['m'], fit_dict['n'], fit_dict['m_error'], fit_dict['n_error'] = myMCMC.stats()['m']['mean'], myMCMC.stats()['n']['mean'], myMCMC.stats()['m']['standard deviation'], myMCMC.stats()['n']['standard deviation']
#         
#         return fit_dict
# 
#     def kellyBces_regression(self):
#         
#         fit_dict                    = OrderedDict()
#             
#         fit_dict['methodology']   = (r'Inferences$_{bces}$')
#         
#         n_tuple, m_tuple, cov       = kelly(x1=self.x_array, x2=self.y_array, x1err=self.x_error, x2err=self.y_error)
#         
#         fit_dict['m'],fit_dict['n'],fit_dict['m_error'],fit_dict['n_error'],fit_dict['cov'] = m_tuple[0], n_tuple[0], m_tuple[1], n_tuple[1], cov
#         
#         return fit_dict 
#         
#     def Outliers_Krough(self):
#         
#         fit_dict                = OrderedDict()
#         
#         fit_dict['methodology'] = r'Outliers Krough'
#         
#         #Initial Guess for fitting
#         Bces_guess              = self.bces_regression()
#         m_0, n_0                = Bces_guess['m'][0], Bces_guess['n'][0]
#         
#         print 'estos son...', m_0, n_0
#         
#         Spread_vector           = ones(len(self.x_array))
#         
#         #Model for outliers detection
#         Outliers_dect_dict      = self.inference_outliers(self.x_array, self.y_array, m_0, n_0, Spread_vector)
#         
#         mcmc = MCMC(Outliers_dect_dict)
#         mcmc.sample(100000, 20000)
#         
#         #Extract the data with the outliers coordinates
#         probability_of_points           = mcmc.trace('inlier')[:].astype(float).mean(0)
#         fit_dict['x_coords_outliers']   = self.x_array[probability_of_points < self.prob_threshold]
#         fit_dict['y_coords_outliers']   = self.y_array[probability_of_points < self.prob_threshold]
#                 
#         return fit_dict
# 
# def Python_linfit(x_true, y, y_err, errors_output = True):
#      
#     Regression_Fit, Uncertainty_Matrix, Red_Chi_Sq, Residuals   = linfit(x_true, y, y_err, cov=True, relsigma=False, chisq=True, residuals=True)
#     m_n_Matrix                                                  = [sqrt(Uncertainty_Matrix[t,t]) for t in range(2)] 
#     R_Factor                                                    = Uncertainty_Matrix[0,1]/(m_n_Matrix[0]*m_n_Matrix[1])
#     m, m_error                                                  = Regression_Fit[0], m_n_Matrix[0]
#     n, n_error                                                  = Regression_Fit[1], m_n_Matrix[1]
#     
#     if errors_output: 
#         return m, m_error, n, n_error
#     
#     else:
#         return m, n    
#     
# class LineMesurer_Log():
#     
#     def __init__(self, Conf_Folder, LinesLogHeader_Name):
#         
#         self.ColumnHeaderVector, LineLog_Width, LineLog_Format = loadtxt(Conf_Folder + LinesLogHeader_Name, dtype=str, usecols = [1, 2, 3], skiprows=1, unpack=True)
#         self.LineLog_FormatDict = dict(zip(self.ColumnHeaderVector, add(LineLog_Width, LineLog_Format)))
#         
#         self.TableHeaderSize                            = 2
#         self.ColumnWidth                                = str(50 + 2)    #len(max(self.ColumnHeaderVector,key=len)) + 2
#     
#     def CleanTableMaker(self, TableAddress, RemakeLog, ColumnsHeaders, ColumnsWidth):
#         
#         if (exists(TableAddress) == False) or (RemakeLog == True):
#             myTextFile = open(TableAddress, "w")   
#             IntroductionLines = ["-File address: " + TableAddress]
#         
#             myTextFile.write(IntroductionLines[0]+'\n')
#         
#             Width = "%" + str(ColumnsWidth) + "s"
#             
#             BlankLine = "".join(Width % i for i in ColumnsHeaders)
#             myTextFile.write(BlankLine + "\n")
#         
#             myTextFile.close() 
#     
#     def GetDataInTable(self, ColumnInput,ColumnOutput,RowOutput,TableAddress,LinesBeforeNumbers):
#         
#         TableFile = open(TableAddress,"r")
#         TableLines = TableFile.readlines()
#         TableFile.close()
#             
#         InPutIndex = "None"
#         OutPutIndex = "None"
#         HeaderRow = TableLines[LinesBeforeNumbers-1].split()
#         
#         for i in range(len(HeaderRow)):
#             Header = HeaderRow[i]
#             if ColumnInput == Header:
#                 InPutIndex = i
#                  
#         if InPutIndex == "None":
#             print "WARNING: Header not found in table (GetData1)"
#             print "Column Input: " +  ColumnInput
#             print "Headers:"
#             print HeaderRow
#           
#         for i in range(len(HeaderRow)):
#             Header = HeaderRow[i]
#             if ColumnOutput == Header:
#                 OutPutIndex = i
#         
#         if OutPutIndex == "None":
#             print "WARNING: Header not found in table (GetData2)"
#             print "ColumnOutput: " +  ColumnOutput
#             print "Headers:"
#             print HeaderRow
#             
#         Parameter = "None"
#         FoundCheck = False
#                 
#         for j in range(LinesBeforeNumbers,len(TableLines)):
#             EmissionRow = TableLines[j].split()
#             if str(RowOutput) in str(EmissionRow[InPutIndex]):                      #This is ineficient: The code should locate the row to interpolate with
#                 Parameter =  EmissionRow[OutPutIndex]
#                 FoundCheck = True
#     
#         if FoundCheck == False:
#             print "WARNING: Damned Parameter not found in table"
#             print "Column Name " + str(ColumnOutput)
#             print "Row output " + str(RowOutput)
#             print "Ouput Index " + str(OutPutIndex)
#             print TableAddress
#         
#         return Parameter
#     
#     def InsertNewLine(self, Elemental_Line_String, Current_Em_El, TableAddress, TableHeaderSize):
# 
#         TableFile = open(TableAddress,"r")
#         TableLines = TableFile.readlines()
#         TableFile.close()
# 
#         #Check if emission line has been previously measured
#         EL_Previously_Measured = False
#         
#         for j in range(TableHeaderSize,len(TableLines)):
#             EmissionRow = TableLines[j].split()
#             if float(Current_Em_El) == float(EmissionRow[2]):                      #We compare the wavelength
#                 EL_Previously_Measured = True
#         
#         if EL_Previously_Measured == False:
#         
#             Wavelengths = []
#             LineLocation = "None"
#             if len(range(TableHeaderSize,len(TableLines))) > 0:
#                 for j in range(TableHeaderSize,len(TableLines)):
#                     Wave = float(TableLines[j].split()[2])
#                     Wavelengths.append(Wave)  
#              
#                 LineLocation = bisect(Wavelengths, Current_Em_El)
#                             
#                 if LineLocation == "None":
#                     print "WARNING: Emission Line cannot be inserted in data list"
# 
#                 TableLines.insert(LineLocation + TableHeaderSize ,Elemental_Line_String+"\n")                               #Insert the line in the right location. Remember to include the size of the header!
#                 
#                 out = open(TableAddress, 'w')
#                 out.writelines(TableLines)
#                 out.close()
#                 
#             else:
#                   
#                 TableLines.append(Elemental_Line_String + "\n")
#                 out = open(TableAddress, 'w')
#                 out.writelines(TableLines)
#                 out.close()
#         return
#     
#     def Replace_Row(self, Index, RowName, DataDict, TableAddress, TableHeaderSize):
# 
#         #Load the table
#         TableFile   = open(TableAddress,"r")
#         TableLines  = TableFile.readlines()
#         TableFile.close()
# 
#         #Find row Index
#         WavelengthColumn = loadtxt(TableAddress, usecols = [2], skiprows=TableHeaderSize)
#         RowIndex = where(WavelengthColumn == RowName)[0] + TableHeaderSize
#         
# 
#         NewRow = ""
#         for i in range(len(self.ColumnHeaderVector)):
#             formatel = "%" + str(self.LineLog_FormatDict[self.ColumnHeaderVector[i]])
#             print 'formatel', formatel, type(formatel)
#             NewRow.join(formatel % DataDict[str(self.ColumnHeaderVector[i])][Index])
# 
#         print 'The final row', NewRow
#         #Replace the row   
#         TableLines[RowIndex] = NewRow + "\n"
#         
#         #Load the column
#         out = open(TableAddress, 'w')
#         out.writelines(TableLines)
#         out.close()
# 
#         return
#     
#     def ReplaceDataInTable(self, ColumnName, RowName, Data, DataFormat, TableAddress, TableHeaderSize, ColumnWidth):
#         
#         #Load the table
#         TableFile   = open(TableAddress,"r")
#         TableLines  = TableFile.readlines()
#         TableFile.close()
#         
#         #Define global columnd width
# #         Width       = "%" + str(ColumnWidth) + "s"
#         DataFormat  = '%' + str(DataFormat)
# 
#         #Find column index        
#         HeaderIndex = where(self.ColumnHeaderVector == ColumnName)[0]
#         
#         #Find row Index
#         LabelsColumn = loadtxt(TableAddress, usecols = [2], skiprows=TableHeaderSize)
#         RowIndex = where(LabelsColumn == RowName)[0] + TableHeaderSize
#         
#         #Change the new value
#         EmissionRow = TableLines[RowIndex].split()
#         if "s" in DataFormat:
#             EmissionRow[HeaderIndex] = DataFormat % (str(Data))    
#         else:
#             EmissionRow[HeaderIndex] = DataFormat % (float(Data))
#         
#         #Set row design                           
# #         NewRow = "".join(Width % z for z in EmissionRow)
#         
#         NewRow = ""
#         print 'ColumnWidth', type(ColumnWidth)
#         for i in range(len(self.ColumnHeaderVector)):
#             formatel = "%" + str(ColumnWidth[str(self.ColumnHeaderVector[i])])
#             print 'formatel', formatel, type(EmissionRow[i])
#             NewRow.join(formatel % EmissionRow[i])
#         
#         #Replace the row   
#         TableLines[RowIndex] = NewRow + "\n"
#         
#         #Load the column
#         out = open(TableAddress, 'w')
#         out.writelines(TableLines)
#         out.close()
#     
#     def CreateEmissionLineReccord(self, Current_Em_El, Current_Ion,  Current_Label, ColumnsHeaders, ColumnWidth):
#             
#         Width = "%" + ColumnWidth + "s"
#         Wave_New, Label, Ion = "None", "None", "None"
#             
#         Elemental_Line_List = [Current_Label, Current_Ion, Current_Em_El]                                                                          #This can be easily changed to a numpy array
#                     
#         for j in range(len(ColumnsHeaders)-3):
#             Elemental_Line_List.append("None")
#             
#         Elemental_Line_String = "".join(Width % z for z in Elemental_Line_List)
#         
#         return Elemental_Line_String
#     
#     def DeleteLine(self,  Current_Em_Wave, TableAddress, TableHeaderSize, Deblend_Check, List_BlendedLines, Current_BlendedGroup):
#         #WARNING CORRECT THIS BLENDED LINES
# 
#         #Open text file and import the lines
#         TableFile = open(TableAddress,"r")
#         TableLines = TableFile.readlines()
#         TableFile.close()
#         
#         #Index for the lines to delete (CHANGE THIS TO THE ENUMERATE SCHEME
#         LocationLineToDelete = None        
#         
#         #FOR THE TIME BEING THIS ONLY WORKS IF WE ARE CURRENTLY MEASURING A LINE
#         if Deblend_Check == False:
#             Measured_Wavelengths = loadtxt(TableAddress, dtype=float, skiprows= TableHeaderSize, usecols = [2])
#             LocationLineToDelete = int(where(Measured_Wavelengths==Current_Em_Wave)[0]) + TableHeaderSize
#             del TableLines[LocationLineToDelete]
#                     
#         else:
#             print 'Current_BlendedGroup', Current_BlendedGroup
#             grouped_lines_wavelengths   = List_BlendedLines[1][Current_BlendedGroup]
#             remove_list = []
# 
#             for Current_Em_Wave in grouped_lines_wavelengths:
#                 for j in range(TableHeaderSize,len(TableLines)):
#                     Wave = float(TableLines[j].split()[2])
#                     if Wave == Current_Em_Wave:
#                         remove_list.append(j)
#                         
#             TableLines = [v for i, v in enumerate(TableLines) if i not in remove_list]
#                                              
#         out = open(TableAddress, 'w')
#         out.writelines(TableLines)
#         out.close()
#         
#         return 
# 
#     def Get_DataLog_Parameter(self, Row_Label, Column_Label, TableAddress, TableHeaderSize):
#         
#         TableFile = open(TableAddress,"r")
#         LogLines = TableFile.readlines()
#         TableFile.close()
#             
#         Parameter = None
#         HeaderIndex = None
#         
#         HeaderLine = LogLines[TableHeaderSize-1].split()
#         
#         for i in range(len(HeaderLine)):
#             item = HeaderLine[i]
#             if item == Column_Label:
#                 HeaderIndex = i
#         
#         if HeaderIndex == None:
#             print "WARNING: Header not found"
#             print Column_Label
#             print HeaderLine
#             
#         for i in range(TableHeaderSize, len(LogLines)):
#             LineElements = LogLines[i].split()
#             if LineElements[0] == Row_Label:
#                 Parameter = LineElements[HeaderIndex]
#     
#         return Parameter
# 
#     def RangesR(self, Selections, RowName, TableAddress, TableHeaderSize):
#         
#         #Use this trick to check we are not measuring a line (and not load new data):       
#         if len(Selections) == 0:
#             
#             #Check if emission line has been previously measured (LOADTXT IS NOT USED PROPERLY GIVES ERROR WHEN NO LINES IN TXT FILE (ANY GIVES ERROR WHEN LOADTXT RETURNS A STRING INSTEAD OF A LIST))
#             TextFile    = open(TableAddress, "r")
#             Filelines   = TextFile.readlines()
#             TextFile.close()
#             
#             if len(Filelines) <= TableHeaderSize:
#                 EL_Previously_Measured                          = False
#             
#             elif len(Filelines) == 3:
#                 Measured_Wavelengths                            = loadtxt(TableAddress, dtype=float, skiprows= TableHeaderSize, usecols = [2])
#                 print Measured_Wavelengths, type(Measured_Wavelengths)
#                 if Measured_Wavelengths == RowName:
#                     EL_Previously_Measured                      = True
#                 else:
#                     EL_Previously_Measured                      = False
# 
#             else:
#                 Measured_Wavelengths                            = loadtxt(TableAddress, dtype=float, skiprows= TableHeaderSize, usecols = [2])
#                 if any(Measured_Wavelengths == RowName):
#                     EL_Previously_Measured                      = True
#                 else:
#                     EL_Previously_Measured                      = False
#         
#             #If it was load the data
#             if EL_Previously_Measured == True:
#                 Wave1 = self.GetDataInTable("TheoWavelength","Wave1",RowName,TableAddress,TableHeaderSize)
#                 Wave2 = self.GetDataInTable("TheoWavelength","Wave2",RowName,TableAddress,TableHeaderSize)
#                 Wave3 = self.GetDataInTable("TheoWavelength","Wave3",RowName,TableAddress,TableHeaderSize)
#                 Wave4 = self.GetDataInTable("TheoWavelength","Wave4",RowName,TableAddress,TableHeaderSize)
#                 Wave5 = self.GetDataInTable("TheoWavelength","Wave5",RowName,TableAddress,TableHeaderSize)
#                 Wave6 = self.GetDataInTable("TheoWavelength","Wave6",RowName,TableAddress,TableHeaderSize)
#                 
#                 if Wave1 != "None":
#                     Selections.append(float(Wave1))
#                 if Wave2 != "None":
#                     Selections.append(float(Wave2))
#                 if Wave3 != "None":
#                     Selections.append(float(Wave3))
#                 if Wave4 != "None":
#                     Selections.append(float(Wave4))
#                 if Wave5 != "None":
#                     Selections.append(float(Wave5))
#                 if Wave6 != "None":
#                     Selections.append(float(Wave6))     
#             
#         return Selections