from collections        import OrderedDict
from numpy              import argmax, append, exp, zeros, pi,argsort, diff, array, sqrt, square, ones, sum, mean, std, linspace, polyfit, vstack, greater, less, searchsorted, sort, empty
from numpy              import float as npfloat, log as nplog, float32, float64, invert, nan as np_nan, copy
from numpy.random       import normal as np_normal_dist
from numpy.linalg       import lstsq
from pymc               import deterministic, stochastic, Normal, Uniform, MCMC, Bernoulli, stochastic_from_dist
from scipy.integrate    import simps
from scipy.signal       import argrelextrema
from uncertainties      import ufloat, unumpy
from bces_script        import bces
from linfit_script      import linfit
from lnr_script         import kelly
from lmfit              import Parameters, minimize as lmfit_minimize, fit_report, Minimizer
from pyneb              import Atom
from lmfit.models       import GaussianModel
from timeit             import default_timer as timer
from scipy.optimize     import curve_fit
from pandas             import Series
from cairo._cairo       import Matrix
from kapteyn            import kmpfit

def lmfit_gaussian_Residual_lmfit(self, params, x, y, zerolev, Ncomps, err):
 
    return (self.gaussian_curve_SingleMixture(params.valuesdict(), x, zerolev, Ncomps) - y) / err

def gaussian_components(params, x, zerolev, Ncomps):
    
    y_vector = empty([Ncomps, len(x)])
    for i in range(Ncomps):
        idx = str(i) 
        A = params['A' + idx] 
        mu = params['mu' + idx] 
        sigma = params['sigma' + idx] 
        y_vector[i,:] = zerolev + A * exp(-(x-mu)*(x-mu)/(2*sigma*sigma))
        
    return y_vector

def gaussian_mixture(params, x, zerolev, components_index):
    
    y_model = 0.0
    for idx in components_index:
        A = params['A' + idx] 
        mu = params['mu' + idx] 
        sigma = params['sigma' + idx] 
        y_model += A * exp(-(x-mu)*(x-mu)/(2*sigma*sigma))        
    
    return y_model + zerolev

def gaussian_curve(A, mu, sigma, x, zerolev):           
    return A * exp(-(x-mu)*(x-mu)/(2 * sigma * sigma)) + zerolev

def gaussian_curveBS(ind_params, A, mu, sigma):
    x, zerolev = ind_params 
    return A * exp(-(x-mu)*(x-mu)/(2 * sigma * sigma)) + zerolev

def gaussian_MixBS(ind_params, *p):
    
    x, zerolev, Ncomps = ind_params
    y = 0.0
    for i in range(Ncomps):
        A, mu, sigma = p[i*3:(i+1)*3]
        y += A * exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))
    return y + zerolev

def residual_gauss(p, x, y, zerolev, err):
    param = p.valuesdict()
    return (gaussian_curve(param['A0'], param['mu0'], param['sigma0'], x, zerolev) - y) / err

def gauss_kmpfit(A, mu, sigma, x, zerolev):
    return A * exp(-(x-mu)*(x-mu)/(2 * sigma * sigma)) + zerolev

def residual_gauss_kmpfit(p, d):
    x, y, zerolev, err = d            # The values for x, y and weights
    A, mu, sigma = p
    return (gauss_kmpfit(A, mu, sigma, x, zerolev) - y) / err

def residual_gaussMix(p, x, y, zerolev, err, components_index):
    params = p.valuesdict()
    return (gaussian_mixture(params, x, zerolev, components_index) - y) / err

def lnprob_gaussCurve(p, x, y, zerolev, err):
    resid = residual_gauss(p, x, y, err, zerolev)
    return -0.5 * sum(((resid - y) / err)**2 + nplog(2 * pi * err**2))

def lnprob_gaussMix(p, x, y, zerolev, err, components_index):
    resid = residual_gaussMix(p, x, y, err, zerolev, components_index)
    return -0.5 * sum(((resid - y) / err)**2 + nplog(2 * pi * err**2))

def bces_regression(x_array, y_array, x_error, y_error, cov = None):
    
    #Rodrigo Nemmen, http://goo.gl/8S1Oo
    
    fit_dict    = OrderedDict()

    if cov == None:
        #This is the covariance between the measurements. If none provided it is assume error is independent between measurments for the 
        cov = zeros(len(x_array))
        
    fit_dict['methodology'] = (r'OLS(Y|X)$_{bces}$', r'OLS(X|Y)$_{bces}$', r'bisector$_{bces}$', r'Orthogonal$_{bces}$')
    
    fit_dict['m'],fit_dict['n'],fit_dict['m_error'],fit_dict['n_error'],fit_dict['cov'] = bces(x_array, x_error, y_array, y_error, cov)
    
    return fit_dict

class Fitting_Gaussians():
    
    def __init__(self):
        
        #Variables included in the series
        self.fitting_parameters =   ['idx0', 'idx1', 'idx2', 'idx3', 'idx4', 'idx5'] 
        self.fitting_parameters +=  ['area_intg', 'area_intg_er', 'flux_gauss', 'flux_gauss_er', 'flux_intg', 'flux_intg_er'] #Additionally there is (A, mu, sigma, Eqw) + idx + _norm + _norm_er
        self.fitting_parameters +=  ['m_zerolev', 'n_zerolev', 'zerolev_mean', 'zerolev_std', 'zerolev_linear', 'zerolev_width', 'continuum_width']
        self.fitting_parameters +=  ['fit_routine', 'MC_iterations', 'blended_check', 'start_treatment', 'line_number', 'add_wide_component', 'wide_component']
        self.fitting_parameters +=  ['params_lmfit', 'params_lmfit_wide', 'parameters_list', 'fit_output']
        self.fitting_parameters +=  ['Wave1', 'Wave2', 'Wave3', 'Wave4', 'Wave5', 'Wave6']
        self.fitting_parameters +=  ['group_label', 'blended_lambdas', 'blended_labels', 'blended_ions']
        self.fitting_parameters +=  ['maxLambdas', 'maxPeaks', 'x_scaler', 'y_scaler', 'x_n', 'y_n', 'zerolev_n', 'sigZerolev_n']

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

    def load_lmfit_parameters(self, x, y, zerolev, err_zerolev, n_comps, wide_component = False, A_limits = 0.30, mu_precission = 2, sigma_limit = 5):
                      
        #Scale parameters
        ind_max = argmax(y)
        self.fit_dict['x_scaler'], self.fit_dict['y_scaler'] = x[ind_max], y[ind_max]

        #Scale the range
        self.fit_dict['x_n'] = x - self.fit_dict.x_scaler
        self.fit_dict['y_n'] = y / self.fit_dict.y_scaler
        self.fit_dict['zerolev_n'] = zerolev / self.fit_dict.y_scaler
        self.fit_dict['sigZerolev_n'] = err_zerolev / self.fit_dict.y_scaler
          
        #Get line maxima and minima
        peak_wave, peak_flux, minima_wave, minima_flux = self.get_lines_peaks(ind_max, n_comps)
        
        #Store peaks location for log        
        self.fit_dict['maxLambdas'] = peak_wave + self.fit_dict['x_scaler']
        self.fit_dict['maxPeaks'] = peak_flux * self.fit_dict['y_scaler']
        self.fit_dict['params_lmfit_wide'] = None
          
        #Lmfit dictionary        
        params = Parameters()
        for i in range(n_comps):  
            index = str(i)
            params.add('A'     + index, value = peak_flux[i] - mean(self.fit_dict.zerolev_n), min = 0.0)
            params.add('mu'    + index, value = peak_wave[i], min = peak_wave[i] - mu_precission, max = peak_wave[i] + mu_precission)
            params.add('sigma' + index, value = 1, min = 0)
            params.add('fwhm'  + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
            params.add('area_G' + index, expr = '{A} * {sigma} * {sqrt2pi}'.format(A = 'A'  + index, sigma = 'sigma' + index, sqrt2pi = self.sqrt2pi))
            
        #For blended components we set the same sigma: #WARNING: We could not just delete this
        if n_comps > 1:
            small_components = range(n_comps)
            Highest_index = argmax(self.fit_dict.maxPeaks)
            
            del small_components[Highest_index]
            
            for indx in small_components: #We set the same sigma               
                expresion = 'sigma{index_big} * ((mu{index_small} + {scaler}) / (mu{index_big} + {scaler}))'.format(
                                index_big = Highest_index, index_small = str(indx), scaler = self.fit_dict['x_scaler'])
                params['sigma' + str(indx)].set(expr = expresion) 
                      
        #Special condition: Wide componentine in Halpha
        wide_params_list = []        
        if self.fit_dict.add_wide_component:
            
            #Additional fitter
            params_W = Parameters()
            
            #TRICK TO ADD AN ADDITIONAL VALUE
            n_nindex = str(n_comps)               
            params_W.add('A'  + n_nindex,       value  = 0.2, min = 0)
            params_W.add('mu' + n_nindex,       value = 0.0)
            params_W.add('sigma' + n_nindex,    value = 6, min = 3, max = 20.0)
            params_W.add('fwhm' + n_nindex,     expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + n_nindex))
            params_W.add('area_G' + n_nindex,    expr = '{A} * {sigma} * {sqrt2pi}'.format(A = 'A'  + n_nindex, sigma = 'sigma' + n_nindex, sqrt2pi = self.sqrt2pi))
            wide_params_list = params_W.keys()
            
            #Update for Nitrogen relation: Mode 1 adjuxt the fluxes
            params['area_G0'].set(expr = 'area_G2 / {N2_ratio}'.format(N2_ratio = 2.94))
            
            self.fit_dict['params_lmfit_wide'] = params_W
            
        #Store the data 
        self.fit_dict['params_lmfit'] = params
        self.fit_dict['parameters_list'] = array(params.keys() + wide_params_list) 
 
        return

    def get_lines_peaks(self, ind_max, Ncomps):
        
        target_wavelengths = None
                    
        #-- Single line  #WARNING: Should this change for fitting absoving lines
        if self.fit_dict.blended_check == False:            
            peak_flux = array([self.fit_dict.y_n[ind_max]]) 
            peak_wave = array([self.fit_dict.x_n[ind_max]])    
            minima_wave, minima_flux = 0.0, 0.0 
            
        #--Blended line
        else:
            max_index, min_index = argrelextrema(self.fit_dict.y_n, greater)[0], argrelextrema(self.fit_dict.y_n, less)[0]
            maxima_wavelengths = sort(self.fit_dict.x_n[max_index])
            minima_wavelengths = sort(self.fit_dict.x_n[min_index])
          
            #With wide component #ONLY WORKS FOR THE BLENDED HALPHA SCHEME
            if self.fit_dict.add_wide_component == False:
                target_wavelengths = array(self.fit_dict.blended_lambdas) - self.fit_dict.x_scaler
            else:
                target_wavelengths = array(self.fit_dict.blended_lambdas + [self.fit_dict.blended_lambdas[1]]) - self.fit_dict.x_scaler
            
            #Determine peak waves and fluxes     
            if len(max_index) == Ncomps:
                peak_flux, minima_flux  = self.fit_dict.y_n[max_index], self.fit_dict.y_n[min_index]
                peak_wave, minima_wave  = maxima_wavelengths, minima_wavelengths
            else:
                closest_indeces = self.search_targets_in_array(maxima_wavelengths, target_wavelengths)
                peak_wave, peak_flux = self.fit_dict.x_n[max_index][closest_indeces], self.fit_dict.y_n[max_index][closest_indeces]
            
            #Append peak waves and fluxes if wide component
            if self.fit_dict.add_wide_component:
                if len(peak_wave) ==  len(target_wavelengths) - 1:
                    peak_wave = append(peak_wave, [0])
                    peak_flux = append(peak_flux, [0.1])
            
            minima_wave, minima_flux =  self.fit_dict.x_n[min_index], self.fit_dict.y_n[min_index]
            
        return peak_wave, peak_flux, minima_wave, minima_flux

    def search_targets_in_array(self, known_array, test_array):
        
        #This function gives the indeces of the closest values within a sorted array
        
        index_sorted        = argsort(known_array)
        known_array_sorted  = known_array[index_sorted]
        known_array_middles = known_array_sorted[1:] - diff(known_array_sorted.astype('f'))/2
        idx1                = searchsorted(known_array_middles, test_array)
        indices             = index_sorted[idx1]
        
        return indices

    def fit_single_line(self, x, y, zero_lev, err_continuum, fitting_parameters, bootstrap_iterations = 1000):

        #Simple fit
        if self.fit_dict['MC_iterations'] == 1:
            fit_output = lmfit_minimize(residual_gauss, fitting_parameters, args=(x, y, zero_lev, err_continuum))
            self.fit_dict['area_intg'] = simps(y, x) - simps(zero_lev, x)
            self.fit_dict['area_intg_err'] = 0.0
             
        #Bootstrap
        else:
            mini_posterior  = Minimizer(lnprob_gaussCurve, fitting_parameters, fcn_args = ([x, y, zero_lev, err_continuum]))
            fit_output      = mini_posterior.emcee(steps=200, params = fitting_parameters)
            
            #Bootstrap for the area of the lines
            area_array = empty(bootstrap_iterations) 
            len_x_array = len(x)
            for i in range(bootstrap_iterations):
                y_new =  y + np_normal_dist(0.0, err_continuum, len_x_array)
                area_array[i] = simps(y_new, x) - simps(zero_lev, x)
            self.fit_dict['area_intg'] = mean(area_array)
            self.fit_dict['area_intg_err'] = std(area_array)           
        
        #Store the fitting parameters
        output_params = fit_output.params
        for key in self.fit_dict['parameters_list']:
            self.fit_dict[key + '_norm'] = output_params[key].value
            self.fit_dict[key + '_norm_er'] = output_params[key].stderr
            
        return

    def fit_single_line_BS(self, x, y, zero_lev, err_continuum, fitting_parameters, bootstrap_iterations = 1000):
        
        #Declare parameters and containers for the fit
        params_dict     = fitting_parameters.valuesdict()
        initial_values  = [params_dict['A0'], params_dict['mu0'], params_dict['sigma0']]
        area_array      = empty(bootstrap_iterations)
        params_array    = empty([3, bootstrap_iterations])
        n_points        = len(x)
                            
#         #Perform the fit using curve fit
#         for i in range(bootstrap_iterations):
#             y_new               = y + np_normal_dist(0, err_continuum, n_points)
#             area_array[i]       = simps(y_new, x) - simps(zero_lev, x)
#             best_vals, covar    = curve_fit(gaussian_curveBS, (x, zero_lev), y_new, p0=initial_values, maxfev = 1600)
#             params_array[:,i]   = best_vals

        #Perform the fit using kapteyn

        for i in range(bootstrap_iterations):
            y_new               = y + np_normal_dist(0, err_continuum, n_points)
            fitobj              = kmpfit.Fitter(residuals=residual_gauss_kmpfit, data=(x, y_new, zero_lev, err_continuum))
            area_array[i]       = simps(y_new, x) - simps(zero_lev, x)
            fitobj.fit(params0  = initial_values)
            params_array[:,i]   = fitobj.params


#         #Perform the fit using kapteyn
#         x_i     = copy(x)
#         y_new   = y
#         fitobj  = kmpfit.Fitter(residuals=residual_gauss_kmpfit, data=(x_i, y_new, zero_lev, err_continuum))
#         for i in range(bootstrap_iterations):
#             print i
#             print y
#             x_i[:]              = x
#             y_new[:]            = y + np_normal_dist(0, err_continuum, n_points)
#             area_array[i]       = simps(y_new, x) - simps(zero_lev, x)
#             fitobj.fit(params0  = initial_values)
#             params_array[:,i]   = fitobj.params
            
        #Compute Bootstrap output
        mean_area, std_area = mean(area_array), std(area_array)
        mean_params_array, stdev_params_array = params_array.mean(1), params_array.std(1)
        
        #Store the data
        self.fit_dict['area_intg'],     self.fit_dict['area_intg_er']   = mean_area, std_area
        self.fit_dict['A0_norm'],       self.fit_dict['A0_norm_er']     = mean_params_array[0], stdev_params_array[0]
        self.fit_dict['mu0_norm'],      self.fit_dict['mu0_norm_er']    = mean_params_array[1], stdev_params_array[1]  
        self.fit_dict['sigma0_norm'],   self.fit_dict['sigma0_norm_er'] = mean_params_array[2], stdev_params_array[2]  
                
        A = ufloat(mean_params_array[0], stdev_params_array[0]) 
        sigma = ufloat(mean_params_array[2], stdev_params_array[2])         
        fwhm0_norm = 2.354820045 * sigma
        areaG0_norm = A * sigma * self.sqrt2pi
      
        self.fit_dict['fwhm0_norm'], self.fit_dict['fwhm0_norm_er'] = fwhm0_norm.nominal_value, fwhm0_norm.std_dev
        self.fit_dict['area_G0_norm'], self.fit_dict['area_G0_norm_er'] = areaG0_norm.nominal_value, areaG0_norm.std_dev
                           
        return
  
    def fit_blended_line_emcee(self, x, y, zero_lev, err_continuum, Ncomps, fitting_parameters, add_wide_component, fitting_parameters_wide, bootstrap_iterations = 1000, MC_iterations = 200):

        #---First we integrate the brute area with all the components
        if self.fit_dict['MC_iterations'] == 1:
            self.fit_dict['area_intg'] = simps(y, x) - simps(zero_lev, x)
            self.fit_dict['area_intg_err'] = 0.0
        
        else:
            area_array  = empty(bootstrap_iterations) 
            len_x_array = len(x)
            for i in range(bootstrap_iterations):
                y_new           = y + np_normal_dist(0.0, err_continuum, len_x_array)
                area_array[i]   = simps(y_new, x) - simps(zero_lev, x)
            self.fit_dict['area_intg'] = mean(area_array)
            self.fit_dict['area_intg_err'] = std(area_array)                   
        
        #---Second we proceed to analyze as gaussian components
        idcs_components = map(str, range(Ncomps))        
        
        mini_posterior  = Minimizer(lnprob_gaussMix, fitting_parameters, fcn_args = ([x, y, zero_lev, err_continuum, idcs_components]), method='powell')
        fit_output      = mini_posterior.emcee(steps=MC_iterations, params = fitting_parameters)
        output_params   = fit_output.params
        
        if add_wide_component: #This currently only valid for Halpha

            sigma_limit = output_params['sigma1'].value
            limit_0, limit_1 = 6548.05 - self.fit_dict['x_scaler'] - sigma_limit * 1.5,     6548.05 - self.fit_dict['x_scaler'] + sigma_limit * 1.5
            limit_2, limit_3 = 0 - sigma_limit * 4,                                         0 + sigma_limit * 4 
            limit_4, limit_5 = 6583.46 - self.fit_dict['x_scaler'] - sigma_limit * 3,       6583.46 - self.fit_dict['x_scaler'] + sigma_limit * 3

            #Get the wide component area
            indeces = ((x >= limit_0) & (x <= limit_1)) + ((x >= limit_2) & (x <= limit_3)) + ((x >= limit_4) & (x <= limit_5))
            mask = invert(indeces) 
            x_wide, y_wide, zero_wide = x[mask], y[mask], zero_lev[mask]
            Ncomps_wide = ['3']
            
            #Fit the wide component without narrow component
            mini_posterior_wide = Minimizer(lnprob_gaussMix, fitting_parameters_wide, fcn_args = ([x_wide, y_wide, zero_wide, err_continuum, Ncomps_wide]), method='powell')
            fit_output_wide     = mini_posterior_wide.emcee(steps=MC_iterations, params = fitting_parameters_wide)
            output_params_wide  = fit_output_wide.params

            #Calculate wide component curve            
            y_wide_fit          = gaussian_mixture(output_params_wide.valuesdict(), x, zero_lev, Ncomps_wide)
            
            #Calculate emission line region again
            y_pure_narrow       = y - y_wide_fit + zero_lev 
            
            #Fit narrow components again
            mini_posterior          = Minimizer(lnprob_gaussMix, fitting_parameters, fcn_args = ([x, y_pure_narrow, zero_lev, err_continuum, idcs_components]), method='powell')
            fit_output_narrow       = mini_posterior.emcee(steps=MC_iterations, params=fitting_parameters)
            output_params_narrow    = fit_output_narrow.params
            
            #Combine the results from both fits
            output_params = output_params_narrow + output_params_wide

            #Add the wide component to the fit we are performing
            self.fit_dict.line_number = self.fit_dict.line_number + 1
        
        for key in self.fit_dict['parameters_list']:
            
            self.fit_dict[key + '_norm'] = output_params[key].value if output_params[key].value is not None else np_nan
            self.fit_dict[key + '_norm_er'] = output_params[key].stderr if output_params[key].stderr is not None else np_nan       
             
        return
    
    def fit_blended_line_BS(self, x, y, zero_lev, err_continuum, Ncomps, fitting_parameters, add_wide_component, fitting_parameters_wide, bootstrap_iterations = 1000, MC_iterations = 200):

        #Declare parameters and containers for the fit
        params_dict     = fitting_parameters.valuesdict()
        n_points        = len(x)
        list_parameters = []
        initial_val     = empty(3 * Ncomps)
        area_array      = empty(bootstrap_iterations)
        params_array    = empty([3 * Ncomps, bootstrap_iterations])
        
#         if add_wide_component:
#             x_scaler        = self.fit_dict['x_scaler']
#             need_to_extract = coso
            
        #Reshape parameters from dict to array (A0, mu0, sigma0, A1, mu1, sigma1...)
        for i in range(Ncomps):
            comp = str(i)
            line_param                  = [param.format(str(comp)) for param in ['A{}', 'mu{}', 'sigma{}']] 
            initial_val[i*3:(i+1)*3]    = params_dict[line_param[0]], params_dict[line_param[1]], params_dict[line_param[2]]
            list_parameters             += line_param
        
        #Perform the fit
        for j in range(bootstrap_iterations):
            y_new               = y + np_normal_dist(0, err_continuum, n_points)
            area_array[j]       = simps(y_new, x) - simps(zero_lev, x)
            best_vals, covar    = curve_fit(gaussian_MixBS, (x, zero_lev, Ncomps), y_new, p0=initial_val) # Need to poot here the new function
            params_array[:,j]   = best_vals
            
#             if add_wide_component:
#                 sigma_limit = best_vals[5]
#                 limit_0, limit_1 = 6548.05 - x_scaler - sigma_limit * 1.5,     6548.05 - x_scaler + sigma_limit * 1.5
#                 limit_2, limit_3 = 0 - sigma_limit * 4,                                         0 + sigma_limit * 4 
#                 limit_4, limit_5 = 6583.46 - x_scaler - sigma_limit * 3,       6583.46 - x_scaler + sigma_limit * 3            
#       
#                 indeces = ((x >= limit_0) & (x <= limit_1)) + ((x >= limit_2) & (x <= limit_3)) + ((x >= limit_4) & (x <= limit_5))
#                 mask = invert(indeces) 
#                 x_wide, y_wide, zero_wide = x[mask], y_new[mask], zero_lev[mask]
#                 
#                 best_vals_wide, covar = curve_fit(gaussian_curveBS, (x_wide, zero_wide), y_wide, p0=initial_wide_values)
#                 
#                 #Calculate wide component curve            
#                 y_wide_fit = gaussian_curveBS((x, zero_lev), *best_vals_wide) 
#                 
#                 #Calculate emission line region again
#                 y_pure_narrow = y_new - y_wide_fit + zero_lev 
#                 
#                 #Recalculate the narrow components
#                 best_vals = curve_fit(gaussian_MixBS, (x, zero_lev, Ncomps), y_pure_narrow, p0=initial_val)        
# 
#         if add_wide_component:
#             self.fit_dict.line_number = self.fit_dict.line_number + 1
#             concantenate Matrix
            
        #Compute Bootstrap output
        mean_area, std_area = mean(area_array), std(area_array)
        mean_params_array, stdev_params_array = params_array.mean(1), params_array.std(1)

        #Store the data
        self.fit_dict[['area_intg', 'area_intg_er']] = mean_area, std_area
        
        #Return to a dictionary format
        for m in range(Ncomps):
            comp            = str(m)
            line_param      = [param.format(str(comp)) for param in ['A{}_norm', 'mu{}_norm', 'sigma{}_norm']]
            line_param_er   = [param.format(str(comp)) for param in ['A{}_norm_er', 'mu{}_norm_er', 'sigma{}_norm_er']]
            for n in range(len(line_param)):
                self.fit_dict[line_param[n]] = mean_params_array[m*3:(m+1)*3][n]
                self.fit_dict[line_param_er[n]] = stdev_params_array[m*3:(m+1)*3][n]
                
        #Calculate the gaussian area
        A_keys, Aerr_keys = ['A{}_norm'.format(str(i)) for i in range(Ncomps)], ['A{}_norm_er'.format(str(i)) for i in range(Ncomps)]
        sigma_keys, sigmaerr_keys = ['sigma{}_norm'.format(str(i)) for i in range(Ncomps)], ['sigma{}_norm_er'.format(str(i)) for i in range(Ncomps)]   
        A_vector = unumpy.uarray(self.fit_dict[A_keys].values, self.fit_dict[Aerr_keys].values) 
        sigma_vector = unumpy.uarray(self.fit_dict[sigma_keys].values, self.fit_dict[sigmaerr_keys].values) 

        fwhm_vector = 2.354820045 * sigma_vector
        areaG_vector = A_vector * sigma_vector * self.sqrt2pi 
    
        #Add areas to dict
        fwhm_keys, fwhmerr_keys = ['fwhm{}'.format(str(i)) for i in range(Ncomps)], ['fwhm{}_er'.format(str(i)) for i in range(Ncomps)]           
        AreaG_keys, AreaGerr_keys = ['area_G{}_norm'.format(str(i)) for i in range(Ncomps)], ['area_G{}_norm_er'.format(str(i)) for i in range(Ncomps)]
        fwhm_nominal, fwhm_std_dev = unumpy.nominal_values(fwhm_vector), unumpy.std_devs(fwhm_vector)
        AreaG_nominal, AreaG_std_dev = unumpy.nominal_values(areaG_vector), unumpy.std_devs(areaG_vector)
        
        for m in range(len(AreaG_keys)):
            self.fit_dict[fwhm_keys[m]]     = fwhm_nominal[m]
            self.fit_dict[fwhmerr_keys[m]]  = fwhm_std_dev[m]
            self.fit_dict[AreaG_keys[m]]    = AreaG_nominal[m]
            self.fit_dict[AreaGerr_keys[m]] = AreaG_std_dev[m]      
                            
        return

    def fit_blended_line_BS_lmfit(self, x, y, zero_lev, err_continuum, Ncomps, fitting_parameters, add_wide_component, fitting_parameters_wide, bootstrap_iterations = 500, MC_iterations = 200):
                        
        #Prepare some data in advance
        n_points        = len(x)
        n_parameters    = len(self.fit_dict['parameters_list'])
        params_list     = self.fit_dict['parameters_list']
        range_param     = range(n_parameters)
        range_boots     = range(bootstrap_iterations)
        area_array      = empty(bootstrap_iterations)
        idcs_components = map(str, range(Ncomps))
        params_matrix   = empty((n_parameters, bootstrap_iterations))
        errors_matrix   = empty((n_parameters, bootstrap_iterations))
        Ncomps_wide     = ['3']
        
        #Loop through the iterations:
        for i in range_boots:
            
            y_new           = y + np_normal_dist(0.0, err_continuum, n_points)
            area_array[i]   = simps(y_new, x) - simps(zero_lev, x)
            fit_Output      = lmfit_minimize(residual_gaussMix, fitting_parameters, args=(x, y_new, zero_lev, err_continuum, idcs_components))
            output_params   = fit_Output.params
            
            #Case with a wide component
            if add_wide_component:
                
                #Only works for Halpha
                sigma_limit = output_params['sigma1'].value
                limit_0, limit_1 = 6548.05 - self.fit_dict['x_scaler'] - sigma_limit * 1.5, 6548.05 - self.fit_dict['x_scaler'] + sigma_limit * 1.5
                limit_2, limit_3 = 0 - sigma_limit * 4,                                     0 + sigma_limit * 4 
                limit_4, limit_5 = 6583.46 - self.fit_dict['x_scaler'] - sigma_limit * 3,   6583.46 - self.fit_dict['x_scaler'] + sigma_limit * 3
                
                indeces = ((x >= limit_0) & (x <= limit_1)) + ((x >= limit_2) & (x <= limit_3)) + ((x >= limit_4) & (x <= limit_5))
                mask = invert(indeces) 
                x_wide, y_wide, zero_wide = x[mask], y[mask], zero_lev[mask]

                #Fit the wide component without narrow component
                fit_output_wide     = lmfit_minimize(residual_gaussMix, fitting_parameters_wide, args=(x_wide, y_wide, zero_wide, err_continuum, Ncomps_wide))
                output_params_wide  = fit_output_wide.params
                
                #Get wide curve
                y_wide = gaussian_mixture(output_params_wide.valuesdict(), x, zero_lev, Ncomps_wide)

                #Calculate emission line region again
                y_pure_narrow       = y - y_wide + zero_lev 
                
                #Fit narrow components again
                fit_Output          = lmfit_minimize(residual_gaussMix, fitting_parameters, args=(x, y_pure_narrow, zero_lev, err_continuum, idcs_components))
                out_params_narrow   = fit_Output.params           

                #Combine the results from both fits
                output_params       = out_params_narrow + output_params_wide
                
            #save parameters to array
            for j in range_param:
                params_matrix[j,i] = output_params[params_list[j]].value
                errors_matrix[j,i] = output_params[params_list[j]].stderr

        #Calculate mean and std deviation values
        mean_area, std_area = mean(area_array), std(area_array)
        mean_params_array, stdev_params_array = params_matrix.mean(1), params_matrix.std(1)
        errorsmean_array = errors_matrix.mean(1)
                
        #Return the data to a dictionary format
        self.fit_dict['area_intg'], self.fit_dict['area_intg_er'] = mean_area, std_area
        
        for j in range(len(self.fit_dict['parameters_list'])):
            param = self.fit_dict['parameters_list'][j]
            self.fit_dict[param+'_norm'] = mean_params_array[j]
            self.fit_dict[param+'_norm_er'] = errorsmean_array[j]

        #Increase the number of components if wide observed
        self.fit_dict.line_number = self.fit_dict.line_number + 1 if add_wide_component else self.fit_dict.line_number

        return
        
    def rescale_lmfit_params(self, line_wave, line_flux, Ncomps, x_scale, y_scale, fitting_method):
        
        #Scale integrated area
        self.fit_dict['flux_intg'], self.fit_dict['flux_intg_er'] = self.fit_dict['area_intg']* y_scale, self.fit_dict['area_intg_er'] * y_scale
                
        for i in range(Ncomps): #WARNING: We have two loops which almost do the same we could forget the previous
            index = str(i)
            self.fit_dict['A'+index]                = self.fit_dict['A'+index+'_norm'] * y_scale
            self.fit_dict['A'+index+'_er']          = self.fit_dict['A'+index+'_norm_er'] * y_scale
            self.fit_dict['mu'+index]               = self.fit_dict['mu'+index+'_norm'] + x_scale
            self.fit_dict['mu'+index+'_er']         = self.fit_dict['mu'+index+'_norm_er']
            self.fit_dict['sigma'+index]            = self.fit_dict['sigma' + index + '_norm']
            self.fit_dict['sigma'+index+'_er']      = self.fit_dict['sigma' + index + '_norm_er']            
            self.fit_dict['flux_gauss'+index]       = self.fit_dict['area_G' + index + '_norm'] * y_scale
            self.fit_dict['flux_gauss'+index+'_er'] = self.fit_dict['area_G' + index + '_norm_er'] * y_scale        
            
        #Calculate the gaussian curves for plotting
        self.fit_dict['x_resample']         = linspace(line_wave[0], line_wave[-1], 50 * Ncomps)
        self.fit_dict['zerolev_resample']   = self.fit_dict['m_zerolev'] * self.fit_dict['x_resample'] + self.fit_dict['n_zerolev']
        
        if self.fit_dict.blended_check == False:
            self.fit_dict['y_resample'] = gaussian_curve(self.fit_dict['A0'], self.fit_dict['mu0'], self.fit_dict['sigma0'], self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'])
        else:
            self.fit_dict['y_resample'] = gaussian_mixture(self.fit_dict, self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'], map(str, range(Ncomps)))
            self.fit_dict['y_comps']    = gaussian_components(self.fit_dict, self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'], Ncomps)
    
        return

class Fitting_Gaussians_v2():
    
    def __init__(self):
        
        #Variables included in the series
        self.fitting_parameters =   ['idx0', 'idx1', 'idx2', 'idx3', 'idx4', 'idx5'] 
        self.fitting_parameters +=  ['area_intg', 'area_intg_er', 'flux_gauss', 'flux_gauss_er', 'flux_intg', 'flux_intg_er'] #Additionally there is (A, mu, sigma, Eqw) + idx + _norm + _norm_er
        self.fitting_parameters +=  ['m_zerolev', 'n_zerolev', 'zerolev_mean', 'zerolev_std', 'zerolev_linear', 'zerolev_width', 'continuum_width']
        self.fitting_parameters +=  ['fit_routine', 'MC_iterations', 'blended_check', 'start_treatment', 'line_number', 'add_wide_component', 'wide_component']
        self.fitting_parameters +=  ['params_lmfit', 'params_lmfit_wide', 'parameters_list', 'fit_output']
        self.fitting_parameters +=  ['Wave1', 'Wave2', 'Wave3', 'Wave4', 'Wave5', 'Wave6']
        self.fitting_parameters +=  ['group_label', 'blended_lambdas', 'blended_labels', 'blended_ions']
        self.fitting_parameters +=  ['maxLambdas', 'maxPeaks', 'x_scaler', 'y_scaler', 'x_n', 'y_n', 'zerolev_n', 'sigZerolev_n']
        self.fitting_parameters +=  ['Blue_wave_zerolev', 'Red_wave_zerolev', 'Blue_flux_zerolev', 'Red_flux_zerolev']
        self.fitting_parameters +=  ['A0_norm', 'A0_norm_er', 'mu0_norm', 'mu0_norm_er', 'sigma0_norm', 'sigma0_norm_er']        
        self.fitting_parameters +=  ['fwhm0_norm', 'fwhm0_norm_er', 'area_G0_norm', 'area_G0_norm_er']
        self.fitting_parameters +=  ['A0', 'mu0', 'sigma0', 'flux_gauss0']        
        self.fitting_parameters +=  ['A0_er', 'mu0_er', 'sigma0_er', 'flux_gauss0_er']        
        self.fitting_parameters +=  ['x_resample', 'y_resample', 'zerolev_resample', 'y_comps']        
        self.fitting_parameters +=  ['eqw0', 'eqw0_er']        
                    
        #Ordered dictionary for lmfit
        self.params_lmfit       = Parameters()

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

    def load_lmfit_parameters(self, x, y, zerolev, err_zerolev, n_comps, wide_component = False, A_limits = 0.30, mu_precission = 2, sigma_limit = 5):
        
        #WARNING: lmfit uses ordered dictionaries... this consumes resources
        
        #Scale parameters
        ind_max = argmax(y)
        self.fit_dict['x_scaler'], self.fit_dict['y_scaler'] = x[ind_max], y[ind_max]
        
        #Scale the range
        self.fit_dict['x_n']                = x - self.fit_dict.x_scaler
        self.fit_dict['y_n']                = y / self.fit_dict.y_scaler
        self.fit_dict['zerolev_n']          = zerolev / self.fit_dict.y_scaler
        self.fit_dict['sigZerolev_n']       = err_zerolev / self.fit_dict.y_scaler

        #Get line maxima and minima
        peak_wave, peak_flux, minima_wave, minima_flux = self.get_lines_peaks(ind_max, n_comps)
        
        #Store peaks location for log        
        self.fit_dict['maxLambdas']         = peak_wave + self.fit_dict['x_scaler']
        self.fit_dict['maxPeaks']           = peak_flux * self.fit_dict['y_scaler']
        self.fit_dict['params_lmfit_wide']  = None
        
        #Clear lmfit dictionary        
        self.params_lmfit.clear()
        
        #Just one component
        if n_comps:
            self.params_lmfit.add('A0',     value = peak_flux[0] - mean(self.fit_dict.zerolev_n), min = 0.0)
            self.params_lmfit.add('mu0',    value = peak_wave[0], min = peak_wave[0] - mu_precission, max = peak_wave[0] + mu_precission)
            self.params_lmfit.add('sigma0', value = 1, min = 0)
            self.params_lmfit.add('fwhm0',  expr = '2.354820045 * {sigma}'.format(sigma = 'sigma0'))
            self.params_lmfit.add('area_G0', expr = 'A0 * sigma0 * {sqrt2pi}'.format(sqrt2pi = self.sqrt2pi))
        
        #Blended structure            
        else:
            for i in range(n_comps):  
                index = str(i)
                self.params_lmfit.add('A'     + index, value = peak_flux[i] - mean(self.fit_dict.zerolev_n), min = 0.0)
                self.params_lmfit.add('mu'    + index, value = peak_wave[i], min = peak_wave[i] - mu_precission, max = peak_wave[i] + mu_precission)
                self.params_lmfit.add('sigma' + index, value = 1, min = 0)
                self.params_lmfit.add('fwhm'  + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
                self.params_lmfit.add('area_G' + index, expr = '{A} * {sigma} * {sqrt2pi}'.format(A = 'A'  + index, sigma = 'sigma' + index, sqrt2pi = self.sqrt2pi))        
                
        #For blended components we set the same sigma: #WARNING: We could not just delete this
        if n_comps > 1:
            small_components = range(n_comps)
            Highest_index = argmax(self.fit_dict.maxPeaks)
            
            del small_components[Highest_index]
            
            for indx in small_components: #We set the same sigma               
                expresion = 'sigma{index_big} * ((mu{index_small} + {scaler}) / (mu{index_big} + {scaler}))'.format(
                                index_big = Highest_index, index_small = str(indx), scaler = self.fit_dict['x_scaler'])
                self.params_lmfit['sigma' + str(indx)].set(expr = expresion) 
                      
        #Special condition: Wide componentine in Halpha
        wide_params_list = []        
        if self.fit_dict.add_wide_component:
            
            #Additional fitter
            params_W = Parameters()
            
            #TRICK TO ADD AN ADDITIONAL VALUE
            n_nindex = str(n_comps)               
            params_W.add('A'  + n_nindex,       value  = 0.2, min = 0)
            params_W.add('mu' + n_nindex,       value = 0.0)
            params_W.add('sigma' + n_nindex,    value = 6, min = 3, max = 20.0)
            params_W.add('fwhm' + n_nindex,     expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + n_nindex))
            params_W.add('area_G' + n_nindex,    expr = '{A} * {sigma} * {sqrt2pi}'.format(A = 'A'  + n_nindex, sigma = 'sigma' + n_nindex, sqrt2pi = self.sqrt2pi))
            wide_params_list = params_W.keys()
            
            #Update for Nitrogen relation: Mode 1 adjuxt the fluxes
            self.params_lmfit['area_G0'].set(expr = 'area_G2 / {N2_ratio}'.format(N2_ratio = 2.94))
            self.fit_dict['params_lmfit_wide'] = params_W
            
        #Store the data 
        self.fit_dict['params_lmfit'] = self.params_lmfit
        self.fit_dict['parameters_list'] = array(self.params_lmfit.keys() + wide_params_list) 

        return

    def get_lines_peaks(self, ind_max, Ncomps):
        
        target_wavelengths = None
                    
        #-- Single line  #WARNING: Should this change for fitting absoving lines
        if self.fit_dict.blended_check == False:            
            peak_flux = array([self.fit_dict.y_n[ind_max]]) 
            peak_wave = array([self.fit_dict.x_n[ind_max]])    
            minima_wave, minima_flux = 0.0, 0.0 
            
        #--Blended line
        else:
            max_index, min_index = argrelextrema(self.fit_dict.y_n, greater)[0], argrelextrema(self.fit_dict.y_n, less)[0]
            maxima_wavelengths = sort(self.fit_dict.x_n[max_index])
            minima_wavelengths = sort(self.fit_dict.x_n[min_index])
          
            #With wide component #ONLY WORKS FOR THE BLENDED HALPHA SCHEME
            if self.fit_dict.add_wide_component == False:
                target_wavelengths = array(self.fit_dict.blended_lambdas) - self.fit_dict.x_scaler
            else:
                target_wavelengths = array(self.fit_dict.blended_lambdas + [self.fit_dict.blended_lambdas[1]]) - self.fit_dict.x_scaler
            
            #Determine peak waves and fluxes     
            if len(max_index) == Ncomps:
                peak_flux, minima_flux  = self.fit_dict.y_n[max_index], self.fit_dict.y_n[min_index]
                peak_wave, minima_wave  = maxima_wavelengths, minima_wavelengths
            else:
                closest_indeces = self.search_targets_in_array(maxima_wavelengths, target_wavelengths)
                peak_wave, peak_flux = self.fit_dict.x_n[max_index][closest_indeces], self.fit_dict.y_n[max_index][closest_indeces]
            
            #Append peak waves and fluxes if wide component
            if self.fit_dict.add_wide_component:
                if len(peak_wave) ==  len(target_wavelengths) - 1:
                    peak_wave = append(peak_wave, [0])
                    peak_flux = append(peak_flux, [0.1])
            
            minima_wave, minima_flux =  self.fit_dict.x_n[min_index], self.fit_dict.y_n[min_index]
            
        return peak_wave, peak_flux, minima_wave, minima_flux

    def search_targets_in_array(self, known_array, test_array):
        
        #This function gives the indeces of the closest values within a sorted array
        
        index_sorted        = argsort(known_array)
        known_array_sorted  = known_array[index_sorted]
        known_array_middles = known_array_sorted[1:] - diff(known_array_sorted.astype('f'))/2
        idx1                = searchsorted(known_array_middles, test_array)
        indices             = index_sorted[idx1]
        
        return indices

    def fit_single_line(self, x, y, zero_lev, err_continuum, fitting_parameters):

        #Declare parameters and containers for the fit
        params_dict         = fitting_parameters.valuesdict()
        initial_values      = (params_dict['A0'], params_dict['mu0'], params_dict['sigma0'])
                            
        #Perform the fit using kapteyn
        fitobj              = kmpfit.Fitter(residuals=residual_gauss_kmpfit, data=(x, y, zero_lev, err_continuum))
        mean_area           = simps(y, x) - simps(zero_lev, x)
        fitobj.fit(params0  = initial_values)
        params_array        = fitobj.params
            
        #Compute Bootstrap output        
        A                   = params_array[0]
        mu                  = params_array[1]
        sigma               = params_array[2]         
        fwhm0_norm          = 2.354820045 * sigma
        areaG0_norm         = A * sigma * self.sqrt2pi

        #Store the data
        self.fit_dict['area_intg'],     self.fit_dict['area_intg_er']       = mean_area, 0.0
        self.fit_dict['A0_norm'],       self.fit_dict['A0_norm_er']         = A, 0.0
        self.fit_dict['mu0_norm'],      self.fit_dict['mu0_norm_er']        = mu, 0.0  
        self.fit_dict['sigma0_norm'],   self.fit_dict['sigma0_norm_er']     = sigma, 0.0  
        self.fit_dict['fwhm0_norm'],    self.fit_dict['fwhm0_norm_er']      = fwhm0_norm, 0.0
        self.fit_dict['area_G0_norm'],  self.fit_dict['area_G0_norm_er']    = areaG0_norm, 0.0
              
        return

    def fit_single_line_BS(self, x, y, zero_lev, err_continuum, fitting_parameters, bootstrap_iterations = 1000):
        
        #Declare parameters and containers for the fit
        params_dict     = fitting_parameters.valuesdict()
        initial_values  = [params_dict['A0'], params_dict['mu0'], params_dict['sigma0']]
        area_array      = empty(bootstrap_iterations)
        params_array    = empty([3, bootstrap_iterations])
        n_points        = len(x)
                            
        #Perform the fit using kapteyn
        for i in range(bootstrap_iterations):
            y_new               = y + np_normal_dist(0, err_continuum, n_points)
            fitobj              = kmpfit.Fitter(residuals=residual_gauss_kmpfit, data=(x, y_new, zero_lev, err_continuum))
            area_array[i]       = simps(y_new, x) - simps(zero_lev, x)
            fitobj.fit(params0  = initial_values)
            params_array[:,i]   = fitobj.params
            
        #Compute Bootstrap output
        mean_area, std_area = mean(area_array), std(area_array)
        mean_params_array, stdev_params_array = params_array.mean(1), params_array.std(1)
        
        A                       = ufloat(mean_params_array[0], stdev_params_array[0]) 
        sigma                   = ufloat(mean_params_array[2], stdev_params_array[2])         
        fwhm0_norm              = 2.354820045 * sigma
        areaG0_norm             = A * sigma * self.sqrt2pi

        #Store the data
        self.fit_dict['area_intg'],     self.fit_dict['area_intg_er']       = mean_area, std_area
        self.fit_dict['A0_norm'],       self.fit_dict['A0_norm_er']         = mean_params_array[0], stdev_params_array[0]
        self.fit_dict['mu0_norm'],      self.fit_dict['mu0_norm_er']        = mean_params_array[1], stdev_params_array[1]  
        self.fit_dict['sigma0_norm'],   self.fit_dict['sigma0_norm_er']     = mean_params_array[2], stdev_params_array[2]  
        self.fit_dict['fwhm0_norm'],    self.fit_dict['fwhm0_norm_er']      = fwhm0_norm.nominal_value, fwhm0_norm.std_dev
        self.fit_dict['area_G0_norm'],  self.fit_dict['area_G0_norm_er']    = areaG0_norm.nominal_value, areaG0_norm.std_dev
                           
        return
      
    def fit_blended_line_BS(self, x, y, zero_lev, err_continuum, Ncomps, fitting_parameters, add_wide_component, fitting_parameters_wide, bootstrap_iterations = 1000, MC_iterations = 200):

        #Declare parameters and containers for the fit
        params_dict     = fitting_parameters.valuesdict()
        n_points        = len(x)
        list_parameters = []
        initial_val     = empty(3 * Ncomps)
        area_array      = empty(bootstrap_iterations)
        params_array    = empty([3 * Ncomps, bootstrap_iterations])
        
#         if add_wide_component:
#             x_scaler        = self.fit_dict['x_scaler']
#             need_to_extract = coso
            
        #Reshape parameters from dict to array (A0, mu0, sigma0, A1, mu1, sigma1...)
        for i in range(Ncomps):
            comp = str(i)
            line_param                  = [param.format(str(comp)) for param in ['A{}', 'mu{}', 'sigma{}']] 
            initial_val[i*3:(i+1)*3]    = params_dict[line_param[0]], params_dict[line_param[1]], params_dict[line_param[2]]
            list_parameters             += line_param
        
        #Perform the fit
        for j in range(bootstrap_iterations):
            y_new               = y + np_normal_dist(0, err_continuum, n_points)
            area_array[j]       = simps(y_new, x) - simps(zero_lev, x)
            best_vals, covar    = curve_fit(gaussian_MixBS, (x, zero_lev, Ncomps), y_new, p0=initial_val) # Need to poot here the new function
            params_array[:,j]   = best_vals
            
#             if add_wide_component:
#                 sigma_limit = best_vals[5]
#                 limit_0, limit_1 = 6548.05 - x_scaler - sigma_limit * 1.5,     6548.05 - x_scaler + sigma_limit * 1.5
#                 limit_2, limit_3 = 0 - sigma_limit * 4,                                         0 + sigma_limit * 4 
#                 limit_4, limit_5 = 6583.46 - x_scaler - sigma_limit * 3,       6583.46 - x_scaler + sigma_limit * 3            
#       
#                 indeces = ((x >= limit_0) & (x <= limit_1)) + ((x >= limit_2) & (x <= limit_3)) + ((x >= limit_4) & (x <= limit_5))
#                 mask = invert(indeces) 
#                 x_wide, y_wide, zero_wide = x[mask], y_new[mask], zero_lev[mask]
#                 
#                 best_vals_wide, covar = curve_fit(gaussian_curveBS, (x_wide, zero_wide), y_wide, p0=initial_wide_values)
#                 
#                 #Calculate wide component curve            
#                 y_wide_fit = gaussian_curveBS((x, zero_lev), *best_vals_wide) 
#                 
#                 #Calculate emission line region again
#                 y_pure_narrow = y_new - y_wide_fit + zero_lev 
#                 
#                 #Recalculate the narrow components
#                 best_vals = curve_fit(gaussian_MixBS, (x, zero_lev, Ncomps), y_pure_narrow, p0=initial_val)        
# 
#         if add_wide_component:
#             self.fit_dict.line_number = self.fit_dict.line_number + 1
#             concantenate Matrix
            
        #Compute Bootstrap output
        mean_area, std_area = mean(area_array), std(area_array)
        mean_params_array, stdev_params_array = params_array.mean(1), params_array.std(1)

        #Store the data
        self.fit_dict[['area_intg', 'area_intg_er']] = mean_area, std_area
        
        #Return to a dictionary format
        for m in range(Ncomps):
            comp            = str(m)
            line_param      = [param.format(str(comp)) for param in ['A{}_norm', 'mu{}_norm', 'sigma{}_norm']]
            line_param_er   = [param.format(str(comp)) for param in ['A{}_norm_er', 'mu{}_norm_er', 'sigma{}_norm_er']]
            for n in range(len(line_param)):
                self.fit_dict[line_param[n]] = mean_params_array[m*3:(m+1)*3][n]
                self.fit_dict[line_param_er[n]] = stdev_params_array[m*3:(m+1)*3][n]
                
        #Calculate the gaussian area
        A_keys, Aerr_keys = ['A{}_norm'.format(str(i)) for i in range(Ncomps)], ['A{}_norm_er'.format(str(i)) for i in range(Ncomps)]
        sigma_keys, sigmaerr_keys = ['sigma{}_norm'.format(str(i)) for i in range(Ncomps)], ['sigma{}_norm_er'.format(str(i)) for i in range(Ncomps)]   
        A_vector = unumpy.uarray(self.fit_dict[A_keys].values, self.fit_dict[Aerr_keys].values) 
        sigma_vector = unumpy.uarray(self.fit_dict[sigma_keys].values, self.fit_dict[sigmaerr_keys].values) 

        fwhm_vector = 2.354820045 * sigma_vector
        areaG_vector = A_vector * sigma_vector * self.sqrt2pi 
    
        #Add areas to dict
        fwhm_keys, fwhmerr_keys = ['fwhm{}'.format(str(i)) for i in range(Ncomps)], ['fwhm{}_er'.format(str(i)) for i in range(Ncomps)]           
        AreaG_keys, AreaGerr_keys = ['area_G{}_norm'.format(str(i)) for i in range(Ncomps)], ['area_G{}_norm_er'.format(str(i)) for i in range(Ncomps)]
        fwhm_nominal, fwhm_std_dev = unumpy.nominal_values(fwhm_vector), unumpy.std_devs(fwhm_vector)
        AreaG_nominal, AreaG_std_dev = unumpy.nominal_values(areaG_vector), unumpy.std_devs(areaG_vector)
        
        for m in range(len(AreaG_keys)):
            self.fit_dict[fwhm_keys[m]]     = fwhm_nominal[m]
            self.fit_dict[fwhmerr_keys[m]]  = fwhm_std_dev[m]
            self.fit_dict[AreaG_keys[m]]    = AreaG_nominal[m]
            self.fit_dict[AreaGerr_keys[m]] = AreaG_std_dev[m]      
                            
        return

    def fit_blended_line_BS_lmfit(self, x, y, zero_lev, err_continuum, Ncomps, fitting_parameters, add_wide_component, fitting_parameters_wide, bootstrap_iterations = 500, MC_iterations = 200):
                        
        #Prepare some data in advance
        n_points        = len(x)
        n_parameters    = len(self.fit_dict['parameters_list'])
        params_list     = self.fit_dict['parameters_list']
        range_param     = range(n_parameters)
        range_boots     = range(bootstrap_iterations)
        area_array      = empty(bootstrap_iterations)
        idcs_components = map(str, range(Ncomps))
        params_matrix   = empty((n_parameters, bootstrap_iterations))
        errors_matrix   = empty((n_parameters, bootstrap_iterations))
        Ncomps_wide     = ['3']
        
        #Loop through the iterations:
        for i in range_boots:
            
            y_new           = y + np_normal_dist(0.0, err_continuum, n_points)
            area_array[i]   = simps(y_new, x) - simps(zero_lev, x)
            fit_Output      = lmfit_minimize(residual_gaussMix, fitting_parameters, args=(x, y_new, zero_lev, err_continuum, idcs_components))
            output_params   = fit_Output.params
            
            #Case with a wide component
            if add_wide_component:
                
                #Only works for Halpha
                sigma_limit = output_params['sigma1'].value
                limit_0, limit_1 = 6548.05 - self.fit_dict['x_scaler'] - sigma_limit * 1.5, 6548.05 - self.fit_dict['x_scaler'] + sigma_limit * 1.5
                limit_2, limit_3 = 0 - sigma_limit * 4,                                     0 + sigma_limit * 4 
                limit_4, limit_5 = 6583.46 - self.fit_dict['x_scaler'] - sigma_limit * 3,   6583.46 - self.fit_dict['x_scaler'] + sigma_limit * 3
                
                indeces = ((x >= limit_0) & (x <= limit_1)) + ((x >= limit_2) & (x <= limit_3)) + ((x >= limit_4) & (x <= limit_5))
                mask = invert(indeces) 
                x_wide, y_wide, zero_wide = x[mask], y[mask], zero_lev[mask]

                #Fit the wide component without narrow component
                fit_output_wide     = lmfit_minimize(residual_gaussMix, fitting_parameters_wide, args=(x_wide, y_wide, zero_wide, err_continuum, Ncomps_wide))
                output_params_wide  = fit_output_wide.params
                
                #Get wide curve
                y_wide = gaussian_mixture(output_params_wide.valuesdict(), x, zero_lev, Ncomps_wide)

                #Calculate emission line region again
                y_pure_narrow       = y - y_wide + zero_lev 
                
                #Fit narrow components again
                fit_Output          = lmfit_minimize(residual_gaussMix, fitting_parameters, args=(x, y_pure_narrow, zero_lev, err_continuum, idcs_components))
                out_params_narrow   = fit_Output.params           

                #Combine the results from both fits
                output_params       = out_params_narrow + output_params_wide
                
            #save parameters to array
            for j in range_param:
                params_matrix[j,i] = output_params[params_list[j]].value
                errors_matrix[j,i] = output_params[params_list[j]].stderr

        #Calculate mean and std deviation values
        mean_area, std_area = mean(area_array), std(area_array)
        mean_params_array, stdev_params_array = params_matrix.mean(1), params_matrix.std(1)
        errorsmean_array = errors_matrix.mean(1)
                
        #Return the data to a dictionary format
        self.fit_dict['area_intg'], self.fit_dict['area_intg_er'] = mean_area, std_area
        
        for j in range(len(self.fit_dict['parameters_list'])):
            param = self.fit_dict['parameters_list'][j]
            self.fit_dict[param+'_norm'] = mean_params_array[j]
            self.fit_dict[param+'_norm_er'] = errorsmean_array[j]

        #Increase the number of components if wide observed
        self.fit_dict.line_number = self.fit_dict.line_number + 1 if add_wide_component else self.fit_dict.line_number

        return
        
    def rescale_lmfit_params(self, line_wave, line_flux, Ncomps, x_scale, y_scale, fitting_method):
                
        #Scale integrated area
        self.fit_dict['flux_intg'], self.fit_dict['flux_intg_er'] = self.fit_dict['area_intg']* y_scale, self.fit_dict['area_intg_er'] * y_scale
                
        for i in range(Ncomps): #WARNING: We have two loops which almost do the same we could forget the previous
            index = str(i)
            self.fit_dict['A'+index]                = self.fit_dict['A'+index+'_norm'] * y_scale
            self.fit_dict['A'+index+'_er']          = self.fit_dict['A'+index+'_norm_er'] * y_scale
            self.fit_dict['mu'+index]               = self.fit_dict['mu'+index+'_norm'] + x_scale
            self.fit_dict['mu'+index+'_er']         = self.fit_dict['mu'+index+'_norm_er']
            self.fit_dict['sigma'+index]            = self.fit_dict['sigma' + index + '_norm']
            self.fit_dict['sigma'+index+'_er']      = self.fit_dict['sigma' + index + '_norm_er']            
            self.fit_dict['flux_gauss'+index]       = self.fit_dict['area_G' + index + '_norm'] * y_scale
            self.fit_dict['flux_gauss'+index+'_er'] = self.fit_dict['area_G' + index + '_norm_er'] * y_scale        
            
        #Calculate the gaussian curves for plotting
        self.fit_dict['x_resample']         = linspace(line_wave[0], line_wave[-1], 50 * Ncomps)
        self.fit_dict['zerolev_resample']   = self.fit_dict['m_zerolev'] * self.fit_dict['x_resample'] + self.fit_dict['n_zerolev']
        
        if self.fit_dict.blended_check == False:
            self.fit_dict['y_resample'] = gaussian_curve(self.fit_dict['A0'], self.fit_dict['mu0'], self.fit_dict['sigma0'], self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'])
        else:
            self.fit_dict['y_resample'] = gaussian_mixture(self.fit_dict, self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'], map(str, range(Ncomps)))
            self.fit_dict['y_comps']    = gaussian_components(self.fit_dict, self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'], Ncomps)
    
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

