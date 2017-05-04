from collections        import OrderedDict
from numpy              import argmax, exp, zeros, pi, asarray, argsort, diff, concatenate, int8, where, array, sqrt, square, ones, power, sum, mean, std, linspace, max, round, median, polyfit, vstack, random, greater, less, searchsorted, sort
from numpy              import float as npfloat, log as nplog, round as np_round, around, float32, invert
from numpy.linalg       import lstsq
from pymc               import deterministic, stochastic, Normal, Uniform, MCMC, Bernoulli, stochastic_from_dist
from scipy.integrate    import simps
from scipy.interpolate  import interp1d
from scipy.optimize     import minimize, curve_fit, leastsq
from scipy.signal       import argrelextrema
from uncertainties      import ufloat
from bces_script        import bces
# from kmpfit_script      import kmpfit_effectivevariance, scipy_ODR
from linfit_script      import linfit
from lnr_script         import kelly
from lmfit              import Parameters, minimize as lmfit_minimize, report_fit,  Minimizer, fit_report
from pyneb              import Atom
from lmfit.models import GaussianModel

from scipy.integrate    import simps
from dazer_methods import Dazer
from numpy import median, random, log as np_ln

class Fitting_Gaussians():
    
    '''This class contains the methods to measure gaussians in a scheme similar to IRAF'''
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
        print self.N2_Ratio
        
        self.s2pi = sqrt(2*pi)

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
#             y_model = y_model + A * exp(-(x-mu)*(x-mu)/(2*sigma*sigma))
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
        
    def Load_lmfit_parameters(self, x, y, zerolev, err_zerolev, n_comps, wide_component = False, A_limits = 0.30, mu_precission = 2, sigma_limit = 5):
        
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
            params.add('fwhm'  + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
            params.add('FluxG' + index, expr = '{A} * {sigma} * {sqrt2pi}'.format(A = 'A'  + index, sigma = 'sigma' + index), sqrt2pi = self.s2pi)
            params.add('FluxG' + index, expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + index, fwhm = 'fwhm' + index))
            
        #For blended components we set the same wavelength:
        if n_comps > 1:
            
            print self.Fitting_dict['blended wavelengths']
            
            Highest_index       = argmax(self.Fitting_dict['peak_Maxima'])
            small_components    = range(n_comps)    
            del small_components[Highest_index]
            
            for indx in small_components:
                
                #We set the same sigma
                index_small = str(indx)
                expresion = 'sigma{index_big} * ( (mu{index_small} + {scaller}) / (mu{index_big} + {scaller}) )'.format(index_big = Highest_index, index_small = index_small, scaller = self.Fitting_dict['x_scaler'])
                params['sigma' + index_small].set(expr = expresion) 
            
#                 #We force the theoretical - biggest mu
#                 expresion = '{mu_small} - mu{index_big}'.format(mu_small = self.Fitting_dict['blended wavelengths'][indx] - self.Fitting_dict['x_scaler'], index_big = Highest_index)
#                 params['mu' + index_small].set(expr = expresion) 
          
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
            params_W.add('fwhm' + n_nindex,     expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + n_nindex))
#             params_W.add('FluxG'+ n_nindex,     expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + n_nindex, fwhm = 'fwhm' + n_nindex))
            params_W.add('FluxG' + n_nindex, expr = '{A} * {sigma} * (2*3.1415)**0.5'.format(A = 'A'  + n_nindex, sigma = 'sigma' + n_nindex))
            Wide_params_list = params_W.keys()
            
            sqrt(2)

            #Update for Nitrogen relation
#             params['FluxG0'].set(expr = 'FluxG2 / {N2_ratio}'.format(N2_ratio = 2.94))
            expression = '(A2*sigma2) / ({N2_ratio}*sigma0) '.format(N2_ratio = 2.94)
            params['A0'].set(expr = expression) 

        #Store the data
        self.Fitting_dict['lmfit_params']       = params
        self.Fitting_dict['lmfit_params_wide']  = params_W
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
                
            minima_wave, minima_flux =  self.Fitting_dict['x_norm'][min_index], self.Fitting_dict['y_norm'][min_index]
        
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
                sigma_limit = fit_Output.params['sigma1'].value
                limit_0 = 6548.05 - self.Fitting_dict['x_scaler']  - sigma_limit * 1.5
                limit_1 = 6548.05 - self.Fitting_dict['x_scaler']  + sigma_limit * 1.5
                limit_2 = 0  - sigma_limit * 4
                limit_3 = 0  + sigma_limit * 4 
                limit_4 = 6583.46 - self.Fitting_dict['x_scaler']  - sigma_limit * 3
                limit_5 = 6583.46 - self.Fitting_dict['x_scaler']  + sigma_limit * 3
                
                indeces = ((x >= limit_0) & (x <= limit_1)) + ((x >= limit_2) & (x <= limit_3)) + ((x >= limit_4) & (x <= limit_5))    
                mask    = invert(indeces)
                self.Fitting_dict['wide mask'] = mask
                
                x_wide      = x[mask]
                y_wide      = y_new[mask]
                zero_wide   = zero_lev[mask]
                Ncomps_wide = ['3']
                
#                 y_wide[(y_wide < (mean(self.Fitting_dict['zerolev_norm']) - self.Fitting_dict['sig_zerolev_norm']))] = mean(self.Fitting_dict['zerolev_norm'])
                
                fit_Output_wide     = lmfit_minimize(self.lmfit_gaussian_Residual_wide, fitting_parameters_wide, args=(x_wide, y_wide, zero_wide, Ncomps_wide, err_continuum))

                #Substract the wide components
                y_wide              = self.gaussian_curve_SingleMixture_wide(fit_Output_wide.params.valuesdict(), x, zero_lev, Ncomps_wide)

                y_pure_emission     = y_new - y_wide + zero_lev 
                self.Fitting_dict['emis_limpio'] = y_pure_emission 
                
                fit_Output_emission = lmfit_minimize(self.lmfit_gaussian_Residual, fitting_parameters, args=(x, y_pure_emission, zero_lev, Ncomps, err_continuum))

                output_params = fit_Output_emission.params + fit_Output_wide.params
                
                print fit_report(output_params)
                print 'flux relation', output_params['FluxG2'] /  output_params['FluxG0']
                print 'A2*sigma2 / (A0*sigma0)', (output_params['A2'] * output_params['sigma2']) /  (output_params['A0'] * output_params['sigma0'])
                
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
                self.Fitting_dict['fwhm' + index]   = ufloat(self.Fitting_dict['fwhm' + index + '_norm'],   self.Fitting_dict['fwhm' + index + '_norm_er'])
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
                self.Fitting_dict['fwhm' + index]   = ufloat(mean(self.Fitting_dict['fwhm' + index + '_norm']),     std(self.Fitting_dict['fwhm' + index + '_norm']))
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

dz = Dazer()
fg = Fitting_Gaussians()

#Define operation
Catalogue_Dic   = dz.import_catalogue()
Pattern         = 'SHOC579_WHT.fits'

#Locate files on hard drive
FilesList = dz.Folder_Explorer(Pattern, Catalogue_Dic['Obj_Folder'], CheckComputer=False)

#Generate plot frame and colors
dz.FigConf(n_colors=10)

for i in range(len(FilesList)):

    CodeName, FileName, FileFolder  = dz.Analyze_Address(FilesList[i])
    Wave, Flux, ExtraData           = dz.File_to_data(FileFolder, FileName)

    Line_label, line_wavelength     = 'N2_6548A', 6548.05
    Line_selections                 = dz.RangesR([], line_wavelength, FileFolder + CodeName + '_WHT_LinesLog_v3.txt')