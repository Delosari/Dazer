'''
Created on Jul 23, 2015

@author: vital
'''
from sys                            import exit
from pyneb                          import RecAtom
from scipy.interpolate              import interp1d
from uncertainties                  import ufloat, UFloat
from uncertainties.unumpy           import uarray, exp as unum_exp, log10 as unum_log10, pow as unum_pow
from uncertainties.umath            import pow as umath_pow, log10 as umath_log10, exp as umath_exp, isnan as un_isnan
from pandas                         import read_csv
from numpy                          import loadtxt, where, array, ndarray, ones, concatenate, power as np_power, dot, arange, empty, isnan

class ReddeningLaws():

    def __init__(self):
                        
        self.SpectraEdges_Limit = 200 
        self.Hbeta_wavelength   = 4861.3316598713955

        #List of hydrogen recombination lines in our EmissionLines coordinates log. WARNING: We have removed the very blended ones
        self.RecombRatios_Ions =  array(['HBal_20_2', 'HBal_19_2','HBal_18_2', 'HBal_17_2', 'HBal_16_2','HBal_15_2','HBal_12_2','HBal_11_2',
                             'HBal_9_2','Hdelta_6_2','Hgamma_5_2','Hbeta_4_2','Halpha_3_2','HPas_20_3','HPas19_3','HPas_18_3','HPas_17_3','HPas_16_3','HPas_15_3','HPas_14_3',
                             'HPas13_3','HPas12_3','HPas_11_3','HPas_10_3','HPas_9_3','HPas_8_3','HPas_7_3'])
        
        #Dictionary with the reddening curves
        self.reddening_curves_calc = {'MM72'            : self.f_Miller_Mathews1972,
                                      'CCM89'           : self.X_x_Cardelli1989,
                                      'G03_bar'         : self.X_x_Gordon2003_bar,
                                      'G03_average'     : self.X_x_Gordon2003_average,
                                      'G03_supershell'  : self.X_x_Gordon2003_supershell
                                      }
    
    def checking_for_ufloat(self, x):
                
        if isinstance(x, UFloat): 
            param = x.nominal_value
               
        else:
            param = x
        
        return param
    
    def compare_RecombCoeffs(self, obj_data, lineslog_frame, spectral_limit = 200):
        
        #Create hidrogen atom object 
        self.H1_atom = RecAtom('H', 1)

        linformat_df = read_csv('/home/vital/workspace/dazer/format/emlines_pyneb_optical_infrared.dz', index_col=0, names=['ion', 'lambda_theo', 'latex_format'], delim_whitespace=True)
        lineslog_frame['latex_format'] = 'none'
        
        for line in lineslog_frame.index:
            if '_w' not in line:    #Structure to avoid wide components
                lineslog_frame.loc[line,'latex_format'] = r'${}$'.format(linformat_df.loc[line, 'latex_format'])
        
        #Load electron temperature and density (if not available it will use Te = 10000K and ne = 100cm^-3)             
        T_e = self.checking_for_ufloat(obj_data.TeSIII) if ~isnan(self.checking_for_ufloat(obj_data.TeSIII)) else 10000.0
        n_e = self.checking_for_ufloat(obj_data.neSII) if ~isnan(self.checking_for_ufloat(obj_data.neSII)) else 100.0
                
        #Get the Hidrogen recombination lines that we have observed in this object
        Obs_Hindx = lineslog_frame.Ion.isin(self.RecombRatios_Ions)
                            
        #Calculate recombination coefficients and reddening curve values
        Obs_H_Emis  = empty(len(Obs_Hindx))
        Obs_Ions    = lineslog_frame.loc[Obs_Hindx, 'Ion'].values
        
        for i in range(len(Obs_Ions)):
            TransitionCode = Obs_Ions[i][Obs_Ions[i].find('_')+1:len(Obs_Ions[i])]
            Obs_H_Emis[i]  = self.H1_atom.getEmissivity(tem = T_e, den = n_e, label = TransitionCode)        
        
        #Normalize by Hbeta (this new constant is necessary to avoid a zero sigma)
        Hbeta_Emis = self.H1_atom.getEmissivity(tem = T_e, den = n_e, label = '4_2')        
        Fbeta_flux = ufloat(lineslog_frame.loc['H1_4861A']['line_Flux'].nominal_value, lineslog_frame.loc['H1_4861A']['line_Flux'].std_dev)
        
        #Load theoretical and observational recombination ratios to data frame
        lineslog_frame.loc[Obs_Hindx,'line_Emissivity'] = Obs_H_Emis
        lineslog_frame.loc[Obs_Hindx,'line_TheoRecombRatio'] = Obs_H_Emis / Hbeta_Emis
        lineslog_frame.loc[Obs_Hindx,'line_ObsRecombRatio'] = lineslog_frame.loc[Obs_Hindx, 'line_Flux'].values / Fbeta_flux
        
        #Get indeces of emissions in each arm
        ObsBlue_Hindx   = Obs_Hindx & (lineslog_frame.lambda_theo < obj_data.join_wavelength)
        ObsRed_Hindx    = Obs_Hindx & (lineslog_frame.lambda_theo > obj_data.join_wavelength)
      
        #Recalculate red arm coefficients so they are normalized by red arm line (the most intense)
        if (ObsRed_Hindx.sum()) > 0: #Can only work if there is at least one line
            idx_Redmax      = lineslog_frame.loc[ObsRed_Hindx]['flux_intg'].idxmax()  #We do not use line_flux because it is a ufloat and it does not work with idxmax 
            H_Redmax_flux   = lineslog_frame.loc[idx_Redmax,'line_Flux']
            H_Redmax_emis   = lineslog_frame.loc[idx_Redmax,'line_Emissivity']
            Flux_Redmax     = ufloat(H_Redmax_flux.nominal_value * Hbeta_Emis / H_Redmax_emis, H_Redmax_flux.std_dev * Hbeta_Emis / H_Redmax_emis)
            lineslog_frame.loc[ObsRed_Hindx,'line_ObsRecombRatio'] = lineslog_frame.loc[ObsRed_Hindx, 'line_Flux'].values / Flux_Redmax
                    
        #Load x axis values: (f_lambda - f_Hbeta) and y axis values: log(F/Fbeta)_theo - log(F/Fbeta)_obs to dataframe
        lineslog_frame.loc[Obs_Hindx,'x axis values'] = lineslog_frame.loc[Obs_Hindx, 'line_f'].values - lineslog_frame.loc['H1_4861A']['line_f']
        lineslog_frame.loc[Obs_Hindx,'y axis values'] = unum_log10(lineslog_frame.loc[Obs_Hindx,'line_TheoRecombRatio'].values) - unum_log10(lineslog_frame.loc[Obs_Hindx,'line_ObsRecombRatio'].values)
            
        #Compute all the possible configuration of points and store them
        output_dict = {}
        
        #--- All points:
        output_dict['all_x']    = lineslog_frame.loc[Obs_Hindx,'x axis values'].values
        output_dict['all_y']    = lineslog_frame.loc[Obs_Hindx,'y axis values'].values
        output_dict['all_ions'] = list(lineslog_frame.loc[Obs_Hindx,'latex_format'].values)
        
        #--- By arm
        output_dict['blue_x']   = lineslog_frame.loc[ObsBlue_Hindx,'x axis values'].values
        output_dict['blue_y']   = lineslog_frame.loc[ObsBlue_Hindx,'y axis values'].values
        output_dict['blue_ions'] = list(lineslog_frame.loc[ObsBlue_Hindx,'latex_format'].values)        
        output_dict['red_x']    = lineslog_frame.loc[ObsRed_Hindx,'x axis values'].values
        output_dict['red_y']    = lineslog_frame.loc[ObsRed_Hindx,'y axis values'].values
        output_dict['red_ions'] = list(lineslog_frame.loc[ObsRed_Hindx,'latex_format'].values)
        
        #--- Store fluxes
        output_dict['Blue_ObsRatio'] = unum_log10(lineslog_frame.loc[ObsBlue_Hindx,'line_ObsRecombRatio'].values)
        output_dict['Red_ObsRatio'] = unum_log10(lineslog_frame.loc[ObsRed_Hindx,'line_ObsRecombRatio'].values)        
        output_dict['line_Flux_Blue'] = lineslog_frame.loc[ObsBlue_Hindx,'line_Flux'].values
        output_dict['line_Flux_Red'] = lineslog_frame.loc[ObsRed_Hindx,'line_Flux'].values
        output_dict['line_wave_Blue'] = lineslog_frame.loc[ObsBlue_Hindx,'lambda_theo'].values
        output_dict['line_wave_Red'] = lineslog_frame.loc[ObsRed_Hindx,'lambda_theo'].values
        
        #--- Inside limits
        if obj_data.h_gamma_valid == 'yes':
            in_idcs = Obs_Hindx
        elif obj_data.h_gamma_valid == 'no':
            wave_idx = ((lineslog_frame.lambda_theo > (obj_data.Wmin_Blue + spectral_limit)) & (lineslog_frame.lambda_theo < (obj_data.join_wavelength - spectral_limit))) \
                    | ((lineslog_frame.lambda_theo > (obj_data.join_wavelength + spectral_limit)) & (lineslog_frame.lambda_theo < (obj_data.Wmax_Red - spectral_limit)))
            in_idcs = Obs_Hindx & wave_idx
           
        output_dict['in_x'] = lineslog_frame.loc[in_idcs,'x axis values'].values
        output_dict['in_y'] = lineslog_frame.loc[in_idcs,'y axis values'].values
        output_dict['in_ions'] = list(lineslog_frame.loc[in_idcs,'latex_format'].values)

        #--- Outside limis
        if obj_data.h_gamma_valid == 'no':
            wave_idx = (lineslog_frame.lambda_theo < (obj_data.Wmin_Blue + spectral_limit)) | (lineslog_frame.lambda_theo > (obj_data.Wmax_Red - spectral_limit))
            out_idcs = Obs_Hindx & wave_idx
            output_dict['out_x'] = lineslog_frame.loc[out_idcs,'x axis values'].values
            output_dict['out_y'] = lineslog_frame.loc[out_idcs,'y axis values'].values
            output_dict['out_ions'] = list(lineslog_frame.loc[out_idcs,'latex_format'].values)
        else:
            output_dict['out_x'] = None
            output_dict['out_y'] = None
            output_dict['out_ions'] = None
           
        return output_dict

    def deredden_lines(self, lines_frame, reddening_curve, cHbeta = None, E_BV = None, R_v = None):
        
        #cHbeta format
        if isinstance(cHbeta, float):
            cHbeta_mag = cHbeta
            
        elif isinstance(cHbeta, UFloat):         
            cHbeta_mag = cHbeta

        #If it is negative we set it to zero
        if cHbeta_mag < 0.0:
            cHbeta_mag = 0.0
        
        #By default we perform the calculation using the colour excess
        E_BV = E_BV if E_BV != None else self.Ebv_from_cHbeta(cHbeta, reddening_curve, R_v)    
                      
        #Get reddening curves f_value
        lines_fluxes                = lines_frame.line_Flux.values
        lines_wavelengths           = lines_frame.lambda_theo.values
        lines_Xx                    = self.reddening_Xx(lines_wavelengths, reddening_curve, R_v)
        lines_frame['line_Xx']      = lines_Xx
        
        #Get line intensities 
        lines_int                   = lines_fluxes * unum_pow(10,  0.4 * lines_Xx * E_BV)
        lines_frame['line_Int']     = lines_int
        
        #We recalculate the equivalent width using the intensity of the lines (instead of the flux)
        continua_flux               = uarray(lines_frame.zerolev_mean.values, lines_frame.zerolev_std.values)
        continua_int                = continua_flux * unum_pow(10,  0.4 * lines_Xx * E_BV)
        lines_frame['con_dered']    = continua_int
        lines_frame['Eqw_dered']    = lines_int / continua_int
        
        #For extra testing we add integrated and gaussian values derredden
        lines_brut_flux                     = uarray(lines_frame['flux_intg'].values,  lines_frame['flux_intg_er'].values)
        lines_gauss_flux                    = uarray(lines_frame['flux_gauss'].values, lines_frame['flux_gauss_er'].values)
        lines_frame['line_IntBrute_dered']  = lines_brut_flux * unum_pow(10, 0.4 * lines_Xx * E_BV)
        lines_frame['line_IntGauss_dered']  = lines_gauss_flux * unum_pow(10, 0.4 * lines_Xx * E_BV)
                
        return
    
    def derreddening_spectrum(self, wave, flux, reddening_curve, cHbeta = None, E_BV = None, R_v = None):

        #cHbeta format
        if isinstance(cHbeta, float):
            cHbeta_mag = cHbeta
            
        elif isinstance(cHbeta, UFloat):         
            cHbeta_mag = cHbeta

        #If it is negative we set it to zero
        if cHbeta_mag < 0.0:
            cHbeta_mag = 0.0

        #By default we perform the calculation using the colour excess
        E_BV = E_BV if E_BV != None else self.Ebv_from_cHbeta(cHbeta, reddening_curve, R_v)
        
        #Perform dereddening
        wavelength_range_Xx = self.reddening_Xx(wave, reddening_curve, R_v)
        flux_range_derred   = flux * np_power(10,  0.4 * wavelength_range_Xx * E_BV)
                               
        return flux_range_derred

    def reddening_spectrum(self, wave, flux, reddening_curve, cHbeta = None, E_BV = None, R_v = None):

        #cHbeta format
        if isinstance(cHbeta, float):
            cHbeta_mag = cHbeta
            
        elif isinstance(cHbeta, UFloat):         
            cHbeta_mag = cHbeta

        #If it is negative we set it to zero
        if cHbeta_mag < 0.0:
            cHbeta_mag = 0.0

        #By default we perform the calculation using the colour excess
        E_BV = E_BV if E_BV != None else self.Ebv_from_cHbeta(cHbeta, reddening_curve, R_v)
        
        #Perform dereddening
        wavelength_range_Xx = self.reddening_Xx(wave, reddening_curve, R_v)
        flux_range_red  = flux * np_power(10,  - 0.4 * wavelength_range_Xx * E_BV)

        return flux_range_red

       
#         if (cHbeta != None) and (cHbeta > 0):
#             f_wave              = self.Reddening_f(Spectrum_WavelengthRange, 3.2)
#             f_Hbeta             = self.Reddening_f(array([4862.683]), 3.2)
#             
#             Power_Coeff         = cHbeta * (1 + f_wave - f_Hbeta)
#             
#             Spectrum_flux       = Spectrum_Int * np_power(10, -Power_Coeff)
#         
#         else:
#             print '- WARNING: c(Hbeta) is less than zero'
#             Spectrum_flux        = Spectrum_Int
                       
    def Ebv_from_cHbeta(self, cHbeta, reddening_curve, R_v):

        if cHbeta == None:
            exit('Warning: no cHbeta or E(B-V) provided to reddening curve, code aborted')
        else:
            if cHbeta != None:
                E_BV = cHbeta * 2.5 / self.reddening_Xx(array([self.Hbeta_wavelength]), reddening_curve, R_v)[0]
                return E_BV
 
    def flambda_from_Xx(self, Xx, reddening_curve, R_v):
        
        X_Hbeta = self.reddening_Xx(array([self.Hbeta_wavelength]), reddening_curve, R_v)[0]
        
        f_lines = Xx/X_Hbeta - 1
        
        return f_lines
 
    def reddening_Xx(self, waves, curve_methodology, R_v):
        
        self.R_v = R_v
        self.wavelength_rc = waves
        return self.reddening_curves_calc[curve_methodology]()
        
    def f_Miller_Mathews1972(self):
                    
        if isinstance(self.wavelength_rc, ndarray):
            y           = 1.0/(self.wavelength_rc/10000.0) 
            y_beta      = 1.0/(4862.683/10000.0)
            
            ind_low     = where(y <= 2.29)[0]
            ind_high    = where(y >  2.29)[0]
            
            dm_lam_low  = 0.74 * y[ind_low]  - 0.34 + 0.341 * self.R_v - 1.014
            dm_lam_high = 0.43 * y[ind_high] + 0.37 + 0.341 * self.R_v - 1.014
            dm_beta     = 0.74 * y_beta      - 0.34 + 0.341 * self.R_v - 1.014
            
            dm_lam = concatenate((dm_lam_low, dm_lam_high))
            
            f = dm_lam/dm_beta - 1
            
        else:
            
            y           = 1.0/(self.wavelength_rc/10000.0) 
            y_beta      = 1.0/(4862.683/10000.0)
            
            if y <= 2.29:
                dm_lam  = 0.74 * y  - 0.34 + 0.341 * self.R_v - 1.014
            else:
                dm_lam = 0.43 * y + 0.37 + 0.341 * self.R_v - 1.014
    
            dm_beta     = 0.74 * y_beta - 0.34 + 0.341 * self.R_v - 1.014
                            
            f = dm_lam/dm_beta - 1
            
        return f
        
    def X_x_Cardelli1989(self):
                        
        x_true = 1.0 / (self.wavelength_rc / 10000.0)
        y = x_true - 1.82
    
        y_coeffs = array([ones(len(y)), y, np_power(y, 2), np_power(y, 3), np_power(y, 4), np_power(y, 5), np_power(y, 6), np_power(y, 7)])
        a_coeffs = array([1, 0.17699,    -0.50447,   -0.02427,   0.72085,    0.01979,    -0.77530,   0.32999])
        b_coeffs = array([0, 1.41338,    2.28305,   1.07233,   -5.38434,    -0.62251,    5.30260,   -2.09002])
        
        a_x =  dot(a_coeffs,y_coeffs)
        b_x =  dot(b_coeffs,y_coeffs)
        
        X_x = a_x + b_x / self.R_v
                        
        return X_x
    
    def X_x_Gordon2003_bar(self):
        
        #Default R_V is 3.4
        R_v = self.R_v if self.R_v != None else 3.4 #This is not very nice   
        x = 1.0 / (self.wavelength_rc / 10000.0)
        
        #This file format has 1/um in column 0 and A_x/A_V in column 1
        file_data = loadtxt('/home/vital/workspace/dazer/bin/lib/Astro_Libraries/gordon_2003_SMC_bar.txt')
        
        #This file has column        
        Xx_interpolator = interp1d(file_data[:, 0], file_data[:, 1])
        X_x = R_v * Xx_interpolator(x)
        return X_x

    def X_x_Gordon2003_average(self):
        
        #Default R_V is 3.4
        R_v = self.R_v if self.R_v != None else 3.4 #This is not very nice   
        x = 1.0 / (self.wavelength_rc / 10000.0)
        
        #This file format has 1/um in column 0 and A_x/A_V in column 1
        file_data = loadtxt('/home/vital/workspace/dazer/bin/lib/Astro_Libraries/gordon_2003_LMC_average.txt')

        #This file has column        
        Xx_interpolator = interp1d(file_data[:, 0], file_data[:, 1])
        X_x = R_v * Xx_interpolator(x)
        return X_x

    def X_x_Gordon2003_supershell(self):
        
        #Default R_V is 3.4
        R_v = self.R_v if self.R_v != None else 3.4 #This is not very nice   
        x = 1.0 / (self.wavelength_rc / 10000.0)
        
        #This file format has 1/um in column 0 and A_x/A_V in column 1
        file_data = loadtxt('/home/vital/workspace/dazer/bin/lib/Astro_Libraries/gordon_2003_LMC2_supershell.txt')
        
        #This file has column        
        Xx_interpolator = interp1d(file_data[:, 0], file_data[:, 1])
        X_x = R_v * Xx_interpolator(x)
        return X_x

    def Epm_ReddeningPoints(self):
         
        x_true      = arange(1.0, 2.8, 0.1) #in microns -1
        X_Angs      = 1 / x_true * 1e4
         
        Xx          = array([1.36, 1.44, 1.84, 2.04, 2.24, 2.44, 2.66, 2.88, 3.14, 3.36, 3.56, 3.77, 3.96, 4.15, 4.26, 4.40, 4.52, 4.64])
        f_lambda    = array([-0.63,-0.61,-0.5, -0.45, -0.39, -0.34, -0.28, -0.22, -0.15, -0.09, -0.03, 0.02, 0.08, 0.13, 0.16, 0.20, 0.23, 0.26])
         
        return x_true, X_Angs, Xx, f_lambda
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        