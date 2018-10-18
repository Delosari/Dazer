from scipy.stats    import truncnorm
from uncertainties  import ufloat, umath, unumpy, UFloat
from pyneb          import RecAtom, Atom, atomicData, Diagnostics
from numpy          import ndarray, zeros, mean, median, where, std, log10, power, random, exp, loadtxt, nan as np_nan, isnan, savetxt, transpose, copy, vstack, nanmean, nanstd, empty, all as np_all, sum as np_sum
from collections    import OrderedDict
from pandas         import Series, isnull, notnull
from lmfit          import Parameters, minimize as lmfit_minimmize, report_fit, models
import numexpr
from timeit import default_timer as timer

def residual_Y_v3(pars, x, y):
    return (y - pars['Y'] * x)

class Chemical_Analysis_pyneb():

    def __init__(self):
        
        self.MC_array_len = 1000
        self.MC_warning_limit = self.MC_array_len * 0.1
        
        self.Hbeta_label  = 'H1_4861A'
    
    def load_elements(self):

        #Set atomic data 
        #atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
        atomicData.setDataFile('s_iii_coll_HRS12.dat')
 
        #Default: 's_iii_atom_PKW09.dat' 
        'S3: All energy and A values: Podobedova, Kelleher, and Wiese 2009, J. Phys. Chem. Ref. Data, Vol.'
        'S3: collision strengths: Tayal & Gupta 1999, ApJ, 526, 544'

        #New Atomic data s_iii_coll_HRS12.dat
        'S3: All energy and A values: Podobedova, Kelleher, and Wiese 2009, J. Phys. Chem. Ref. Data, Vol.'
        'S3: collision strengths: Hudson, Ramsbottom & Scott 2012, ApJ, 750, 65'
        
        #Declare ions
        self.S2_atom            = Atom('S', 2)
        self.S3_atom            = Atom('S', 3)
        self.Ar3_atom           = Atom('Ar', 3)
        self.Ar4_atom           = Atom('Ar', 4)
        self.N2_atom            = Atom('N', 2)
        self.O2_atom            = Atom('O', 2)
        self.O3_atom            = Atom('O', 3)
        self.H1_atom            = RecAtom('H', 1)
        self.He1_atom           = RecAtom('He', 1)
        self.He2_atom           = RecAtom('He', 2)
        
        #Pyneb objects
        self.diags              = Diagnostics()
        
        #Ohrs 2016 relation for the OI_SI gradient
        self.logSI_OI_Gradient  = random.normal(-1.53,  0.05, size = self.MC_array_len) # random.normal(-1.78,  0.03, size = self.MC_array_len)
        self.OI_SI              = power(10, -self.logSI_OI_Gradient)
        
        #Theoretical ratios
        self.S3_ratio           = self.S3_atom.getEmissivity(10000, 100, wave = 9531) / self.S3_atom.getEmissivity(10000, 100, wave = 9069)
        self.S3_9000_ratio      = random.normal(self.S3_atom.getEmissivity(10000, 100, wave = 9531) / self.S3_atom.getEmissivity(10000, 100, wave = 9069),  0.01, size = self.MC_array_len)
        self.N2_6000_ratio      = self.N2_atom.getEmissivity(10000, 100, wave = 6584) / self.N2_atom.getEmissivity(10000, 100, wave = 6548)
        self.O3_5000_ratio      = self.O3_atom.getEmissivity(10000, 100, wave = 5007) / self.O3_atom.getEmissivity(10000, 100, wave = 4959)
        
        #Factors to speed calculations
        self.lines_factors = {}
        self.lines_factors['S3_9069A'] = 1 + self.S3_ratio
        self.lines_factors['S3_9531A'] = 1 + 1/self.S3_ratio
        
        #Cloudy models for the SIV contribution
        
        self.m_SIV_correction   = random.normal(1.1628,  0.00559, size = self.MC_array_len)
        self.n_SIV_correction   = random.normal(0.0470,  0.0097, size = self.MC_array_len)
        #self.m_SIV_correction   = random.normal(1.109,  0.01, size = self.MC_array_len)
        #self.n_SIV_correction   = random.normal(0.135,  0.0173, size = self.MC_array_len)
 
        #CHAOS relation TNII-TSIII
        #T[SIII]  = 1.312(+-0.075)T[NII]-0.313(+-0.058)
                #TNII     = (0.762+-0.044)*TSIII  + 0.239+-0.046
        self.m_TNII_correction   = random.normal(0.762,  0.044, size = self.MC_array_len)
        self.n_TNII_correction   = random.normal(0.239,  0.046, size = self.MC_array_len)       
        
        
        #Truncated gaussian for the density
        lower_trunc, upper_trunc = (1.0 - 50.0) / 25.0, (100 - 50) / 25.0
        self.Truncated_gaussian  = truncnorm(lower_trunc,  upper_trunc, loc = 50, scale = 25)
       
        print '-Elements loaded\n'
       
        return
            
    def declare_object(self, lines_log_frame):
        
        #List of all parameters
#         lineRatios      = ['R_SII', 'R_SII_prime', 'R_SIII', 'R_NII', 'R_OII', 'R_OII_prime', 'R_OIII']
#         elecProperties  = ['neSII', 'neOII', 'TeOII', 'TeSII', 'TeNII', 'TeOIII', 'TeSIII', 'TeOII_from_TeOIII', 'TeNII_from_TeOIII', 'TeSIII_from_TeOIII', 'TeOIII_from_TeSIII']        
#         ionicAbund      = ['SII_HII', 'SIII_HII', 'SIV_HII', 'OII_HII', 'OII_HII_3279A', 'OII_HII_7319A', 'NII_HII', 'ArIII_HII', 'ArIV_HII', 'HeII_HII_from_O',
#                             'HeIII_HII_from_O', 'HeII_HII_from_S', 'HeIII_HII_from_S']
#         elemAbund       = ['SI_HI', 'OI_HI', 'NI_OI', 'NI_HI', 'HeI_HI_from_O', 'HeI_HI_from_S', 'Ymass_O', 'Ymass_S']
    
        self.abunData = Series()
        
        self.Hbeta_flux = random.normal(lines_log_frame.loc['H1_4861A', 'line_Int'].nominal_value, lines_log_frame.loc['H1_4861A', 'line_Int'].std_dev, size = self.MC_array_len)
                
        self.low_density_dist = self.Truncated_gaussian.rvs(self.MC_array_len)

        #Generate a dictionary to store the random array for all lines
        self.lines_dict = OrderedDict()
        
        #Dictionary with lines which my need special treatements
        Blended_lines = {}
        Blended_lines['O2_3726A'] = ('O2_3726A', 'O2_3729A')
        Blended_lines['O2_7319A'] = ('O2_7319A', 'O2_7330A')
        NoError_lines = {}
        NoError_lines['N2_6548A'] = ('N2_6548A')

        #Generate
        for line in lines_log_frame.index.values:
            
            #Start with the particular cases: Blended lines
            if line in Blended_lines:
                blended_lines = Blended_lines[line]

                if set(lines_log_frame.index) >= set(blended_lines):
                    label_line = line + '+'
                                        
                    #Lines are blended we use integrated flux else we add the individual integrated
                    if lines_log_frame.loc[blended_lines[0], 'flux_intg'] == lines_log_frame.loc[blended_lines[1]]['flux_intg']:
                        flux_line   = lines_log_frame.loc[blended_lines[0], 'line_IntBrute_dered'].nominal_value
                        error_line  = lines_log_frame.loc[blended_lines[0], 'line_IntBrute_dered'].std_dev
           
                    else:
                        line_sum    = lines_log_frame.loc[blended_lines[0], 'line_IntBrute_dered'] + lines_log_frame.loc[blended_lines[1], 'line_IntBrute_dered']
                        flux_line   = line_sum.nominal_value
                        error_line  = line_sum.std_dev
                                                        
                #Case only one of the lines was measured
                else:
                    label_line  = line
                    flux_line   = lines_log_frame.loc[line, 'line_Int'].nominal_value
                    error_line  = lines_log_frame.loc[line, 'line_Int'].std_dev                                           
            
            #Lines with not error
            elif (line in NoError_lines) and (lines_log_frame.loc[line, 'line_Int'].std_dev == 0.0):
                label_line  = line
                flux_line   = lines_log_frame.loc[line, 'line_Int'].nominal_value 
                error_line  = lines_log_frame.loc['N2_6584A', 'line_Int'].std_dev / self.N2_6000_ratio
            
            #None-issue lines
            else:
                label_line  = line
                flux_line   = lines_log_frame.loc[line, 'line_Int'].nominal_value
                error_line  = lines_log_frame.loc[line, 'line_Int'].std_dev 

            #Generate line gaussian shaped array
            self.lines_dict[label_line] = random.normal(flux_line, error_line, size = self.MC_array_len)
                                
        return
                               
    def den_temp_diagnostic_pair(self,  diagProperties, den_distribution = None, atom_temp = None):
        
        #Check if all necessary lines are there
        if self.lines_dict.viewkeys() >= set(diagProperties['required_denlines'] + diagProperties['required_temlines']):
            
            if den_distribution is None:
                den_ratio = numexpr.evaluate(diagProperties['den_ratio'], self.lines_dict)
                tem_ratio = numexpr.evaluate(diagProperties['tem_ratio'], self.lines_dict)
                
                Te, ne = self.diags.getCrossTemDen(diag_tem = diagProperties['diag_tem'], diag_den = diagProperties['diag_den'], 
                                                   value_tem = tem_ratio, value_den = den_ratio) 
                
            else:
                tem_ratio = numexpr.evaluate(diagProperties['tem_ratio'], self.lines_dict)
                Te  = atom_temp.getTemDen(tem_ratio, den = den_distribution, to_eval = diagProperties['atom_temdiag'])
                ne  = den_distribution
                
        #else empty (nan) arrays
        else:
            Te, ne          = empty(self.MC_array_len), empty(self.MC_array_len)
            Te[:], ne[:]    = np_nan, np_nan
                
        return Te, ne

    def determine_electron_parameters(self, objectData):
        
        #----------To start make sure we are not in the very low density regimes, 
        low_density_dist = None
        if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:
            S2_ratio = mean(self.lines_dict['S2_6716A']) / mean(self.lines_dict['S2_6731A'])
            if S2_ratio > 1.35:
                print '--Low density object'
                lower, upper, mu, sigma = 1.0, 100.0, 50.0, 25.0
                X_func                  = truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)    
                low_density_dist        = X_func.rvs(self.MC_array_len)
                self.abunData['neSII']  = low_density_dist
                if low_density_dist is None:
                    print 'WARNING: QUE PASA!!!!!!'
                    print lower, upper, mu, sigma
                    print xrange
                
        #-----------Sulfur
        diagProperties = {}
        diagProperties['required_denlines'] = ['S2_6716A', 'S2_6731A']
        diagProperties['required_temlines'] = ['S3_9069A', 'S3_9531A', 'S3_6312A'] if objectData.SIII_lines == 'BOTH' else [objectData.SIII_lines] + ['S3_6312A']
        diagProperties['diag_den']          = '[SII] 6731/6716'
        diagProperties['diag_tem']          = '[SIII] 6312/9200+' 
        diagProperties['atom_temdiag']      = 'L(6312)/(L(9069)+L(9531))'
        diagProperties['den_ratio']         = 'S2_6731A/S2_6716A'
        diagProperties['tem_ratio']         = 'S3_6312A/(S3_9069A+S3_9531A)' if objectData.SIII_lines == 'BOTH' else 'S3_6312A/({valid_line} * {line_factor})'.format(valid_line = objectData.SIII_lines, line_factor=self.lines_factors[objectData.SIII_lines])
        
        if '*' in diagProperties['tem_ratio']:
            print '-- Using factor', diagProperties['tem_ratio']
                
        S3_lines = ['S3_9069A', 'S3_9531A', 'S3_6312A'] if objectData.SIII_lines == 'BOTH' else [objectData.SIII_lines] + ['S3_6312A']

        #--Calculate NeSII and TeSIII
        self.abunData['TeSIII'],  neSII_TSIII = self.den_temp_diagnostic_pair(diagProperties, low_density_dist, atom_temp=self.S3_atom)
        
        #--Determine empirical TOIII from TSIII  #Epm & Diaz 2005
        self.abunData['TeOIII_from_TeSIII'] = (0.95 * self.abunData.TeSIII/10000 + 0.076) * 10000
         
        diagProperties = {}
        diagProperties['required_denlines'] = ['S2_6716A', 'S2_6731A']
        diagProperties['required_temlines'] = ['S2_4069A', 'S2_4076A']
        diagProperties['diag_den']          = '[SII] 6731/6716'
        diagProperties['diag_tem']          = '[SII] 4069/4076'
        diagProperties['atom_temdiag']      = 'L(4069)/L(4076)'
        diagProperties['den_ratio']         = 'S2_6731A/S2_6716A'
        diagProperties['tem_ratio']         = 'S2_4069A/S2_4076A'        
        
#       #--Calculate NeSII and TeSII
        self.abunData['TeSII'],  neSII_TSII = self.den_temp_diagnostic_pair(diagProperties, low_density_dist, atom_temp=self.S2_atom)

        #-----------Oxygen
        diagProperties = {}
        diagProperties['required_denlines'] = ['S2_6716A', 'S2_6731A']
        diagProperties['required_temlines'] = ['O3_4363A', 'O3_4959A', 'O3_5007A']
        diagProperties['diag_den']          = '[SII] 6731/6716'
        diagProperties['diag_tem']          = '[OIII] 4363/5007+' 
        diagProperties['atom_temdiag']      = 'L(4363)/(L(5007)+L(4959))'
        diagProperties['den_ratio']         = 'S2_6731A/S2_6716A'
        diagProperties['tem_ratio']         = 'O3_4363A/(O3_4959A+O3_5007A)'

        #--Calculate NeSII and TeOIII
        self.abunData['TeOIII'], neSII_OIII = self.den_temp_diagnostic_pair(diagProperties, low_density_dist, atom_temp=self.O3_atom)

        #--Determine empirical TOIII from TSIII #Epm & Diaz 2005
        self.abunData['TeSIII_from_TeOIII'] = (1.05 * self.abunData.TeOIII/10000 - 0.08) * 10000

        #--Determine empirical TOII from TOIII #Dors Jr 2006
        self.abunData['TeOII_from_TeOIII']  = (1.397 / ((1/(self.abunData.TeOIII/10000)) + 0.385) ) * 10000

        #--Determine empirical TNII from TOIII #Epm 2014
        self.abunData['TeNII_from_TeOIII']  = (1.452 / ((1/(self.abunData.TeOIII/10000)) + 0.479) ) * 10000

        #-----------Nitrogen
        diagProperties = {}
        diagProperties['required_denlines'] = ['S2_6716A', 'S2_6731A']
        diagProperties['required_temlines'] = ['N2_5755A', 'N2_6548A', 'N2_6584A']
        diagProperties['diag_den']          = '[SII] 6731/6716'
        diagProperties['diag_tem']          = '[NII] 5755/6584+'
        diagProperties['atom_temdiag']      = '(L(6584) + L(6548)) / L(5755)'
        diagProperties['den_ratio']         = 'S2_6731A/S2_6716A'
        diagProperties['tem_ratio']         = '(N2_6548A+N2_6584A)/N2_5755A'

        #--Calculate Ne_SII and Te_NII
        self.abunData['TeNII'],  neSII_TNII = self.den_temp_diagnostic_pair(diagProperties, low_density_dist, atom_temp=self.N2_atom)

        #Assign object density from lines or from the low density distribution
        #--This code favors the neSII calculated from SIII-SII line pai+
        if 'neSII' not in self.abunData: 
            if np_sum(isnan(neSII_TSIII)) < self.MC_array_len:
                self.abunData['neSII'] = neSII_TSIII
            elif np_sum(isnan(neSII_OIII)) < self.MC_array_len:
                self.abunData['neSII'] = neSII_OIII
            else:
                self.abunData['neSII'] = neSII_TSIII
                
        #--Check if some results contain nan entries
        nanCount = OrderedDict()
        for electron_parameter in self.abunData.index:
            variable_array  = self.abunData[electron_parameter]
            nan_count       = np_sum(isnan(variable_array))
            if nan_count > self.MC_array_len * 0.90:
                self.abunData.drop(electron_parameter, inplace=True)            
            elif nan_count > 0:
                mag, error = nanmean(self.abunData[electron_parameter]), nanstd(self.abunData[electron_parameter])
                self.abunData[electron_parameter] = random.normal(mag, error, size = self.MC_array_len)
                if nan_count > self.MC_warning_limit:
                    nanCount[electron_parameter] = nan_count
        
        #Display calculations with issues
        if len(nanCount) > 0:
            print '-Issues calculating:'
            for element in nanCount:
                print '--', element, nanCount[element] 
        
        return

    def determine_ionic_abundance(self, abund_code, atom, diagnos_eval, diagnos_mag, tem, den):
                
        try:
            hbeta_flux = self.Hbeta_flux
        except AttributeError:
            hbeta_flux = self.H1_atom.getEmissivity(tem=tem, den=den, label = '4_2', product = False)
            print '--Warning using theoretical Hbeta emissivity'
            
        #Ionic abundance calculation using pyneb
        ionic_abund = atom.getIonAbundance(int_ratio = diagnos_mag, tem=tem, den=den, to_eval = diagnos_eval, Hbeta = hbeta_flux)

        #Evaluate the nan array
        nan_idcs    = isnan(ionic_abund)
        nan_count   = np_sum(nan_idcs)
        
        #Directly save if not nan
        if nan_count == 0:
            self.abunData[abund_code] = ionic_abund
            
        #Remove the nan entries performing a normal distribution
        elif nan_count < 0.90 * self.MC_array_len:
            mag, error = nanmean(ionic_abund), nanstd(ionic_abund)
            
            #Generate truncated array to store the data
            a, b = (0 - mag) / error, (1000 * mag - mag) / error
            new_samples = truncnorm(a, b, loc = mag, scale=error).rvs(size = nan_count)
            
            #Replace nan entries
            ionic_abund[nan_idcs] = new_samples
            self.abunData[abund_code] = ionic_abund
        
            if nan_count > self.MC_warning_limit:
                print '-- {} calculated with {}'.format(abund_code, nan_count)
                        
        return
        
    def check_obsLines(self, lines_list, just_one_line = False):
        
        #WARNING it would be better something that reads a standard preference over some.
        eval_lines = map(lambda x : 'L({})'.format(x[x.find('_') + 1 :len(x)-1]), lines_list) #Right format for pyneb eval: Ar3_7751A -> L(7751)                
        diagnos_eval = None
        
        #Case all lines are there
        if self.lines_dict.viewkeys() >= set(lines_list):
            diagnos_mag = zeros(self.MC_array_len)           
            for i in range(len(lines_list)):
                diagnos_mag += self.lines_dict[lines_list[i]]
            diagnos_eval = '+'.join(eval_lines)
        
        #Case we can use any line: #WARNING last line is favoured
        elif just_one_line:
            diagnos_mag = zeros(self.MC_array_len)           
            for i in range(len(lines_list)):
                if lines_list[i] in self.lines_dict:
                    diagnos_mag = self.lines_dict[lines_list[i]]
                    diagnos_eval  = eval_lines[i]
                   
        #Case none of the lines
        if  diagnos_eval is None:
            diagnos_mag = self.generate_nan_array()
            diagnos_eval = '+'.join(eval_lines)
                      
        return diagnos_eval, diagnos_mag

    def argon_abundance_scheme(self, Tlow, Thigh,  ne):

        #Calculate the Ar_+2 abundance according to the lines observed
        Ar3_lines = ['Ar3_7136A', 'Ar3_7751A']
        diagnos_eval, diagnos_mag   = self.check_obsLines(Ar3_lines, just_one_line = True)
        self.determine_ionic_abundance('ArIII_HII', self.Ar3_atom, diagnos_eval, diagnos_mag, Tlow, ne)
                
        #Calculate the Ar_+3 abundance according to the lines observed
        Ar4_lines = ['Ar4_4740A', 'Ar4_4711A']
        diagnos_eval, diagnos_mag   = self.check_obsLines(Ar4_lines, just_one_line = True)
        self.determine_ionic_abundance('ArIV_HII', self.Ar4_atom, diagnos_eval, diagnos_mag, Thigh, ne)
        
    def oxygen_abundance_scheme(self,  Tlow, Thigh, ne):                       

        #Calculate the O_+1 abundances from 3200+ lines
        O2_lines = ['O2_3726A+']
        diagnos_eval, diagnos_mag = self.check_obsLines(O2_lines)
        diagnos_eval = 'L(3726)+L(3729)'
        self.determine_ionic_abundance('OII_HII_3279A', self.O2_atom, diagnos_eval, diagnos_mag, Tlow, ne)
        
        #Calculate the O_+1 abundances from 7300+ lines
        O2_lines = ['O2_7319A+']
        diagnos_eval, diagnos_mag = self.check_obsLines(O2_lines)
        diagnos_eval = 'L(7319)+L(7330)'
        self.determine_ionic_abundance('OII_HII_7319A', self.O2_atom, diagnos_eval, diagnos_mag, Tlow, ne)
 
        #--Correction for recombination contribution Liu2000
        if 'OII_HII_7319A' in self.abunData:
            
            try:
                hbeta_flux = self.Hbeta_flux
            except AttributeError:
                hbeta_flux = self.H1_atom.getEmissivity(tem=Tlow, den=ne, label = '4_2', product = False)
                print '--Warning using theoretical Hbeta emissivity'
                        
            Lines_Correction = (9.36 * power((Tlow/10000), 0.44) * self.abunData.OII_HII_7319A) * hbeta_flux            
            ratio = self.lines_dict['O2_7319A+'] - Lines_Correction
            self.determine_ionic_abundance('OII_HII_7319A', self.O2_atom, diagnos_eval, ratio, Tlow, ne)
            
        #Get the ratios for empirical relation between OII lines
        if 'O2_3726A+' in self.lines_dict:
            self.abunData['O_R3200']    = self.lines_dict['O2_3726A+'] / self.Hbeta_flux
            print 'O_R3200', mean(self.abunData['O_R3200'])
            print 'OII_HII_3279A', mean(self.abunData['OII_HII_3279A'])
            print 'Original flux', mean(self.lines_dict['O2_3726A+'])
            
        if 'O2_7319A+' in self.lines_dict:
            self.abunData['O_R7300']    = self.lines_dict['O2_7319A+'] / self.Hbeta_flux
            print 'OII_HII_7319A', mean(self.abunData['OII_HII_7319A'])
        if self.lines_dict.viewkeys()   >= set(['O3_5007A']):
            self.abunData['O_R3']       = self.lines_dict['O3_5007A'] / self.Hbeta_flux       
        
        #Calculate the abundance from the empirical O_R3200_ffO2
        if  set(self.abunData.index) >=  {'O_R7300', 'O_R3'}:
            logRO2  = 1.913 + log10(self.abunData['O_R7300']) - 0.374*log10(self.abunData['O_R3']) / 0.806
            print 'logRO2', mean(logRO2)
            RO2 = power(10, logRO2)
            self.abunData['O_R3200_ffO2'] = RO2
            print 'O_R3200_ffO2', mean(self.abunData['O_R3200_ffO2'])
            print 'RO2*Hbeta', mean(RO2*self.Hbeta_flux)            
            diagnos_eval = 'L(3726)+L(3729)'
            self.determine_ionic_abundance('OII_HII_ffO2', self.O2_atom, diagnos_eval, RO2*self.Hbeta_flux, Tlow, ne)            
            print 'OII_HII_ffO2', mean(self.abunData['OII_HII_ffO2'])
            
        #Calculate the O_+2 abundance
        O3_lines = ['O3_4959A', 'O3_5007A']
        diagnos_eval, diagnos_mag = self.check_obsLines(O3_lines)    
        self.determine_ionic_abundance('OIII_HII', self.O3_atom, diagnos_eval, diagnos_mag, Thigh, ne)
         
        #Determine the O/H abundance (favoring the value from OII_HII
        if set(self.abunData.index) >= {'OII_HII_3279A', 'OIII_HII'}:
            self.abunData['OII_HII']    = self.abunData['OII_HII_3279A']
            self.abunData['OI_HI']      = self.abunData['OII_HII_3279A'] + self.abunData['OIII_HII']
        elif set(self.abunData.index) >= {'OII_HII_7319A', 'OIII_HII'}:
            self.abunData['OII_HII']    = self.abunData['OII_HII_7319A']
            self.abunData['OI_HI']      = self.abunData['OII_HII_7319A'] + self.abunData['OIII_HII']         
    
        if set(self.abunData.index) >= {'OII_HII_ffO2', 'OIII_HII'}:
            if set(self.abunData.index) >= {'OII_HII_3279A'}:
                self.abunData['OI_HI_ff02'] = self.abunData['OII_HII_3279A'] + self.abunData['OIII_HII']         
            else:
                self.abunData['OI_HI_ff02'] = self.abunData['OII_HII_ffO2'] + self.abunData['OIII_HII']         
              
        return
 
    def nitrogen_abundance_scheme(self, Tlow, ne):
        
        #Calculate TNII temperature from the CHAOS relation
        T_NII = Tlow #self.m_TNII_correction * Tlow + self.n_TNII_correction
        
        #Calculate the N+1 abundance
        N2_lines = ['N2_6548A', 'N2_6584A']
        diagnos_eval, diagnos_mag = self.check_obsLines(N2_lines)
        self.determine_ionic_abundance('NII_HII', self.N2_atom, diagnos_eval, diagnos_mag, T_NII, ne)
        
        #Calculate NI_HI using the OI_HI
        if set(self.abunData.index) >= {'NII_HII', 'OI_HI'}:
 
            #Compute  NI_OI 
            self.abunData['NI_OI']  = self.abunData['NII_HII'] / self.abunData['OII_HII']
            self.abunData['NI_HI']  = self.abunData['NI_OI'] * self.abunData['OI_HI']
                       
#             #Repeat calculation if 5755 line was observed to include the recombination contribution
#             if self.lines_dict.viewkeys() >= {'N2_5755A'}:         
#                  
#                 NIII_HI             = self.abunData.NI_HI - self.abunData['NII_HII']
#                 Lines_Correction    = 3.19 * power((Thigh/10000), 0.30) * NIII_HI * self.Hbeta_flux      
#                 self.abunData['TNII'], nSII = self.diags.getCrossTemDen(diag_tem = '[NII] 5755/6584+',
#                                                                         diag_den  = '[SII] 6731/6716',
#                                                                         value_tem = (self.lines_dict['N2_5755A'] - Lines_Correction)/(self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']),
#                                                                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#                  
#                 Ratio = self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']
#                 self.determine_ionic_abundance('NII_HII', self.N2_atom, Ratio, diagnos_mag, self.abunData['TNII'], ne)
#                                                  
#                 self.abunData['NI_OI'] = self.abunData['NII_HII'] / self.abunData['OII_HII']
#                 self.abunData['NI_HI'] = self.abunData['NI_OI'] * self.abunData['OI_HI']
                 
        return
        
    def sulfur_abundance_scheme(self, Tlow, ne, SIII_lines_to_use):

        print 'Metiendo esto', SIII_lines_to_use

        #Calculate the S+1 abundance
        S2_lines = ['S2_6716A', 'S2_6731A']
        diagnos_eval, diagnos_mag   = self.check_obsLines(S2_lines)
        self.determine_ionic_abundance('SII_HII', self.S2_atom, diagnos_eval, diagnos_mag, Tlow, ne)
        
        
        #Calculate the S+2 abundance
        S3_lines = ['S3_9069A', 'S3_9531A'] if SIII_lines_to_use == 'BOTH' else [SIII_lines_to_use]
        diagnos_eval, diagnos_mag   = self.check_obsLines(S3_lines)
        if set(S3_lines) != set(['S3_9069A', 'S3_9531A']):
            print '-- Using SIII lines', diagnos_eval
        
        self.determine_ionic_abundance('SIII_HII', self.S3_atom, diagnos_eval, diagnos_mag, Tlow, ne)

        #Calculate the total sulfur abundance
        if set(self.abunData.index) >= {'SII_HII', 'SIII_HII'}:
            
            self.abunData['SI_HI'] = self.abunData['SII_HII'] + self.abunData['SIII_HII']
            
            #Add the S+3 component if the argon correction is found
            if set(self.abunData.index) >= {'ArIII_HII', 'ArIV_HII'}:

                logAr2Ar3   = log10(self.abunData['ArIII_HII'] / self.abunData['ArIV_HII'])
                logSIV      = log10(self.abunData['SIII_HII']) - (logAr2Ar3 - self.n_SIV_correction) / self.m_SIV_correction
                SIV_HII     = power(10, logSIV)

                # Evaluate the nan array
                nan_idcs = isnan(SIV_HII)
                nan_count = np_sum(nan_idcs)

                # Directly save if not nan
                if nan_count == 0:
                    self.abunData['SIV_HII'] = SIV_HII

                # Remove the nan entries performing a normal distribution
                elif nan_count < 0.90 * self.MC_array_len:
                    mag, error = nanmean(SIV_HII), nanstd(SIV_HII)

                    # Generate truncated array to store the data
                    a, b = (0 - mag) / error, (1000 * mag - mag) / error
                    new_samples = truncnorm(a, b, loc=mag, scale=error).rvs(size=nan_count)

                    # Replace nan entries
                    SIV_HII[nan_idcs] = new_samples
                    self.abunData['SIV_HII'] = SIV_HII

                    if nan_count > self.MC_warning_limit:
                        print '-- {} calculated with {}'.format('SIV_HII', nan_count)

                self.abunData['SI_HI']   = self.abunData['SII_HII'] + self.abunData['SIII_HII'] + self.abunData['SIV_HII']
                self.abunData['ICF_SIV'] = self.abunData['SI_HI'] / (self.abunData['SII_HII'] + self.abunData['SIII_HII'])
        
        return
                              
    def helium_abundance_elementalScheme(self, Te, ne, lineslog_frame, metal_ext = ''):

        #Check temperatures are not nan before starting the treatment
        if (not isinstance(Te, float)) and (not isinstance(ne, float)): 
                        
            #HeI_indices = (lineslog_frame.Ion.str.contains('HeI_')) & (lineslog_frame.index != 'He1_8446A')  & (lineslog_frame.index != 'He1_7818A') & (lineslog_frame.index != 'He1_5016A')
            HeI_indices = (lineslog_frame.Ion.str.contains('HeI_')) & (lineslog_frame.index.isin(['He1_4472A','He1_5876A','He1_6678A']))
            HeI_labels  = lineslog_frame.loc[HeI_indices].index.values
            HeI_ions    = lineslog_frame.loc[HeI_indices].Ion.values

            Emis_Hbeta  = self.H1_atom.getEmissivity(tem=Te, den=ne, label = '4_2', product = False)
            
            #--Generating matrices with fluxes and emissivities
            for i in range(len(HeI_labels)):
                
                pyneb_code                  = float(HeI_ions[i][HeI_ions[i].find('_')+1:len(HeI_ions[i])]) 
                line_relative_Flux          = self.lines_dict[HeI_labels[i]] / self.Hbeta_flux
                line_relative_emissivity    = self.He1_atom.getEmissivity(tem = Te, den = ne, wave = pyneb_code, product = False) / Emis_Hbeta
                line_relative_emissivity    = self.check_nan_entries(line_relative_emissivity)

                if i == 0:
                    matrix_HeI_fluxes       = copy(line_relative_Flux)
                    matrix_HeI_emis         = copy(line_relative_emissivity)
                else:
                    matrix_HeI_fluxes       = vstack((matrix_HeI_fluxes, line_relative_Flux)) 
                    matrix_HeI_emis         = vstack((matrix_HeI_emis,   line_relative_emissivity))   

            matrix_HeI_fluxes = transpose(matrix_HeI_fluxes)
            matrix_HeI_emis   = transpose(matrix_HeI_emis)           
            
            #Perform the fit
            params = Parameters()
            params.add('Y', value = 0.01)
            HeII_HII_array = zeros(len(matrix_HeI_fluxes))
            HeII_HII_error = zeros(len(matrix_HeI_fluxes))
            for i in range(len(matrix_HeI_fluxes)):
                fit_Output = lmfit_minimmize(residual_Y_v3, params, args=(matrix_HeI_emis[i], matrix_HeI_fluxes[i]))
                HeII_HII_array[i] = fit_Output.params['Y'].value
                HeII_HII_error[i] = fit_Output.params['Y'].stderr
            
            #NO SUMANDO LOS ERRORES CORRECTOS?
            #self.abunData['HeII_HII_from_' + metal_ext] = random.normal(mean(HeII_HII_array), mean(HeII_HII_error), size = self.MC_array_len)
            ionic_abund = random.normal(mean(HeII_HII_array), mean(HeII_HII_error), size = self.MC_array_len)
            
            #Evaluate the nan array
            nan_count = np_sum(isnan(ionic_abund))
            if nan_count == 0:
                self.abunData['HeII_HII_from_' + metal_ext] = ionic_abund
            #Remove the nan entries performing a normal distribution
            elif nan_count < 0.90 * self.MC_array_len:
                mag, error = nanmean(ionic_abund), nanstd(ionic_abund)
                self.abunData['HeII_HII_from_' + metal_ext] = random.normal(mag, error, size = self.MC_array_len)
                if nan_count >  self.MC_warning_limit:
                    print '-- {} calculated with {}'.format('HeII_HII_from_' + metal_ext, nan_count)

            #Calculate the He+2 abundance 
            if self.lines_dict.viewkeys() >= {'He2_4686A'}:
                #self.abunData['HeIII_HII_from_' + metal_ext] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=Te, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
                self.determine_ionic_abundance('HeIII_HII_from_' + metal_ext, self.He2_atom, 'L(4685)', self.lines_dict['He2_4686A'], Te, ne)

            
            #Calculate elemental abundance
            Helium_element_keys = ['HeII_HII_from_' + metal_ext, 'HeIII_HII_from_' + metal_ext]
            if set(self.abunData.index) >= set(Helium_element_keys):
                self.abunData['HeI_HI_from_' + metal_ext] = self.abunData[Helium_element_keys[0]] + self.abunData[Helium_element_keys[1]]
            else:
                self.abunData['HeI_HI_from_' + metal_ext] = self.abunData[Helium_element_keys[0]]

            #Proceed to get the Helium mass fraction Y
            Element_abund = metal_ext + 'I_HI'
            Y_fraction, Helium_abund = 'Ymass_' + metal_ext, 'HeI_HI_from_' + metal_ext
            if set(self.abunData.index) >= {Helium_abund, Element_abund}:
                self.abunData[Y_fraction] = (4 * self.abunData[Helium_abund] * (1 - 20 * self.abunData[Element_abund])) / (1 + 4 * self.abunData[Helium_abund])
    
    def store_abundances_excel(self, objCode, catalogue_df, extension = ''):
                        
        #Store the values using the mean and the std from the array
        for parameter in self.abunData.index:

            mean_value, std_value = mean(self.abunData[parameter]), std(self.abunData[parameter])

            if (~isnan(mean_value)) & (~isnan(std_value)):
                catalogue_df.loc[objCode, parameter + extension] = ufloat(mean_value, std_value)
            else:
                catalogue_df.loc[objCode, parameter + extension] = np_nan

            print parameter, mean_value, std_value

        return

    def generate_nan_array(self):
        
        nan_array = empty(self.MC_array_len)
        nan_array[:] = np_nan    
        
        return nan_array

    def check_nan_entries(self, input_array):
        
        nan_count = np_sum(isnan(input_array))
        
        if  nan_count > 0:
            mag, error = nanmean(input_array), nanstd(input_array)
            new_distr = random.normal(mag, error, size = self.MC_array_len)
            if nan_count > 0.1 *self.MC_array_len:
                print '--Helium issue with {} nans'.format(nan_count)
        else:
            new_distr = input_array
        
        return new_distr
            
