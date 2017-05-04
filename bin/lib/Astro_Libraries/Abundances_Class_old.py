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
        atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
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
        self.S3_ratio           = self.S3_atom.getEmissivity(10000, 1000, wave = 9531) / self.S3_atom.getEmissivity(10000, 1000, wave = 9069)
        self.S3_9000_ratio      = random.normal(self.S3_atom.getEmissivity(10000, 1000, wave = 9531) / self.S3_atom.getEmissivity(10000, 1000, wave = 9069),  0.01, size = self.MC_array_len)
        self.N2_6000_ratio      = self.N2_atom.getEmissivity(10000, 1000, wave = 6584) / self.N2_atom.getEmissivity(10000, 1000, wave = 6548)
        
        #Factors to speed calculations
        self.lines_factors = {}
        self.lines_factors['S3_9069A'] = 1 + self.S3_ratio
        self.lines_factors['S3_9531A'] = 1 + 1/self.S3_ratio
        
        #Cloudy models for the SIV contribution
        
        self.m_SIV_correction   = random.normal(1.1628,  0.00559, size = self.MC_array_len)
        self.n_SIV_correction   = random.normal(0.0470,  0.0097, size = self.MC_array_len)
        #self.m_SIV_correction   = random.normal(1.109,  0.01, size = self.MC_array_len)
        #self.n_SIV_correction   = random.normal(0.135,  0.0173, size = self.MC_array_len)
 
        
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
 
        #Get the ratios for empirical relation between OII lines
        if 'O2_3726A+' in self.lines_dict:
            self.abunData['O_R3200']    = self.lines_dic / self.Hbeta_flux
        if 'O2_7319A+' in self.lines_dict:
            self.abunData['O_R7300']    = self.lines_dic / self.Hbeta_flux
        if self.lines_dict.viewkeys()   >= ('O3_4959A', 'O3_5007A'):
            self.abunData['O_R3']       = self.lines_dic / self.Hbeta_flux       
        
        
        
                Blended_lines['O2_3726A'] = ('O2_3726A', 'O2_3729A')
        Blended_lines['O2_7319A'] = ('O2_7319A', 'O2_7330A')
        
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
                                   
        return
 
    def nitrogen_abundance_scheme(self, Thigh, ne):
 
        #Calculate the N+1 abundance
        N2_lines = ['N2_6548A', 'N2_6584A']
        diagnos_eval, diagnos_mag = self.check_obsLines(N2_lines)
        self.determine_ionic_abundance('NII_HII', self.N2_atom, diagnos_eval, diagnos_mag, Thigh, ne)
        
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
                
                self.abunData['SIV_HII'] = power(10, logSIV)
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
        
        return

    def generate_nan_array(self):
        
        nan_array = empty(self.MC_array_len)
        nan_array[:] = np_nan    
        
        return nan_array



#     def store_dict_ObjLog(self, FileFolder, CodeName, verbose = False):
#         
#         keys, values = self.abunData.keys(), self.abunData.values()
#         
#         for i in range(len(keys)):
#             if values[i] is not None:
#                 
#                 #This removes the nan elements from the treatment values[i][~isnan(values[i])]
#                 Median_Value    = median(values[i][~isnan(values[i])])
#                 error           = std(values[i][~isnan(values[i])])  
#                 self.SaveParameter_ObjLog(CodeName, FileFolder, keys[i], ufloat(Median_Value, error))
# 
#                 if verbose == True:
#                     print i, keys[i], Median_Value, error
#                 
#             else:
#                 self.SaveParameter_ObjLog(CodeName, FileFolder, keys[i], values[i])
#                 if verbose == True:
#                     print i, keys[i], values[i]  
#         return
#
#     def abundances_to_catalogueDf(self, objCode, catalogue_df):
#         
#         for parameter in self.abunData.index:
#             
#             mean_value, std_value = mean(self.abunData[parameter]), std(self.abunData[parameter])
#             if (~isnan(mean_value)) & (~isnan(std_value)):
#                 catalogue_df.loc[objCode, parameter] = ufloat(mean_value, std_value)
#             else:
#                 catalogue_df.loc[objCode, parameter] =  np_nan        
#         
#         return
#     
#     def store_abundances(self, objCode, catalogue_df, catalogue_address):
#                 
#         for parameter in self.abunData.index:
#             
#             mean_value, std_value = mean(self.abunData[parameter]), std(self.abunData[parameter])
#             if (~isnan(mean_value)) & (~isnan(std_value)):
#                 catalogue_df.loc[objCode, parameter] = ufloat(mean_value, std_value)
#             else:
#                 catalogue_df.loc[objCode, parameter] =  np_nan
#         
#         catalogue_df.to_pickle(catalogue_address)
#         
#         return
#
#     def low_density_function(self, density_array):
#         
#         inds_nan    = where(isnan(density_array))
#         
#         if len(inds_nan[0]) > 0:
#             
#             if len(inds_nan[0]) > (0.1 * self.MC_array_len):
#                 print '--WARNING 10% boundary: ', len(inds_nan[0]), 'values with nan'            
#             elif len(inds_nan[0]) > (0.05 * self.MC_array_len):
#                 print '--WARNING 5% boundary: ', len(inds_nan[0]), 'values with nan'            
#                                 
#             density_array[inds_nan]= self.low_density_dist[inds_nan]
#             
#         mean_density    = mean(density_array) 
#         std_density     = std(density_array)   
#         
#         if isnan(mean_density): 
#             print '--WARNING: Density array cannot be treated numerically'
#             
#         if mean_density < 75:
#             #print '--WARNING: Low density', mean_density, '+/-', std_density, 'Points below limit', len(where(density_array < 75)[0])
#             return  self.low_density_dist
#         else:                              
#             return density_array
# 
#     def high_temp_function(self, temperature_array):
# 
#         inds_nan    = where(isnan(temperature_array))
#         n_nan = len(inds_nan[0])
#         if n_nan > 0:
# 
#             Te_mean = mean(temperature_array)
#             Te_std = std(temperature_array)
#             
#             
#             if n_nan > (0.1 * self.MC_array_len):
#                 print '--WARNING 10% boundary: ', n_nan, 'values with nan'            
#             elif n_nan > (0.05 * self.MC_array_len):
#                 print '--WARNING 5% boundary: ', n_nan, 'values with nan'   
#                         
#             temperature_array[inds_nan] = random.normal(Te_mean, Te_std, size = n_nan)
#                 
#         return temperature_array
#
#     def argon_abundance_scheme_old(self):
#         
#         #Physical parameters
#         TSIII, TOIII, ne = None, None, None
#         
#         #Select the physical parameters
#         if 'TeSIII' in self.abunData:
#             TSIII = self.abunData.TeSIII
#             TOIII = self.abunData.TeOIII_from_TeSIII
#         elif 'TeOIII' in self.abunData:
#             TSIII = self.abunData.TeSIII_from_TeOIII
#             TOIII = self.abunData.TeOIII
#          
#         #Select the density
#         if 'neSII' in self.abunData:
#             ne = self.abunData.neSII
#         elif 'neOII' in self.abunData:
#             ne = self.abunData.neOII
#                  
#         if ne is not None:
#         
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)            
#         
#             #Calculate the Ar_+2 abundance
#             if TSIII is not None:
#                 if self.lines_dict.viewkeys() >= {'Ar3_7136A', 'Ar3_7751'}:    
#                     Ratio = self.lines_dict['Ar3_7751A'] + self.lines_dict['Ar3_7136A']
#                     self.abunData['ArIII_HII'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, to_eval = 'L(7751) + L(7136)', Hbeta = self.Hbeta_flux)     
#                 
#                 elif self.lines_dict.viewkeys() >= {'Ar3_7136A'}:
#                     Ratio = self.lines_dict['Ar3_7136A']
#                     self.abunData['ArIII_HII'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, wave = 7136, Hbeta = self.Hbeta_flux)
#                                  
#                 elif self.lines_dict.viewkeys() >= {'Ar3_7751A'}:
#                     Ratio = self.lines_dict['Ar3_7751A']
#                     self.abunData['ArIII_HII'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, wave = 7751, Hbeta = self.Hbeta_flux)   
# 
#             #Calculate the Ar_+3 abundance
#             if TOIII is not None:
#                 if self.lines_dict.viewkeys() >= {'Ar4_4740A', 'Ar4_4711A'}:
#                     Ratio = self.lines_dict['Ar4_4740A'] + self.lines_dict['Ar4_4711A']
#                     self.abunData['ArIV_HII'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, to_eval = 'L(4740) + L(4711)', Hbeta = self.Hbeta_flux)
#                                            
#                 elif self.lines_dict.viewkeys() >= {'Ar4_4740A'}:
#                     Ratio = self.lines_dict['Ar4_4740A']
#                     self.abunData['ArIV_HII'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, wave = 4740, Hbeta = self.Hbeta_flux)
#                     
#                 elif self.lines_dict.viewkeys() >= {'Ar4_4711A'}:
#                     Ratio = self.lines_dict['Ar4_4711A']
#                     self.abunData['ArIV_HII'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, wave = 4711, Hbeta = self.Hbeta_flux)

#     def sulfur_abundance_scheme_old(self):
# 
#         #Physical parameters
#         Te, ne = None, None
#                 
#         #Select the temperature
#         if 'TeSIII' in self.abunData:
#             Te = self.abunData.TeSIII
#         elif 'TeSIII_from_TeOIII' in self.abunData:
#             Te = self.abunData.TeSIII_from_TeOIII
# 
#         #Select the density
#         if 'neSII' in self.abunData:
#             ne = self.abunData.neSII
#         
#         #Calculating the R_S2 ratio
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:
#             self.abunData['R_SII_abund'] = self.lines_dict['S2_6716A'] + self.lines_dict['S2_6731A']
#                 
#         #Calculating the R_S3 ratio
#         if self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_9531A'}:
#             self.abunData['R_SIII_abund'] = self.lines_dict['S3_9069A'] + self.lines_dict['S3_9531A']
#         
#         elif self.lines_dict.viewkeys() >= {'S3_9069A'}:
#             self.abunData['R_SIII_abund'] = self.lines_dict['S3_9069A'] * (1 + self.S3_9000_ratio)
#             
#         elif self.lines_dict.viewkeys() >= {'S3_9531A'}:
#             self.abunData['R_SIII_abund'] = self.lines_dict['S3_9531A'] * (1 + 1/self.S3_9000_ratio)          
#                 
#         #Determine ionic abundances
#         if (Te is not None) and (ne is not None):
#                         
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             
#             #Calculate the S+1 abundance
#             if 'R_SII_abund' in self.abunData:
#                 self.abunData['SII_HII'] = self.S2_atom.getIonAbundance(int_ratio=self.abunData.R_SII_abund, tem=Te, den=ne, to_eval = 'L(6731)+L(6716)', Hbeta = self.Hbeta_flux)        
#             
#             #Calculate the S+2 abundance
#             if 'R_SIII_abund' in self.abunData:
#                 self.abunData['SIII_HII'] = self.S3_atom.getIonAbundance(int_ratio=self.abunData.R_SIII_abund, tem=Te, den=ne, to_eval = 'L(9069)+L(9531)', Hbeta = self.Hbeta_flux)
#     
#             #Calculating Sulfur abundance
#             if ('SII_HII' in self.abunData) and ('SIII_HII' in self.abunData):
#                 
#                 self.abunData['SI_HI'] = self.abunData['SII_HII'] + self.abunData['SIII_HII']
#                 
#                 #Calculating the S+3 abundance from Ar - Sulfur from Popstar - Cloudy models if Argon metallicity previously measured                 
#                 if ('ArIII_HII' in self.abunData) and ('ArIV_HII' in self.abunData):
#                     logAr2Ar3   = log10(self.abunData['ArIII_HII'] / self.abunData['ArIV_HII'])
#                     logSIV      = log10(self.abunData['SIII_HII']) - (logAr2Ar3 - self.n_SIV_correction) / self.m_SIV_correction
#                     
#                     #Saving S+3 value
#                     self.abunData['SIV_HII'] = power(10, logSIV)
#                                         
#                     #Adding S+3 value to S total value
#                     self.abunData['SI_HI'] = self.abunData['SI_HI'] + self.abunData['SIV_HII']
#                     
#         return

#     def oxygen_abundance_scheme_old(self):
# 
#         #Physical parameters
#         TOII, TOIII, ne = None, None, None
#         
#         #Select the physical parameters
#         if 'TeOIII' in self.abunData:
#             TOIII = self.abunData.TeOIII
#         elif 'TeOIII_from_TeSIII' in self.abunData:
#             TOIII = self.abunData.TeOIII_from_TeSIII
#         
#         if 'TeSIII' in self.abunData:
#             TOII = self.abunData.TeSIII
#         elif 'TeSIII_from_TeOIII' in self.abunData:
#             TOII = self.abunData.TeSIII_from_TeOIII
#         
#         #Select the density
#         if 'neSII' in self.abunData:
#             ne = self.abunData.neSII
#         elif 'neOII' in self.abunData:
#             ne = self.abunData.neOII
#      
#         if ne is not None:
#               
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)   
#             
#             #Calculate the O+1 ion abundance (If O2_3726A, O2_3729A are available we prioritize them)
#             if TOII is not None:
#                 
#                 if self.lines_dict.viewkeys() >= {'O2_3726A+'}: #WE ASSUME THEY ARE BLENDED
#                     ratio = self.lines_dict['O2_3726A+']
#                     self.abunData['OII_HII_3279A'] = self.O2_atom.getIonAbundance(int_ratio = ratio, tem=TOII, den=ne, to_eval = 'L(3726)+L(3729)', Hbeta = self.Hbeta_flux)
#                             
#                 if self.lines_dict.viewkeys() >= {'O2_7319A+'}: #WE ASSUME THEY ARE BLENDED
#                     ratio = self.lines_dict['O2_7319A+'] 
#                     OII_HII_7319_initial = self.O2_atom.getIonAbundance(int_ratio = ratio, tem=TOII, den=ne, to_eval = 'L(7319)+L(7330)', Hbeta = self.Hbeta_flux)     
#                     
#                     #Correction for recombination contribution Liu2000
#                     Lines_Correction = (9.36 * power((TOII/10000), 0.44) * OII_HII_7319_initial) * self.Hbeta_flux      
#                     ratio = self.lines_dict['O2_7319A+'] - Lines_Correction
#                     self.abunData['OII_HII_7319A'] = self.O2_atom.getIonAbundance(int_ratio = ratio, tem=TOII, den=ne, to_eval = 'L(7319)+L(7330)', Hbeta = self.Hbeta_flux) 
#    
#             #Calculate the O+2 ion abundance
#             if TOIII is not None:
#                 if self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A'}:
#                     ratio = self.lines_dict['O3_4959A'] + self.lines_dict['O3_5007A']
#                     self.abunData['OIII_HII'] = self.O3_atom.getIonAbundance(int_ratio = ratio, tem=TOIII, den=ne, to_eval = 'L(4959)+L(5007)', Hbeta = self.Hbeta_flux)                           
#                           
#         #Determine the O/H abundance
#         if ('OII_HII_3279A' in self.abunData) and ('OIII_HII' in self.abunData):
#             self.abunData['OI_HI'] = self.abunData['OII_HII_3279A'] + self.abunData['OIII_HII'] 
#             self.abunData['OII_HII'] = self.abunData['OII_HII_3279A']
#             
#         elif ('OII_HII_7319A' in self.abunData) and ('OIII_HII' in self.abunData):
#             self.abunData['OI_HI'] = self.abunData['OII_HII_7319A'] + self.abunData['OIII_HII']   
#             self.abunData['OII_HII'] = self.abunData['OII_HII_7319A']  
#         
#         return
   
# class pyneb_tools():
#     
#     def __init__(self):
#  
#         #Set atomic data 
#         atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
#         atomicData.setDataFile('s_iii_coll_HRS12.dat')
#  
#         #Declare hydrogen atom so it becomes readily available:
#         self.H1_atom    = RecAtom('H', 1)
#         self.S2_atom    = Atom('S', 2)
#         self.S3_atom    = Atom('S', 3)
#         self.He1_atom   = RecAtom('He', 1)
#         self.He2_atom   = RecAtom('He', 2)
#         
#     def density_pyneb(self, Ion, Temp):
#         
#         element     = Ion[0:-1]
#         ion_n       = int(Ion[-1])
#                 
#         if Ion == 'S2':
#             
#             Ion_Atom    = Atom(element, ion_n)
#             
#             emLine      = ["S2_6716A","S2_6731A"]
#             Flux_dict   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#             
#             FluxRatio   = Flux_dict['S2_6716A'] / Flux_dict['S2_6731A']
#             
#             Den         = Ion_Atom.getTemDen(FluxRatio, tem = Temp, wave1=6716, wave2=6731)
#             
#             return Den
#         
#     def temperature_pyneb(self, Ion, Den):
#         
#         element     = Ion[0:-1]
#         ion_n       = int(Ion[-1])
#         
#         if Ion == 'O3':        
#         
#             Ion_Atom    = Atom(element, ion_n)
#             
#             emLine      = ["S2_6716A","S2_6731A"]
#             Flux_dict   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)        
#         
#             FluxRatio   = (Flux_dict['O3_4959A'] + Flux_dict['O3_5007A']) / Flux_dict['O3_4363A']
# 
#             Temp        = Ion_Atom.getTemDen(FluxRatio, den = Den, to_eval = "(L(4959) + L(5007)) / L(4363)")
#         
#             return Temp
#    
# class direct_Method(ReddeningLaws, pyneb_tools):
#     
#     def __init__(self):
#         
#         pyneb_tools.__init__(self)
#         ReddeningLaws.__init__(self)
#     
#         self.den_Methodology        = None
#         self.tem_Methodology        = None
#         
#         self.m_SIV_correction       = ufloat(1.008,0.037)
#         self.n_SIV_correction       = ufloat(1.576,0.066)
#                         
#     def DeclareObject(self, FileFolder, FileName, CodeName, AbundancesFileExtension = '_WHT_LinesLog_v3.txt', cHbeta_type = 'cHBeta_red'):
#         #WARNING: This method should be define at a lower level
#                 
#         self.FileFolder             = FileFolder
#         self.FileName               = FileName
#         self.CodeName               = CodeName
#     
#         EmLine_List                 = ['R_SII', 'R_SIIprime', 'R_SIII', 'R_NII', 'R_OII', 'R_OIII']
#         Den_List                    = ['nSII']
#         Temp_List                   = ['TOIII', 'TOII', 'TSII', 'TSIII','TNII', 'TOII_approx_TOIII', 'TSIII_approx_TOIII', 'TOIII_approx_TSIII', 'TNII_approx_TOIII']
#         IonicAbund_List             = ['SII_HII', 'SIII_HII', 'SIV_HII', 'OII_HII', 'OII_HII_3279A', 'OII_HII_7319A', 'NII_HII', 'HeII_HII', 'HeIII_HII', 'ArIII_HII', 'ArIV_HII']
#         Element_List                = ['SI_HI', 'SI_HI_ArCorr', 'OI_HI', 'NI_OI', 'NI_HI', 'HeI_HI']
#             
#         self.Properties_dict        = dict.fromkeys(Element_List + IonicAbund_List + Den_List + Temp_List + EmLine_List)
#         
#         self.cHbeta                 = self.GetParameter_ObjLog(self.CodeName, self.FileFolder, Parameter = cHbeta_type, Assumption = 'float')
#         
#         self.AbundancesExtension    = AbundancesFileExtension
#         
#     def density_determination(self, Ion, Temp = None, methodology = None):
#                 
#         #Paremetric approach
#         if methodology == 'Epm2014':
#             if Ion == 'S2':
#                 
#                 emLine                      = ["S2_6716A", "S2_6731A"]
#                 Flux_dict                   = self.getLinesFlux_dered(emLine, self.cHbeta)
#                 Flux_dict['Temp']           = Temp
#                 den                         = self.empiric_formulae(Flux_dict, methodology, Ion, 'nSII')
#                 
#                 #Check for non physical output 
#                 den                         = self.check_issues(magnitude = den, parameter_type = 'density')
#                 return den
#         
#         #Using pyneb 
#         if methodology == 'pyneb':
#             
#             #Calculate density using pyneb
#             den = self.density_pyneb(Ion, Temp)
#             
#             #Check for non physical output #WARNING: For pyneb this is too general
#             den = self.check_issues(magnitude = den, parameter_type = 'density')
#             
#             return den
#     
#     def temperature_determination(self, Ion, dens = None, methodology = None, Temp = None, RecCorr = None):
#                 
#         #Direct methods
#         if methodology == 'Epm2014':
#     
#             if Ion == 'O3':
#                 
#                 emLine                          = ["O3_4959A","O3_5007A","O3_4363A"]
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TOIII')
#  
#             if Ion == 'O2':               
#                 data_dict                       = {'Temp': Temp}                 
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TOII_approx_TOIII')
#                                 
#             if Ion == 'S3':
#                 
#                 emLine                          = ["S3_9069A","S3_9531A","S3_6312A"]
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TSIII')
#                             
#             if Ion == 'N2':
#                 
#                 emLine                          = ["N2_6548A","N2_6584A", "N2_5755A", 'H1_4861A']
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 if RecCorr != None:
#                     if data_dict["N2_5755A"] != None:
#                         data_dict["N2_5755A"] = data_dict["N2_5755A"] - RecCorr * data_dict['H1_4861A']
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TNII')
# 
#         if methodology == 'Haegele2008':
# 
#             if Ion == 'S2':
#                 
#                 emLine                          = ["S2_6716A","S2_6731A", "S2_4069A", 'S2_4076A']
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta =  self.cHbeta)
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TSII')
# 
#             if Ion == 'O2':
# 
#                 emLine                          = ['O2_3726A',"O2_3729A","O2_7319A",'O2_7330A']
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TOII')
# 
#         
#         #Using pyneb 
#         if methodology == 'pyneb':
#             
#             temp = self.temperature_pyneb(Ion, dens)
# 
#         return temp
# 
#     def argon_IonAbundance(self, Temp, Den, Ion, methodology = 'Haegele2008'):
#         
#         if methodology == 'Haegele2008':
#             
#             if Ion == 'Ar3':
#                 
#                 emLine                              = ["Ar3_7136A","H1_4861A"]
#                 data_dict                           = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                   = Temp
#                 self.Properties_dict['ArIII_HII']   = self.empiric_formulae(data_dict, methodology, Ion, 'ArIII_HII')
#                 print 'The argon thing III', self.Properties_dict['ArIII_HII']
#                 return
# 
#             if Ion == 'Ar4':
#                 
#                 emLine                              = ["Ar4_4740A","H1_4861A"]
#                 data_dict                           = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                   = Temp             
#                 self.Properties_dict['ArIV_HII']    = self.empiric_formulae(data_dict, methodology, Ion, 'ArIV_HII')
#                 print 'The argon thing IV', self.Properties_dict['ArIV_HII']
#                 return
#                        
#     def oxygen_IonAbundance(self, Temp, Den, Ion, methodology = 'Epm2014', RecCorr = None):
# 
#         if methodology == 'Epm2014':
#     
#             if Ion == 'O3':
#                 
#                 emLine                                      = ["O3_4959A","O3_5007A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                          = Temp
#                 self.Properties_dict['OIII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'OIII_HII')
#                 
#                 return
# 
#             if Ion == 'O2':
#                                                             #This assumes both lines are blended
#                 emLine                                      = ["O2_3726A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta, Mode='Integrated')
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den 
#                 self.Properties_dict['OII_HII_3279A']       = self.empiric_formulae(data_dict, methodology, Ion, 'OII_HII')
#                 
#         if methodology == 'Fabian2006':
# 
#             if Ion == 'O2':
#                 
#                 emLine                                      = ["O2_7319A","O2_7330A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']        = Temp, Den 
#                 
#                 if RecCorr == None:
#                     self.Properties_dict['OII_HII_7319A']   = self.empiric_formulae(data_dict, methodology, Ion, 'OII_HII_7319A')
#                 else:
#                     data_dict["O2_7319A"]                   = data_dict["O2_7319A"] - RecCorr
#                     self.Properties_dict['OII_HII_7319A']   = self.empiric_formulae(data_dict, methodology, Ion, 'OII_HII_7319A')
#                     
#                 return
#             
#     def nitrogen_IonAbundance(self, Temp, Den, Ion, methodology = 'Epm2014'):
#             
#         if methodology == 'Epm2014':
#             
#             if Ion == 'N2':
#             
#                 emLine                                      = ["N2_6548A","N2_6584A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)            
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['NII_HII']             = self.empiric_formulae(data_dict, methodology, Ion, 'NII_HII')
#                 
#                 return
#             
#     def sulfur_IonAbundance(self,  Temp, Den, Ion, methodology = 'Haegele2008'):
#         
#         if methodology == 'Haegele2008':
#             
#             if Ion == 'S3':
#                 
#                 emLine                                      = ["S3_9069A","S3_9531A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['SIII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'SIII_HII')
#                 
#                 return
#       
#             if Ion == 'S2':
#                 
#                 emLine                                      = ["S2_6716A","S2_6731A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['SII_HII']             = self.empiric_formulae(data_dict, methodology, Ion, 'SII_HII')
# 
#         if methodology == 'Vital2015':
#             if Ion == 'S3':                
#                 emLine                                      = ["S3_9069A","S3_9531A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['SIII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'SIII_HII')
# 
#                 return
#       
#             if Ion == 'S2':
#                 methodology = 'Haegele2008'
#                 emLine                                      = ["S2_6716A","S2_6731A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['SII_HII']             = self.empiric_formulae(data_dict, methodology, Ion, 'SII_HII')
# 
#                 
#                 return
#     
#     def helium_IonAbundance(self, Temp, Den, Ion, methodology = 'Fabian2006'):
#         
#         if methodology == 'Fabian2006':
#             
#             if Ion == 'He1':                
#                 emLine                                      = ["He1_4472A", "He1_5876A", "He1_6678A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['HeII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'HeII_HII')
#                 
#                 return
#             
#             if Ion == 'He2':
#                 emLine                                      = ["He2_4686A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['HeIII_HII']           = self.empiric_formulae(data_dict, methodology, Ion, 'HeIII_HII')
#                  
#                 return
#                  
#         if methodology == 'Vital2015':
#             
#             if Ion == 'He1':                
#                 emLine                                      = ["He1_4472A", "He1_5876A", "He1_6678A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['HeII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'HeII_HII')
#                             
#             if Ion == 'He2':
#                 emLine                                      = ["He2_4686A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['HeIII_HII']           = self.empiric_formulae(data_dict, methodology, Ion, 'HeIII_HII')
#         
#         if methodology == 'PyNeb':
# 
#             if Ion == 'He1':                
#                 emLine                                      = ["He1_4472A", "He1_5876A", "He1_6678A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['HeII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'HeII_HII')
# 
#                 return
# 
#             if Ion == 'He2':
#         
#                 return
# 
#     def empiric_formulae(self, data_dict, methodology, Ion, Parameter):
#     
#         Magnitude = None
#         if data_dict != None:
#             if (None not in data_dict.values()):
#                 if methodology == 'Epm2014':
#                     if Parameter == 'nSII':
#                     
#                         self.Properties_dict['R_SII']               = data_dict['S2_6716A'] / data_dict['S2_6731A']
#                         t4                                          = data_dict['Temp'] / 10000.0       
#                         a_0                                         = 16.054  - 7.79  / t4 - 11.32  * t4
#                         a_1                                         = -22.66  + 11.08 / t4 + 16.02  * t4
#                         b_0                                         = -21.61  + 11.89 / t4 + 14.59  * t4
#                         b_1                                         = 9.17    - 5.09  / t4 - 6.18   * t4
#                         Magnitude                                   = 1000.0 * (self.Properties_dict['R_SII']*a_0 + a_1) / (self.Properties_dict['R_SII']*b_0 + b_1)
#                         #self.Properties_dict['nSII']               = 1000.0 * (self.Properties_dict['R_SII']*a_0 + a_1) / (self.Properties_dict['R_SII']*b_0 + b_1)
#                         
#                     elif Parameter == 'TOIII':
#                         self.Properties_dict['R_OIII']              = (data_dict['O3_4959A'] + data_dict['O3_5007A']) / data_dict['O3_4363A']
#                         self.Properties_dict['TOIII']               = (0.7840 - 0.0001357 * self.Properties_dict['R_OIII'] + (48.44 / self.Properties_dict['R_OIII'])) * 10000
#                         Magnitude                                   = self.Properties_dict['TOIII']
#                         
#                     elif Parameter == 'TOII_approx_TOIII':
#                         t4 = data_dict['Temp'] / 10000
#                         self.Properties_dict['TOII_approx_TOIII']   = (1.397 /( (1/t4) + 0.385) ) * 10000
#                         Magnitude                                   = self.Properties_dict['TOII_approx_TOIII']
#                         
#                     elif Parameter == 'TSIII_approx_TOIII':
#                         #From equation TOIII = 1.0807 * TSIII - 0.0846
#                         t4 = data_dict['Temp'] / 10000
#                         self.Properties_dict['TSIII_approx_TOIII']  = ((t4 + 0.0846) / 1.0807)*10000
# 
#                     elif Parameter == 'TOIII_approx_TSIII':    
#                         t4 = data_dict['Temp'] / 10000
#                         self.Properties_dict['TOIII_approx_TSIII']  = (1.0807 * t4 - 0.0846)*10000
#                         Magnitude                                   = self.Properties_dict['TOIII_approx_TSIII']
# 
#                     elif Parameter == 'TSIII':
#                         self.Properties_dict['R_SIII']              = (data_dict['S3_9069A'] + data_dict['S3_9531A']) / data_dict['S3_6312A']
#                         self.Properties_dict['TSIII']               = (0.5147 + 0.0003187 * self.Properties_dict['R_SIII'] + (23.6404 / self.Properties_dict['R_SIII'])) * 10000                         
#                         Magnitude                                   = self.Properties_dict['TSIII']
#     
#                     elif Parameter == 'TNII':
#                         self.Properties_dict['R_NII']               = (data_dict['N2_6548A'] + data_dict['N2_6584A']) / data_dict['N2_5755A']
#                         self.Properties_dict['TNII']                = (0.6153 - 0.0001529 * self.Properties_dict['R_NII'] + (35.3641 / self.Properties_dict['R_NII'])) * 10000
#                         Magnitude                                   = self.Properties_dict['TNII']
# 
#                     elif Parameter == 'TNII_approx_TOIII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         self.Properties_dict['TNII_approx_TOIII']   = (1.452 /( (1/t4) + 0.479))*10000
#                         Magnitude                                   = self.Properties_dict['TNII_approx_TOIII']
# 
#                     elif Parameter == 'OIII_HII':
#                         t4 = data_dict['Temp'] / 10000
#                         logOIII_logHII                              = -12 + umath.log10((data_dict['O3_4959A'] + data_dict['O3_5007A']) / data_dict['H1_4861A']) + 6.1868 + 1.2491 / t4 - 0.5816 * umath.log10(t4)
#                         Magnitude                                   = umath.pow(10, logOIII_logHII)
#                     
#                     elif Parameter == 'OII_HII_3279A':
#                         t4 = data_dict['Temp'] / 10000
#                         ne = data_dict['Den']
#                         logOII_logHII                               = -12 + umath.log10((data_dict['O2_3726A']) / data_dict['H1_4861A']) + 5.887 + 1.641 / t4 - 0.543 * umath.log10(t4) + 0.000114 * ne
#                         Magnitude                                   = umath.pow(10, logOII_logHII) 
#                         
#                     elif Parameter == 'NII_HII':
#                         t4  = data_dict['Temp'] / 10000                       
#                         logNII_logHII                               = -12 + umath.log10((data_dict['N2_6548A'] + data_dict['N2_6584A']) / data_dict['H1_4861A']) + 6.291 + 0.90221 / t4 - 0.5511 * umath.log10(t4)
#                         Magnitude                                   = umath.pow(10, logNII_logHII)
# 
#                 elif methodology == 'Angeles2015':
# 
#                     if Parameter == 'SIV_HII':
#  
#                         ArIII_HII   = data_dict['ArIII_HII']
#                         ArIV_HII    = data_dict['ArIV_HII']
#                         SIII_HII    = data_dict['SIII_HII']
#                                                  
#                         if (ArIII_HII != None) and (ArIV_HII != None) and (ArIV_HII != 0.0):   #Somehow ArIV should be saved as None 
#                             logAr2Ar3   = umath.log10(ArIII_HII/ArIV_HII)
# 
#                         else:
#                             logAr2Ar3 = ufloat(0.0, 0.0)
#                       
#                         logSIV   = umath.log10(SIII_HII) - (logAr2Ar3 - self.n_SIV_correction) / self.m_SIV_correction
# 
#                                   
#                         self.Properties_dict['SIV_HII'] = umath.pow(10, logSIV)
#                         Magnitude = self.Properties_dict['SIV_HII']
# 
# #                         Old formulation
# #                         ArIII_HII   = data_dict['ArIII_HII']
# #                         ArIV_HII    = data_dict['ArIV_HII']
# #                         SIII_HII    = data_dict['SIII_HII']
# #                                                 
# #                         if (ArIII_HII != None) and (ArIV_HII != 0.0):    
# #                             y_Ar = umath.log10(ArIII_HII/ArIV_HII)
# #                         else:
# #                             y_Ar = ufloat(0.0, 0.0)
# #                      
# #                         x_S = (y_Ar - 0.126) / 1.192
# #                 
# #                         self.Properties_dict['SIV_HII'] = SIII_HII / umath.pow(10, x_S)
# #                         Magnitude = self.Properties_dict['SIV_HII']
#                         
#                 elif methodology == 'Haegele2008':
#                     
#                     if Parameter == 'TSII':
#                         self.Properties_dict['R_SIIprime']          = (data_dict['S2_6716A'] + data_dict['S2_6731A']) / (data_dict['S2_4069A'] + data_dict['S2_4076A'])
#                         self.Properties_dict['TSII']               = (1.92 - 0.0375 *self.Properties_dict['R_SIIprime'] - 14.5 / self.Properties_dict['R_SIIprime'] + 105.64 / (self.Properties_dict['R_SIIprime'] * self.Properties_dict['R_SIIprime'])) * 10000
#                         Magnitude                                   = self.Properties_dict['TSII']
# 
#                     if Parameter == 'TOII':
#                         ne = data_dict['Den']
#                         self.Properties_dict['R_OII']               = (data_dict['O2_3726A'] + data_dict['O2_3729A']) / (data_dict['O2_7319A'] + data_dict['O2_7330A'])
#                         self.Properties_dict['TOII']               = (0.23 + 0.0017 * self.Properties_dict['R_OII'] + 38.3 / self.Properties_dict['R_OII'] + umath.log10(1 + 0.0001 * ne)) * 10000
#                         Magnitude                                   = self.Properties_dict['TOII']
# 
#                     elif Parameter == 'ArIII_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         logArIII_logHII                             = -12 + umath.log10(data_dict['Ar3_7136A'] / data_dict['H1_4861A']) +  6.157 + 0.808/t4 - 0.508 * umath.log10(t4) 
#                         Magnitude                                   = umath.pow(10, logArIII_logHII)
# 
#                     elif Parameter == 'ArIV_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         logArIV_logHII                              = -12 + umath.log10(data_dict['Ar4_4740A'] / data_dict['H1_4861A']) +  5.705 + 1.246/t4 - 0.156 * umath.log10(t4) 
#                         Magnitude                                   = umath.pow(10, logArIV_logHII)
# 
#                     elif Parameter == 'SII_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         ne = data_dict['Den']
#                         logSII_logHII                               = -12 + umath.log10((data_dict['S2_6716A'] + data_dict['S2_6731A']) / data_dict['H1_4861A']) + 5.423 + 0.929/t4 - 0.280 * umath.log10(t4) + umath.log10(1 + 0.0001 * ne)
#                         Magnitude                                   = umath.pow(10, logSII_logHII)
# 
#                     elif Parameter == 'SIII_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         logSIII_logHII                              = -12 + umath.log10((data_dict['S3_9069A'] + data_dict['S3_9531A']) / data_dict['H1_4861A']) + 5.8 + 0.77 / t4 - 0.22 * umath.log10(t4) 
#                         Magnitude                                   = umath.pow(10, logSIII_logHII)
#                 
#                 elif methodology == 'Vital2015':
# 
#                     if Parameter == 'SIII_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         logSIII_logHII                              = -12 + umath.log10((data_dict['S3_9069A'] + data_dict['S3_9531A']) / data_dict['H1_4861A']) + 6.012 + 0.6309 / t4 - 0.5722 * umath.log10(t4) 
#                         Magnitude                                   = umath.pow(10, logSIII_logHII)   
#                         
#                     if Parameter == 'HeIII_HII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         y_PPI4686                                   = 0.081 * umath.pow(t4, 0.15) * data_dict['He2_4686A']  /  data_dict['H1_4861A']
#                         Magnitude                                   = y_PPI4686                                            
#                 
#                     if Parameter == 'HeII_HII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         ne                                          = data_dict['Den']
#                                         
#                         y0_I4471                                    = 2.039   * umath.pow(t4,0.102) * data_dict['He1_4472A']  / data_dict['H1_4861A']
#                         y0_I5876                                    = 0.732  * umath.pow(t4,0.171) * data_dict['He1_5876A']  / data_dict['H1_4861A']
#                         y0_I6678                                    = 2.573   * umath.pow(t4,0.18) * data_dict['He1_6678A']  / data_dict['H1_4861A']
#                         
#                         D                                           = 1.0 + 3110.0 * umath.pow(t4,-0.51) * (1.0/ne)
#                         g_I4471                                     = 6.11 * umath.pow(t4,0.02)  * (umath.exp(-4.544)) / D
#                         g_I5876                                     = (7.12   * umath.pow(t4,0.14)  * (umath.exp(-3.776/t4)) + 1.47 * umath.pow(t4,-0.28) * (umath.exp(-4.544/t4))) / D
#                         g_I6678                                     = (3.27   * umath.pow(t4,-0.41) * (umath.exp(-3.777/t4)) + 0.49 * umath.pow(t4,-0.52) * (umath.exp(-4.544/t4))) / D
# 
#                         y_I4471                                     = y0_I4471 / (1.0 + g_I4471)
#                         y_I5876                                     = y0_I5876 / (1.0 + g_I5876)
#                         y_I6678                                     = y0_I6678 / (1.0 + g_I6678)
#                     
#                         Magnitude                                   = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)                
#                 
#                 elif methodology == 'Liu2000':
#                     
#                     if Parameter == 'O2_RecombCorr':
#                         t4 = data_dict['Temp'] / 10000  
#                         OII_HII                                     = data_dict['OII_HII_7319A']
#                         OII_RecombCor                               = 9.36 * umath.pow(t4, 0.44) * OII_HII
#                         Magnitude                                   = OII_RecombCor
#  
#                     if Parameter == 'N2_RecombCorr':
#                         t4 = data_dict['Temp'] / 10000  
#                         NIII_HII                                    = data_dict['NIII_HII']
#                         NII_RecombCor                               = 3.19 * umath.pow(t4, 0.30) * NIII_HII
#                         Magnitude                                   = NII_RecombCor                        
#                 
#                 elif methodology == 'Epm2007':
#                 
#                     if Parameter == 'O2':
#                         ne                                          = data_dict['Den']
#                         t4                                          = data_dict['Temp'] / 10000
#                         O2plus_Hplus                                = data_dict['OIII_HII']
#                         I7320_IBeta_R                               = 9.36 * umath.pow(t4, 0.44) * O2plus_Hplus
#                         logOII_logHII                               = -12 + umath.log10((data_dict['O2_7319A'] + data_dict['O2_7330A']) / data_dict['H1_4861A'] - I7320_IBeta_R) + 6.895 + 2.44/t4 - 0.58*umath.log10(t4) - umath.log10(1.0 + 0.0047 * ne)
#                         Magnitude                                   = umath.pow(10, logOII_logHII)
#                                 
#                 elif methodology == 'Fabian2006':
#                 
#                     if Parameter == 'HeII_HII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         ne                                          = data_dict['Den']                
#                 
#                         y0_I4471                                    = 2.04   * umath.pow(t4,0.13) * data_dict['He1_4472A']  / data_dict['H1_4861A']
#                         y0_I5876                                    = 0.738  * umath.pow(t4,0.23) * data_dict['He1_5876A']  / data_dict['H1_4861A']
#                         y0_I6678                                    = 2.58   * umath.pow(t4,0.25) * data_dict['He1_6678A']  / data_dict['H1_4861A']
#                         
#                         D                                           = 1.0 + 3110.0 * umath.pow(t4,-0.51) * (1.0/ne)
#                         g_I4471                                     = 6.11 * umath.pow(t4,0.02)  * (umath.exp(-4.544)) / D
#                         g_I5876                                     = (7.12   * umath.pow(t4,0.14)  * (umath.exp(-3.776/t4)) + 1.47 * umath.pow(t4,-0.28) * (umath.exp(-4.544/t4))) / D
#                         g_I6678                                     = (3.27   * umath.pow(t4,-0.41) * (umath.exp(-3.777/t4)) + 0.49 * umath.pow(t4,-0.52) * (umath.exp(-4.544/t4))) / D
# 
#                                 
#                         y_I4471                                     = y0_I4471 / (1.0 + g_I4471)
#                         y_I5876                                     = y0_I5876 / (1.0 + g_I5876)
#                         y_I6678                                     = y0_I6678 / (1.0 + g_I6678)
#                     
#                         Magnitude                                   = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)                
#                 
#                     elif Parameter == 'HeIII_HII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         y_PPI4686                                   = 0.084 * umath.pow(t4, 0.14) * data_dict['He2_4686A']  /  data_dict['H1_4861A']
#                         Magnitude                                   = y_PPI4686
#                 
#                     elif Parameter == 'OII_HII_7319A':
#                         t4 = data_dict['Temp'] / 10000
#                         ne = data_dict['Den']                        
#                         logOII_logHII                               = -12 + umath.log10((data_dict['O2_7319A'] + data_dict['O2_7330A']) / data_dict['H1_4861A']) + 6.895 + 2.44 / t4 - 0.58 * umath.log10(t4) - umath.log10(1 + 0.0047 * ne)
#                         Magnitude                                   = umath.pow(10, logOII_logHII)
#                         
#                 elif methodology == 'PyNeb':
#                         
#                     if Parameter == 'HeII_HII':
#                         
#                         te                                          = data_dict['Temp']
#                         ne                                          = data_dict['Den']
#                                         
#                         self.He1.getIonAbundance()
# 
#                         y0_I4471                                    = self.He1_atom.getIonAbundance(int_ratio = 100 * data_dict['He1_4472A']/data_dict['H1_4861A'] , tem=te, den=ne, label='He1_4472A')
#                         y0_I5876                                    = self.He1_atom.getIonAbundance(int_ratio = 100 * data_dict['He1_5876A']/data_dict['H1_4861A'] , tem=te, den=ne, label='He1_5876A')
#                         y0_I6678                                    = self.He1_atom.getIonAbundance(int_ratio = 100 * data_dict['He1_6678A']/data_dict['H1_4861A'] , tem=te, den=ne, label='He1_6678A')
#                                                 
# #                         y0_I4471                                    = 2.04   * umath.pow(t4,0.13) * data_dict['He1_4472A']  / data_dict['H1_4861A']
# #                         y0_I5876                                    = 0.738  * umath.pow(t4,0.23) * data_dict['He1_5876A']  / data_dict['H1_4861A']
# #                         y0_I6678                                    = 2.58   * umath.pow(t4,0.25) * data_dict['He1_6678A']  / data_dict['H1_4861A']
#                         
#                         D                                           = 1.0 + 3110.0 * umath.pow(t4,-0.51) * (1.0/ne)
#                         g_I4471                                     = 6.11 * umath.pow(t4,0.02)  * (umath.exp(-4.544)) / D
#                         g_I5876                                     = (7.12   * umath.pow(t4,0.14)  * (umath.exp(-3.776/t4)) + 1.47 * umath.pow(t4,-0.28) * (umath.exp(-4.544/t4))) / D
#                         g_I6678                                     = (3.27   * umath.pow(t4,-0.41) * (umath.exp(-3.777/t4)) + 0.49 * umath.pow(t4,-0.52) * (umath.exp(-4.544/t4))) / D
# 
#                                 
#                         y_I4471                                     = y0_I4471 / (1.0 + g_I4471)
#                         y_I5876                                     = y0_I5876 / (1.0 + g_I5876)
#                         y_I6678                                     = y0_I6678 / (1.0 + g_I6678)
#                     
#                         Magnitude                                   = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)
#                         
#                         
#         return Magnitude
#                       
#     def check_issues(self, magnitude, parameter_type):
#         
#         #Security check in the case we are getting a negative density:
#         if parameter_type == 'density':
#             if magnitude != None:
#                 
#                 if magnitude <= 0:
#                     return ufloat(50.0, 50.0)            
#                 else:
#                     return magnitude
#             return magnitude
#         
#         #Security check in case some flux is missing
#         elif parameter_type == 'EmFlux':
#             if  magnitude[1] == None:
#                 return None
#             else:
#                 return ufloat(magnitude[1], magnitude[2]) 
#     
# class Chemical_Analysis(direct_Method):
# 
#     def __init__(self):
#         
#         direct_Method.__init__(self)
#         
#         #logSI_OI_Gradient = ufloat(-1.53, 0.05)
#         logSI_OI_Gradient = ufloat(-1.78, 0.03)
#                 
#         self.OI_SI      = umath.pow(10, -logSI_OI_Gradient)
#         self.He1        = RecAtom('He', 1)
#         atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
#         
#     def set_element(self, element):
#         
#         if element == 'Argon':
#             
#             self.argon_abundance_scheme()
#             
#         elif element == 'Oxygen':
#             
#             self.oxygen_abundance_scheme()
#             
#         elif element == 'Nitrogen':
#         
#             self.nitrogen_abundance_scheme()
#         
#         elif element == 'Sulfur': 
# 
#             self.sulfur_abundance_scheme()
# 
#         elif element == 'Helium':
# 
#             self.helium_abundance_scheme()
# 
#         return
# 
#     def argon_abundance_scheme(self):
# 
#         #Determine the sulfur density and temperature    
#         T_SII,  T_SIII, ne_SII  = self.sulfur_density_scheme()
#         T_OIII, T_OII           = self.oxygen_temperature_scheme()
#         
#         #We try to calculate the T_ArIII from the sulfur lines
#         T_ArIII                 = T_SIII
#         T_ArIV                  = T_OIII
#         
#         #We try to calculate the T_ArIV from the sulfur lines, if not we use the Oxygen ones
#         if  (T_OIII == None) and (T_SIII != None):
#             data_dict           = {'Temp' : T_SIII}
#             T_OIII_approx       = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion='O3', Parameter = 'TOIII_approx_TSIII')
#             T_ArIV              = T_OIII_approx
# 
#         elif (T_SIII == None) and (T_OIII != None):
#             data_dict           = {'Temp' : T_OIII}
#             T_SIII_approx       = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion='S3', Parameter = 'TSIII_approx_TOIII')
#             T_ArIV              = T_SIII_approx
#             
#         #Calculate the ArIII abundance
#         self.argon_IonAbundance(T_ArIII, ne_SII, Ion = 'Ar3', methodology = 'Haegele2008')
#   
#         #Calculate the ArIV abundance
#         self.argon_IonAbundance(T_ArIV, ne_SII, Ion = 'Ar4', methodology = 'Haegele2008')  
#    
#     def sulfur_density_scheme(self, TempIn = None):
#         
#         #We use this set up mechanism to calculate the electron density for other elements 
#         if TempIn == None:
#             #Determine the Te[SII] and Te[SIII]
#             T_SII   = self.temperature_determination(Ion = 'S2', methodology = 'Haegele2008')
#             T_SIII  = self.temperature_determination(Ion = 'S3', methodology = 'Epm2014')
#         else:
#             T_SII   = None
#             T_SIII  = TempIn
#                 
#         #Determine ne_SII using the TSIII. If not available use TOIII to calculate TSIII
#         #If a TempIn is not available it will use the standard procedure 
#         if T_SIII  != None:
#             ne_SII = self.density_determination(Ion = 'S2', Temp = T_SIII, methodology='Epm2014')
#         else:
#             T_OIII, T_OII       = self.oxygen_temperature_scheme()
#             data_dict           = {}
#             data_dict['Temp']   = T_OIII
#             T_SIII              = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion = 'S3', Parameter='TSIII_approx_TOIII')
#             ne_SII              = self.density_determination(Ion = 'S2', Temp = T_SIII, methodology='Epm2014')
#         
#         return T_SII, T_SIII, ne_SII
#     
#     def sulfur_abundance_scheme(self):
#         
#         #Get the sulfur electron density and temperature
#         T_SII, T_SIII, ne_SII = self.sulfur_density_scheme(TempIn = None)
#         
#         self.Properties_dict['nSII'] = ne_SII
#                            
#         #Determine the SII/HII abundance. 
#         #If the T_SII temperature is available use it. If not T_SII = T_SIII:
#         if T_SII != None:
#             self.sulfur_IonAbundance(T_SII, ne_SII, Ion = 'S2', methodology = 'Haegele2008')
#         else:
#             self.sulfur_IonAbundance(T_SIII, ne_SII, Ion = 'S2', methodology = 'Haegele2008')
# 
#         #Determine the SIII/HII abundance
# #         self.sulfur_IonAbundance(T_SIII, ne_SII, Ion = 'S3', methodology = 'Haegele2008')
#         self.sulfur_IonAbundance(T_SIII, ne_SII, Ion = 'S3', methodology = 'Vital2015')
# 
#         #Determine the S/H abundance
#         #Include the Argon correction if the ArIV lines were observed
#         if (self.Properties_dict['SII_HII'] != None) and (self.Properties_dict['SIII_HII'] != None):
#             self.Properties_dict['SI_HI'] = self.Properties_dict['SII_HII'] + self.Properties_dict['SIII_HII']
#                         
#             #Calculate the [SIV] correction factor from the Argon abundance
#             #This structures enforces the calculation of the Argon apriori
#             data_dict               = {}
#             data_dict['ArIII_HII']  = self.Properties_dict['ArIII_HII']   
#             data_dict['ArIV_HII']   = self.Properties_dict['ArIV_HII']   
#             data_dict['SIII_HII']   = self.Properties_dict['SIII_HII']               
#             
#             #The code does not work if we introduce a None entry. However, in this correction there is still a quantity even if the ArIV is no there
#             if data_dict['ArIV_HII'] == None:
#                 data_dict['ArIV_HII'] = ufloat(0.0,0.0) 
#             
#             #Compute the SIV_HII component
#             SIV_HII = self.empiric_formulae(data_dict, methodology = 'Angeles2015', Ion = 'S4', Parameter='SIV_HII')         
#             
#             self.Properties_dict['SI_HI_ArCorr'] = self.Properties_dict['SII_HII'] + self.Properties_dict['SIII_HII'] + SIV_HII
#              
#     def oxygen_temperature_scheme(self):
# 
#         #Determine the Te[OIII]
#         T_OIII          = self.temperature_determination(Ion = 'O3', methodology = 'Epm2014')
#         
#         #Determine the Te[OII] using all the [OII] lines
#         T_OII           = self.temperature_determination(Ion = 'O2', methodology = 'Epm2014')
#         
#         #If lines not observed use approximation from T[OIII]
#         if T_OII == None:
#             T_OII       = self.temperature_determination(Ion = 'O2', methodology = 'Epm2014', Temp = T_OIII)
# 
#         return T_OIII, T_OII
# 
#     def oxygen_abundance_scheme(self):
#         
#         #Determine the oxygen temperatures
#         T_OIII, T_OII = self.oxygen_temperature_scheme()
#         
#         #Get the electron density from the sulfur lines using the oxigen temperature if observed. If not used sulfur lines
#         T_SII, T_SIII, ne_SII = self.sulfur_density_scheme()
#         
#         #Calculate the OIII ion abundance
#         self.oxygen_IonAbundance(Temp = T_OIII, Den = ne_SII, Ion = 'O3', methodology = 'Epm2014')
#         
#         #Calculate the OII ion abundance
#         self.oxygen_IonAbundance(Temp = T_OII, Den = ne_SII, Ion = 'O2', methodology = 'Epm2014')
#         self.oxygen_IonAbundance(Temp = T_OII, Den = ne_SII, Ion = 'O2', methodology = 'Fabian2006')
#         
#         #Correct the OII from 7319A, 7330A lines from the recombination contribution
#         if self.Properties_dict['OII_HII_7319A'] != None:
#             OII_HII_7319A   = self.Properties_dict['OII_HII_7319A']
#             data_dict       = {'Temp'           : T_OII, 
#                                'OII_HII_7319A'  : OII_HII_7319A}
#             
#             RecCorr_O2_7319 = self.empiric_formulae(data_dict, methodology = 'Liu2000', Ion = 'O2', Parameter='O2_RecombCorr')
#             self.oxygen_IonAbundance(Temp = T_OII, Den = ne_SII, Ion = 'O2', methodology = 'Epm2014', RecCorr = RecCorr_O2_7319)    
#         
#         #If O2_3726A, O2_3729A lines are not available we use O2_7319A, O2_7330A lines
#         if self.Properties_dict['OII_HII_3279A'] != None:
#             self.Properties_dict['OII_HII'] = self.Properties_dict['OII_HII_3279A']
#             
#         elif (self.Properties_dict['OII_HII_7319A'] != None): 
#             self.Properties_dict['OII_HII']     = self.Properties_dict['OII_HII_7319A']
#         
#         #Determine the S/H abundance
#         #Include the Argon correction if the ArIV lines were observed
#         if (self.Properties_dict['OII_HII'] != None) and (self.Properties_dict['OIII_HII'] != None):
#             self.Properties_dict['OI_HI'] = self.Properties_dict['OII_HII'] + self.Properties_dict['OIII_HII']
#             
#     def nitrogen_abundance_scheme(self):
#         #WARNING: THIS COULD BE FASTER WITH A MECHANISM WITH PRELOADS ALL THE EMISSION LINES
#         
#         #Determine the Te[NSII]
#         T_NII                   = self.temperature_determination(Ion = 'N2', methodology = 'Epm2014')
#                 
#         if T_NII != None:
#             #YOU CAN TEST THIS WITH OBJECT 51959-092
#             T_SII, T_SIII, ne_SII   = self.sulfur_density_scheme(TempIn = T_NII)
#             self.nitrogen_IonAbundance(Temp = T_NII, Den = ne_SII, Ion = 'N2', methodology = 'Epm2014')
#             
#         #If Te[NSII] cannot be calculated try to determine the Te[OII]
#         elif T_NII == None:
#             T_OIII, T_OII       = self.oxygen_temperature_scheme()            
#             data_dict           = {'Temp' : T_OIII}
#             TNII_approx         = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion='N2', Parameter = 'TNII_approx_TOIII')
#         
#             #Get the sulfur electron density and temperature
#             T_SII, T_SIII, ne_SII   = self.sulfur_density_scheme(TempIn = TNII_approx)
#             self.nitrogen_IonAbundance(Temp = TNII_approx, Den = ne_SII, Ion = 'N2', methodology = 'Epm2014')
# 
#         #Determine the N/H abundance
#         #We need O/H abundance in order to calculate N/H abundance
#         if (self.Properties_dict['OI_HI'] != None) and (self.Properties_dict['NII_HII'] != None):
#             self.Properties_dict['NI_OI'] = self.Properties_dict['NII_HII'] / self.Properties_dict['OII_HII']
#             self.Properties_dict['NI_HI'] = self.Properties_dict['NI_OI'] * self.Properties_dict['OI_HI']
#             
#             #Take into account the recombination contribution from the 
#             if T_NII != None:
#                 NIII_HI = self.Properties_dict['NI_HI'] - self.Properties_dict['NII_HII']
#                 
#                 data_dict       = {'Temp'           : T_NII, 
#                                    'NIII_HII'       : NIII_HI}
#             
#                 RecCorr_N2_5755 = self.empiric_formulae(data_dict, methodology = 'Liu2000', Ion = 'N2', Parameter = 'N2_RecombCorr')
#                 
#                 #Repeat the previous steps starting with the calculation of TNII
#                 T_NII   = self.temperature_determination(Ion = 'N2', methodology = 'Epm2014', RecCorr = RecCorr_N2_5755)
#                 
#                 T_SII, T_SIII, ne_SII   = self.sulfur_density_scheme(TempIn = T_NII)
#                 self.nitrogen_IonAbundance(Temp = T_NII, Den = ne_SII, Ion = 'N2', methodology = 'Epm2014')
#                 self.Properties_dict['NI_OI'] = self.Properties_dict['NII_HII'] / self.Properties_dict['OII_HII']
#                 self.Properties_dict['NI_HI'] = self.Properties_dict['NI_OI'] * self.Properties_dict['OI_HI']
#                                        
#     def helium_abundance_scheme(self):
#         
#         #Determine the sulfur density and temperature    
#         T_SII, T_SIII, ne_SII = self.sulfur_density_scheme()
#     
#         #In case we do not have the sulfur temperature we use the oxygen one
#         if T_SIII != None:
#             data_dict   = {'Temp' : T_SIII}
#             T_OIII      = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion = 'O3', Parameter = 'TOIII_approx_TSIII')
#         
#         else:
#             T_OIII, T_OII = self.oxygen_temperature_scheme()
#         
#         #Calculate the HeII ion  abundance
#         self.helium_IonAbundance(Temp = T_OIII, Den = ne_SII , Ion = 'He1', methodology = 'Vital2015')
# 
#         #Calculate the HeIII ion abundance
#         self.helium_IonAbundance(Temp = T_OIII, Den = ne_SII , Ion = 'He2', methodology = 'Vital2015')
# 
#         if self.Properties_dict['HeII_HII'] != None:
#             if self.Properties_dict['HeIII_HII'] != None:
#                 self.Properties_dict['HeI_HI'] = self.Properties_dict['HeII_HII'] + self.Properties_dict['HeIII_HII']
#             elif self.Properties_dict['HeIII_HII'] == None:
#                 self.Properties_dict['HeI_HI'] = self.Properties_dict['HeII_HII']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.Properties_dict['HeI_HI'] != None):
#             
#             #From Oxygen:
#             if self.Properties_dict['OI_HI'] != None:
#                 self.Properties_dict['Y_Mass_O'] = (4 * self.Properties_dict['HeI_HI'] * (1 - 20 * self.Properties_dict['OI_HI'])) / (1 + 4 * self.Properties_dict['HeI_HI']) 
#             
#             #From Sulfur
#             if self.Properties_dict['SI_HI'] != None:
#                 self.Properties_dict['Y_Mass_S'] = (4 * self.Properties_dict['HeI_HI'] * (1 - 20 * self.OI_SI * self.Properties_dict['SI_HI_ArCorr'])) / (1 + 4 * self.Properties_dict['HeI_HI']) 
# 
#     def heliumPyneb_abundance_scheme(self):
#         
#         #Determine the sulfur density and temperature    
#         T_SII, T_SIII, ne_SII = self.sulfur_density_scheme()
#     
#         #In case we do not have the sulfur temperature we use the oxygen one
#         if T_SIII != None:
#             data_dict   = {'Temp' : T_SIII}
#             T_OIII      = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion = 'O3', Parameter = 'TOIII_approx_TSIII')
#         
#         else:
#             T_OIII, T_OII = self.oxygen_temperature_scheme()
#         
#         #Calculate the HeII ion  abundance
#         self.helium_IonAbundance(Temp = T_OIII, Den = ne_SII , Ion = 'He1', methodology = 'PyNeb')        
#         
#         #Calculate the HeIII ion abundance
#         self.helium_IonAbundance(Temp = T_OIII, Den = ne_SII , Ion = 'He2', methodology = 'PyNeb')
#                 
#         return 
#            
#     def store_dict_ObjLog(self, varaiables_dict, verbose = False):
#         
#         keys, values = varaiables_dict.keys(), varaiables_dict.values()
#         
#         for i in range(len(keys)):
#             if verbose:
#                 print i, keys[i], values[i]
#             self.SaveParameter_ObjLog(self.CodeName, self.FileFolder, keys[i], values[i])
#         
#         return    
#     
# class Parametrized_Emissivities(Txt_Files_Manager):
#     
#     def __init__(self):
# 
#         self.LiteratureParametrization = None
#         Txt_Files_Manager.__init__(self)
# 
#     def AOS2010(self, Te, ne):
#         
#         #phi has the shape : phi_lambda = Em(Hbeta) / E(He_lambda) * (1) / (1 + C/R(lambda)) 
#         #Hence Em(He_lambda)/Em(H\beta) = 1 / phi_lambda
#                 
#         Phi3889 = power(0.8779 * Te, -0.128 -0.00041 * ne )
#         Phi4026 = power(4.233 * Te, 0.085 - 0.00012 * ne )
#         Phi4471 = power(2.021 * Te, 0.121 - 0.00020 * ne )
#         Phi5876 = power(0.754 * Te, 0.212 - 0.00051 * ne )
#         Phi6678 = power(2.639 * Te, 0.244 - 0.00054 * ne )
#         Phi7065 = power(5.903 * Te, -0.519) / (1.462 - Te * (0.127 - 0.00076 * ne * + 0.000000255 * ne*ne))
#         
#         Em_He3889 = 1/Phi3889
#         Em_He4026 = 1/Phi4026
#         Em_He4471 = 1/Phi4471
#         Em_He5876 = 1/Phi5876
#         Em_He6678 = 1/Phi6678
#         Em_He7065 = 1/Phi7065
#         
#         return Em_He3889, Em_He4026, Em_He4471, Em_He4471, Em_He5876, Em_He5876, Em_He6678, Em_He7065
#     
#     def Haegele2008_Fluxes(self):
#         #SDSS J1729
#         Flux_dict = {}
#         
#         Flux_dict["H1_4861A"]   = ufloat(10000, 83)
#         
#         Flux_dict["O2_3726A"]   = ufloat(17622, 243)
#         Flux_dict["O2_3729A"]   = ufloat(0.0, 0.0)
#         Flux_dict["O2_7319A"]   = ufloat(233,8)
#         Flux_dict['O2_7330A']   = ufloat(194,9)        
#         
#         Flux_dict["O3_4363A"]   = ufloat(660, 30)
#         Flux_dict["O3_4959A"]   = ufloat(17097, 146)
#         Flux_dict["O3_5007A"]   = ufloat(51541, 418)        
#     
#         Flux_dict["N2_5755A"]   = ufloat(60, 6)
#         Flux_dict["N2_6548A"]   = ufloat(771, 21)
#         Flux_dict["N2_6584A"]   = ufloat(2189, 52)        
#         
#         Flux_dict['Ar3_7136A']  = ufloat(864,41)
#         Flux_dict['Ar4_4740A']  = ufloat(40,0.06)  
#         
#         Flux_dict["S2_6716A"]   = ufloat(1286, 32)
#         Flux_dict["S2_6731A"]   = ufloat(1009, 30)
#         Flux_dict["S2_4069A"]   = ufloat(122, 11)
#         Flux_dict['S2_4076A']   = ufloat(0.0,0.0)
#         
#         Flux_dict["S3_6312A"]   = ufloat(176, 11)
#         Flux_dict["S3_9069A"]   = ufloat(2092, 95)
#         Flux_dict["S3_9531A"]   = ufloat(4718, 250)
#          
#         Flux_dict["He1_4472A"]   = ufloat(516, 24)
#         Flux_dict["He1_5876A"]   = ufloat(1261, 37)
#         Flux_dict["He1_6678A"]   = ufloat(355, 2)
#         Flux_dict['He2_4686A']   = ufloat(329,43)         
# 
#     def Abundances_FabianSample(self):
#         
#         Table_Address   = '/home/vital/Workspace/Dazer/Astro_Libraries/Fabian2006_CatalogueAbundances.txt'
# 
#         Headers_List    = ['O/H_10^5_Flux', 'O/H_0^5_Error', 'N/H_10^6_Flux', 'N/H_10^6_Error', 'Y_magnitude', 'Y_error', 'y_magnitude', 'y_error']
#         
#         Oxygen_Flux, Oxygen_Error, Nitrogen_Flux, Nitrogen_Error, Y, Yerror, y_magnitude, y_error = self.get_ColumnData(Headers_List, Table_Address, HeaderSize=1, unpack_check = True)
#         
#         return Oxygen_Flux * 1e-5, Oxygen_Error*1e-5, Nitrogen_Flux*1e-6, Nitrogen_Error*1e-6, Y, Yerror,  y_magnitude, y_error     
    
# def residual_Y_v3(pars, x, y):
#     return (y - pars['Y'] * x)
# 
# class Chemical_Analysis_pyneb():
# 
#     def __init__(self):
#         
#         self.MC_array_len           = 100
#         
#         self.Hbeta_label            = 'H1_4861A'
#     
#     def load_elements(self):
# 
#         #Set atomic data 
#         atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
#         atomicData.setDataFile('s_iii_coll_HRS12.dat')
#  
#         #Declare hydrogen atom so it becomes readily available:
#         self.S2_atom            = Atom('S', 2)
#         self.S3_atom            = Atom('S', 3)
#         self.Ar3_atom           = Atom('Ar', 3)
#         self.Ar4_atom           = Atom('Ar', 4)
#         self.N2_atom            = Atom('N', 2)
#         self.O2_atom            = Atom('O', 2)
#         self.O3_atom            = Atom('O', 3)
#                     
#         self.H1_atom            = RecAtom('H', 1)
#         self.He1_atom           = RecAtom('He', 1)
#         self.He2_atom           = RecAtom('He', 2)
# 
#         self.diags              = Diagnostics()
#         
#         self.logSI_OI_Gradient  = random.normal(-1.53,  0.05, size = self.MC_array_len) # random.normal(-1.78,  0.03, size = self.MC_array_len)
# 
#         self.OI_SI              = power(10, -self.logSI_OI_Gradient)
#         
#         self.S3_9000_ratio      = random.normal(self.S3_atom.getEmissivity(10000, 1000, wave = 9531) / self.S3_atom.getEmissivity(10000, 1000, wave = 9069),  0.01, size = self.MC_array_len)
#         
#         self.m_SIV_correction   = random.normal(1.109,  0.1109, size = self.MC_array_len)
#         self.n_SIV_correction   = random.normal(0.135,  0.013, size = self.MC_array_len)
#         
#         lower_trunc, upper_trunc = (0.0 - 50.0) / 25.0, (100 - 50) / 25.0
#         self.Truncated_gaussian  = truncnorm(lower_trunc,  upper_trunc, loc = 50, scale = 25)
#        
#         return
#             
#     def DeclareObject(self, lines_log_frame):
#         
#         #List of all parameters
#         lineRatios      = ['R_SII', 'R_SII_prime', 'R_SIII', 'R_NII', 'R_OII', 'R_OII_prime', 'R_OIII']
#         elecProperties  = ['neSII', 'neOII', 'TeOII', 'TeSII', 'TeNII', 'TeOIII', 'TeSIII', 'TeOII_from_TeOIII', 'TeNII_from_TeOIII', 'TeSIII_from_TeOIII', 'TeOIII_from_TeSIII']        
#         ionicAbund      = ['SII_HII', 'SIII_HII', 'SIV_HII', 'OII_HII', 'OII_HII_3279A', 'OII_HII_7319A', 'NII_HII', 'ArIII_HII', 'ArIV_HII', 'HeII_HII_from_O', 'HeIII_HII_from_O', 'HeII_HII_from_S', 'HeIII_HII_from_S']
#         elemAbund       = ['SI_HI', 'OI_HI', 'NI_OI', 'NI_HI', 'HeI_HI_from_O', 'HeI_HI_from_S', 'Ymass_O', 'Ymass_S']
#     
#         self.data    = OrderedDict.fromkeys(lineRatios + elecProperties + ionicAbund + elemAbund)
#                 
#         self.Hbeta_flux         = random.normal(lines_log_frame.loc['H1_4861A']['line_Int'].nominal_value, lines_log_frame.loc['H1_4861A']['line_Int'].std_dev, size = self.MC_array_len)
#               
#         self.low_density_dist   = self.Truncated_gaussian.rvs(self.MC_array_len)
# 
#         #Generate a dictionary to store the random array for all lines
#         self.lines_dict = OrderedDict()
# 
#         #Generate
#         for line in lines_log_frame.index.values:
#             if line == 'O2_3726A': 
#                 if lines_log_frame.loc['O2_3726A']['flux_intg'] == lines_log_frame.loc['O2_3726A']['flux_intg']:
#                     flux, error = lines_log_frame.loc[line]['line_IntBrute_dered'].nominal_value, lines_log_frame.loc[line]['line_IntBrute_dered'].std_dev
# 
#             elif line == 'O2_7319A': 
#                 if lines_log_frame.loc['O2_7319A']['flux_intg'] == lines_log_frame.loc['O2_7330A']['flux_intg']:
#                     flux, error = lines_log_frame.loc[line]['line_IntBrute_dered'].nominal_value, lines_log_frame.loc[line]['line_IntBrute_dered'].std_dev
#                 else:
#                     line_sum = lines_log_frame.loc['O2_7319A']['line_IntBrute_dered'] + lines_log_frame.loc['O2_7330A']['line_IntBrute_dered']
#                     flux, error = line_sum.nominal_value, line_sum.std_dev
#                 
#             else:
#                 flux, error         = lines_log_frame.loc[line]['line_Int'].nominal_value, lines_log_frame.loc[line]['line_Int'].std_dev
#             
#             self.lines_dict[line]   = random.normal(flux, error, size = self.MC_array_len)
#             
#         return
#                               
#     def low_density_function(self, density_array):
#         
#         inds_nan    = where(isnan(density_array))
#         
#         if len(inds_nan[0]) > 0:
#             
#             if len(inds_nan[0]) > (0.1 * self.MC_array_len):
#                 print '--WARNING 10% boundary: ', len(inds_nan[0]), 'values with nan'            
#             elif len(inds_nan[0]) > (0.05 * self.MC_array_len):
#                 print '--WARNING 5% boundary: ', len(inds_nan[0]), 'values with nan'            
#                                 
#             density_array[inds_nan]= self.low_density_dist[inds_nan]
#             
#         mean_density    = mean(density_array) 
#         std_density     = std(density_array)   
#         
#         if isnan(mean_density): 
#             print '--WARNING: Density array cannot be treated numerically'
#             
#         if mean_density < 75:
#             print '--WARNING: Low density', mean_density, '+/-', std_density, 'Points below limit', len(where(density_array < 75)[0])
#             return  self.low_density_dist
#         else:                              
#             return density_array
#             
#     def determine_electron_parameters(self):
#                 
# #-------Starting with the Sulfur
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:    
#               
#             #Check TSIII lines
#             if self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_9531A', 'S3_6312A'}:
#         
#                 self.Properties_dict['TSIII_pn'], self.Properties_dict['nSII_pn'] = self.diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+',
#                                                                                     diag_den = '[SII] 6731/6716',
#                                                                                     value_tem = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] + self.lines_dict['S3_9531A']),
#                                                                                     value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']) 
#                                 
#             #Missing S3_9531A line
#             elif self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_6312A'}:
#                 
#                 self.Properties_dict['TSIII_pn'], self.Properties_dict['nSII_pn'] = self.diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+',
#                                         diag_den = '[SII] 6731/6716',
#                                         value_tem = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] * (1 + self.S3_9000_ratio)),
#                                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']) 
# 
#             #Missing S3_9069A line
#             elif self.lines_dict.viewkeys() >= {'S3_9531A', 'S3_6312A'}:
#                 self.Properties_dict['TSIII_pn'], self.Properties_dict['nSII_pn'] = self.diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+',
#                                         diag_den = '[SII] 6731/6716',
#                                         value_tem = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9531A'] * (1 + 1/self.S3_9000_ratio)),
#                                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#             
#             #Missing the sulfur lines we use the ones from oxygen
#             elif self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A', 'O3_4363A'}:         
#                 
#                 self.Properties_dict['TOIII_pn'], nSII = self.diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+',
#                         diag_den = '[SII] 6731/6716',
#                         value_tem = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A']),
#                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#             
#                 self.Properties_dict['TSIII_approxfrom_TOIII_pn'] = ((self.Properties_dict['TOIII_pn']/10000 + 0.0846) / 1.0807) * 10000
#                 self.Properties_dict['nSII_pn'] = self.S2_atom.getTemDen((self.lines_dict['S2_6731A'])/(self.lines_dict['S2_6716A']), tem = self.Properties_dict['TSIII_approxfrom_TOIII_pn'], to_eval = 'L(6731) / L(6716)') 
#                 
#             
#             #Try to measure the TSII temperature
#             if self.lines_dict.viewkeys() >= {'S2_4069A', 'S2_4076A'}:
# 
#                 self.Properties_dict['TSII_pn'], nSII = self.diags.getCrossTemDen(diag_tem = '[SII] 4069/4076',
#                         diag_den  = '[SII] 6731/6716',
#                         value_tem = self.lines_dict['S2_4069A']/self.lines_dict['S2_4076A'],
#                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#             
# #-------Following with the oxygen
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:
#                 
#             #Missing the sulfur lines we use the ones from oxygen
#             if self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A', 'O3_4363A'}:         
#                 self.Properties_dict['TOIII_pn'], nSII_x = self.diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+',
#                         diag_den  = '[SII] 6731/6716',
#                         value_tem = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A']),
#                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#                         
#             #Determine the TOIII from the sulfur temperature as a proxy
#             if self.Properties_dict['TSIII_pn'] is not None:
#                 
#                 self.Properties_dict['TOIII_approxfrom_TSIII_pn'] = (1.0807 * self.Properties_dict['TSIII_pn']/10000 - 0.0846) * 10000
#             
#             #Determine the TOII value from the TOIII proxy
#             if self.Properties_dict['TOIII_pn'] is not None:
#                 self.Properties_dict['TOII_approxfrom_TOIII_pn'] = (1.397 /( (1/(self.Properties_dict['TOIII_pn']/10000)) + 0.385) ) * 10000
# 
# #             #Calculate the oxygen density if the 3721A+ lines can be measured        
# #             if self.lines_dict.viewkeys() >= {'O2_3726A', 'O2_3729A', 'O2_7319A', 'O2_7330A', 'O3_4959A', 'O3_5007A', 'O3_4363A'}:
# # 
# #                 self.Properties_dict['TOIII_pn'], self.Properties_dict['nOII_pn'] = self.diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+',
# #                                                                                         diag_den = '[OII] 3726/3729',
# #                                                                                         value_tem = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A']),
# #                                                                                         value_den = (self.lines_dict['O2_3726A']/self.lines_dict['O2_3729A']))
# 
#             #Calculate the oxygen temperature cannot be measure use the SIII approximation       
#             elif self.lines_dict.viewkeys() >= {'O2_3726A', 'O2_3729A', 'O2_7319A', 'O2_7330A'}:
#                 if self.Properties_dict['TOIII_approxfrom_TSIII_pn'] is not None:
#                     self.Properties_dict['nOII_pn'] = self.O2_atom.getTemDen((self.lines_dict['O2_3726A']/self.lines_dict['O2_3729A']), 
#                                                                              tem = self.Properties_dict['TOIII_approxfrom_TSIII_pn'], 
#                                                                              to_eval = '(L(3726))/(L(3729))')
#                 
# #-------Following with the Nitrogen
#         if self.lines_dict.viewkeys() >= {'N2_5755A', 'N2_6548A', 'N2_6584A'}:         
# 
#             self.Properties_dict['TNII_pn'], nSII = self.diags.getCrossTemDen(diag_tem  = '[NII] 5755/6584+',
#                                                                               diag_den  = '[SII] 6731/6716',
#                                                                               value_tem = self.lines_dict['N2_5755A']/(self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']),
#                                                                               value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
# 
#         if self.Properties_dict['TOIII_pn'] is not None:
#             self.Properties_dict['TNII_approxfrom_TOIII_pn'] =  (1.452 /( (1/(self.Properties_dict['TOIII_pn']/10000)) + 0.479)) * 10000
#              
#         return
# 
#     def argon_abundance_scheme(self):
#  
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TSIII_pn'] is not None:
#             TSIII = self.Properties_dict['TSIII_pn']
#             TOIII = self.Properties_dict['TOIII_approxfrom_TSIII_pn']
# 
#         else:
#             TSIII = self.Properties_dict['TSIII_approxfrom_TOIII_pn']
#             TOIII = self.Properties_dict['TOIII_pn']
#              
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = self.Properties_dict['nOII_pn']
#                  
#         if ne is not None:
#         
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)            
#         
#             #Calculate the Ar_+2 abundance
#             if (TSIII is not None):
#                 if self.lines_dict.viewkeys() >= {'Ar3_7136A', 'Ar3_7751'}:
#                     Ratio = self.lines_dict['Ar3_7751A'] + self.lines_dict['Ar3_7136A']
#                     self.Properties_dict['ArIII_HII_pn'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, to_eval = 'L(7751) + L(7136)', Hbeta = self.Hbeta_flux)     
#                 elif self.lines_dict.viewkeys() >= {'Ar3_7136A'}:
#                     Ratio = self.lines_dict['Ar3_7136A']
#                     self.Properties_dict['ArIII_HII_pn'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, wave = 7136, Hbeta = self.Hbeta_flux)             
#                 elif self.lines_dict.viewkeys() >= {'Ar3_7751A'}:
#                     Ratio = self.lines_dict['Ar3_7751A']
#                     self.Properties_dict['ArIII_HII_pn'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, wave = 7751, Hbeta = self.Hbeta_flux)   
# 
#             #Calculate the Ar_+3 abundance
#             if (TOIII is not None):
#                 if self.lines_dict.viewkeys() >= {'Ar4_4740A', 'Ar4_4711A'}:
#                     Ratio = self.lines_dict['Ar4_4740A'] + self.lines_dict['Ar4_4711A']
#                     self.Properties_dict['ArIV_HII_pn'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, to_eval = 'L(4740) + L(4711)', Hbeta = self.Hbeta_flux)                             
#                 elif self.lines_dict.viewkeys() >= {'Ar4_4740A'}:
#                     Ratio = self.lines_dict['Ar4_4740A']
#                     self.Properties_dict['ArIV_HII_pn'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, wave = 4740, Hbeta = self.Hbeta_flux)
#                 elif self.lines_dict.viewkeys() >= {'Ar4_4711A'}:
#                     Ratio = self.lines_dict['Ar4_4711A']
#                     self.Properties_dict['ArIV_HII_pn'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, wave = 4711, Hbeta = self.Hbeta_flux)
#                                           
#     def sulfur_abundance_scheme(self):
#         
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TSIII_pn'] is not None:
#             Te = self.Properties_dict['TSIII_pn']
#         else:
#             Te = self.Properties_dict['TSIII_approxfrom_TOIII_pn']
# 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = self.Properties_dict['nOII_pn']
# 
#         #Calculating the R_S2 ratio
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:
#             self.Properties_dict['R_SII_pn'] = self.lines_dict['S2_6716A'] + self.lines_dict['S2_6731A']
#                 
#         #Calculating the R_S3 ratio
#         if self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_9531A'}:
#             self.Properties_dict['R_SIII_pn'] = self.lines_dict['S3_9069A'] + self.lines_dict['S3_9531A']
#         
#         elif self.lines_dict.viewkeys() >= {'S3_9069A'}:
#             self.Properties_dict['R_SIII_pn'] = self.lines_dict['S3_9069A'] * (1 + self.S3_9000_ratio)
#             
#         elif self.lines_dict.viewkeys() >= {'S3_9531A'}:
#             self.Properties_dict['R_SIII_pn'] = self.lines_dict['S3_9531A'] * (1 + 1/self.S3_9000_ratio)          
#                 
#         #Determine ionic abundances
#         if (Te is not None) and (ne is not None):
#                         
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             
#             #Calculate the S+1 abundance
#             if self.Properties_dict['R_SIII_pn'] is not None:
#                 self.Properties_dict['SII_HII_pn'] = self.S2_atom.getIonAbundance(int_ratio=self.Properties_dict['R_SII_pn'], tem=Te, den=ne, to_eval = 'L(6731)+L(6716)', Hbeta = self.Hbeta_flux)        
#             
#             #Calculate the S_+2 abundance
#             if self.Properties_dict['R_SIII_pn'] is not None:
#                 self.Properties_dict['SIII_HII_pn'] = self.S3_atom.getIonAbundance(int_ratio=self.Properties_dict['R_SIII_pn'], tem=Te, den=ne, to_eval = 'L(9069)+L(9531)', Hbeta = self.Hbeta_flux)
#     
#             #Calculating the S_+3 abundance
#             if (self.Properties_dict['SII_HII_pn'] is not None) and (self.Properties_dict['SIII_HII_pn'] is not None):
#                                 
#                 #Apply correction from Ar - Sulfur from Popstar - Cloudy models if Argon metallicity previously measured     
#                 if (self.Properties_dict['ArIII_HII_pn'] is not None) and (self.Properties_dict['ArIV_HII_pn'] is not None):
#                     logAr2Ar3   = log10(self.Properties_dict['ArIII_HII_pn']/self.Properties_dict['ArIV_HII_pn'])
#                     logSIV  = log10(self.Properties_dict['SIII_HII_pn']) - (logAr2Ar3 - self.n_SIV_correction) / self.m_SIV_correction
#                     self.Properties_dict['SIV_HII_pn'] = power(10, logSIV)
# 
#                 else:
#                     self.Properties_dict['SIV_HII_pn'] = zeros(self.MC_array_len)
#         
#         #Determine elemental abundance
#         if (self.Properties_dict['SII_HII_pn'] is not None) and (self.Properties_dict['SIII_HII_pn'] is not None) and (self.Properties_dict['SIV_HII_pn'] is not None):
#             self.Properties_dict['SI_HI_pn'] = self.Properties_dict['SII_HII_pn'] + self.Properties_dict['SIII_HII_pn'] + self.Properties_dict['SIV_HII_pn']
# 
#         if (self.Properties_dict['SII_HII_pn'] is not None) and (self.Properties_dict['SIII_HII_pn'] is not None): 
#             self.Properties_dict['SI_HI_pn'] = self.Properties_dict['SII_HII_pn'] + self.Properties_dict['SIII_HII_pn']
# 
#         return
#     
#     def oxygen_abundance_scheme(self):
#   
#         #Select the temperature and density for the procedure49
#         if self.Properties_dict['TOIII_pn'] is not None:
#             TOIII   = self.Properties_dict['TOIII_pn']
#             TOII    = self.Properties_dict['TOII_approxfrom_TOIII_pn']   
#         else:
#             TOIII   = self.Properties_dict['TOIII_approxfrom_TSIII_pn']
#             TOII    = self.Properties_dict['TSIII_pn']
#                 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = self.Properties_dict['nOII_pn']    
#      
#         if ne is not None:
#               
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)   
#             
#             #Calculate the O+1 ion abundance (If O2_3726A, O2_3729A are available we prioritize them)
#             if TOII is not None:
#                 
#                 if self.lines_dict.viewkeys() >= {'O2_3726A'}: #WE ASSUME THEY ARE BLENDED
#                     Ratio                                       = self.lines_dict['O2_3726A']
#                     self.Properties_dict['OII_HII_3279A_pn']    = self.O2_atom.getIonAbundance(int_ratio = Ratio, tem=TOII, den=ne, to_eval = 'L(3726)+L(3729)', Hbeta = self.Hbeta_flux) 
#                 elif  self.lines_dict.viewkeys() >= {'O2_3726A'}:
#                     Ratio                                       = self.lines_dict['O2_3726A']
#                     self.Properties_dict['OII_HII_3279A_pn']    = self.O2_atom.getIonAbundance(int_ratio = Ratio, tem=TOII, den=ne, to_eval = 'L(3726)+L(3729)', Hbeta = self.Hbeta_flux) 
#         
#                 if self.lines_dict.viewkeys() >= {'O2_7319A'}: 
#                     #Initial value
#                     Ratio                                       = self.lines_dict['O2_7319A'] #WE ASSUME THEY ARE BLENDED
#                     OII_HII_pn_7319_initial                     = self.O2_atom.getIonAbundance(int_ratio = Ratio, tem=TOII, den=ne, to_eval = 'L(7319)+L(7330)', Hbeta = self.Hbeta_flux)     
#                     #Correction for recombination contribution Liu2000
#                     Lines_Correction                            = (9.36 * power((TOII/10000), 0.44) * OII_HII_pn_7319_initial) * self.Hbeta_flux      
#                     Ratio                                       = self.lines_dict['O2_7319A'] - Lines_Correction
#                     self.Properties_dict['OII_HII_7319A_pn']    = self.O2_atom.getIonAbundance(int_ratio = Ratio, tem=TOII, den=ne, to_eval = 'L(7319)+L(7330)', Hbeta = self.Hbeta_flux) 
#    
#             #Calculate the O+2 ion abundance
#             if (TOIII is not None):         
#                 if self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A'}:
#                     Ratio = self.lines_dict['O3_4959A'] + self.lines_dict['O3_5007A']
#                     self.Properties_dict['OIII_HII_pn'] = self.O3_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, to_eval = 'L(4959)+L(5007)', Hbeta = self.Hbeta_flux)                           
#                           
#         #Determine the O/H abundance
#         if (self.Properties_dict['OII_HII_3279A_pn'] is not None) and (self.Properties_dict['OIII_HII_pn'] is not None):
#             self.Properties_dict['OI_HI_pn'] = self.Properties_dict['OII_HII_3279A_pn'] + self.Properties_dict['OIII_HII_pn'] 
#             self.Properties_dict['OII_HII_pn'] = self.Properties_dict['OII_HII_3279A_pn']  
#         elif (self.Properties_dict['OII_HII_7319A_pn'] is not None) and (self.Properties_dict['OIII_HII_pn'] is not None):
#             self.Properties_dict['OI_HI_pn'] = self.Properties_dict['OII_HII_7319A_pn'] + self.Properties_dict['OIII_HII_pn']   
#             self.Properties_dict['OII_HII_pn'] = self.Properties_dict['OII_HII_7319A_pn']  
#         
#         return
#         
#     def nitrogen_abundance_scheme(self):
#         
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TSIII_pn'] is not None: 
#             Te  = self.Properties_dict['TSIII_pn'] #YOU CAN TEST THIS WITH OBJECT 51959-092
#         else:
#             Te  = self.Properties_dict['TNII_pn']
#               
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne  = self.Properties_dict['nSII_pn']
#         else:
#             ne  = self.Properties_dict['nOII_pn']
# 
#         if ne is not None:
# 
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
# 
#             #Calculate the N+1 ion abundance (If O2_3726A, O2_3729A are available we prioritize them)
#             if Te is not None:
#                                 
#                 if self.lines_dict.viewkeys() >= {'N2_6548A', 'N2_6584A'}:         
#                     Ratio                               = self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']
#                     self.Properties_dict['NII_HII_pn']  = self.N2_atom.getIonAbundance(int_ratio = Ratio, tem=Te, den=ne, to_eval = 'L(6548)+L(6584)', Hbeta = self.Hbeta_flux)                           
#         
#             
#         #Calculate the elemental abundance
#         if (self.Properties_dict['OI_HI_pn'] is not None) and (self.Properties_dict['NII_HII_pn'] is not None):
#             
#             #Standard method 
#             self.Properties_dict['NI_OI_pn'] = self.Properties_dict['NII_HII_pn'] / self.Properties_dict['OII_HII_pn']
#             NI_HI_initial = self.Properties_dict['NI_OI_pn'] * self.Properties_dict['OI_HI_pn']
# 
#             #Repeat calculation if 5755 line was observed to include the recombination contribution
#             if self.lines_dict.viewkeys() >= {'N2_5755A'}:         
#                 
#                 NIII_HI             = NI_HI_initial - self.Properties_dict['NII_HII_pn']
#                 Lines_Correction    = 3.19 * power((Te/10000), 0.30) * NIII_HI * self.Hbeta_flux      
#                 
#                 self.Properties_dict['TNII_pn'], nSII = self.diags.getCrossTemDen(diag_tem      = '[NII] 5755/6584+',
#                                                                                     diag_den    = '[SII] 6731/6716',
#                                                                                     value_tem   = (self.lines_dict['N2_5755A'] - Lines_Correction)/(self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']),
#                                                                                     value_den   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#                 
#                 Ratio                               = self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']
#                 self.Properties_dict['NII_HII_pn']  = self.N2_atom.getIonAbundance(int_ratio = Ratio, tem=self.Properties_dict['TNII_pn'], den=ne, to_eval = 'L(6548)+L(6584)', Hbeta = self.Hbeta_flux)                 
#                 self.Properties_dict['NI_OI_pn'] = self.Properties_dict['NII_HII_pn'] / self.Properties_dict['OII_HII_pn']
#                 self.Properties_dict['NI_HI_pn'] = self.Properties_dict['NI_OI_pn'] * self.Properties_dict['OI_HI_pn']
#             
#             else:
#                 self.Properties_dict['NI_HI_pn'] = NI_HI_initial
#     
#         return
#     
#     def helium_abundance_scheme_fab(self):
#         
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TOIII_approxfrom_TSIII_pn'] is not None:
#             TOIII   = self.Properties_dict['TOIII_approxfrom_TSIII_pn']
#         else:
#             TOIII   = self.Properties_dict['TOIII_pn']
#                 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn'] 
#         else:
#             ne = self.Properties_dict['nOII_pn']       
#                 
#         if (TOIII is not None) and (ne is not None):
#             
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)    
#             
#             #Calculate the He+1 abundance
#             if self.lines_dict.viewkeys() >= {'He1_4472A', 'He1_5876A', 'He1_6678A'}:
# 
#                 y0_I4471    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_4472A'], tem=TOIII, den=ne, wave = 4471, Hbeta = self.Hbeta_flux)
#                 y0_I5876    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_5876A'], tem=TOIII, den=ne, wave = 5876, Hbeta = self.Hbeta_flux)
#                 y0_I6678    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_6678A'], tem=TOIII, den=ne, wave = 6678, Hbeta = self.Hbeta_flux)
#                
#                 t4          = TOIII / 10000
# 
#                 D           = 1.0 + 3110.0 * power(t4,-0.51) * (1.0/ne)
#                 g_I4471     = 6.11 * power(t4,0.02)  * (exp(-4.544)) / D
#                 g_I5876     = (7.12   * power(t4,0.14)  * (exp(-3.776/t4)) + 1.47 * power(t4,-0.28) * (exp(-4.544/t4))) / D
#                 g_I6678     = (3.27   * power(t4,-0.41) * (exp(-3.777/t4)) + 0.49 * power(t4,-0.52) * (exp(-4.544/t4))) / D
#                 
#                 y_I4471     = y0_I4471 / (1.0 + g_I4471)
#                 y_I5876     = y0_I5876 / (1.0 + g_I5876)
#                 y_I6678     = y0_I6678 / (1.0 + g_I6678)
#                      
#                 self.Properties_dict['HeII_HII_pn'] = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)     
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.Properties_dict['HeIII_HII_pn'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TOIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
#         #Calculate elemental abundance
#         if (self.Properties_dict['HeII_HII_pn'] is not None) and (self.Properties_dict['HeIII_HII_pn'] is not None):
#             self.Properties_dict['HeI_HI_pn'] = self.Properties_dict['HeII_HII_pn'] + self.Properties_dict['HeIII_HII_pn']
#         elif (self.Properties_dict['HeII_HII_pn'] is not None):
#             self.Properties_dict['HeI_HI_pn'] = self.Properties_dict['HeII_HII_pn']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.Properties_dict['HeI_HI_pn'] is not None):
#             #From Oxygen:
#             if self.Properties_dict['OI_HI_pn'] is not None:
#                 self.Properties_dict['Y_Mass_O_pn'] = (4 * self.Properties_dict['HeI_HI_pn'] * (1 - 20 * self.Properties_dict['OI_HI_pn'])) / (1 + 4 * self.Properties_dict['HeI_HI_pn']) 
#             #From Sulfur
#             if self.Properties_dict['SI_HI_pn'] is not None:
#                 self.Properties_dict['Y_Mass_S_pn'] = (4 * self.Properties_dict['HeI_HI_pn'] * (1 - 20 * self.OI_SI * self.Properties_dict['SI_HI_pn'])) / (1 + 4 * self.Properties_dict['HeI_HI_pn']) 
#     
#         return
# 
#     def helium_abundance_scheme_Oxygen(self, lineslog_frame):
# 
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TOIII_pn'] is not None:
#             TOIII   = self.Properties_dict['TOIII_pn']
#         else:
#             TOIII   = None
#                 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = None
#                      
#         if (TOIII is not None) and (ne is not None):
# 
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             
#             HeI_indices = (lineslog_frame.Ion.str.contains('HeI_')) & (lineslog_frame.index != 'He1_8446A')  & (lineslog_frame.index != 'He1_7818A')
#             HeI_labels  = lineslog_frame.loc[HeI_indices].index.values
#             HeI_waves   = lineslog_frame.loc[HeI_indices].TheoWavelength.values
#             
#             Flux_Hbeta  = self.lines_dict['H1_4861A']
#             Emis_Hbeta  = self.H1_atom.getEmissivity(tem=TOIII, den=ne, label = '4_2', product = False)
#     
#             for i in range(len(HeI_labels)):
#                 
#                 print 'Esta es', HeI_labels[i]
#                 print 'HeI_waves[i]', round(HeI_waves[i])
#                 
#                 line_relative_Flux          = self.lines_dict[HeI_labels[i]] / Flux_Hbeta
#                 line_relative_emissivity    = self.He1_atom.getEmissivity(tem = TOIII, den = ne, wave = round(HeI_waves[i]), product = False) / Emis_Hbeta
#                  
#                 if i == 0:
#                     matrix_HeI_fluxes       = copy(line_relative_Flux)
#                     matrix_HeI_emis         = copy(line_relative_emissivity)
#          
#                 else:
#                     matrix_HeI_fluxes       = vstack((matrix_HeI_fluxes, line_relative_Flux)) 
#                     matrix_HeI_emis         = vstack((matrix_HeI_emis,   line_relative_emissivity))    
#     
#             matrix_HeI_fluxes = transpose(matrix_HeI_fluxes)
#             matrix_HeI_emis   = transpose(matrix_HeI_emis)           
#             
#             #Perform the fits
#             params = Parameters()
#             params.add('Y', value = 0.01)
#             HeII_HII_array = zeros(len(matrix_HeI_fluxes))
#             HeII_HII_error = zeros(len(matrix_HeI_fluxes))
#            
#             for i in range(len(matrix_HeI_fluxes)):
#                 fit_Output = lmfit_minimmize(residual_Y_v3, params, args=(matrix_HeI_emis[i], matrix_HeI_fluxes[i]))
#                 HeII_HII_array[i] = fit_Output.params['Y'].value
#                 HeII_HII_error[i] = fit_Output.params['Y'].stderr
#             
#             #NO SUMANDO LOS ERRORES CORRECTOS?
# #             self.Properties_dict['HeII_HII_from_O_pn'] = unumpy.uarray(HeII_HII_array, mean(HeII_HII_error))
# #             self.Properties_dict['HeII_HII_from_O_pn'] = HeII_HII_array
#             self.Properties_dict['HeII_HII_from_O_pn'] = random.normal(mean(HeII_HII_array), mean(HeII_HII_error), size = self.MC_array_len)
# 
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.Properties_dict['HeIII_HII_from_O_pn'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TOIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
# 
#         #Calculate elemental abundance
#         if (self.Properties_dict['HeII_HII_from_O_pn'] is not None) and (self.Properties_dict['HeIII_HII_from_O_pn'] is not None):
#             self.Properties_dict['HeI_HI_from_O_pn'] = self.Properties_dict['HeII_HII_from_O_pn'] + self.Properties_dict['HeIII_HII_from_O_pn']
#         elif (self.Properties_dict['HeII_HII_from_O_pn'] is not None):
#             self.Properties_dict['HeI_HI_from_O_pn'] = self.Properties_dict['HeII_HII_from_O_pn']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.Properties_dict['HeI_HI_from_O_pn'] is not None):
#             #From Oxygen:
#             if self.Properties_dict['OI_HI_pn'] is not None:
#                 self.Properties_dict['Y_Mass_O_pn'] = (4 * self.Properties_dict['HeI_HI_from_O_pn'] * (1 - 20 * self.Properties_dict['OI_HI_pn'])) / (1 + 4 * self.Properties_dict['HeI_HI_from_O_pn']) 
# 
#         return
# 
#     def helium_abundance_scheme_Sulfur(self, lineslog_frame):
# 
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TOIII_approxfrom_TSIII_pn'] is not None:
#             TOIII   = self.Properties_dict['TOIII_approxfrom_TSIII_pn']
#         else:
#             TOIII   = None
#                 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = None
#                      
#         if (TOIII is not None) and (ne is not None):
# 
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             
#             HeI_indices = lineslog_frame.Ion.str.contains('HeI_') & (lineslog_frame.index != 'He1_8446A') & (lineslog_frame.index != 'He1_7818A')
#             HeI_labels  = lineslog_frame.loc[HeI_indices].index.values
#             HeI_waves   = lineslog_frame.loc[HeI_indices].TheoWavelength.values
#             
#             Flux_Hbeta  = self.lines_dict['H1_4861A']
#             Emis_Hbeta  = self.H1_atom.getEmissivity(tem=TOIII, den=ne, label = '4_2', product = False)
#     
#             for i in range(len(HeI_labels)):
#                  
#                 line_relative_Flux          = self.lines_dict[HeI_labels[i]] / Flux_Hbeta
#                 line_relative_emissivity    = self.He1_atom.getEmissivity(tem = TOIII, den = ne, wave = round(HeI_waves[i]), product = False) / Emis_Hbeta
#                  
#                 if i == 0:
#                     matrix_HeI_fluxes       = copy(line_relative_Flux)
#                     matrix_HeI_emis         = copy(line_relative_emissivity)
#          
#                 else:
#                     matrix_HeI_fluxes       = vstack((matrix_HeI_fluxes, line_relative_Flux)) 
#                     matrix_HeI_emis         = vstack((matrix_HeI_emis,   line_relative_emissivity))    
#     
#             matrix_HeI_fluxes = transpose(matrix_HeI_fluxes)
#             matrix_HeI_emis   = transpose(matrix_HeI_emis)           
#             
#             #Perform the fits
#             params = Parameters()
#             params.add('Y', value = 0.01)
#             HeII_HII_array = zeros(len(matrix_HeI_fluxes))
#             HeII_HII_error = zeros(len(matrix_HeI_fluxes))
#            
#             for i in range(len(matrix_HeI_fluxes)):
#                 fit_Output = lmfit_minimmize(residual_Y_v3, params, args=(matrix_HeI_emis[i], matrix_HeI_fluxes[i]))
#                 HeII_HII_array[i] = fit_Output.params['Y'].value
#                 HeII_HII_error[i] = fit_Output.params['Y'].stderr
#                 
#             #NO SUMANDO LOS ERRORES CORRECTOS?
# #             self.Properties_dict['HeII_HII_from_S_pn'] = unumpy.uarray(HeII_HII_array, mean(HeII_HII_error))
# #             self.Properties_dict['HeII_HII_from_S_pn'] = HeII_HII_array
#             self.Properties_dict['HeII_HII_from_S_pn'] = random.normal(mean(HeII_HII_array), mean(HeII_HII_error), size = self.MC_array_len)
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.Properties_dict['HeIII_HII_from_S_pn'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TOIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
# 
#         #Calculate elemental abundance
#         if (self.Properties_dict['HeII_HII_from_S_pn'] is not None) and (self.Properties_dict['HeIII_HII_from_S_pn'] is not None):
#             self.Properties_dict['HeI_HI_from_S_pn'] = self.Properties_dict['HeII_HII_from_S_pn'] + self.Properties_dict['HeIII_HII_from_S_pn']
#         elif (self.Properties_dict['HeII_HII_from_S_pn'] is not None):
#             self.Properties_dict['HeI_HI_from_S_pn'] = self.Properties_dict['HeII_HII_from_S_pn']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.Properties_dict['HeI_HI_from_S_pn'] is not None):
#             #From Oxygen:
#             if self.Properties_dict['SI_HI_pn'] is not None:
#                 self.Properties_dict['Y_Mass_S_pn'] = (4 * self.Properties_dict['HeI_HI_from_S_pn'] * (1 - 20 * self.OI_SI * self.Properties_dict['SI_HI_pn'])) / (1 + 4 * self.Properties_dict['HeI_HI_from_S_pn']) 
# 
#         return
# 
#     def store_dict_ObjLog(self, FileFolder, CodeName, verbose = False):
#         
#         keys, values = self.Properties_dict.keys(), self.Properties_dict.values()
#         
#         for i in range(len(keys)):
#             if values[i] is not None:
#                 #This removes the nan elements from the treatment values[i][~isnan(values[i])]
#                 Median_Value    = median(values[i][~isnan(values[i])])
#                 error           = std(values[i][~isnan(values[i])])  
#                 self.SaveParameter_ObjLog(CodeName, FileFolder, keys[i], ufloat(Median_Value, error))
# 
#                 if verbose == True:
#                     print i, keys[i], Median_Value, error
#                 
#             else:
#                 self.SaveParameter_ObjLog(CodeName, FileFolder, keys[i], values[i])
#                 if verbose == True:
#                     print i, keys[i], values[i]  
#         return
#     
#     def SaveParameter_ObjLog(self, CodeName, FileFolder, Parameter, Magnitude, Error = None, Assumption = None, Log_extension = None):
#                 
#         if Log_extension == None:
#             Log_extension = self.ObjectLog_extension
#             ObjLog_Address                                                          = FileFolder + CodeName + Log_extension
#             ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments    = loadtxt(ObjLog_Address, dtype=str, usecols = [0,1,2,3], unpack = True)
#         
#         else:
#             ObjLog_Address                                                          = FileFolder + CodeName
#             ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments    = loadtxt(ObjLog_Address, dtype=str, usecols = [0,1,2,3], unpack = True)
# 
#         #Saving the new value for the given parameter
#         Parameter_Index                                                         = where(ObjLog_Parameters == Parameter)[0][0]
# 
#         #Save the error in the case of a nparray
#         if isinstance(Magnitude, UFloat) and (Error == None):
#             ObjLog_Magnitudes[Parameter_Index]                                  = Magnitude.nominal_value
#             ObjLog_Errors[Parameter_Index]                                      = Magnitude.std_dev
#         elif Error != None:
#             ObjLog_Magnitudes[Parameter_Index]                                  = str(Magnitude)
#             ObjLog_Errors[Parameter_Index]                                      = str(Error)
#         elif Magnitude != None:
#             ObjLog_Magnitudes[Parameter_Index]                                  = str(Magnitude)
#         elif Magnitude == None:
#             ObjLog_Magnitudes[Parameter_Index]                                  = 'None'
#             ObjLog_Errors[Parameter_Index]                                      = '-'
#             
#         #Saving the text file
#         savetxt(ObjLog_Address, transpose((ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments)), fmt='%s')
#     
#         return    
# class pyneb_tools():
#     
#     def __init__(self):
#  
#         #Set atomic data 
#         atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
#         atomicData.setDataFile('s_iii_coll_HRS12.dat')
#  
#         #Declare hydrogen atom so it becomes readily available:
#         self.H1_atom    = RecAtom('H', 1)
#         self.S2_atom    = Atom('S', 2)
#         self.S3_atom    = Atom('S', 3)
#         self.He1_atom   = RecAtom('He', 1)
#         self.He2_atom   = RecAtom('He', 2)
#         
#     def density_pyneb(self, Ion, Temp):
#         
#         element     = Ion[0:-1]
#         ion_n       = int(Ion[-1])
#                 
#         if Ion == 'S2':
#             
#             Ion_Atom    = Atom(element, ion_n)
#             
#             emLine      = ["S2_6716A","S2_6731A"]
#             Flux_dict   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#             
#             FluxRatio   = Flux_dict['S2_6716A'] / Flux_dict['S2_6731A']
#             
#             Den         = Ion_Atom.getTemDen(FluxRatio, tem = Temp, wave1=6716, wave2=6731)
#             
#             return Den
#         
#     def temperature_pyneb(self, Ion, Den):
#         
#         element     = Ion[0:-1]
#         ion_n       = int(Ion[-1])
#         
#         if Ion == 'O3':        
#         
#             Ion_Atom    = Atom(element, ion_n)
#             
#             emLine      = ["S2_6716A","S2_6731A"]
#             Flux_dict   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)        
#         
#             FluxRatio   = (Flux_dict['O3_4959A'] + Flux_dict['O3_5007A']) / Flux_dict['O3_4363A']
# 
#             Temp        = Ion_Atom.getTemDen(FluxRatio, den = Den, to_eval = "(L(4959) + L(5007)) / L(4363)")
#         
#             return Temp
#    
# class direct_Method(ReddeningLaws, pyneb_tools):
#     
#     def __init__(self):
#         
#         pyneb_tools.__init__(self)
#         ReddeningLaws.__init__(self)
#     
#         self.den_Methodology        = None
#         self.tem_Methodology        = None
#         
#         self.m_SIV_correction       = ufloat(1.008,0.037)
#         self.n_SIV_correction       = ufloat(1.576,0.066)
#                         
#     def DeclareObject(self, FileFolder, FileName, CodeName, AbundancesFileExtension = '_WHT_LinesLog_v3.txt', cHbeta_type = 'cHBeta_red'):
#         #WARNING: This method should be define at a lower level
#                 
#         self.FileFolder             = FileFolder
#         self.FileName               = FileName
#         self.CodeName               = CodeName
#     
#         EmLine_List                 = ['R_SII', 'R_SIIprime', 'R_SIII', 'R_NII', 'R_OII', 'R_OIII']
#         Den_List                    = ['nSII']
#         Temp_List                   = ['TOIII', 'TOII', 'TSII', 'TSIII','TNII', 'TOII_approx_TOIII', 'TSIII_approx_TOIII', 'TOIII_approx_TSIII', 'TNII_approx_TOIII']
#         IonicAbund_List             = ['SII_HII', 'SIII_HII', 'SIV_HII', 'OII_HII', 'OII_HII_3279A', 'OII_HII_7319A', 'NII_HII', 'HeII_HII', 'HeIII_HII', 'ArIII_HII', 'ArIV_HII']
#         Element_List                = ['SI_HI', 'SI_HI_ArCorr', 'OI_HI', 'NI_OI', 'NI_HI', 'HeI_HI']
#             
#         self.Properties_dict        = dict.fromkeys(Element_List + IonicAbund_List + Den_List + Temp_List + EmLine_List)
#         
#         self.cHbeta                 = self.GetParameter_ObjLog(self.CodeName, self.FileFolder, Parameter = cHbeta_type, Assumption = 'float')
#         
#         self.AbundancesExtension    = AbundancesFileExtension
#         
#     def density_determination(self, Ion, Temp = None, methodology = None):
#                 
#         #Paremetric approach
#         if methodology == 'Epm2014':
#             if Ion == 'S2':
#                 
#                 emLine                      = ["S2_6716A", "S2_6731A"]
#                 Flux_dict                   = self.getLinesFlux_dered(emLine, self.cHbeta)
#                 Flux_dict['Temp']           = Temp
#                 den                         = self.empiric_formulae(Flux_dict, methodology, Ion, 'nSII')
#                 
#                 #Check for non physical output 
#                 den                         = self.check_issues(magnitude = den, parameter_type = 'density')
#                 return den
#         
#         #Using pyneb 
#         if methodology == 'pyneb':
#             
#             #Calculate density using pyneb
#             den = self.density_pyneb(Ion, Temp)
#             
#             #Check for non physical output #WARNING: For pyneb this is too general
#             den = self.check_issues(magnitude = den, parameter_type = 'density')
#             
#             return den
#     
#     def temperature_determination(self, Ion, dens = None, methodology = None, Temp = None, RecCorr = None):
#                 
#         #Direct methods
#         if methodology == 'Epm2014':
#     
#             if Ion == 'O3':
#                 
#                 emLine                          = ["O3_4959A","O3_5007A","O3_4363A"]
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TOIII')
#  
#             if Ion == 'O2':               
#                 data_dict                       = {'Temp': Temp}                 
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TOII_approx_TOIII')
#                                 
#             if Ion == 'S3':
#                 
#                 emLine                          = ["S3_9069A","S3_9531A","S3_6312A"]
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TSIII')
#                             
#             if Ion == 'N2':
#                 
#                 emLine                          = ["N2_6548A","N2_6584A", "N2_5755A", 'H1_4861A']
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 if RecCorr != None:
#                     if data_dict["N2_5755A"] != None:
#                         data_dict["N2_5755A"] = data_dict["N2_5755A"] - RecCorr * data_dict['H1_4861A']
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TNII')
# 
#         if methodology == 'Haegele2008':
# 
#             if Ion == 'S2':
#                 
#                 emLine                          = ["S2_6716A","S2_6731A", "S2_4069A", 'S2_4076A']
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta =  self.cHbeta)
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TSII')
# 
#             if Ion == 'O2':
# 
#                 emLine                          = ['O2_3726A',"O2_3729A","O2_7319A",'O2_7330A']
#                 data_dict                       = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 temp                            = self.empiric_formulae(data_dict, methodology, Ion, 'TOII')
# 
#         
#         #Using pyneb 
#         if methodology == 'pyneb':
#             
#             temp = self.temperature_pyneb(Ion, dens)
# 
#         return temp
# 
#     def argon_IonAbundance(self, Temp, Den, Ion, methodology = 'Haegele2008'):
#         
#         if methodology == 'Haegele2008':
#             
#             if Ion == 'Ar3':
#                 
#                 emLine                              = ["Ar3_7136A","H1_4861A"]
#                 data_dict                           = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                   = Temp
#                 self.Properties_dict['ArIII_HII']   = self.empiric_formulae(data_dict, methodology, Ion, 'ArIII_HII')
#                 print 'The argon thing III', self.Properties_dict['ArIII_HII']
#                 return
# 
#             if Ion == 'Ar4':
#                 
#                 emLine                              = ["Ar4_4740A","H1_4861A"]
#                 data_dict                           = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                   = Temp             
#                 self.Properties_dict['ArIV_HII']    = self.empiric_formulae(data_dict, methodology, Ion, 'ArIV_HII')
#                 print 'The argon thing IV', self.Properties_dict['ArIV_HII']
#                 return
#                        
#     def oxygen_IonAbundance(self, Temp, Den, Ion, methodology = 'Epm2014', RecCorr = None):
# 
#         if methodology == 'Epm2014':
#     
#             if Ion == 'O3':
#                 
#                 emLine                                      = ["O3_4959A","O3_5007A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                          = Temp
#                 self.Properties_dict['OIII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'OIII_HII')
#                 
#                 return
# 
#             if Ion == 'O2':
#                                                             #This assumes both lines are blended
#                 emLine                                      = ["O2_3726A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta, Mode='Integrated')
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den 
#                 self.Properties_dict['OII_HII_3279A']       = self.empiric_formulae(data_dict, methodology, Ion, 'OII_HII')
#                 
#         if methodology == 'Fabian2006':
# 
#             if Ion == 'O2':
#                 
#                 emLine                                      = ["O2_7319A","O2_7330A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']        = Temp, Den 
#                 
#                 if RecCorr == None:
#                     self.Properties_dict['OII_HII_7319A']   = self.empiric_formulae(data_dict, methodology, Ion, 'OII_HII_7319A')
#                 else:
#                     data_dict["O2_7319A"]                   = data_dict["O2_7319A"] - RecCorr
#                     self.Properties_dict['OII_HII_7319A']   = self.empiric_formulae(data_dict, methodology, Ion, 'OII_HII_7319A')
#                     
#                 return
#             
#     def nitrogen_IonAbundance(self, Temp, Den, Ion, methodology = 'Epm2014'):
#             
#         if methodology == 'Epm2014':
#             
#             if Ion == 'N2':
#             
#                 emLine                                      = ["N2_6548A","N2_6584A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)            
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['NII_HII']             = self.empiric_formulae(data_dict, methodology, Ion, 'NII_HII')
#                 
#                 return
#             
#     def sulfur_IonAbundance(self,  Temp, Den, Ion, methodology = 'Haegele2008'):
#         
#         if methodology == 'Haegele2008':
#             
#             if Ion == 'S3':
#                 
#                 emLine                                      = ["S3_9069A","S3_9531A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['SIII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'SIII_HII')
#                 
#                 return
#       
#             if Ion == 'S2':
#                 
#                 emLine                                      = ["S2_6716A","S2_6731A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['SII_HII']             = self.empiric_formulae(data_dict, methodology, Ion, 'SII_HII')
# 
#         if methodology == 'Vital2015':
#             if Ion == 'S3':                
#                 emLine                                      = ["S3_9069A","S3_9531A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['SIII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'SIII_HII')
# 
#                 return
#       
#             if Ion == 'S2':
#                 methodology = 'Haegele2008'
#                 emLine                                      = ["S2_6716A","S2_6731A","H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['SII_HII']             = self.empiric_formulae(data_dict, methodology, Ion, 'SII_HII')
# 
#                 
#                 return
#     
#     def helium_IonAbundance(self, Temp, Den, Ion, methodology = 'Fabian2006'):
#         
#         if methodology == 'Fabian2006':
#             
#             if Ion == 'He1':                
#                 emLine                                      = ["He1_4472A", "He1_5876A", "He1_6678A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['HeII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'HeII_HII')
#                 
#                 return
#             
#             if Ion == 'He2':
#                 emLine                                      = ["He2_4686A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['HeIII_HII']           = self.empiric_formulae(data_dict, methodology, Ion, 'HeIII_HII')
#                  
#                 return
#                  
#         if methodology == 'Vital2015':
#             
#             if Ion == 'He1':                
#                 emLine                                      = ["He1_4472A", "He1_5876A", "He1_6678A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['HeII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'HeII_HII')
#                             
#             if Ion == 'He2':
#                 emLine                                      = ["He2_4686A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp']                           = Temp
#                 self.Properties_dict['HeIII_HII']           = self.empiric_formulae(data_dict, methodology, Ion, 'HeIII_HII')
#         
#         if methodology == 'PyNeb':
# 
#             if Ion == 'He1':                
#                 emLine                                      = ["He1_4472A", "He1_5876A", "He1_6678A", "H1_4861A"]
#                 data_dict                                   = self.getLinesFlux_dered(Lines_List = emLine, cHbeta = self.cHbeta)
#                 data_dict['Temp'], data_dict['Den']         = Temp, Den
#                 self.Properties_dict['HeII_HII']            = self.empiric_formulae(data_dict, methodology, Ion, 'HeII_HII')
# 
#                 return
# 
#             if Ion == 'He2':
#         
#                 return
# 
#     def empiric_formulae(self, data_dict, methodology, Ion, Parameter):
#     
#         Magnitude = None
#         if data_dict != None:
#             if (None not in data_dict.values()):
#                 if methodology == 'Epm2014':
#                     if Parameter == 'nSII':
#                     
#                         self.Properties_dict['R_SII']               = data_dict['S2_6716A'] / data_dict['S2_6731A']
#                         t4                                          = data_dict['Temp'] / 10000.0       
#                         a_0                                         = 16.054  - 7.79  / t4 - 11.32  * t4
#                         a_1                                         = -22.66  + 11.08 / t4 + 16.02  * t4
#                         b_0                                         = -21.61  + 11.89 / t4 + 14.59  * t4
#                         b_1                                         = 9.17    - 5.09  / t4 - 6.18   * t4
#                         Magnitude                                   = 1000.0 * (self.Properties_dict['R_SII']*a_0 + a_1) / (self.Properties_dict['R_SII']*b_0 + b_1)
#                         #self.Properties_dict['nSII']               = 1000.0 * (self.Properties_dict['R_SII']*a_0 + a_1) / (self.Properties_dict['R_SII']*b_0 + b_1)
#                         
#                     elif Parameter == 'TOIII':
#                         self.Properties_dict['R_OIII']              = (data_dict['O3_4959A'] + data_dict['O3_5007A']) / data_dict['O3_4363A']
#                         self.Properties_dict['TOIII']               = (0.7840 - 0.0001357 * self.Properties_dict['R_OIII'] + (48.44 / self.Properties_dict['R_OIII'])) * 10000
#                         Magnitude                                   = self.Properties_dict['TOIII']
#                         
#                     elif Parameter == 'TOII_approx_TOIII':
#                         t4 = data_dict['Temp'] / 10000
#                         self.Properties_dict['TOII_approx_TOIII']   = (1.397 /( (1/t4) + 0.385) ) * 10000
#                         Magnitude                                   = self.Properties_dict['TOII_approx_TOIII']
#                         
#                     elif Parameter == 'TSIII_approx_TOIII':
#                         #From equation TOIII = 1.0807 * TSIII - 0.0846
#                         t4 = data_dict['Temp'] / 10000
#                         self.Properties_dict['TSIII_approx_TOIII']  = ((t4 + 0.0846) / 1.0807)*10000
# 
#                     elif Parameter == 'TOIII_approx_TSIII':    
#                         t4 = data_dict['Temp'] / 10000
#                         self.Properties_dict['TOIII_approx_TSIII']  = (1.0807 * t4 - 0.0846)*10000
#                         Magnitude                                   = self.Properties_dict['TOIII_approx_TSIII']
# 
#                     elif Parameter == 'TSIII':
#                         self.Properties_dict['R_SIII']              = (data_dict['S3_9069A'] + data_dict['S3_9531A']) / data_dict['S3_6312A']
#                         self.Properties_dict['TSIII']               = (0.5147 + 0.0003187 * self.Properties_dict['R_SIII'] + (23.6404 / self.Properties_dict['R_SIII'])) * 10000                         
#                         Magnitude                                   = self.Properties_dict['TSIII']
#     
#                     elif Parameter == 'TNII':
#                         self.Properties_dict['R_NII']               = (data_dict['N2_6548A'] + data_dict['N2_6584A']) / data_dict['N2_5755A']
#                         self.Properties_dict['TNII']                = (0.6153 - 0.0001529 * self.Properties_dict['R_NII'] + (35.3641 / self.Properties_dict['R_NII'])) * 10000
#                         Magnitude                                   = self.Properties_dict['TNII']
# 
#                     elif Parameter == 'TNII_approx_TOIII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         self.Properties_dict['TNII_approx_TOIII']   = (1.452 /( (1/t4) + 0.479))*10000
#                         Magnitude                                   = self.Properties_dict['TNII_approx_TOIII']
# 
#                     elif Parameter == 'OIII_HII':
#                         t4 = data_dict['Temp'] / 10000
#                         logOIII_logHII                              = -12 + umath.log10((data_dict['O3_4959A'] + data_dict['O3_5007A']) / data_dict['H1_4861A']) + 6.1868 + 1.2491 / t4 - 0.5816 * umath.log10(t4)
#                         Magnitude                                   = umath.pow(10, logOIII_logHII)
#                     
#                     elif Parameter == 'OII_HII_3279A':
#                         t4 = data_dict['Temp'] / 10000
#                         ne = data_dict['Den']
#                         logOII_logHII                               = -12 + umath.log10((data_dict['O2_3726A']) / data_dict['H1_4861A']) + 5.887 + 1.641 / t4 - 0.543 * umath.log10(t4) + 0.000114 * ne
#                         Magnitude                                   = umath.pow(10, logOII_logHII) 
#                         
#                     elif Parameter == 'NII_HII':
#                         t4  = data_dict['Temp'] / 10000                       
#                         logNII_logHII                               = -12 + umath.log10((data_dict['N2_6548A'] + data_dict['N2_6584A']) / data_dict['H1_4861A']) + 6.291 + 0.90221 / t4 - 0.5511 * umath.log10(t4)
#                         Magnitude                                   = umath.pow(10, logNII_logHII)
# 
#                 elif methodology == 'Angeles2015':
# 
#                     if Parameter == 'SIV_HII':
#  
#                         ArIII_HII   = data_dict['ArIII_HII']
#                         ArIV_HII    = data_dict['ArIV_HII']
#                         SIII_HII    = data_dict['SIII_HII']
#                                                  
#                         if (ArIII_HII != None) and (ArIV_HII != None) and (ArIV_HII != 0.0):   #Somehow ArIV should be saved as None 
#                             logAr2Ar3   = umath.log10(ArIII_HII/ArIV_HII)
# 
#                         else:
#                             logAr2Ar3 = ufloat(0.0, 0.0)
#                       
#                         logSIV   = umath.log10(SIII_HII) - (logAr2Ar3 - self.n_SIV_correction) / self.m_SIV_correction
# 
#                                   
#                         self.Properties_dict['SIV_HII'] = umath.pow(10, logSIV)
#                         Magnitude = self.Properties_dict['SIV_HII']
# 
# #                         Old formulation
# #                         ArIII_HII   = data_dict['ArIII_HII']
# #                         ArIV_HII    = data_dict['ArIV_HII']
# #                         SIII_HII    = data_dict['SIII_HII']
# #                                                 
# #                         if (ArIII_HII != None) and (ArIV_HII != 0.0):    
# #                             y_Ar = umath.log10(ArIII_HII/ArIV_HII)
# #                         else:
# #                             y_Ar = ufloat(0.0, 0.0)
# #                      
# #                         x_S = (y_Ar - 0.126) / 1.192
# #                 
# #                         self.Properties_dict['SIV_HII'] = SIII_HII / umath.pow(10, x_S)
# #                         Magnitude = self.Properties_dict['SIV_HII']
#                         
#                 elif methodology == 'Haegele2008':
#                     
#                     if Parameter == 'TSII':
#                         self.Properties_dict['R_SIIprime']          = (data_dict['S2_6716A'] + data_dict['S2_6731A']) / (data_dict['S2_4069A'] + data_dict['S2_4076A'])
#                         self.Properties_dict['TSII']               = (1.92 - 0.0375 *self.Properties_dict['R_SIIprime'] - 14.5 / self.Properties_dict['R_SIIprime'] + 105.64 / (self.Properties_dict['R_SIIprime'] * self.Properties_dict['R_SIIprime'])) * 10000
#                         Magnitude                                   = self.Properties_dict['TSII']
# 
#                     if Parameter == 'TOII':
#                         ne = data_dict['Den']
#                         self.Properties_dict['R_OII']               = (data_dict['O2_3726A'] + data_dict['O2_3729A']) / (data_dict['O2_7319A'] + data_dict['O2_7330A'])
#                         self.Properties_dict['TOII']               = (0.23 + 0.0017 * self.Properties_dict['R_OII'] + 38.3 / self.Properties_dict['R_OII'] + umath.log10(1 + 0.0001 * ne)) * 10000
#                         Magnitude                                   = self.Properties_dict['TOII']
# 
#                     elif Parameter == 'ArIII_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         logArIII_logHII                             = -12 + umath.log10(data_dict['Ar3_7136A'] / data_dict['H1_4861A']) +  6.157 + 0.808/t4 - 0.508 * umath.log10(t4) 
#                         Magnitude                                   = umath.pow(10, logArIII_logHII)
# 
#                     elif Parameter == 'ArIV_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         logArIV_logHII                              = -12 + umath.log10(data_dict['Ar4_4740A'] / data_dict['H1_4861A']) +  5.705 + 1.246/t4 - 0.156 * umath.log10(t4) 
#                         Magnitude                                   = umath.pow(10, logArIV_logHII)
# 
#                     elif Parameter == 'SII_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         ne = data_dict['Den']
#                         logSII_logHII                               = -12 + umath.log10((data_dict['S2_6716A'] + data_dict['S2_6731A']) / data_dict['H1_4861A']) + 5.423 + 0.929/t4 - 0.280 * umath.log10(t4) + umath.log10(1 + 0.0001 * ne)
#                         Magnitude                                   = umath.pow(10, logSII_logHII)
# 
#                     elif Parameter == 'SIII_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         logSIII_logHII                              = -12 + umath.log10((data_dict['S3_9069A'] + data_dict['S3_9531A']) / data_dict['H1_4861A']) + 5.8 + 0.77 / t4 - 0.22 * umath.log10(t4) 
#                         Magnitude                                   = umath.pow(10, logSIII_logHII)
#                 
#                 elif methodology == 'Vital2015':
# 
#                     if Parameter == 'SIII_HII':
#                         t4  = data_dict['Temp'] / 10000
#                         logSIII_logHII                              = -12 + umath.log10((data_dict['S3_9069A'] + data_dict['S3_9531A']) / data_dict['H1_4861A']) + 6.012 + 0.6309 / t4 - 0.5722 * umath.log10(t4) 
#                         Magnitude                                   = umath.pow(10, logSIII_logHII)   
#                         
#                     if Parameter == 'HeIII_HII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         y_PPI4686                                   = 0.081 * umath.pow(t4, 0.15) * data_dict['He2_4686A']  /  data_dict['H1_4861A']
#                         Magnitude                                   = y_PPI4686                                            
#                 
#                     if Parameter == 'HeII_HII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         ne                                          = data_dict['Den']
#                                         
#                         y0_I4471                                    = 2.039   * umath.pow(t4,0.102) * data_dict['He1_4472A']  / data_dict['H1_4861A']
#                         y0_I5876                                    = 0.732  * umath.pow(t4,0.171) * data_dict['He1_5876A']  / data_dict['H1_4861A']
#                         y0_I6678                                    = 2.573   * umath.pow(t4,0.18) * data_dict['He1_6678A']  / data_dict['H1_4861A']
#                         
#                         D                                           = 1.0 + 3110.0 * umath.pow(t4,-0.51) * (1.0/ne)
#                         g_I4471                                     = 6.11 * umath.pow(t4,0.02)  * (umath.exp(-4.544)) / D
#                         g_I5876                                     = (7.12   * umath.pow(t4,0.14)  * (umath.exp(-3.776/t4)) + 1.47 * umath.pow(t4,-0.28) * (umath.exp(-4.544/t4))) / D
#                         g_I6678                                     = (3.27   * umath.pow(t4,-0.41) * (umath.exp(-3.777/t4)) + 0.49 * umath.pow(t4,-0.52) * (umath.exp(-4.544/t4))) / D
# 
#                         y_I4471                                     = y0_I4471 / (1.0 + g_I4471)
#                         y_I5876                                     = y0_I5876 / (1.0 + g_I5876)
#                         y_I6678                                     = y0_I6678 / (1.0 + g_I6678)
#                     
#                         Magnitude                                   = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)                
#                 
#                 elif methodology == 'Liu2000':
#                     
#                     if Parameter == 'O2_RecombCorr':
#                         t4 = data_dict['Temp'] / 10000  
#                         OII_HII                                     = data_dict['OII_HII_7319A']
#                         OII_RecombCor                               = 9.36 * umath.pow(t4, 0.44) * OII_HII
#                         Magnitude                                   = OII_RecombCor
#  
#                     if Parameter == 'N2_RecombCorr':
#                         t4 = data_dict['Temp'] / 10000  
#                         NIII_HII                                    = data_dict['NIII_HII']
#                         NII_RecombCor                               = 3.19 * umath.pow(t4, 0.30) * NIII_HII
#                         Magnitude                                   = NII_RecombCor                        
#                 
#                 elif methodology == 'Epm2007':
#                 
#                     if Parameter == 'O2':
#                         ne                                          = data_dict['Den']
#                         t4                                          = data_dict['Temp'] / 10000
#                         O2plus_Hplus                                = data_dict['OIII_HII']
#                         I7320_IBeta_R                               = 9.36 * umath.pow(t4, 0.44) * O2plus_Hplus
#                         logOII_logHII                               = -12 + umath.log10((data_dict['O2_7319A'] + data_dict['O2_7330A']) / data_dict['H1_4861A'] - I7320_IBeta_R) + 6.895 + 2.44/t4 - 0.58*umath.log10(t4) - umath.log10(1.0 + 0.0047 * ne)
#                         Magnitude                                   = umath.pow(10, logOII_logHII)
#                                 
#                 elif methodology == 'Fabian2006':
#                 
#                     if Parameter == 'HeII_HII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         ne                                          = data_dict['Den']                
#                 
#                         y0_I4471                                    = 2.04   * umath.pow(t4,0.13) * data_dict['He1_4472A']  / data_dict['H1_4861A']
#                         y0_I5876                                    = 0.738  * umath.pow(t4,0.23) * data_dict['He1_5876A']  / data_dict['H1_4861A']
#                         y0_I6678                                    = 2.58   * umath.pow(t4,0.25) * data_dict['He1_6678A']  / data_dict['H1_4861A']
#                         
#                         D                                           = 1.0 + 3110.0 * umath.pow(t4,-0.51) * (1.0/ne)
#                         g_I4471                                     = 6.11 * umath.pow(t4,0.02)  * (umath.exp(-4.544)) / D
#                         g_I5876                                     = (7.12   * umath.pow(t4,0.14)  * (umath.exp(-3.776/t4)) + 1.47 * umath.pow(t4,-0.28) * (umath.exp(-4.544/t4))) / D
#                         g_I6678                                     = (3.27   * umath.pow(t4,-0.41) * (umath.exp(-3.777/t4)) + 0.49 * umath.pow(t4,-0.52) * (umath.exp(-4.544/t4))) / D
# 
#                                 
#                         y_I4471                                     = y0_I4471 / (1.0 + g_I4471)
#                         y_I5876                                     = y0_I5876 / (1.0 + g_I5876)
#                         y_I6678                                     = y0_I6678 / (1.0 + g_I6678)
#                     
#                         Magnitude                                   = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)                
#                 
#                     elif Parameter == 'HeIII_HII':
#                         t4                                          = data_dict['Temp'] / 10000
#                         y_PPI4686                                   = 0.084 * umath.pow(t4, 0.14) * data_dict['He2_4686A']  /  data_dict['H1_4861A']
#                         Magnitude                                   = y_PPI4686
#                 
#                     elif Parameter == 'OII_HII_7319A':
#                         t4 = data_dict['Temp'] / 10000
#                         ne = data_dict['Den']                        
#                         logOII_logHII                               = -12 + umath.log10((data_dict['O2_7319A'] + data_dict['O2_7330A']) / data_dict['H1_4861A']) + 6.895 + 2.44 / t4 - 0.58 * umath.log10(t4) - umath.log10(1 + 0.0047 * ne)
#                         Magnitude                                   = umath.pow(10, logOII_logHII)
#                         
#                 elif methodology == 'PyNeb':
#                         
#                     if Parameter == 'HeII_HII':
#                         
#                         te                                          = data_dict['Temp']
#                         ne                                          = data_dict['Den']
#                                         
#                         self.He1.getIonAbundance()
# 
#                         y0_I4471                                    = self.He1_atom.getIonAbundance(int_ratio = 100 * data_dict['He1_4472A']/data_dict['H1_4861A'] , tem=te, den=ne, label='He1_4472A')
#                         y0_I5876                                    = self.He1_atom.getIonAbundance(int_ratio = 100 * data_dict['He1_5876A']/data_dict['H1_4861A'] , tem=te, den=ne, label='He1_5876A')
#                         y0_I6678                                    = self.He1_atom.getIonAbundance(int_ratio = 100 * data_dict['He1_6678A']/data_dict['H1_4861A'] , tem=te, den=ne, label='He1_6678A')
#                                                 
# #                         y0_I4471                                    = 2.04   * umath.pow(t4,0.13) * data_dict['He1_4472A']  / data_dict['H1_4861A']
# #                         y0_I5876                                    = 0.738  * umath.pow(t4,0.23) * data_dict['He1_5876A']  / data_dict['H1_4861A']
# #                         y0_I6678                                    = 2.58   * umath.pow(t4,0.25) * data_dict['He1_6678A']  / data_dict['H1_4861A']
#                         
#                         D                                           = 1.0 + 3110.0 * umath.pow(t4,-0.51) * (1.0/ne)
#                         g_I4471                                     = 6.11 * umath.pow(t4,0.02)  * (umath.exp(-4.544)) / D
#                         g_I5876                                     = (7.12   * umath.pow(t4,0.14)  * (umath.exp(-3.776/t4)) + 1.47 * umath.pow(t4,-0.28) * (umath.exp(-4.544/t4))) / D
#                         g_I6678                                     = (3.27   * umath.pow(t4,-0.41) * (umath.exp(-3.777/t4)) + 0.49 * umath.pow(t4,-0.52) * (umath.exp(-4.544/t4))) / D
# 
#                                 
#                         y_I4471                                     = y0_I4471 / (1.0 + g_I4471)
#                         y_I5876                                     = y0_I5876 / (1.0 + g_I5876)
#                         y_I6678                                     = y0_I6678 / (1.0 + g_I6678)
#                     
#                         Magnitude                                   = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)
#                         
#                         
#         return Magnitude
#                       
#     def check_issues(self, magnitude, parameter_type):
#         
#         #Security check in the case we are getting a negative density:
#         if parameter_type == 'density':
#             if magnitude != None:
#                 
#                 if magnitude <= 0:
#                     return ufloat(50.0, 50.0)            
#                 else:
#                     return magnitude
#             return magnitude
#         
#         #Security check in case some flux is missing
#         elif parameter_type == 'EmFlux':
#             if  magnitude[1] == None:
#                 return None
#             else:
#                 return ufloat(magnitude[1], magnitude[2]) 
#     
# class Chemical_Analysis(direct_Method):
# 
#     def __init__(self):
#         
#         direct_Method.__init__(self)
#         
#         #logSI_OI_Gradient = ufloat(-1.53, 0.05)
#         logSI_OI_Gradient = ufloat(-1.78, 0.03)
#                 
#         self.OI_SI      = umath.pow(10, -logSI_OI_Gradient)
#         self.He1        = RecAtom('He', 1)
#         atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
#         
#     def set_element(self, element):
#         
#         if element == 'Argon':
#             
#             self.argon_abundance_scheme()
#             
#         elif element == 'Oxygen':
#             
#             self.oxygen_abundance_scheme()
#             
#         elif element == 'Nitrogen':
#         
#             self.nitrogen_abundance_scheme()
#         
#         elif element == 'Sulfur': 
# 
#             self.sulfur_abundance_scheme()
# 
#         elif element == 'Helium':
# 
#             self.helium_abundance_scheme()
# 
#         return
# 
#     def argon_abundance_scheme(self):
# 
#         #Determine the sulfur density and temperature    
#         T_SII,  T_SIII, ne_SII  = self.sulfur_density_scheme()
#         T_OIII, T_OII           = self.oxygen_temperature_scheme()
#         
#         #We try to calculate the T_ArIII from the sulfur lines
#         T_ArIII                 = T_SIII
#         T_ArIV                  = T_OIII
#         
#         #We try to calculate the T_ArIV from the sulfur lines, if not we use the Oxygen ones
#         if  (T_OIII == None) and (T_SIII != None):
#             data_dict           = {'Temp' : T_SIII}
#             T_OIII_approx       = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion='O3', Parameter = 'TOIII_approx_TSIII')
#             T_ArIV              = T_OIII_approx
# 
#         elif (T_SIII == None) and (T_OIII != None):
#             data_dict           = {'Temp' : T_OIII}
#             T_SIII_approx       = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion='S3', Parameter = 'TSIII_approx_TOIII')
#             T_ArIV              = T_SIII_approx
#             
#         #Calculate the ArIII abundance
#         self.argon_IonAbundance(T_ArIII, ne_SII, Ion = 'Ar3', methodology = 'Haegele2008')
#   
#         #Calculate the ArIV abundance
#         self.argon_IonAbundance(T_ArIV, ne_SII, Ion = 'Ar4', methodology = 'Haegele2008')  
#    
#     def sulfur_density_scheme(self, TempIn = None):
#         
#         #We use this set up mechanism to calculate the electron density for other elements 
#         if TempIn == None:
#             #Determine the Te[SII] and Te[SIII]
#             T_SII   = self.temperature_determination(Ion = 'S2', methodology = 'Haegele2008')
#             T_SIII  = self.temperature_determination(Ion = 'S3', methodology = 'Epm2014')
#         else:
#             T_SII   = None
#             T_SIII  = TempIn
#                 
#         #Determine ne_SII using the TSIII. If not available use TOIII to calculate TSIII
#         #If a TempIn is not available it will use the standard procedure 
#         if T_SIII  != None:
#             ne_SII = self.density_determination(Ion = 'S2', Temp = T_SIII, methodology='Epm2014')
#         else:
#             T_OIII, T_OII       = self.oxygen_temperature_scheme()
#             data_dict           = {}
#             data_dict['Temp']   = T_OIII
#             T_SIII              = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion = 'S3', Parameter='TSIII_approx_TOIII')
#             ne_SII              = self.density_determination(Ion = 'S2', Temp = T_SIII, methodology='Epm2014')
#         
#         return T_SII, T_SIII, ne_SII
#     
#     def sulfur_abundance_scheme(self):
#         
#         #Get the sulfur electron density and temperature
#         T_SII, T_SIII, ne_SII = self.sulfur_density_scheme(TempIn = None)
#         
#         self.Properties_dict['nSII'] = ne_SII
#                            
#         #Determine the SII/HII abundance. 
#         #If the T_SII temperature is available use it. If not T_SII = T_SIII:
#         if T_SII != None:
#             self.sulfur_IonAbundance(T_SII, ne_SII, Ion = 'S2', methodology = 'Haegele2008')
#         else:
#             self.sulfur_IonAbundance(T_SIII, ne_SII, Ion = 'S2', methodology = 'Haegele2008')
# 
#         #Determine the SIII/HII abundance
# #         self.sulfur_IonAbundance(T_SIII, ne_SII, Ion = 'S3', methodology = 'Haegele2008')
#         self.sulfur_IonAbundance(T_SIII, ne_SII, Ion = 'S3', methodology = 'Vital2015')
# 
#         #Determine the S/H abundance
#         #Include the Argon correction if the ArIV lines were observed
#         if (self.Properties_dict['SII_HII'] != None) and (self.Properties_dict['SIII_HII'] != None):
#             self.Properties_dict['SI_HI'] = self.Properties_dict['SII_HII'] + self.Properties_dict['SIII_HII']
#                         
#             #Calculate the [SIV] correction factor from the Argon abundance
#             #This structures enforces the calculation of the Argon apriori
#             data_dict               = {}
#             data_dict['ArIII_HII']  = self.Properties_dict['ArIII_HII']   
#             data_dict['ArIV_HII']   = self.Properties_dict['ArIV_HII']   
#             data_dict['SIII_HII']   = self.Properties_dict['SIII_HII']               
#             
#             #The code does not work if we introduce a None entry. However, in this correction there is still a quantity even if the ArIV is no there
#             if data_dict['ArIV_HII'] == None:
#                 data_dict['ArIV_HII'] = ufloat(0.0,0.0) 
#             
#             #Compute the SIV_HII component
#             SIV_HII = self.empiric_formulae(data_dict, methodology = 'Angeles2015', Ion = 'S4', Parameter='SIV_HII')         
#             
#             self.Properties_dict['SI_HI_ArCorr'] = self.Properties_dict['SII_HII'] + self.Properties_dict['SIII_HII'] + SIV_HII
#              
#     def oxygen_temperature_scheme(self):
# 
#         #Determine the Te[OIII]
#         T_OIII          = self.temperature_determination(Ion = 'O3', methodology = 'Epm2014')
#         
#         #Determine the Te[OII] using all the [OII] lines
#         T_OII           = self.temperature_determination(Ion = 'O2', methodology = 'Epm2014')
#         
#         #If lines not observed use approximation from T[OIII]
#         if T_OII == None:
#             T_OII       = self.temperature_determination(Ion = 'O2', methodology = 'Epm2014', Temp = T_OIII)
# 
#         return T_OIII, T_OII
# 
#     def oxygen_abundance_scheme(self):
#         
#         #Determine the oxygen temperatures
#         T_OIII, T_OII = self.oxygen_temperature_scheme()
#         
#         #Get the electron density from the sulfur lines using the oxigen temperature if observed. If not used sulfur lines
#         T_SII, T_SIII, ne_SII = self.sulfur_density_scheme()
#         
#         #Calculate the OIII ion abundance
#         self.oxygen_IonAbundance(Temp = T_OIII, Den = ne_SII, Ion = 'O3', methodology = 'Epm2014')
#         
#         #Calculate the OII ion abundance
#         self.oxygen_IonAbundance(Temp = T_OII, Den = ne_SII, Ion = 'O2', methodology = 'Epm2014')
#         self.oxygen_IonAbundance(Temp = T_OII, Den = ne_SII, Ion = 'O2', methodology = 'Fabian2006')
#         
#         #Correct the OII from 7319A, 7330A lines from the recombination contribution
#         if self.Properties_dict['OII_HII_7319A'] != None:
#             OII_HII_7319A   = self.Properties_dict['OII_HII_7319A']
#             data_dict       = {'Temp'           : T_OII, 
#                                'OII_HII_7319A'  : OII_HII_7319A}
#             
#             RecCorr_O2_7319 = self.empiric_formulae(data_dict, methodology = 'Liu2000', Ion = 'O2', Parameter='O2_RecombCorr')
#             self.oxygen_IonAbundance(Temp = T_OII, Den = ne_SII, Ion = 'O2', methodology = 'Epm2014', RecCorr = RecCorr_O2_7319)    
#         
#         #If O2_3726A, O2_3729A lines are not available we use O2_7319A, O2_7330A lines
#         if self.Properties_dict['OII_HII_3279A'] != None:
#             self.Properties_dict['OII_HII'] = self.Properties_dict['OII_HII_3279A']
#             
#         elif (self.Properties_dict['OII_HII_7319A'] != None): 
#             self.Properties_dict['OII_HII']     = self.Properties_dict['OII_HII_7319A']
#         
#         #Determine the S/H abundance
#         #Include the Argon correction if the ArIV lines were observed
#         if (self.Properties_dict['OII_HII'] != None) and (self.Properties_dict['OIII_HII'] != None):
#             self.Properties_dict['OI_HI'] = self.Properties_dict['OII_HII'] + self.Properties_dict['OIII_HII']
#             
#     def nitrogen_abundance_scheme(self):
#         #WARNING: THIS COULD BE FASTER WITH A MECHANISM WITH PRELOADS ALL THE EMISSION LINES
#         
#         #Determine the Te[NSII]
#         T_NII                   = self.temperature_determination(Ion = 'N2', methodology = 'Epm2014')
#                 
#         if T_NII != None:
#             #YOU CAN TEST THIS WITH OBJECT 51959-092
#             T_SII, T_SIII, ne_SII   = self.sulfur_density_scheme(TempIn = T_NII)
#             self.nitrogen_IonAbundance(Temp = T_NII, Den = ne_SII, Ion = 'N2', methodology = 'Epm2014')
#             
#         #If Te[NSII] cannot be calculated try to determine the Te[OII]
#         elif T_NII == None:
#             T_OIII, T_OII       = self.oxygen_temperature_scheme()            
#             data_dict           = {'Temp' : T_OIII}
#             TNII_approx         = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion='N2', Parameter = 'TNII_approx_TOIII')
#         
#             #Get the sulfur electron density and temperature
#             T_SII, T_SIII, ne_SII   = self.sulfur_density_scheme(TempIn = TNII_approx)
#             self.nitrogen_IonAbundance(Temp = TNII_approx, Den = ne_SII, Ion = 'N2', methodology = 'Epm2014')
# 
#         #Determine the N/H abundance
#         #We need O/H abundance in order to calculate N/H abundance
#         if (self.Properties_dict['OI_HI'] != None) and (self.Properties_dict['NII_HII'] != None):
#             self.Properties_dict['NI_OI'] = self.Properties_dict['NII_HII'] / self.Properties_dict['OII_HII']
#             self.Properties_dict['NI_HI'] = self.Properties_dict['NI_OI'] * self.Properties_dict['OI_HI']
#             
#             #Take into account the recombination contribution from the 
#             if T_NII != None:
#                 NIII_HI = self.Properties_dict['NI_HI'] - self.Properties_dict['NII_HII']
#                 
#                 data_dict       = {'Temp'           : T_NII, 
#                                    'NIII_HII'       : NIII_HI}
#             
#                 RecCorr_N2_5755 = self.empiric_formulae(data_dict, methodology = 'Liu2000', Ion = 'N2', Parameter = 'N2_RecombCorr')
#                 
#                 #Repeat the previous steps starting with the calculation of TNII
#                 T_NII   = self.temperature_determination(Ion = 'N2', methodology = 'Epm2014', RecCorr = RecCorr_N2_5755)
#                 
#                 T_SII, T_SIII, ne_SII   = self.sulfur_density_scheme(TempIn = T_NII)
#                 self.nitrogen_IonAbundance(Temp = T_NII, Den = ne_SII, Ion = 'N2', methodology = 'Epm2014')
#                 self.Properties_dict['NI_OI'] = self.Properties_dict['NII_HII'] / self.Properties_dict['OII_HII']
#                 self.Properties_dict['NI_HI'] = self.Properties_dict['NI_OI'] * self.Properties_dict['OI_HI']
#                                        
#     def helium_abundance_scheme(self):
#         
#         #Determine the sulfur density and temperature    
#         T_SII, T_SIII, ne_SII = self.sulfur_density_scheme()
#     
#         #In case we do not have the sulfur temperature we use the oxygen one
#         if T_SIII != None:
#             data_dict   = {'Temp' : T_SIII}
#             T_OIII      = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion = 'O3', Parameter = 'TOIII_approx_TSIII')
#         
#         else:
#             T_OIII, T_OII = self.oxygen_temperature_scheme()
#         
#         #Calculate the HeII ion  abundance
#         self.helium_IonAbundance(Temp = T_OIII, Den = ne_SII , Ion = 'He1', methodology = 'Vital2015')
# 
#         #Calculate the HeIII ion abundance
#         self.helium_IonAbundance(Temp = T_OIII, Den = ne_SII , Ion = 'He2', methodology = 'Vital2015')
# 
#         if self.Properties_dict['HeII_HII'] != None:
#             if self.Properties_dict['HeIII_HII'] != None:
#                 self.Properties_dict['HeI_HI'] = self.Properties_dict['HeII_HII'] + self.Properties_dict['HeIII_HII']
#             elif self.Properties_dict['HeIII_HII'] == None:
#                 self.Properties_dict['HeI_HI'] = self.Properties_dict['HeII_HII']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.Properties_dict['HeI_HI'] != None):
#             
#             #From Oxygen:
#             if self.Properties_dict['OI_HI'] != None:
#                 self.Properties_dict['Y_Mass_O'] = (4 * self.Properties_dict['HeI_HI'] * (1 - 20 * self.Properties_dict['OI_HI'])) / (1 + 4 * self.Properties_dict['HeI_HI']) 
#             
#             #From Sulfur
#             if self.Properties_dict['SI_HI'] != None:
#                 self.Properties_dict['Y_Mass_S'] = (4 * self.Properties_dict['HeI_HI'] * (1 - 20 * self.OI_SI * self.Properties_dict['SI_HI_ArCorr'])) / (1 + 4 * self.Properties_dict['HeI_HI']) 
# 
#     def heliumPyneb_abundance_scheme(self):
#         
#         #Determine the sulfur density and temperature    
#         T_SII, T_SIII, ne_SII = self.sulfur_density_scheme()
#     
#         #In case we do not have the sulfur temperature we use the oxygen one
#         if T_SIII != None:
#             data_dict   = {'Temp' : T_SIII}
#             T_OIII      = self.empiric_formulae(data_dict, methodology = 'Epm2014', Ion = 'O3', Parameter = 'TOIII_approx_TSIII')
#         
#         else:
#             T_OIII, T_OII = self.oxygen_temperature_scheme()
#         
#         #Calculate the HeII ion  abundance
#         self.helium_IonAbundance(Temp = T_OIII, Den = ne_SII , Ion = 'He1', methodology = 'PyNeb')        
#         
#         #Calculate the HeIII ion abundance
#         self.helium_IonAbundance(Temp = T_OIII, Den = ne_SII , Ion = 'He2', methodology = 'PyNeb')
#                 
#         return 
#            
#     def store_dict_ObjLog(self, varaiables_dict, verbose = False):
#         
#         keys, values = varaiables_dict.keys(), varaiables_dict.values()
#         
#         for i in range(len(keys)):
#             if verbose:
#                 print i, keys[i], values[i]
#             self.SaveParameter_ObjLog(self.CodeName, self.FileFolder, keys[i], values[i])
#         
#         return    
#     
# class Parametrized_Emissivities(Txt_Files_Manager):
#     
#     def __init__(self):
# 
#         self.LiteratureParametrization = None
#         Txt_Files_Manager.__init__(self)
# 
#     def AOS2010(self, Te, ne):
#         
#         #phi has the shape : phi_lambda = Em(Hbeta) / E(He_lambda) * (1) / (1 + C/R(lambda)) 
#         #Hence Em(He_lambda)/Em(H\beta) = 1 / phi_lambda
#                 
#         Phi3889 = power(0.8779 * Te, -0.128 -0.00041 * ne )
#         Phi4026 = power(4.233 * Te, 0.085 - 0.00012 * ne )
#         Phi4471 = power(2.021 * Te, 0.121 - 0.00020 * ne )
#         Phi5876 = power(0.754 * Te, 0.212 - 0.00051 * ne )
#         Phi6678 = power(2.639 * Te, 0.244 - 0.00054 * ne )
#         Phi7065 = power(5.903 * Te, -0.519) / (1.462 - Te * (0.127 - 0.00076 * ne * + 0.000000255 * ne*ne))
#         
#         Em_He3889 = 1/Phi3889
#         Em_He4026 = 1/Phi4026
#         Em_He4471 = 1/Phi4471
#         Em_He5876 = 1/Phi5876
#         Em_He6678 = 1/Phi6678
#         Em_He7065 = 1/Phi7065
#         
#         return Em_He3889, Em_He4026, Em_He4471, Em_He4471, Em_He5876, Em_He5876, Em_He6678, Em_He7065
#     
#     def Haegele2008_Fluxes(self):
#         #SDSS J1729
#         Flux_dict = {}
#         
#         Flux_dict["H1_4861A"]   = ufloat(10000, 83)
#         
#         Flux_dict["O2_3726A"]   = ufloat(17622, 243)
#         Flux_dict["O2_3729A"]   = ufloat(0.0, 0.0)
#         Flux_dict["O2_7319A"]   = ufloat(233,8)
#         Flux_dict['O2_7330A']   = ufloat(194,9)        
#         
#         Flux_dict["O3_4363A"]   = ufloat(660, 30)
#         Flux_dict["O3_4959A"]   = ufloat(17097, 146)
#         Flux_dict["O3_5007A"]   = ufloat(51541, 418)        
#     
#         Flux_dict["N2_5755A"]   = ufloat(60, 6)
#         Flux_dict["N2_6548A"]   = ufloat(771, 21)
#         Flux_dict["N2_6584A"]   = ufloat(2189, 52)        
#         
#         Flux_dict['Ar3_7136A']  = ufloat(864,41)
#         Flux_dict['Ar4_4740A']  = ufloat(40,0.06)  
#         
#         Flux_dict["S2_6716A"]   = ufloat(1286, 32)
#         Flux_dict["S2_6731A"]   = ufloat(1009, 30)
#         Flux_dict["S2_4069A"]   = ufloat(122, 11)
#         Flux_dict['S2_4076A']   = ufloat(0.0,0.0)
#         
#         Flux_dict["S3_6312A"]   = ufloat(176, 11)
#         Flux_dict["S3_9069A"]   = ufloat(2092, 95)
#         Flux_dict["S3_9531A"]   = ufloat(4718, 250)
#          
#         Flux_dict["He1_4472A"]   = ufloat(516, 24)
#         Flux_dict["He1_5876A"]   = ufloat(1261, 37)
#         Flux_dict["He1_6678A"]   = ufloat(355, 2)
#         Flux_dict['He2_4686A']   = ufloat(329,43)         
# 
#     def Abundances_FabianSample(self):
#         
#         Table_Address   = '/home/vital/Workspace/Dazer/Astro_Libraries/Fabian2006_CatalogueAbundances.txt'
# 
#         Headers_List    = ['O/H_10^5_Flux', 'O/H_0^5_Error', 'N/H_10^6_Flux', 'N/H_10^6_Error', 'Y_magnitude', 'Y_error', 'y_magnitude', 'y_error']
#         
#         Oxygen_Flux, Oxygen_Error, Nitrogen_Flux, Nitrogen_Error, Y, Yerror, y_magnitude, y_error = self.get_ColumnData(Headers_List, Table_Address, HeaderSize=1, unpack_check = True)
#         
#         return Oxygen_Flux * 1e-5, Oxygen_Error*1e-5, Nitrogen_Flux*1e-6, Nitrogen_Error*1e-6, Y, Yerror,  y_magnitude, y_erro
#     def helium_abundance_scheme_fab(self):
#         
#         #Select the temperature and density for the procedure
#         if self.abunData['TOIII_approxfrom_TSIII'] is not None:
#             TOIII   = self.abunData['TOIII_approxfrom_TSIII']
#         else:
#             TOIII   = self.abunData['TOIII']
#                 
#         if self.abunData['nSII'] is not None:
#             ne = self.abunData['nSII'] 
#         else:
#             ne = self.abunData['nOII']       
#                 
#         if (TOIII is not None) and (ne is not None):
#             
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)    
#             
#             #Calculate the He+1 abundance
#             if self.lines_dict.viewkeys() >= {'He1_4472A', 'He1_5876A', 'He1_6678A'}:
# 
#                 y0_I4471    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_4472A'], tem=TOIII, den=ne, wave = 4471, Hbeta = self.Hbeta_flux)
#                 y0_I5876    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_5876A'], tem=TOIII, den=ne, wave = 5876, Hbeta = self.Hbeta_flux)
#                 y0_I6678    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_6678A'], tem=TOIII, den=ne, wave = 6678, Hbeta = self.Hbeta_flux)
#                
#                 t4          = TOIII / 10000
# 
#                 D           = 1.0 + 3110.0 * power(t4,-0.51) * (1.0/ne)
#                 g_I4471     = 6.11 * power(t4,0.02)  * (exp(-4.544)) / D
#                 g_I5876     = (7.12   * power(t4,0.14)  * (exp(-3.776/t4)) + 1.47 * power(t4,-0.28) * (exp(-4.544/t4))) / D
#                 g_I6678     = (3.27   * power(t4,-0.41) * (exp(-3.777/t4)) + 0.49 * power(t4,-0.52) * (exp(-4.544/t4))) / D
#                 
#                 y_I4471     = y0_I4471 / (1.0 + g_I4471)
#                 y_I5876     = y0_I5876 / (1.0 + g_I5876)
#                 y_I6678     = y0_I6678 / (1.0 + g_I6678)
#                      
#                 self.abunData['HeII_HII'] = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)     
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.abunData['HeIII_HII'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TOIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
#         #Calculate elemental abundance
#         if (self.abunData['HeII_HII'] is not None) and (self.abunData['HeIII_HII'] is not None):
#             self.abunData['HeI_HI'] = self.abunData['HeII_HII'] + self.abunData['HeIII_HII']
#         elif (self.abunData['HeII_HII'] is not None):
#             self.abunData['HeI_HI'] = self.abunData['HeII_HII']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.abunData['HeI_HI'] is not None):
#             #From Oxygen:
#             if self.abunData['OI_HI'] is not None:
#                 self.abunData['Y_Mass_O'] = (4 * self.abunData['HeI_HI'] * (1 - 20 * self.abunData['OI_HI'])) / (1 + 4 * self.abunData['HeI_HI']) 
#             #From Sulfur
#             if self.abunData['SI_HI'] is not None:
#                 self.abunData['Y_Mass_S'] = (4 * self.abunData['HeI_HI'] * (1 - 20 * self.OI_SI * self.abunData['SI_HI'])) / (1 + 4 * self.abunData['HeI_HI']) 
#     
#         return

# class Chemical_Analysis_pyneb_old():
# 
#     def __init__(self):
#         
#         self.MC_array_len           = 100
#         
#         self.Hbeta_label            = 'H1_4861A'
#     
#     def load_elements(self):
# 
#         #Set atomic data 
#         atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
#         atomicData.setDataFile('s_iii_coll_HRS12.dat')
#  
#         #Declare hydrogen atom so it becomes readily available:
#         self.S2_atom            = Atom('S', 2)
#         self.S3_atom            = Atom('S', 3)
#         self.Ar3_atom           = Atom('Ar', 3)
#         self.Ar4_atom           = Atom('Ar', 4)
#         self.N2_atom            = Atom('N', 2)
#         self.O2_atom            = Atom('O', 2)
#         self.O3_atom            = Atom('O', 3)
#                     
#         self.H1_atom            = RecAtom('H', 1)
#         self.He1_atom           = RecAtom('He', 1)
#         self.He2_atom           = RecAtom('He', 2)
# 
#         self.diags              = Diagnostics()
#         
#         self.logSI_OI_Gradient  = random.normal(-1.53,  0.05, size = self.MC_array_len) # random.normal(-1.78,  0.03, size = self.MC_array_len)
# 
#         self.OI_SI              = power(10, -self.logSI_OI_Gradient)
#         
#         self.S3_9000_ratio      = random.normal(self.S3_atom.getEmissivity(10000, 1000, wave = 9531) / self.S3_atom.getEmissivity(10000, 1000, wave = 9069),  0.01, size = self.MC_array_len)
#         
#         self.m_SIV_correction   = random.normal(1.109,  0.1109, size = self.MC_array_len)
#         self.n_SIV_correction   = random.normal(0.135,  0.013, size = self.MC_array_len)
#         
#         return
#             
#     def DeclareObject(self, lines_log_frame):
#         
#         #List of all parameters
#         EmLine_List             = ['R_SII_pn', 'R_SIIprime_pn', 'R_SIII_pn', 'R_NII_pn', 'R_OII_pn', 'R_OIII_pn']
#         Den_List                = ['nSII_pn', 'nOII_pn']
#         Temp_List               = ['TOIII_pn', 'TOII_pn', 'TSII_pn', 'TSIII_pn','TNII_pn', 'TOII_approxfrom_TOIII_pn', 'TSIII_approxfrom_TOIII_pn', 'TOIII_approxfrom_TSIII_pn', 'TNII_approxfrom_TOIII_pn']
#         IonicAbund_List         = ['SII_HII_pn', 'SIII_HII_pn', 'SIV_HII_pn', 'OII_HII_pn', 'OII_HII_3279A_pn', 'OII_HII_7319A_pn', 'NII_HII_pn', 'HeII_HII_from_O_pn', 'HeIII_HII_from_O_pn', 'HeII_HII_from_S_pn', 'HeIII_HII_from_S_pn' ,'ArIII_HII_pn', 'ArIV_HII_pn']
#         Element_List            = ['SI_HI_pn', 'OI_HI_pn', 'NI_OI_pn', 'NI_HI_pn', 'HeI_HI_from_O_pn', 'HeI_HI_from_S_pn', 'Y_Mass_O_pn', 'Y_Mass_S_pn']
#             
#         self.Properties_dict    = OrderedDict.fromkeys(Element_List + IonicAbund_List + Den_List + Temp_List + EmLine_List)
#                 
#         self.Hbeta_flux         = random.normal(lines_log_frame.loc['H1_4861A']['line_Int'].nominal_value, lines_log_frame.loc['H1_4861A']['line_Int'].std_dev, size = self.MC_array_len)
#               
#         self.low_density_dist   = self.Truncated_gaussian.rvs(self.MC_array_len)
# 
#         #Generate a dictionary to store the random array for all lines
#         self.lines_dict = OrderedDict()
# 
#         #Generate
#         for line in lines_log_frame.index.values:
#             if line == 'O2_3726A': 
#                 if lines_log_frame.loc['O2_3726A']['flux_intg'] == lines_log_frame.loc['O2_3726A']['flux_intg']:
#                     flux, error = lines_log_frame.loc[line]['line_IntBrute_dered'].nominal_value, lines_log_frame.loc[line]['line_IntBrute_dered'].std_dev
# 
#             elif line == 'O2_7319A': 
#                 if lines_log_frame.loc['O2_7319A']['flux_intg'] == lines_log_frame.loc['O2_7330A']['flux_intg']:
#                     flux, error = lines_log_frame.loc[line]['line_IntBrute_dered'].nominal_value, lines_log_frame.loc[line]['line_IntBrute_dered'].std_dev
#                 else:
#                     line_sum = lines_log_frame.loc['O2_7319A']['line_IntBrute_dered'] + lines_log_frame.loc['O2_7330A']['line_IntBrute_dered']
#                     flux, error = line_sum.nominal_value, line_sum.std_dev
#                 
#             else:
#                 flux, error         = lines_log_frame.loc[line]['line_Int'].nominal_value, lines_log_frame.loc[line]['line_Int'].std_dev
#             
# #             self.lines_dict[line]   = random.normal(flux, error, size = self.MC_array_len)
#             
#             self.lines_dict[line]   = random.normal(flux, error, size = self.MC_array_len)
#             
#         return
#                               
#     def low_density_function(self, density_array):
#         
#         inds_nan    = where(isnan(density_array))
#         
#         if len(inds_nan[0]) > 0:
#             
#             if len(inds_nan[0]) > (0.1 * self.MC_array_len):
#                 print '--WARNING 10% boundary: ', len(inds_nan[0]), 'values with nan'            
#             elif len(inds_nan[0]) > (0.05 * self.MC_array_len):
#                 print '--WARNING 5% boundary: ', len(inds_nan[0]), 'values with nan'            
#                                 
#             density_array[inds_nan]= self.low_density_dist[inds_nan]
#             
#         mean_density    = mean(density_array) 
#         std_density     = std(density_array)   
#         
#         if isnan(mean_density): 
#             print '--WARNING: Density array cannot be treated numerically'
#             
#         if mean_density < 75:
#             #print '--WARNING: Low density', mean_density, '+/-', std_density, 'Points below limit', len(where(density_array < 75)[0])
#             return  self.low_density_dist
#         else:                              
#             return density_array
#             
#     def determine_electron_parameters(self):
#                 
# #-------Starting with the Sulfur
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:    
#               
#             #Check TSIII lines
#             if self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_9531A', 'S3_6312A'}:
#         
#                 self.Properties_dict['TSIII_pn'], self.Properties_dict['nSII_pn'] = self.diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+',
#                                                                                     diag_den = '[SII] 6731/6716',
#                                                                                     value_tem = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] + self.lines_dict['S3_9531A']),
#                                                                                     value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']) 
#                                 
#             #Missing S3_9531A line
#             elif self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_6312A'}:
#                 
#                 self.Properties_dict['TSIII_pn'], self.Properties_dict['nSII_pn'] = self.diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+',
#                                         diag_den = '[SII] 6731/6716',
#                                         value_tem = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] * (1 + self.S3_9000_ratio)),
#                                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']) 
# 
#             #Missing S3_9069A line
#             elif self.lines_dict.viewkeys() >= {'S3_9531A', 'S3_6312A'}:
#                 self.Properties_dict['TSIII_pn'], self.Properties_dict['nSII_pn'] = self.diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+',
#                                         diag_den = '[SII] 6731/6716',
#                                         value_tem = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9531A'] * (1 + 1/self.S3_9000_ratio)),
#                                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#             
#             #Missing the sulfur lines we use the ones from oxygen
#             elif self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A', 'O3_4363A'}:         
#                 
#                 self.Properties_dict['TOIII_pn'], nSII = self.diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+',
#                         diag_den = '[SII] 6731/6716',
#                         value_tem = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A']),
#                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#             
#                 self.Properties_dict['TSIII_approxfrom_TOIII_pn'] = ((self.Properties_dict['TOIII_pn']/10000 + 0.0846) / 1.0807) * 10000
#                 self.Properties_dict['nSII_pn'] = self.S2_atom.getTemDen((self.lines_dict['S2_6731A'])/(self.lines_dict['S2_6716A']), tem = self.Properties_dict['TSIII_approxfrom_TOIII_pn'], to_eval = 'L(6731) / L(6716)') 
#                 
#             
#             #Try to measure the TSII temperature
#             if self.lines_dict.viewkeys() >= {'S2_4069A', 'S2_4076A'}:
# 
#                 self.Properties_dict['TSII_pn'], nSII = self.diags.getCrossTemDen(diag_tem = '[SII] 4069/4076',
#                         diag_den  = '[SII] 6731/6716',
#                         value_tem = self.lines_dict['S2_4069A']/self.lines_dict['S2_4076A'],
#                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#             
# #-------Following with the oxygen
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:
#                 
#             #Missing the sulfur lines we use the ones from oxygen
#             if self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A', 'O3_4363A'}:         
#                 self.Properties_dict['TOIII_pn'], nSII_x = self.diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+',
#                         diag_den  = '[SII] 6731/6716',
#                         value_tem = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A']),
#                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#                         
#             #Determine the TOIII from the sulfur temperature as a proxy
#             if self.Properties_dict['TSIII_pn'] is not None:
#                 
#                 self.Properties_dict['TOIII_approxfrom_TSIII_pn'] = (1.0807 * self.Properties_dict['TSIII_pn']/10000 - 0.0846) * 10000
#             
#             #Determine the TOII value from the TOIII proxy
#             if self.Properties_dict['TOIII_pn'] is not None:
#                 self.Properties_dict['TOII_approxfrom_TOIII_pn'] = (1.397 /( (1/(self.Properties_dict['TOIII_pn']/10000)) + 0.385) ) * 10000
# 
# #             #Calculate the oxygen density if the 3721A+ lines can be measured        
# #             if self.lines_dict.viewkeys() >= {'O2_3726A', 'O2_3729A', 'O2_7319A', 'O2_7330A', 'O3_4959A', 'O3_5007A', 'O3_4363A'}:
# # 
# #                 self.Properties_dict['TOIII_pn'], self.Properties_dict['nOII_pn'] = self.diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+',
# #                                                                                         diag_den = '[OII] 3726/3729',
# #                                                                                         value_tem = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A']),
# #                                                                                         value_den = (self.lines_dict['O2_3726A']/self.lines_dict['O2_3729A']))
# 
#             #Calculate the oxygen temperature cannot be measure use the SIII approximation       
#             elif self.lines_dict.viewkeys() >= {'O2_3726A', 'O2_3729A', 'O2_7319A', 'O2_7330A'}:
#                 if self.Properties_dict['TOIII_approxfrom_TSIII_pn'] is not None:
#                     self.Properties_dict['nOII_pn'] = self.O2_atom.getTemDen((self.lines_dict['O2_3726A']/self.lines_dict['O2_3729A']), 
#                                                                              tem = self.Properties_dict['TOIII_approxfrom_TSIII_pn'], 
#                                                                              to_eval = '(L(3726))/(L(3729))')
#                 
# #-------Following with the Nitrogen
#         if self.lines_dict.viewkeys() >= {'N2_5755A', 'N2_6548A', 'N2_6584A'}:         
# 
#             self.Properties_dict['TNII_pn'], nSII = self.diags.getCrossTemDen(diag_tem  = '[NII] 5755/6584+',
#                                                                               diag_den  = '[SII] 6731/6716',
#                                                                               value_tem = self.lines_dict['N2_5755A']/(self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']),
#                                                                               value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
# 
#         if self.Properties_dict['TOIII_pn'] is not None:
#             self.Properties_dict['TNII_approxfrom_TOIII_pn'] =  (1.452 /( (1/(self.Properties_dict['TOIII_pn']/10000)) + 0.479)) * 10000
#              
#         return
# 
#     def argon_abundance_scheme(self):
#  
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TSIII_pn'] is not None:
#             TSIII = self.Properties_dict['TSIII_pn']
#             TOIII = self.Properties_dict['TOIII_approxfrom_TSIII_pn']
# 
#         else:
#             TSIII = self.Properties_dict['TSIII_approxfrom_TOIII_pn']
#             TOIII = self.Properties_dict['TOIII_pn']
#              
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = self.Properties_dict['nOII_pn']
#                  
#         if ne is not None:
#         
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)            
#         
#             #Calculate the Ar_+2 abundance
#             if (TSIII is not None):
#                 if self.lines_dict.viewkeys() >= {'Ar3_7136A', 'Ar3_7751'}:
#                     Ratio = self.lines_dict['Ar3_7751A'] + self.lines_dict['Ar3_7136A']
#                     self.Properties_dict['ArIII_HII_pn'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, to_eval = 'L(7751) + L(7136)', Hbeta = self.Hbeta_flux)     
#                 elif self.lines_dict.viewkeys() >= {'Ar3_7136A'}:
#                     Ratio = self.lines_dict['Ar3_7136A']
#                     self.Properties_dict['ArIII_HII_pn'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, wave = 7136, Hbeta = self.Hbeta_flux)             
#                 elif self.lines_dict.viewkeys() >= {'Ar3_7751A'}:
#                     Ratio = self.lines_dict['Ar3_7751A']
#                     self.Properties_dict['ArIII_HII_pn'] = self.Ar3_atom.getIonAbundance(int_ratio = Ratio, tem=TSIII, den=ne, wave = 7751, Hbeta = self.Hbeta_flux)   
# 
#             #Calculate the Ar_+3 abundance
#             if (TOIII is not None):
#                 if self.lines_dict.viewkeys() >= {'Ar4_4740A', 'Ar4_4711A'}:
#                     Ratio = self.lines_dict['Ar4_4740A'] + self.lines_dict['Ar4_4711A']
#                     self.Properties_dict['ArIV_HII_pn'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, to_eval = 'L(4740) + L(4711)', Hbeta = self.Hbeta_flux)                             
#                 elif self.lines_dict.viewkeys() >= {'Ar4_4740A'}:
#                     Ratio = self.lines_dict['Ar4_4740A']
#                     self.Properties_dict['ArIV_HII_pn'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, wave = 4740, Hbeta = self.Hbeta_flux)
#                 elif self.lines_dict.viewkeys() >= {'Ar4_4711A'}:
#                     Ratio = self.lines_dict['Ar4_4711A']
#                     self.Properties_dict['ArIV_HII_pn'] = self.Ar4_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, wave = 4711, Hbeta = self.Hbeta_flux)
#                                           
#     def sulfur_abundance_scheme(self):
#         
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TSIII_pn'] is not None:
#             Te = self.Properties_dict['TSIII_pn']
#         else:
#             Te = self.Properties_dict['TSIII_approxfrom_TOIII_pn']
# 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = self.Properties_dict['nOII_pn']
# 
#         #Calculating the R_S2 ratio
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:
#             self.Properties_dict['R_SII_pn'] = self.lines_dict['S2_6716A'] + self.lines_dict['S2_6731A']
#                 
#         #Calculating the R_S3 ratio
#         if self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_9531A'}:
#             self.Properties_dict['R_SIII_pn'] = self.lines_dict['S3_9069A'] + self.lines_dict['S3_9531A']
#         
#         elif self.lines_dict.viewkeys() >= {'S3_9069A'}:
#             self.Properties_dict['R_SIII_pn'] = self.lines_dict['S3_9069A'] * (1 + self.S3_9000_ratio)
#             
#         elif self.lines_dict.viewkeys() >= {'S3_9531A'}:
#             self.Properties_dict['R_SIII_pn'] = self.lines_dict['S3_9531A'] * (1 + 1/self.S3_9000_ratio)          
#                 
#         #Determine ionic abundances
#         if (Te is not None) and (ne is not None):
#                         
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             
#             #Calculate the S+1 abundance
#             if self.Properties_dict['R_SIII_pn'] is not None:
#                 self.Properties_dict['SII_HII_pn'] = self.S2_atom.getIonAbundance(int_ratio=self.Properties_dict['R_SII_pn'], tem=Te, den=ne, to_eval = 'L(6731)+L(6716)', Hbeta = self.Hbeta_flux)        
#             
#             #Calculate the S_+2 abundance
#             if self.Properties_dict['R_SIII_pn'] is not None:
#                 self.Properties_dict['SIII_HII_pn'] = self.S3_atom.getIonAbundance(int_ratio=self.Properties_dict['R_SIII_pn'], tem=Te, den=ne, to_eval = 'L(9069)+L(9531)', Hbeta = self.Hbeta_flux)
#     
#             #Calculating the S_+3 abundance
#             if (self.Properties_dict['SII_HII_pn'] is not None) and (self.Properties_dict['SIII_HII_pn'] is not None):
#                                 
#                 #Apply correction from Ar - Sulfur from Popstar - Cloudy models if Argon metallicity previously measured     
#                 if (self.Properties_dict['ArIII_HII_pn'] is not None) and (self.Properties_dict['ArIV_HII_pn'] is not None):
#                     logAr2Ar3   = log10(self.Properties_dict['ArIII_HII_pn']/self.Properties_dict['ArIV_HII_pn'])
#                     logSIV  = log10(self.Properties_dict['SIII_HII_pn']) - (logAr2Ar3 - self.n_SIV_correction) / self.m_SIV_correction
#                     self.Properties_dict['SIV_HII_pn'] = power(10, logSIV)
# 
#                 else:
#                     self.Properties_dict['SIV_HII_pn'] = zeros(self.MC_array_len)
#         
#         #Determine elemental abundance
#         if (self.Properties_dict['SII_HII_pn'] is not None) and (self.Properties_dict['SIII_HII_pn'] is not None) and (self.Properties_dict['SIV_HII_pn'] is not None):
#             self.Properties_dict['SI_HI_pn'] = self.Properties_dict['SII_HII_pn'] + self.Properties_dict['SIII_HII_pn'] + self.Properties_dict['SIV_HII_pn']
# 
#         if (self.Properties_dict['SII_HII_pn'] is not None) and (self.Properties_dict['SIII_HII_pn'] is not None): 
#             self.Properties_dict['SI_HI_pn'] = self.Properties_dict['SII_HII_pn'] + self.Properties_dict['SIII_HII_pn']
# 
#         return
#     
#     def oxygen_abundance_scheme(self):
#   
#         #Select the temperature and density for the procedure49
#         if self.Properties_dict['TOIII_pn'] is not None:
#             TOIII   = self.Properties_dict['TOIII_pn']
#             TOII    = self.Properties_dict['TOII_approxfrom_TOIII_pn']   
#         else:
#             TOIII   = self.Properties_dict['TOIII_approxfrom_TSIII_pn']
#             TOII    = self.Properties_dict['TSIII_pn']
#                 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = self.Properties_dict['nOII_pn']    
#      
#         if ne is not None:
#               
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)   
#             
#             #Calculate the O+1 ion abundance (If O2_3726A, O2_3729A are available we prioritize them)
#             if TOII is not None:
#                 
#                 if self.lines_dict.viewkeys() >= {'O2_3726A'}: #WE ASSUME THEY ARE BLENDED
#                     Ratio                                       = self.lines_dict['O2_3726A']
#                     self.Properties_dict['OII_HII_3279A_pn']    = self.O2_atom.getIonAbundance(int_ratio = Ratio, tem=TOII, den=ne, to_eval = 'L(3726)+L(3729)', Hbeta = self.Hbeta_flux) 
#                 elif  self.lines_dict.viewkeys() >= {'O2_3726A'}:
#                     Ratio                                       = self.lines_dict['O2_3726A']
#                     self.Properties_dict['OII_HII_3279A_pn']    = self.O2_atom.getIonAbundance(int_ratio = Ratio, tem=TOII, den=ne, to_eval = 'L(3726)+L(3729)', Hbeta = self.Hbeta_flux) 
#         
#                 if self.lines_dict.viewkeys() >= {'O2_7319A'}: 
#                     #Initial value
#                     Ratio                                       = self.lines_dict['O2_7319A'] #WE ASSUME THEY ARE BLENDED
#                     OII_HII_pn_7319_initial                     = self.O2_atom.getIonAbundance(int_ratio = Ratio, tem=TOII, den=ne, to_eval = 'L(7319)+L(7330)', Hbeta = self.Hbeta_flux)     
#                     #Correction for recombination contribution Liu2000
#                     Lines_Correction                            = (9.36 * power((TOII/10000), 0.44) * OII_HII_pn_7319_initial) * self.Hbeta_flux      
#                     Ratio                                       = self.lines_dict['O2_7319A'] - Lines_Correction
#                     self.Properties_dict['OII_HII_7319A_pn']    = self.O2_atom.getIonAbundance(int_ratio = Ratio, tem=TOII, den=ne, to_eval = 'L(7319)+L(7330)', Hbeta = self.Hbeta_flux) 
#    
#             #Calculate the O+2 ion abundance
#             if (TOIII is not None):         
#                 if self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A'}:
#                     Ratio = self.lines_dict['O3_4959A'] + self.lines_dict['O3_5007A']
#                     self.Properties_dict['OIII_HII_pn'] = self.O3_atom.getIonAbundance(int_ratio = Ratio, tem=TOIII, den=ne, to_eval = 'L(4959)+L(5007)', Hbeta = self.Hbeta_flux)                           
#                           
#         #Determine the O/H abundance
#         if (self.Properties_dict['OII_HII_3279A_pn'] is not None) and (self.Properties_dict['OIII_HII_pn'] is not None):
#             self.Properties_dict['OI_HI_pn'] = self.Properties_dict['OII_HII_3279A_pn'] + self.Properties_dict['OIII_HII_pn'] 
#             self.Properties_dict['OII_HII_pn'] = self.Properties_dict['OII_HII_3279A_pn']  
#         elif (self.Properties_dict['OII_HII_7319A_pn'] is not None) and (self.Properties_dict['OIII_HII_pn'] is not None):
#             self.Properties_dict['OI_HI_pn'] = self.Properties_dict['OII_HII_7319A_pn'] + self.Properties_dict['OIII_HII_pn']   
#             self.Properties_dict['OII_HII_pn'] = self.Properties_dict['OII_HII_7319A_pn']  
#         
#         return
#         
#     def nitrogen_abundance_scheme(self):
#         
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TSIII_pn'] is not None: 
#             Te  = self.Properties_dict['TSIII_pn'] #YOU CAN TEST THIS WITH OBJECT 51959-092
#         else:
#             Te  = self.Properties_dict['TNII_pn']
#               
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne  = self.Properties_dict['nSII_pn']
#         else:
#             ne  = self.Properties_dict['nOII_pn']
# 
#         if ne is not None:
# 
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
# 
#             #Calculate the N+1 ion abundance (If O2_3726A, O2_3729A are available we prioritize them)
#             if Te is not None:
#                                 
#                 if self.lines_dict.viewkeys() >= {'N2_6548A', 'N2_6584A'}:         
#                     Ratio                               = self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']
#                     self.Properties_dict['NII_HII_pn']  = self.N2_atom.getIonAbundance(int_ratio = Ratio, tem=Te, den=ne, to_eval = 'L(6548)+L(6584)', Hbeta = self.Hbeta_flux)                           
#         
#             
#         #Calculate the elemental abundance
#         if (self.Properties_dict['OI_HI_pn'] is not None) and (self.Properties_dict['NII_HII_pn'] is not None):
#             
#             #Standard method 
#             self.Properties_dict['NI_OI_pn'] = self.Properties_dict['NII_HII_pn'] / self.Properties_dict['OII_HII_pn']
#             NI_HI_initial = self.Properties_dict['NI_OI_pn'] * self.Properties_dict['OI_HI_pn']
# 
#             #Repeat calculation if 5755 line was observed to include the recombination contribution
#             if self.lines_dict.viewkeys() >= {'N2_5755A'}:         
#                 
#                 NIII_HI             = NI_HI_initial - self.Properties_dict['NII_HII_pn']
#                 Lines_Correction    = 3.19 * power((Te/10000), 0.30) * NIII_HI * self.Hbeta_flux      
#                 
#                 self.Properties_dict['TNII_pn'], nSII = self.diags.getCrossTemDen(diag_tem      = '[NII] 5755/6584+',
#                                                                                     diag_den    = '[SII] 6731/6716',
#                                                                                     value_tem   = (self.lines_dict['N2_5755A'] - Lines_Correction)/(self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']),
#                                                                                     value_den   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#                 
#                 Ratio                               = self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']
#                 self.Properties_dict['NII_HII_pn']  = self.N2_atom.getIonAbundance(int_ratio = Ratio, tem=self.Properties_dict['TNII_pn'], den=ne, to_eval = 'L(6548)+L(6584)', Hbeta = self.Hbeta_flux)                 
#                 self.Properties_dict['NI_OI_pn'] = self.Properties_dict['NII_HII_pn'] / self.Properties_dict['OII_HII_pn']
#                 self.Properties_dict['NI_HI_pn'] = self.Properties_dict['NI_OI_pn'] * self.Properties_dict['OI_HI_pn']
#             
#             else:
#                 self.Properties_dict['NI_HI_pn'] = NI_HI_initial
#     
#         return
#     
#     def helium_abundance_scheme_fab(self):
#         
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TOIII_approxfrom_TSIII_pn'] is not None:
#             TOIII   = self.Properties_dict['TOIII_approxfrom_TSIII_pn']
#         else:
#             TOIII   = self.Properties_dict['TOIII_pn']
#                 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn'] 
#         else:
#             ne = self.Properties_dict['nOII_pn']       
#                 
#         if (TOIII is not None) and (ne is not None):
#             
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)    
#             
#             #Calculate the He+1 abundance
#             if self.lines_dict.viewkeys() >= {'He1_4472A', 'He1_5876A', 'He1_6678A'}:
# 
#                 y0_I4471    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_4472A'], tem=TOIII, den=ne, wave = 4471, Hbeta = self.Hbeta_flux)
#                 y0_I5876    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_5876A'], tem=TOIII, den=ne, wave = 5876, Hbeta = self.Hbeta_flux)
#                 y0_I6678    = self.He1_atom.getIonAbundance(int_ratio = self.lines_dict['He1_6678A'], tem=TOIII, den=ne, wave = 6678, Hbeta = self.Hbeta_flux)
#                
#                 t4          = TOIII / 10000
# 
#                 D           = 1.0 + 3110.0 * power(t4,-0.51) * (1.0/ne)
#                 g_I4471     = 6.11 * power(t4,0.02)  * (exp(-4.544)) / D
#                 g_I5876     = (7.12   * power(t4,0.14)  * (exp(-3.776/t4)) + 1.47 * power(t4,-0.28) * (exp(-4.544/t4))) / D
#                 g_I6678     = (3.27   * power(t4,-0.41) * (exp(-3.777/t4)) + 0.49 * power(t4,-0.52) * (exp(-4.544/t4))) / D
#                 
#                 y_I4471     = y0_I4471 / (1.0 + g_I4471)
#                 y_I5876     = y0_I5876 / (1.0 + g_I5876)
#                 y_I6678     = y0_I6678 / (1.0 + g_I6678)
#                      
#                 self.Properties_dict['HeII_HII_pn'] = (3.0/5.0) * (y_I5876 + (1.0/3.0) * y_I4471 + (1.0/3.0) * y_I6678)     
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.Properties_dict['HeIII_HII_pn'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TOIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
#         #Calculate elemental abundance
#         if (self.Properties_dict['HeII_HII_pn'] is not None) and (self.Properties_dict['HeIII_HII_pn'] is not None):
#             self.Properties_dict['HeI_HI_pn'] = self.Properties_dict['HeII_HII_pn'] + self.Properties_dict['HeIII_HII_pn']
#         elif (self.Properties_dict['HeII_HII_pn'] is not None):
#             self.Properties_dict['HeI_HI_pn'] = self.Properties_dict['HeII_HII_pn']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.Properties_dict['HeI_HI_pn'] is not None):
#             #From Oxygen:
#             if self.Properties_dict['OI_HI_pn'] is not None:
#                 self.Properties_dict['Y_Mass_O_pn'] = (4 * self.Properties_dict['HeI_HI_pn'] * (1 - 20 * self.Properties_dict['OI_HI_pn'])) / (1 + 4 * self.Properties_dict['HeI_HI_pn']) 
#             #From Sulfur
#             if self.Properties_dict['SI_HI_pn'] is not None:
#                 self.Properties_dict['Y_Mass_S_pn'] = (4 * self.Properties_dict['HeI_HI_pn'] * (1 - 20 * self.OI_SI * self.Properties_dict['SI_HI_pn'])) / (1 + 4 * self.Properties_dict['HeI_HI_pn']) 
#     
#         return
# 
#     def helium_abundance_scheme_Oxygen(self, lineslog_frame):
# 
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TOIII_pn'] is not None:
#             TOIII   = self.Properties_dict['TOIII_pn']
#         else:
#             TOIII   = None
#                 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = None
#                      
#         if (TOIII is not None) and (ne is not None):
# 
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             
#             HeI_indices = (lineslog_frame.Ion.str.contains('HeI_')) & (lineslog_frame.index != 'He1_8446A')  & (lineslog_frame.index != 'He1_7818A')
#             HeI_labels  = lineslog_frame.loc[HeI_indices].index.values
#             HeI_waves   = lineslog_frame.loc[HeI_indices].TheoWavelength.values
#             
#             Flux_Hbeta  = self.lines_dict['H1_4861A']
#             Emis_Hbeta  = self.H1_atom.getEmissivity(tem=TOIII, den=ne, label = '4_2', product = False)
#     
#             for i in range(len(HeI_labels)):
#                 
#                 print 'Esta es', HeI_labels[i]
#                 print 'HeI_waves[i]', round(HeI_waves[i])
#                 
#                 line_relative_Flux          = self.lines_dict[HeI_labels[i]] / Flux_Hbeta
#                 line_relative_emissivity    = self.He1_atom.getEmissivity(tem = TOIII, den = ne, wave = round(HeI_waves[i]), product = False) / Emis_Hbeta
#                  
#                 if i == 0:
#                     matrix_HeI_fluxes       = copy(line_relative_Flux)
#                     matrix_HeI_emis         = copy(line_relative_emissivity)
#          
#                 else:
#                     matrix_HeI_fluxes       = vstack((matrix_HeI_fluxes, line_relative_Flux)) 
#                     matrix_HeI_emis         = vstack((matrix_HeI_emis,   line_relative_emissivity))    
#     
#             matrix_HeI_fluxes = transpose(matrix_HeI_fluxes)
#             matrix_HeI_emis   = transpose(matrix_HeI_emis)           
#             
#             #Perform the fits
#             params = Parameters()
#             params.add('Y', value = 0.01)
#             HeII_HII_array = zeros(len(matrix_HeI_fluxes))
#             HeII_HII_error = zeros(len(matrix_HeI_fluxes))
#            
#             for i in range(len(matrix_HeI_fluxes)):
#                 fit_Output = lmfit_minimmize(residual_Y_v3, params, args=(matrix_HeI_emis[i], matrix_HeI_fluxes[i]))
#                 HeII_HII_array[i] = fit_Output.params['Y'].value
#                 HeII_HII_error[i] = fit_Output.params['Y'].stderr
#             
#             #NO SUMANDO LOS ERRORES CORRECTOS?
# #             self.Properties_dict['HeII_HII_from_O_pn'] = unumpy.uarray(HeII_HII_array, mean(HeII_HII_error))
# #             self.Properties_dict['HeII_HII_from_O_pn'] = HeII_HII_array
#             self.Properties_dict['HeII_HII_from_O_pn'] = random.normal(mean(HeII_HII_array), mean(HeII_HII_error), size = self.MC_array_len)
# 
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.Properties_dict['HeIII_HII_from_O_pn'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TOIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
# 
#         #Calculate elemental abundance
#         if (self.Properties_dict['HeII_HII_from_O_pn'] is not None) and (self.Properties_dict['HeIII_HII_from_O_pn'] is not None):
#             self.Properties_dict['HeI_HI_from_O_pn'] = self.Properties_dict['HeII_HII_from_O_pn'] + self.Properties_dict['HeIII_HII_from_O_pn']
#         elif (self.Properties_dict['HeII_HII_from_O_pn'] is not None):
#             self.Properties_dict['HeI_HI_from_O_pn'] = self.Properties_dict['HeII_HII_from_O_pn']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.Properties_dict['HeI_HI_from_O_pn'] is not None):
#             #From Oxygen:
#             if self.Properties_dict['OI_HI_pn'] is not None:
#                 self.Properties_dict['Y_Mass_O_pn'] = (4 * self.Properties_dict['HeI_HI_from_O_pn'] * (1 - 20 * self.Properties_dict['OI_HI_pn'])) / (1 + 4 * self.Properties_dict['HeI_HI_from_O_pn']) 
# 
#         return
# 
#     def helium_abundance_scheme_Sulfur(self, lineslog_frame):
# 
#         #Select the temperature and density for the procedure
#         if self.Properties_dict['TOIII_approxfrom_TSIII_pn'] is not None:
#             TOIII   = self.Properties_dict['TOIII_approxfrom_TSIII_pn']
#         else:
#             TOIII   = None
#                 
#         if self.Properties_dict['nSII_pn'] is not None:
#             ne = self.Properties_dict['nSII_pn']
#         else:
#             ne = None
#                      
#         if (TOIII is not None) and (ne is not None):
# 
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             
#             HeI_indices = lineslog_frame.Ion.str.contains('HeI_') & (lineslog_frame.index != 'He1_8446A') & (lineslog_frame.index != 'He1_7818A')
#             HeI_labels  = lineslog_frame.loc[HeI_indices].index.values
#             HeI_waves   = lineslog_frame.loc[HeI_indices].TheoWavelength.values
#             
#             Flux_Hbeta  = self.lines_dict['H1_4861A']
#             Emis_Hbeta  = self.H1_atom.getEmissivity(tem=TOIII, den=ne, label = '4_2', product = False)
#     
#             for i in range(len(HeI_labels)):
#                  
#                 line_relative_Flux          = self.lines_dict[HeI_labels[i]] / Flux_Hbeta
#                 line_relative_emissivity    = self.He1_atom.getEmissivity(tem = TOIII, den = ne, wave = round(HeI_waves[i]), product = False) / Emis_Hbeta
#                  
#                 if i == 0:
#                     matrix_HeI_fluxes       = copy(line_relative_Flux)
#                     matrix_HeI_emis         = copy(line_relative_emissivity)
#          
#                 else:
#                     matrix_HeI_fluxes       = vstack((matrix_HeI_fluxes, line_relative_Flux)) 
#                     matrix_HeI_emis         = vstack((matrix_HeI_emis,   line_relative_emissivity))    
#     
#             matrix_HeI_fluxes = transpose(matrix_HeI_fluxes)
#             matrix_HeI_emis   = transpose(matrix_HeI_emis)           
#             
#             #Perform the fits
#             params = Parameters()
#             params.add('Y', value = 0.01)
#             HeII_HII_array = zeros(len(matrix_HeI_fluxes))
#             HeII_HII_error = zeros(len(matrix_HeI_fluxes))
#            
#             for i in range(len(matrix_HeI_fluxes)):
#                 fit_Output = lmfit_minimmize(residual_Y_v3, params, args=(matrix_HeI_emis[i], matrix_HeI_fluxes[i]))
#                 HeII_HII_array[i] = fit_Output.params['Y'].value
#                 HeII_HII_error[i] = fit_Output.params['Y'].stderr
#                 
#             #NO SUMANDO LOS ERRORES CORRECTOS?
# #             self.Properties_dict['HeII_HII_from_S_pn'] = unumpy.uarray(HeII_HII_array, mean(HeII_HII_error))
# #             self.Properties_dict['HeII_HII_from_S_pn'] = HeII_HII_array
#             self.Properties_dict['HeII_HII_from_S_pn'] = random.normal(mean(HeII_HII_array), mean(HeII_HII_error), size = self.MC_array_len)
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.Properties_dict['HeIII_HII_from_S_pn'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TOIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
# 
#         #Calculate elemental abundance
#         if (self.Properties_dict['HeII_HII_from_S_pn'] is not None) and (self.Properties_dict['HeIII_HII_from_S_pn'] is not None):
#             self.Properties_dict['HeI_HI_from_S_pn'] = self.Properties_dict['HeII_HII_from_S_pn'] + self.Properties_dict['HeIII_HII_from_S_pn']
#         elif (self.Properties_dict['HeII_HII_from_S_pn'] is not None):
#             self.Properties_dict['HeI_HI_from_S_pn'] = self.Properties_dict['HeII_HII_from_S_pn']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if (self.Properties_dict['HeI_HI_from_S_pn'] is not None):
#             #From Oxygen:
#             if self.Properties_dict['SI_HI_pn'] is not None:
#                 self.Properties_dict['Y_Mass_S_pn'] = (4 * self.Properties_dict['HeI_HI_from_S_pn'] * (1 - 20 * self.OI_SI * self.Properties_dict['SI_HI_pn'])) / (1 + 4 * self.Properties_dict['HeI_HI_from_S_pn']) 
# 
#         return
# 
#     def store_dict_ObjLog(self, FileFolder, CodeName, verbose = False):
#         
#         keys, values = self.Properties_dict.keys(), self.Properties_dict.values()
#         
#         for i in range(len(keys)):
#             if values[i] is not None:
#                 #This removes the nan elements from the treatment values[i][~isnan(values[i])]
#                 Median_Value    = median(values[i][~isnan(values[i])])
#                 error           = std(values[i][~isnan(values[i])])  
#                 self.SaveParameter_ObjLog(CodeName, FileFolder, keys[i], ufloat(Median_Value, error))
# 
#                 if verbose == True:
#                     print i, keys[i], Median_Value, error
#                 
#             else:
#                 self.SaveParameter_ObjLog(CodeName, FileFolder, keys[i], values[i])
#                 if verbose == True:
#                     print i, keys[i], values[i]  
#         return
#     
#     def SaveParameter_ObjLog(self, CodeName, FileFolder, Parameter, Magnitude, Error = None, Assumption = None, Log_extension = None):
#                 
#         if Log_extension == None:
#             Log_extension = self.ObjectLog_extension
#             ObjLog_Address                                                          = FileFolder + CodeName + Log_extension
#             ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments    = loadtxt(ObjLog_Address, dtype=str, usecols = [0,1,2,3], unpack = True)
#         
#         else:
#             ObjLog_Address                                                          = FileFolder + CodeName
#             ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments    = loadtxt(ObjLog_Address, dtype=str, usecols = [0,1,2,3], unpack = True)
# 
#         #Saving the new value for the given parameter
#         Parameter_Index                                                         = where(ObjLog_Parameters == Parameter)[0][0]
# 
#         #Save the error in the case of a nparray
#         if isinstance(Magnitude, UFloat) and (Error == None):
#             ObjLog_Magnitudes[Parameter_Index]                                  = Magnitude.nominal_value
#             ObjLog_Errors[Parameter_Index]                                      = Magnitude.std_dev
#         elif Error != None:
#             ObjLog_Magnitudes[Parameter_Index]                                  = str(Magnitude)
#             ObjLog_Errors[Parameter_Index]                                      = str(Error)
#         elif Magnitude != None:
#             ObjLog_Magnitudes[Parameter_Index]                                  = str(Magnitude)
#         elif Magnitude == None:
#             ObjLog_Magnitudes[Parameter_Index]                                  = 'None'
#             ObjLog_Errors[Parameter_Index]                                      = '-'
#             
#         #Saving the text file
#         savetxt(ObjLog_Address, transpose((ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments)), fmt='%s')
#     
#         return   

#     def determine_electron_parameters_old(self):
#                 
# #-------Starting with the Sulfur
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:    
#             
#             self.abunData['R_SII'] = self.lines_dict['S2_6716A']/self.lines_dict['S2_6731A']
#                       
#             #Check TSIII lines
#             if self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_9531A', 'S3_6312A'}:
#                 
#                 if self.low_density_obj:
#                     temp_diag   = 'L(6312)/(L(9069)+L(9531))' 
#                     temp_ratio  = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] + self.lines_dict['S3_9531A'])                    
#                     Te          = self.S3_atom.getTemDen(temp_ratio, den = self.den_dist, to_eval = temp_diag)
#                     ne          = self.den_dist
#                 else:
#                     den_diag    = '[SII] 6731/6716'
#                     den_ratio   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']
#                     temp_diag   = '[SIII] 6312/9200+' 
#                     temp_ratio  = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] + self.lines_dict['S3_9531A'])
#                     Te, ne      = self.diags.getCrossTemDen(diag_tem = temp_diag, diag_den = den_diag, value_tem = temp_ratio, value_den = den_ratio) 
# 
#                 self.abunData['TeSIII'],  self.abunData['neSII'] = Te, ne
#                     
#             #Missing S3_9531A line
#             elif self.lines_dict.viewkeys() >= {'S3_9069A', 'S3_6312A'}:
#                 
#                 if self.low_density_obj:
#                     temp_diag   = 'L(6312)/(L(9069)+L(9531))' 
#                     temp_ratio  = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] * (1 + self.S3_9000_ratio))                 
#                     Te          = self.S3_atom.getTemDen(temp_ratio, den = self.den_dist, to_eval = temp_diag)
#                     ne          =  self.den_dist
#                 else:         
#                     den_diag    = '[SII] 6731/6716'
#                     den_ratio   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']
#                     temp_diag   = '[SIII] 6312/9200+' 
#                     temp_ratio  = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] * (1 + self.S3_9000_ratio))
#                     Te, ne      = self.diags.getCrossTemDen(diag_tem = temp_diag, diag_den = den_diag, value_tem = temp_ratio, value_den = den_ratio) 
#                 
#                 self.abunData['TeSIII'],  self.abunData['neSII'] = Te, ne
#                 
#             #Missing S3_9069A line
#             elif self.lines_dict.viewkeys() >= {'S3_9531A', 'S3_6312A'}:
#                 
#                 if self.low_density_obj:
#                     temp_diag   = 'L(6312)/(L(9069)+L(9531))' 
#                     temp_ratio  = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9069A'] * (1 + self.S3_9000_ratio))                 
#                     Te          = self.S3_atom.getTemDen(temp_ratio, den = self.den_dist, to_eval = temp_diag)
#                     ne          =  self.den_dist
#                     
#                 else:         
#                     den_diag    = '[SII] 6731/6716'
#                     den_ratio   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']
#                     temp_diag   = '[SIII] 6312/9200+' 
#                     temp_ratio  = self.lines_dict['S3_6312A']/(self.lines_dict['S3_9531A'] * (1 + 1/self.S3_9000_ratio))
#                     Te, ne      = self.diags.getCrossTemDen(diag_tem = temp_diag, diag_den = den_diag, value_tem = temp_ratio, value_den = den_ratio)                
#                 
#                 self.abunData['TeSIII'],  self.abunData['neSII'] = Te, ne
#             
#             #Missing the sulfur lines we use the ones from oxygen
#             elif self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A', 'O3_4363A'}:     
#                 
#                 if self.low_density_obj:
#                     temp_diag   = 'L(4363)/(L(5007)+L(4959))'
#                     temp_ratio  = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A'])               
#                     Te          = self.O3_atom.getTemDen(temp_ratio, den = self.den_dist, to_eval = temp_diag)
#                     ne          =  self.den_dist
#                     
#                 else:         
#                     den_diag    = '[SII] 6731/6716'
#                     den_ratio   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']
#                     temp_diag   = '[OIII] 4363/5007+' 
#                     temp_ratio  = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A'])
#                     Te, ne      = self.diags.getCrossTemDen(diag_tem = temp_diag, diag_den = den_diag, value_tem = temp_ratio, value_den = den_ratio)  
# 
#                 self.abunData['TeSIII_from_TeOIII'] = ((Te/10000 + 0.0846) / 1.0807) * 10000
#                 self.abunData['neSII'] = ne
#                                 
# #             #Try to measure the TSII temperature
# #             if self.lines_dict.viewkeys() >= {'S2_4069A', 'S2_4076A'}:
# #                 
# #                 if self.low_density_obj:
# #                     temp_diag   = 'L(4069)/L(4076)'
# #                     temp_ratio  = self.lines_dict['S2_4069A']/self.lines_dict['S2_4076A']           
# #                     Te          = self.S2_atom.getTemDen(temp_ratio, den = self.den_dist, to_eval = temp_diag)
# #                     ne          = self.den_dist
# #                     
# #                 else:        
# #                     den_diag    = '[SII] 6731/6716'
# #                     den_ratio   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']
# #                     temp_diag   = '[SII] 4069/4076' 
# #                     temp_ratio  = self.lines_dict['S2_4069A']/self.lines_dict['S2_4076A']
# #                     Te, ne      = self.diags.getCrossTemDen(diag_tem = temp_diag, diag_den = den_diag, value_tem = temp_ratio, value_den = den_ratio)                 
# #                 
# #                 self.abunData['TeSII'] = Te
#             
# #-------Following with the oxygen
#         if self.lines_dict.viewkeys() >= {'S2_6716A', 'S2_6731A'}:
#                 
#             #Missing the sulfur lines we use the ones from oxygen
#             if self.lines_dict.viewkeys() >= {'O3_4959A', 'O3_5007A', 'O3_4363A'}:
#                 
#                 if self.low_density_obj:
#                     
#                     temp_diag   = 'L(4363)/(L(5007)+L(4959))'
#                     temp_ratio  = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A'])               
#                     Te          = self.O3_atom.getTemDen(temp_ratio, den = self.den_dist, to_eval = temp_diag)
#                     ne          = self.den_dist
#                 else:         
#                     den_diag    = '[SII] 6731/6716'
#                     den_ratio   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']
#                     temp_diag   = '[OIII] 4363/5007+' 
#                     temp_ratio  = self.lines_dict['O3_4363A']/(self.lines_dict['O3_5007A'] + self.lines_dict['O3_4959A'])
#                     Te, ne      = self.diags.getCrossTemDen(diag_tem = temp_diag, diag_den = den_diag, value_tem = temp_ratio, value_den = den_ratio)  
# 
#                 self.abunData['TeOIII'] = Te
#                 self.abunData['TeSIII_from_TeOIII'] = ((Te/10000 + 0.0846) / 1.0807) * 10000
#                 
#                                               
#             #Determine the TOIII from the sulfur temperature as a proxy
#             if 'TeSIII' in self.abunData:
#                 self.abunData['TeOIII_from_TeSIII'] = (1.0807 * self.abunData.TeSIII/10000 - 0.0846) * 10000
#             
#             #Determine the TOII value from the TOIII proxy
#             if 'TeOIII' in self.abunData:
#                 self.abunData['TeOII_from_TeOIII']  = (1.397 /( (1/(self.abunData.TeOIII/10000)) + 0.385) ) * 10000
# 
#             #Calculate the oxygen temperature cannot be measure use the SIII approximation       
#             elif self.lines_dict.viewkeys() >= {'O2_3726A', 'O2_3729A'}:
#                 if ('TeOII_from_TeOIII' in self.abunData):
#                     den_diag = '(L(3726))/(L(3729))'
#                     den_ratio = self.lines_dict['O2_3726A'] / self.lines_dict['O2_3729A']                          
#                     self.abunData['neOII'] = self.O2_atom.getTemDen(den_ratio, tem = self.abunData['TeOII_from_TeOIII'], to_eval = den_diag)
#                 
# #-------Following with the Nitrogen
#         if self.lines_dict.viewkeys() >= {'N2_5755A', 'N2_6548A', 'N2_6584A'}:         
# 
#             if self.low_density_obj:
#                 temp_diag   = '(L(6584) + L(6548)) / L(5755)'
#                 temp_ratio  = self.lines_dict['N2_5755A']/(self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A'])          
#                 Te          = self.N2_atom.getTemDen(temp_ratio, den = self.den_dist, to_eval = temp_diag)
#         
#             else:         
#                 den_diag    = '[SII] 6731/6716'
#                 den_ratio   = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A']
#                 temp_diag   = '[NII] 5755/6584+'
#                 temp_ratio  = self.lines_dict['N2_5755A']/(self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A'])
#                 Te, ne      = self.diags.getCrossTemDen(diag_tem = temp_diag, diag_den = den_diag, value_tem = temp_ratio, value_den = den_ratio) 
# 
#             self.abunData['TNII']=  Te
# 
#         if 'TeOIII' in self.abunData:
#             self.abunData['TeNII_from_TeOIII'] = (1.452 /( (1/(self.abunData.TeOIII/10000)) + 0.479)) * 10000
#            
#         #Avoid nan values
#         for variable in ['TeOII', 'TeSII', 'TeNII', 'TeOIII', 'TeSIII', 'TeOII_from_TeOIII', 'TeNII_from_TeOIII', 'TeSIII_from_TeOIII', 'TeOIII_from_TeSIII']:
#             if variable in self.abunData.index:
#                 if isnan(self.abunData[variable]).any():
#                     print variable, isnan(self.abunData[variable]).sum()
#                     mag, error = nanmean(self.abunData[variable]), nanstd(self.abunData[variable])
#                     print variable, mag, error
#                     self.abunData[variable] = random.normal(mag, error, size = self.MC_array_len)
#         return
#     def nitrogen_abundance_scheme_old(self):
# 
#         #Physical parameters
#         Te, ne = None, None
# 
#         #Select the temperature
#         if 'TeOIII' in self.abunData:
#             Te = self.abunData.TeOIII
#         elif 'TeOIII_from_TeSIII' in self.abunData:
#             Te = self.abunData.TeOIII_from_TeSIII
# 
#         #Select the density
#         if 'neSII' in self.abunData:
#             ne = self.abunData.neSII
#         elif 'neOII' in self.abunData:
#             ne = self.abunData.neOII
# 
#         if ne is not None:
# 
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
# 
#             #Calculate the N+1 ion abundance (If O2_3726A, O2_3729A are available we prioritize them)
#             if (Te is not None):
#                                 
#                 if self.lines_dict.viewkeys() >= {'N2_6548A', 'N2_6584A'}:         
#                     Ratio                               = self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']
#                     self.abunData['NII_HII']  = self.N2_atom.getIonAbundance(int_ratio = Ratio, tem=Te, den=ne, to_eval = 'L(6548)+L(6584)', Hbeta = self.Hbeta_flux)                           
#         
#             
#         #Calculate the elemental abundance
#         if ('OI_HI' in self.abunData) and ('NII_HII' in self.abunData):
#             
#             #Standard method 
#             self.abunData['NI_OI'] = self.abunData['NII_HII'] / self.abunData['OII_HII']
#             NI_HI_initial = self.abunData['NI_OI'] * self.abunData['OI_HI']
# 
#             #Repeat calculation if 5755 line was observed to include the recombination contribution
#             if self.lines_dict.viewkeys() >= {'N2_5755A'}:         
#                 
#                 NIII_HI             = NI_HI_initial - self.abunData['NII_HII']
#                 Lines_Correction    = 3.19 * power((Te/10000), 0.30) * NIII_HI * self.Hbeta_flux      
#                 
#                 self.abunData['TNII'], nSII = self.diags.getCrossTemDen(diag_tem = '[NII] 5755/6584+',
#                                                                         diag_den  = '[SII] 6731/6716',
#                                                                         value_tem = (self.lines_dict['N2_5755A'] - Lines_Correction)/(self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']),
#                                                                         value_den = self.lines_dict['S2_6731A']/self.lines_dict['S2_6716A'])
#                 
#                 Ratio = self.lines_dict['N2_6548A'] + self.lines_dict['N2_6584A']
#                 self.abunData['NII_HII']  = self.N2_atom.getIonAbundance(int_ratio = Ratio, tem=self.abunData['TNII'], den=ne, to_eval = 'L(6548)+L(6584)', Hbeta = self.Hbeta_flux)                 
#                 self.abunData['NI_OI'] = self.abunData['NII_HII'] / self.abunData['OII_HII']
#                 self.abunData['NI_HI'] = self.abunData['NI_OI'] * self.abunData['OI_HI']
#             
#             else:
#                 self.abunData['NI_HI'] = NI_HI_initial
#     
#         return

#     def helium_abundance_scheme_Oxygen_old(self, lineslog_frame):
# 
#         #Physical parameters
#         TOIII, ne = None, None
# 
#         #Select the temperature
#         if 'TeOIII' in self.abunData:
#             TOIII = self.abunData.TeOIII
# 
#         #Select the density
#         if 'neSII' in self.abunData:
#             ne = self.abunData.neSII
#                          
#         if (TOIII is not None) and (ne is not None):
#                         
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             TOIII = self.high_temp_function(TOIII)
#                     
#             HeI_indices = (lineslog_frame.Ion.str.contains('HeI_')) & (lineslog_frame.index != 'He1_8446A')  & (lineslog_frame.index != 'He1_7818A') & (lineslog_frame.index != 'He1_5016A')
#             HeI_labels  = lineslog_frame.loc[HeI_indices].index.values
#             HeI_ions    = lineslog_frame.loc[HeI_indices].Ion.values
#             
#             Flux_Hbeta  = self.lines_dict['H1_4861A']
#             Emis_Hbeta  = self.H1_atom.getEmissivity(tem=TOIII, den=ne, label = '4_2', product = False)
#     
#             for i in range(len(HeI_labels)):
#             
#                 pyneb_code                  = float(HeI_ions[i][HeI_ions[i].find('_')+1:len(HeI_ions[i])]) 
#                 line_relative_Flux          = self.lines_dict[HeI_labels[i]] / Flux_Hbeta
#                 line_relative_emissivity    = self.He1_atom.getEmissivity(tem = TOIII, den = ne, wave = pyneb_code, product = False) / Emis_Hbeta
# 
#                 if isnan(line_relative_emissivity).any():
#                     print HeI_labels[i], isnan(line_relative_emissivity).sum()
#                     mag, error = nanmean(line_relative_emissivity), nanstd(line_relative_emissivity)
#                     line_relative_emissivity = random.normal(mag, error, size = self.MC_array_len)
# 
#                 if i == 0:
#                     matrix_HeI_fluxes       = copy(line_relative_Flux)
#                     matrix_HeI_emis         = copy(line_relative_emissivity)
#          
#                 else:
#                     matrix_HeI_fluxes       = vstack((matrix_HeI_fluxes, line_relative_Flux)) 
#                     matrix_HeI_emis         = vstack((matrix_HeI_emis,   line_relative_emissivity))    
# 
# 
# 
# 
# 
# 
# 
#               
# #                 print HeI_labels[i]
# #                 print 'fluxes', line_relative_Flux
# #                 print 'Emissivities', line_relative_emissivity
# #                 print
# #                 print
#      
#             matrix_HeI_fluxes = transpose(matrix_HeI_fluxes)
#             matrix_HeI_emis   = transpose(matrix_HeI_emis)           
#                         
#             #Perform the fits
#             params = Parameters()
#             params.add('Y', value = 0.01)
#             HeII_HII_array = zeros(len(matrix_HeI_fluxes))
#             HeII_HII_error = zeros(len(matrix_HeI_fluxes))
#            
#             for i in range(len(matrix_HeI_fluxes)):
#                 fit_Output = lmfit_minimmize(residual_Y_v3, params, args=(matrix_HeI_emis[i], matrix_HeI_fluxes[i]))
#                 HeII_HII_array[i] = fit_Output.params['Y'].value
#                 HeII_HII_error[i] = fit_Output.params['Y'].stderr
#             
#             #NO SUMANDO LOS ERRORES CORRECTOS?
# #             self.abunData['HeII_HII_from_O'] = unumpy.uarray(HeII_HII_array, mean(HeII_HII_error))
# #             self.abunData['HeII_HII_from_O'] = HeII_HII_array
#             self.abunData['HeII_HII_from_O'] = random.normal(mean(HeII_HII_array), mean(HeII_HII_error), size = self.MC_array_len)
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.abunData['HeIII_HII_from_O'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TOIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
#         #Calculate elemental abundance
#         if ('HeII_HII_from_O' in self.abunData) and ('HeIII_HII_from_O' in self.abunData):
#             self.abunData['HeI_HI_from_O'] = self.abunData['HeII_HII_from_O'] + self.abunData['HeIII_HII_from_O']
#         elif ('HeII_HII_from_O' in self.abunData):
#             self.abunData['HeI_HI_from_O'] = self.abunData['HeII_HII_from_O']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if ('HeI_HI_from_O' in self.abunData) and ('OI_HI' in self.abunData):
#             self.abunData['Ymass_O'] = (4 * self.abunData['HeI_HI_from_O'] * (1 - 20 * self.abunData['OI_HI'])) / (1 + 4 * self.abunData['HeI_HI_from_O']) 
# 
#         return
# 
#     def helium_abundance_scheme_Sulfur_old(self, lineslog_frame):
# 
#         #Physical parameters
#         TSIII, ne = None, None
# 
#         #Select the temperature
#         if 'TeSIII' in self.abunData:
#             TSIII = self.abunData.TeSIII
# 
#         #Select the density
#         if 'neSII' in self.abunData:
#             ne = self.abunData.neSII
#                      
#         if (TSIII is not None) and (ne is not None):
# 
#             #Security check for objects with density below 75 cm-3
#             ne = self.low_density_function(ne)
#             TSIII = self.high_temp_function(TSIII)
#            
#             HeI_indices = lineslog_frame.Ion.str.contains('HeI_') & (lineslog_frame.index != 'He1_8446A') & (lineslog_frame.index != 'He1_7818A') & (lineslog_frame.index != 'He1_5016A')
#             HeI_labels  = lineslog_frame.loc[HeI_indices].index.values
#             HeI_ions    = lineslog_frame.loc[HeI_indices].Ion.values
#             
#             Flux_Hbeta  = self.lines_dict['H1_4861A']
#             Emis_Hbeta  = self.H1_atom.getEmissivity(tem=TSIII, den=ne, label = '4_2', product = False)
#     
#             for i in range(len(HeI_labels)):
# 
#                 pyneb_code                  = float(HeI_ions[i][HeI_ions[i].find('_')+1:len(HeI_ions[i])]) 
#                 line_relative_Flux          = self.lines_dict[HeI_labels[i]] / Flux_Hbeta
#                 line_relative_emissivity    = self.He1_atom.getEmissivity(tem = TSIII, den = ne, wave = pyneb_code, product = False) / Emis_Hbeta
# 
#                 if isnan(line_relative_emissivity).any():
#                     print HeI_labels[i], isnan(line_relative_emissivity).sum()
#                     mag, error = nanmean(line_relative_emissivity), nanstd(line_relative_emissivity)
#                     line_relative_emissivity = random.normal(mag, error, size = self.MC_array_len)
#                
#                 if i == 0:
#                     matrix_HeI_fluxes       = copy(line_relative_Flux)
#                     matrix_HeI_emis         = copy(line_relative_emissivity)
#          
#                 else:
#                     matrix_HeI_fluxes       = vstack((matrix_HeI_fluxes, line_relative_Flux)) 
#                     matrix_HeI_emis         = vstack((matrix_HeI_emis,   line_relative_emissivity))    
#     
#             matrix_HeI_fluxes = transpose(matrix_HeI_fluxes)
#             matrix_HeI_emis   = transpose(matrix_HeI_emis)           
#             
#             #Perform the fits
#             params = Parameters()
#             params.add('Y', value = 0.01)
#             HeII_HII_array = zeros(len(matrix_HeI_fluxes))
#             HeII_HII_error = zeros(len(matrix_HeI_fluxes))
#            
#             for i in range(len(matrix_HeI_fluxes)):
#                 fit_Output = lmfit_minimmize(residual_Y_v3, params, args=(matrix_HeI_emis[i], matrix_HeI_fluxes[i]))
#                 HeII_HII_array[i] = fit_Output.params['Y'].value
#                 HeII_HII_error[i] = fit_Output.params['Y'].stderr
#                 
#             #NO SUMANDO LOS ERRORES CORRECTOS?
# #             self.abunData['HeII_HII_from_S'] = unumpy.uarray(HeII_HII_array, mean(HeII_HII_error))
# #             self.abunData['HeII_HII_from_S'] = HeII_HII_array
#             self.abunData['HeII_HII_from_S'] = random.normal(mean(HeII_HII_array), mean(HeII_HII_error), size = self.MC_array_len)
# 
#             #Calculate the He+2 abundance 
#             if self.lines_dict.viewkeys() >= {'He2_4686A'}:
#                 self.abunData['HeIII_HII_from_S'] = self.He2_atom.getIonAbundance(int_ratio = self.lines_dict['He2_4686A'], tem=TSIII, den=ne, wave = 4685.6, Hbeta = self.Hbeta_flux)
# 
#         #Calculate elemental abundance
#         if ('HeII_HII_from_S' in self.abunData) and ('HeIII_HII_from_S' in self.abunData):
#             self.abunData['HeI_HI_from_S'] = self.abunData['HeII_HII_from_S'] + self.abunData['HeIII_HII_from_S']
#         
#         elif 'HeII_HII_from_S' in  self.abunData:
#             self.abunData['HeI_HI_from_S'] = self.abunData['HeII_HII_from_S']
#                         
#         #Proceed to get the Helium mass fraction Y
#         if ('HeI_HI_from_S' in self.abunData) and ('SI_HI' in self.abunData):
#             self.abunData['Ymass_S'] = (4 * self.abunData['HeI_HI_from_S'] * (1 - 20 * self.OI_SI * self.abunData['SI_HI'])) / (1 + 4 * self.abunData['HeI_HI_from_S']) 
# 
#         return
#

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
            
