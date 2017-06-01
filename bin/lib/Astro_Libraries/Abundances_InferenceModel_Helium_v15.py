from collections            import OrderedDict
from pyneb                  import atomicData, RecAtom
from uncertainties.unumpy   import nominal_values, std_devs
from numpy                  import array, loadtxt, genfromtxt, isnan, arange, insert, concatenate, power, exp, zeros, square
from Reddening_Corrections  import ReddeningLaws 
import pymc as pymc2

class Import_model_data():

    def __init__(self):
        
        self.Hydrogen_CollCoeff_TableAddress    = '/home/vital/workspace/dazer/bin/lib/Astro_Libraries/Neutral_Hydrogen_Collisional_Correction_coef.txt'
        self.Helium_CollCoeff_TableAddress      = '/home/vital/workspace/dazer/bin/lib/Astro_Libraries/Neutral_Helium_Collisional_Correction_coef.txt'
        self.Helium_OpticalDepth_TableAddress   = '/home/vital/workspace/dazer/bin/lib/Astro_Libraries/Helium_OpticalDepthFunction_Coefficients.txt'
                
        #Declare Hydrogen and Helium lines for the analysis
        self.posHydrogen_Lines      = ['H1_4102A',      'H1_4340A', 'H1_6563A']
        self.Hydrogen_Wavelengths   = array([4101.742,  4340.471,   6562.819])
        
        self.posHelium_Lines        = ['He1_3889A',  'He1_4026A',    'He1_4387A',    'He1_4471A',    'He1_4686A',    'He1_4714A',    'He1_4922A',   'He1_5876A',    'He1_6678A',   'He1_7065A',    'He1_7281A',    'He1_10830A']
        self.Helium_Wavelengths     = array([3889.0,  4026.0,         4387.0,         4471.0,         4686.0,         4714.0,         4922.0,         5876.0,         6678.0,       7065.0,         7281.0,        10830.0])

    def import_table_data(self, Address, Columns):
        
        Imported_Array  = genfromtxt(Address, dtype=float, usecols = Columns, skip_header = 2).T
        Datarray = Imported_Array[:,~isnan(Imported_Array).all(0)]
                
        return Datarray
     
    def import_coll_coeff_table(self, HydrogenLines, HeliumLines):
        
        Data_dict = OrderedDict()
                
        for i in range(len(HydrogenLines)):
            
            Line            = HydrogenLines[i]
            Data_Columns    = [0 + 3*i, 1 + 3*i, 2 + 3*i]
            Data_dict[Line] = self.import_table_data(self.Hydrogen_CollCoeff_TableAddress, Data_Columns)
                    
        return Data_dict
     
    def import_optical_depth_coeff_table(self, HeliumLines):
        
        Data_dict = OrderedDict()
        
        for i in range(len(HeliumLines)):
            
            Line = HeliumLines[i]
            
            if Line != self.H13889A_label:
                Data_dict[Line] = loadtxt(self.Helium_OpticalDepth_TableAddress, dtype = float, skiprows = 2, usecols = (i,))
                
            else:
                Data_dict[self.He3889_label] = loadtxt(self.Helium_OpticalDepth_TableAddress, dtype = float, skiprows = 2, usecols = (i,))
            
        return Data_dict

    def load_synthetic_data(self, model):

        self.obj_data = dict()
                 
        if model == 'Model3':
            
            self.Hbeta_Eqw                  = 250.0
            self.Hbeta_EqwErr               = 12.5
            self.Hbeta_hlambda              = 1.0
           
            self.obj_data['H_labels']       = ['H1_4102A',          'H1_4340A',         'H1_6563A']
            #self.obj_data['H_Flux']         = array([0.249665005739, 0.457121205187,    2.9720213093])
            #self.obj_data['H_error']        = self.obj_data['H_Flux'] * 0.01
            self.obj_data['H_wave']         = array([4101.742,    4340.471,       6562.819])
            self.obj_data['H_Eqw']          = array([50.0,        50.0,           300.0])
            self.obj_data['H_EqwErr']       = array([1.0,         1.0,            14])
            self.obj_data['H_hlambda']      = array([1.0,         1.0,            1.0])
            self.obj_data['H_pyneb_code']   = array(['6_2',       '5_2',          '3_2'])

            self.obj_data['He_labels']      = ['He1_3889A',  'He1_4026A',   'He1_4471A',   'He1_5876A',  'He1_6678A',    'He1_7065A',    'He1_10830A']
            #self.obj_data['He_Flux']        = array([0.102807348487,    0.0188761305051,     0.0411173796552,     0.128607653966,  0.0373393280888,    0.04339571902,      0.761144585562])
            #self.obj_data['He_error']       = self.obj_data['He_Flux'] * 0.02
            self.obj_data['He_wave']        = array([ 3889.0,       4026.0,         4471.0,         5876.0,      6678.0,        7065.0,         10830.0])
            self.obj_data['He_Eqw']         = array([10.0,          10.0,           10.0,           10.0,        10.0,          10.0,           10.0])
            self.obj_data['He_EqwErr']      = array([1.0,                 1.0,            1.0,            1.0,         1.0,           1.0,            1.0])
            self.obj_data['He_hlambda']     = array([1.0,                 1.0,            1.0,            1.0,         1.0,           1.0,            1.0])
            self.obj_data['He_pyneb_code']  = array(['3889.0',            '4026.0',       '4471.0',       '5876.0',    '6678.0',      '7065.0',       '10830.0'])

            self.obj_data['y_plus']         = 0.085
            self.obj_data['n_e']            = 500.0
            self.obj_data['tau']            = 1.0
            self.obj_data['Te_0']           = 16000.0
            self.obj_data['cHbeta']         = 0.1
            self.obj_data['xi']             = 1.0  
            self.obj_data['TOIII']          = 17000.0
            self.obj_data['TOIII_error']    = self.obj_data['TOIII'] * 0.2
            self.obj_data['nSII']           = 500.0
            self.obj_data['nSII_error']     = 200.0        
                                 
        return

    def load_obs_data(self, lines_df, obj_series, extension_treat = '', Deblend_Check = True):
        
        #Empty dictionary to store the data
        self.obj_data = {}

        #--Define HBeta values in advance
        idx_Hbeta                       = lines_df.index == 'H1_4861A'
        self.Hbeta_Flux                 = nominal_values(lines_df.loc[idx_Hbeta, 'line_Flux'].values)
        self.Hbeta_error                = std_devs(lines_df.loc[idx_Hbeta, 'line_Flux'].values)
        self.Hbeta_Eqw                  = nominal_values(lines_df.loc[idx_Hbeta, 'line_Eqw'].values) #WARNING: in cases with removed stellar this eqw gives trouble
        self.Hbeta_EqwErr               = std_devs(lines_df.loc[idx_Hbeta, 'line_Eqw'].values)
        self.Hbeta_hlambda              = nominal_values(lines_df.loc[idx_Hbeta, 'Continuum_Median'].values)

        #Load Hydrogen lines data in the dictionary
        idcs_H                          = lines_df.index.isin(self.posHydrogen_Lines)   #Hydrogen lines
        self.obj_data['H_labels']       = lines_df.loc[idcs_H].index.values
        self.obj_data['H_Ions']         = lines_df.loc[idcs_H, 'Ion'].values
        self.obj_data['H_Flux']         = nominal_values(lines_df.loc[idcs_H, 'line_Flux'].values)
        self.obj_data['H_error']        = std_devs(lines_df.loc[idcs_H, 'line_Flux'].values)
        self.obj_data['H_wave']         = lines_df.loc[idcs_H, 'TheoWavelength'].values
        self.obj_data['H_Eqw']          = nominal_values(lines_df.loc[idcs_H, 'line_Eqw'].values) #WARNING: in cases with removed stellar this eqw gives trouble
        self.obj_data['H_EqwErr']       = std_devs(lines_df.loc[idcs_H, 'line_Eqw'].values)
        self.obj_data['H_hlambda']      = nominal_values(lines_df.loc[idcs_H, 'Continuum_Median'].values)
        self.obj_data['H_pyneb_code']   = lines_df.loc[idcs_H, 'Ion'].str[lines_df.loc[idcs_H, 'Ion'].str.find('_'):]

        #Load Helium lines data
        idcs_He                         = lines_df.index.isin(self.posHelium_Lines)
        self.obj_data['He_labels']      = lines_df.loc[idcs_He].index.values
        self.obj_data['He_Ions']        = lines_df.loc[idcs_He, 'Ion'].values
        self.obj_data['He_Flux']        = nominal_values(lines_df.loc[idcs_He, 'line_Flux'].values)
        self.obj_data['He_error']       = std_devs(lines_df.loc[idcs_He, 'line_Flux'].values)
        self.obj_data['He_wave']        = lines_df.loc[idcs_He, 'TheoWavelength'].values
        self.obj_data['He_Eqw']         = nominal_values(lines_df.loc[idcs_He, 'line_Eqw'].values) #WARNING: in cases with removed stellar this eqw gives trouble
        self.obj_data['He_EqwErr']      = std_devs(lines_df.loc[idcs_He, 'line_Eqw'].values)
        self.obj_data['He_hlambda']     = nominal_values(lines_df.loc[idcs_He, 'Continuum_Median'].values)
        self.obj_data['He_pyneb_code']  = lines_df.loc[idcs_He].index.str[lines_df.loc[idcs_He].index.str.find('_')+1:-1]

        #Load physical parameters data       
        Thigh_key                       = obj_series['T_high']
        self.obj_data['TOIII']          = obj_series[Thigh_key + extension_treat].nominal_value
        self.obj_data['TOIII_error']    = obj_series[Thigh_key + extension_treat].std_dev
        self.obj_data['nSII']           = obj_series['neSII' + extension_treat].nominal_value
        self.obj_data['nSII_error']     = obj_series['neSII' + extension_treat].std_dev
        self.obj_data['cHbeta_obs']     = obj_series['cHbeta' + extension_treat].nominal_value

        #Check if He1_3889A is on the observations (In the current analysis it appears as 'H1_3889A')
        if 'H1_3889A' in lines_df.index:
            idx_He3889A = (lines_df.index == 'H1_4861A')
            
            #This should be calculated using the error
            He3889A_debFlux = self.deblend_He3889A(Te = self.obj_data['TOIII'], ne = self.obj_data['nSII'], cHbeta = self.obj_data['cHbeta_obs'])
    
            insert(self.obj_data['He_labels'],      0, lines_df.loc[idx_He3889A].index.values)
            insert(self.obj_data['He_Ions'],        0, lines_df.loc[idx_He3889A, 'Ion'].values)
            insert(self.obj_data['He_Flux'],        0, He3889A_debFlux.nominal_value)
            insert(self.obj_data['He_error'],       0, He3889A_debFlux.std_dev)
            insert(self.obj_data['He_wave'],        0, lines_df.loc[idx_He3889A, 'TheoWavelength'].values)
            insert(self.obj_data['He_Eqw'],         0, nominal_values(lines_df.loc[idx_He3889A, 'line_Eqw'].values)) #WARNING: in cases with removed stellar this eqw gives trouble
            insert(self.obj_data['He_EqwErr'],      0, std_devs(lines_df.loc[idx_He3889A, 'line_Eqw'].values))
            insert(self.obj_data['He_hlambda'],     0, nominal_values(lines_df.loc[idx_He3889A, 'Continuum_Median'].values))
            insert(self.obj_data['He_pyneb_code'],  0, lines_df.loc[idx_He3889A].index.str[lines_df.loc[idx_He3889A].index.str.find('_')+1:-1])
        
        return

    def prepare_run_data(self, norm_by_Hbeta = False, deblend_Check = True, red_curve = 'G03', Rv = 3.4):

        #Variables to speed up the code
        self.nHydrogen          = len(self.obj_data['H_labels'])
        self.nHydrogen_range    = arange(self.nHydrogen)
        
        self.nHelium            = len(self.obj_data['He_labels'])
        self.nHelium_range      = arange(self.nHelium)
        self.idx_meet           = self.nHydrogen + 1
       
        #Load the reddening parameters
        H_xX        = self.reddening_Xx(self.obj_data['H_wave'], red_curve, Rv)
        He_xX       = self.reddening_Xx(self.obj_data['He_wave'], red_curve, Rv)
        Hbeta_xX    = self.reddening_Xx(array([self.Hbeta_wave]), red_curve, Rv)[0]
        
        #Dictionary with the vectors
        self.data_dic = {'Emissivity_H_vector'  : zeros(self.nHydrogen),
                         'Emissivity_He_vector' : zeros(self.nHelium),
                         'Kalpha_vector'        : zeros(self.nHydrogen),
                         'ftau_He_vector'       : zeros(self.nHelium),
                         'flambda_H_vector'     : H_xX/Hbeta_xX - 1.0,
                         'flambda_He_vector'    : He_xX/Hbeta_xX - 1.0}    
                
        #Normalize the fluxes by Hbeta if necessary
        if norm_by_Hbeta: 
            self.obj_data['H_Flux']   = self.obj_data['H_Flux']     / self.Hbeta_Flux
            self.obj_data['He_Flux']  = self.obj_data['He_Flux']    / self.Hbeta_Flux
            self.obj_data['H_error']  = self.obj_data['H_error']    / self.Hbeta_Flux
            self.obj_data['He_error'] = self.obj_data['He_error']   / self.Hbeta_Flux
        
        return

class Recombination_FluxCalculation(ReddeningLaws, Import_model_data):
    
    def __init__(self):
        
        #Import sub classes
        ReddeningLaws.__init__(self)
        Import_model_data.__init__(self)
    
        #Define indexes and labels to speed up the code
        self.Hbeta_label            = 'H1_4861A'
        self.Hbeta_wave             = 4862.683
        self.Hbeta_pynebCode        = '4_2'
        self.H13889A_label          = 'H1_3889A'          
        self.He3889_label           = 'He1_3889A'
        self.He3889_Check           = None
        
        #Set up the right emissivities
        atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
                
        #Declare pyneb Hydrogen and Helium atoms to calculate emissivities
        self.H1                     = RecAtom('H', 1)
        self.He1                    = RecAtom('He', 1) 
    
        #Make sure we are using the right helium emissivities
        print '--Helium emissivities: '
        self.He1.printSources()

        #Import collisional coefficients table #We add the Hbeta to get its coefficients
        self.Coef_Kalpha_dict       = self.import_coll_coeff_table(self.posHydrogen_Lines + [self.Hbeta_label], None)

        #Import Optical depth function
        self.Coef_ftau_dict         = self.import_optical_depth_coeff_table(self.posHelium_Lines) 
                                                                               
    def calculate_synthFluxes(self, model):
        
        #Get physical parameters from model HII region
        self.load_synthetic_data(model = model)
        
        #Prepare the data for the run
        self.prepare_run_data(norm_by_Hbeta = False, deblend_Check=False)
             
        if '_abs' not in model:
            
            self.obj_data['H_Flux']     = self.hydrogen_theo_flux(self.obj_data['Te_0'], self.obj_data['n_e'], self.obj_data['xi'], self.obj_data['cHbeta'])
            self.obj_data['He_Flux']    = self.helium_theo_flux(self.obj_data['Te_0'], self.obj_data['n_e'], self.obj_data['tau'], self.obj_data['xi'], self.obj_data['cHbeta'], self.obj_data['y_plus'])
            
            self.obj_data['H_error']    = self.obj_data['H_Flux'] * 0.01
            self.obj_data['He_error']   = self.obj_data['He_Flux'] * 0.02
    
            print '\nInput Parameters'
            print '-y_plus',   self.obj_data['y_plus']
            print '-n_e',      self.obj_data['n_e']
            print '-Te_0',     self.obj_data['Te_0']
            print '-cHbeta',   self.obj_data['cHbeta']
            print '-tau',      self.obj_data['tau']
            print '-xi',       self.obj_data['xi']
            print '-TOIII',    self.obj_data['TOIII'],'+/-',self.obj_data['TOIII_error']
            print '-nSII',     self.obj_data['nSII'],'+/-',self.obj_data['nSII']
            
        print '\nHydrogen emissivities'
        for i in self.nHydrogen_range:
            print '-{} {}'.format(self.obj_data['H_labels'][i], self.obj_data['H_Flux'][i])

        print '\nHelium emissivities'
        for i in self.nHelium_range:
            print '-{} {}'.format(self.obj_data['He_labels'][i], self.obj_data['He_Flux'][i])
        
        #Fluxes combined to speed the process
        self.H_He_Obs_Flux  = concatenate([self.obj_data['H_Flux'], self.obj_data['He_Flux']])                        
        self.H_He_Obs_Error = concatenate([self.obj_data['H_error'], self.obj_data['He_error']])  
        self.H_He_Theo_Flux = zeros(self.nHydrogen_range + self.nHelium_range)
        
        return
  
    def calculate_H_params(self, Te, ne):    
        
        #Calculate in advance the T_4 parameter
        T_4 = Te / 10000.0
        
        #Calculate the hbeta parameters
        self.Hbeta_Kalpha   = self.Kalpha_Ratio_H(T_4 = T_4, H_label = self.Hbeta_label)
        self.Hbeta_emis     = self.H1.getEmissivity(Te, ne, label = self.Hbeta_pynebCode)

        #Calculate physical parameters for Hydrogen lines
        for i in self.nHydrogen_range:
            self.data_dic['Emissivity_H_vector'][i] = self.H1.getEmissivity(Te, ne, label = self.obj_data['H_pyneb_code'][i])
            self.data_dic['Kalpha_vector'][i]       = self.Kalpha_Ratio_H(T_4 = T_4, H_label = self.obj_data['H_labels'][i])
                        
        return 
        
    def calculate_He_params(self, Te, ne, tau):
    
        #Calculate in advance the T_4 parameter
        T_4 = Te / 10000.0
            
        #Calculate physical parameters for Helium lines
        for i in self.nHelium_range:
            self.data_dic['Emissivity_He_vector'][i]    = self.He1.getEmissivity(Te, ne, label = self.obj_data['He_pyneb_code'][i])
            self.data_dic['ftau_He_vector'][i]          = self.OpticalDepth_He(tau = tau, T_4 = T_4, ne = ne, He_label = self.obj_data['He_labels'][i])    
                    
        return

    def hydrogen_theo_flux(self, Te, ne, xi, cHbeta):
        
        #Hydrogen calculation parameters
        self.calculate_H_params(Te, ne)

        #Calculate the emissivities for each of the lines for the given temperature and density
        emis_module = self.data_dic['Emissivity_H_vector'] / self.Hbeta_emis
                
        #Calculate the Collisional excitation fraction
        CR_Module   = (1.0 + 0.0001* xi * self.data_dic['Kalpha_vector']) / (1.0 + 0.0001* xi * self.Hbeta_Kalpha)
                
        #Calculate the reddening component
        f_module    = power(10, -1 * self.data_dic['flambda_H_vector'] * cHbeta)

        #Calculate theoretical Hydrogen flux for each emission line
        H_Flux      = emis_module * CR_Module * f_module
                     
        return H_Flux

    def helium_theo_flux(self, Te, ne, tau, xi, cHbeta, y_plus):
        
        #Helium calculation parameters
        self.calculate_He_params(Te, ne, tau)
             
        #Calculate the emissivities for each of the lines for the given temperature and density
        emis_module = self.data_dic['Emissivity_He_vector'] / self.Hbeta_emis
                
        #Calculate the collisional excitation fraction
        CR_Module   = 1 / (1.0 + 0.0001* xi * self.Hbeta_Kalpha)

        #Calculate the reddening component
        f_module    = power(10, self.data_dic['flambda_He_vector'] * cHbeta)
             
        #Calculate theoretical Hydrogen flux for each emission line
        He_Flux     = y_plus * emis_module * self.data_dic['ftau_He_vector'] * CR_Module * f_module

        return He_Flux
                          
    def Kalpha_Ratio_H(self, T_4, H_label):
                      
        K_alpha_Ratio   = sum(self.Coef_Kalpha_dict[H_label][0] * exp(self.Coef_Kalpha_dict[H_label][1] / T_4) * power(T_4, self.Coef_Kalpha_dict[H_label][2]))
        
        return K_alpha_Ratio

    def OpticalDepth_He(self, tau, T_4, ne, He_label):
                        
        f_tau   = 1 + (tau/2) * (self.Coef_ftau_dict[He_label][0] + (self.Coef_ftau_dict[He_label][1] + self.Coef_ftau_dict[He_label][2]*ne + self.Coef_ftau_dict[He_label][3]*ne*ne) * T_4)

        return f_tau
        
class Inference_AbundanceModel(Recombination_FluxCalculation):

    def __init__(self):

        Recombination_FluxCalculation.__init__(self)

    def helium_abundance_model(self):
        
        y_plus      =   pymc2.Uniform(           'y_plus',  0.050,                       0.15)
        ne          =   pymc2.TruncatedNormal(   'ne',      self.obj_data['nSII'],       self.obj_data['nSII_error']**-2,    a = 0.0 ,   b = 1000.0)
        Te          =   pymc2.Normal(            'Te',      self.obj_data['TOIII'],      self.obj_data['TOIII_error']**-2)
        tau         =   pymc2.TruncatedNormal(   'tau',      0.75,                        0.5**-2,    a = 0.0,    b = 7.0)
        cHbeta      =   pymc2.TruncatedNormal(   'cHbeta',  0.15,                        0.05**-2,   a = 0.0,    b = 3.0)
        xi          =   pymc2.TruncatedNormal(   'xi',       1,                           200**-2,    a = 0.0,    b = 1000.0)

        @pymc2.deterministic #Calculate Hydrogen theoretical flux
        def det_HFlux_theo(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta):
            return self.hydrogen_theo_flux(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta)
           
        @pymc2.deterministic #Calculate Helium theoretical flux 
        def det_He_Flux_theo_nof(Te=Te, ne=ne, tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus):
            return self.helium_theo_flux(Te=Te, ne=ne,  tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus)
        
        @pymc2.deterministic #Combine theoretical fluxes into a single array
        def H_He_TheoFlux(HFlux_theo=det_HFlux_theo, HeFlux_theo=det_He_Flux_theo_nof):
            self.H_He_Theo_Flux[:self.nHydrogen_range] = HFlux_theo[:]
            self.H_He_Theo_Flux[self.nHydrogen_range:] = HeFlux_theo[:]
            return self.H_He_Theo_Flux
                    
        @pymc2.stochastic(observed=True) #Likelihood
        def Likelihood_model(value = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2
    
        @pymc2.deterministic() #Deterministic method to track the evolution of the chi:
        def ChiSq(H_He_ObsFlux = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
            return - chi_F / 2
      
        return locals()
   
    def helium_abundance_model_abs(self):

        print 'Physical parameters:'
        print self.obj_data['nSII'], self.obj_data['nSII_error']
        print self.obj_data['TOIII'], self.obj_data['TOIII_error']
        
        y_plus      =   pymc2.Uniform(           'He_abud',      0.050,                              0.15)
        ne          =   pymc2.TruncatedNormal(   'n_e',          self.obj_data['nSII'],             self.obj_data['nSII_error']**-2,    a = 0.0 ,   b = 1000.0)
        Te          =   pymc2.Normal(            'T_e',          self.obj_data['TOIII'],            self.obj_data['TOIII_error']**-2)
        tau         =   pymc2.TruncatedNormal(   'Tau',          0.75,                               0.5**-2,    a = 0.0,    b = 7.0)
        cHbeta      =   pymc2.TruncatedNormal(   'c_Hbeta',      0.15,                               0.05**-2,   a = 0.0,    b = 3.0)
        xi          =   pymc2.TruncatedNormal(   'Xi',           1,                                  200**-2,    a = 0.0,    b = 1000.0)
        a_H         =   pymc2.TruncatedNormal(   'abs_H',        1.0,                                1**-2,      a = 0.0,    b = 14.0)
        a_He        =   pymc2.TruncatedNormal(   'abs_He',       1.0,                                0.5**-2,    a = 0.0,    b = 5.0)
        Eqw_hbeta   =   pymc2.Normal(            'Eqw_hbeta',    self.EqW_dict['H1_4861A'],          self.EqWerror_dict['H1_4861A']**-2)

        #Calculate Hydrogen theoretical flux
        @pymc2.deterministic
        def det_HFlux_theo_abs(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta, a_H=a_H, Eqw_hbeta=Eqw_hbeta):
            return self.H_Flux_theo_abs(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta, a_H=a_H, Eqw_hbeta=Eqw_hbeta)
        
        #Calculate Helium theoretical flux    
        @pymc2.deterministic
        def det_He_Flux_theo_abs(Te=Te, ne=ne, tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus, a_H=a_H, a_He=a_He, Eqw_hbeta=Eqw_hbeta):
            return self.He_Flux_theo_abs(Te=Te, ne=ne,  tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus, a_H=a_H, a_He=a_He, Eqw_hbeta=Eqw_hbeta)
        
        #Combine theoretical fluxes into a single array
        @pymc2.deterministic
        def H_He_TheoFlux(HFlux_theo=det_HFlux_theo_abs, HeFlux_theo=det_He_Flux_theo_abs):
            return concatenate([HFlux_theo, HeFlux_theo])
                
        #Likelihood
        @pymc2.stochastic(observed=True)
        def Likelihood_model(value = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2
        
        #Deterministic method to track the evolution of the chi:
        
        @pymc2.deterministic()
        def ChiSq(H_He_ObsFlux = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
            return - chi_F / 2
       
        return locals()

class Run_MCMC(Inference_AbundanceModel):

    def __init__(self):
        
        #Import parent classes
        Inference_AbundanceModel.__init__(self)
            
    def select_inference_model(self, model_code):
        
        #Declare inference model
        if '_abs' not in model_code:
            self.inf_dict = self.helium_abundance_model()
        else:
            self.inf_dict = self.helium_abundance_model_abs()
                    
    def run_pymc2(self, iterations = 10000, burn_phase = 1500, thining=2, ext='', default_db_name = 'HeAbund_pymc2', default_db_folder = '/home/vital/Astrodata/Inferece_output/', variables_list = None):
        
        #Ideal cond1: iter=30000, burn=5000, thin=10
        #Ideal cond2: iter=10000, burn=1500, thin=2
        
        #Prefit:
        print '\n--Starting Fmin_powell prefit'
        self.MAP_Model  = pymc2.MAP(self.inf_dict)
        self.MAP_Model.fit(method = 'fmin_powell')
        
        #Print prefit data
        self.display_run_data(self.MAP_Model, variables_list)
        
        #Declare database name
        self.db_address = '{}{}_it{}_burn{}_thin{}_{}'.format(default_db_folder, default_db_name, iterations, burn_phase, thining, ext)
        
        #Launch sample
        self.pymc2_M    = pymc2.MCMC(self.MAP_Model.variables, db = 'pickle', dbname =  self.db_address)
        self.pymc2_M.sample(iter=iterations, burn=burn_phase, thin=thining)
        
        #Save the output csv mean data
        if variables_list != None:
            print '--Saving results in csv'
            self.csv_address = self.db_address + '_Parameters'
            self.pymc2_M.write_csv(self.csv_address, variables=variables_list)
            
        #Print again the             
        self.display_run_data(self.MAP_Model, variables_list)

        #Close the database
        self.pymc2_M.db.close()            
                
    def display_run_data(self, database, variables_list):
        
        for param in variables_list:
            param_entry = getattr(database, param, None)
            if param_entry is not None:
                print '-{} {}'.format(param, param_entry.value)
                
#         #Display prefit results:
#         if variables_list != None:
#                         
#             #Loop through the variales we want to display
#             for param in variables_list:
#                 print param, database.trace(param).stats()['mean']

        
                                                                                                        
