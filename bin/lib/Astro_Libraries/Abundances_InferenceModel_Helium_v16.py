from collections            import OrderedDict
from pyneb                  import atomicData, RecAtom, Atom
from uncertainties.unumpy   import nominal_values, std_devs
from Reddening_Corrections  import ReddeningLaws 
from numpy                  import array, loadtxt, genfromtxt, isnan, arange, insert, concatenate, power, exp, zeros, square, empty, percentile
import pymc as pymc2
import pymc3

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

        self.posSII_Lines           = ['S2_6716A', 'S2_6731A'] 
        self.SII_Wavelenths         = array([6716.44, 6730.81])
        
        self.posSIII_Lines          = ['S3_6312A', 'S3_9069A', 'S3_9531A']
        self.SIII_Wavelenths        = array([6312.06, 9068.6, 9531.1])      

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

        if model == 'Model3_sulfur':
            
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

            self.obj_data['S2_labels']     = ['S2_6716A', 'S2_6731A']
            self.obj_data['S2_wave']       = array([6716.44, 6730.81])
            self.obj_data['S2_pyneb_code'] = array([6716, 6730])
            self.obj_data['S3_labels']     = ['S3_6312A', 'S3_9069A', 'S3_9531A']
            self.obj_data['S3_wave']       = array([6312.06, 9068.6, 9531.1])
            self.obj_data['S3_pyneb_code'] = array([6312, 9069, 9531])
            self.obj_data['S2_abund']      = 0.000015
            self.obj_data['S3_abund']      = 0.000055
             
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
        self.nHydrogen              = len(self.obj_data['H_labels'])
        self.nHydrogen_range        = arange(self.nHydrogen)
        
        self.nHelium                = len(self.obj_data['He_labels'])
        self.nHelium_range          = arange(self.nHelium)

        self.obj_data['nS2']        = len(self.obj_data['S2_labels'])
        self.obj_data['nS2_range']  = arange(self.obj_data['nS2'])

        self.obj_data['nS3']        = len(self.obj_data['S3_labels'])
        self.obj_data['nS3_range']  = arange(self.obj_data['nS3'])
       
        self.nTotal                 = self.nHydrogen + self.nHelium
       
        #Load the reddening parameters
        H_xX        = self.reddening_Xx(self.obj_data['H_wave'], red_curve, Rv)
        He_xX       = self.reddening_Xx(self.obj_data['He_wave'], red_curve, Rv)
        Hbeta_xX    = self.reddening_Xx(array([self.Hbeta_wave]), red_curve, Rv)[0]
        S2_xX       = self.reddening_Xx(self.obj_data['S2_wave'], red_curve, Rv)
        S3_xX       = self.reddening_Xx(self.obj_data['S3_wave'], red_curve, Rv)
        
        #Dictionary with the vectors
        self.data_dic = {'Emissivity_H_vector'  : zeros(self.nHydrogen),
                         'Emissivity_He_vector' : zeros(self.nHelium),
                         'Emissivity_S2_vector' : zeros(self.obj_data['nS2']),
                         'Emissivity_S3_vector' : zeros(self.obj_data['nS3']),
                         'Kalpha_vector'        : zeros(self.nHydrogen),
                         'ftau_He_vector'       : zeros(self.nHelium),
                         'flambda_H_vector'     : H_xX/Hbeta_xX - 1.0,
                         'flambda_He_vector'    : He_xX/Hbeta_xX - 1.0,    
                         'flambda_S2_vector'    : S2_xX/Hbeta_xX - 1.0,
                         'flambda_S3_vector'    : S3_xX/Hbeta_xX - 1.0}    
                
        #Normalize the fluxes by Hbeta if necessary
        if norm_by_Hbeta: 
            self.obj_data['H_Flux']   = self.obj_data['H_Flux']     / self.Hbeta_Flux
            self.obj_data['He_Flux']  = self.obj_data['He_Flux']    / self.Hbeta_Flux
            self.obj_data['H_error']  = self.obj_data['H_error']    / self.Hbeta_Flux
            self.obj_data['He_error'] = self.obj_data['He_error']   / self.Hbeta_Flux

            self.obj_data['S2_Flux']  = self.obj_data['S2_Flux']    / self.Hbeta_Flux
            self.obj_data['S3_Flux']  = self.obj_data['S3_Flux']    / self.Hbeta_Flux
            self.obj_data['S2_error'] = self.obj_data['S2_error']   / self.Hbeta_Flux
            self.obj_data['S3_error'] = self.obj_data['S3_error']   / self.Hbeta_Flux
       
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
        atomicData.setDataFile('s_iii_coll_HRS12.dat')
        
        #Declare pyneb Hydrogen and Helium atoms to calculate emissivities
        self.H1                     = RecAtom('H', 1)
        self.He1                    = RecAtom('He', 1) 
        self.ionDict                = {}
        self.ionDict['S2']          = Atom('S', 2)
        self.ionDict['S3']          = Atom('S', 3) 
    
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
            self.obj_data['S2_Flux']    = self.metal_theo_flux('S2', self.obj_data['Te_0'], self.obj_data['n_e'], self.obj_data['cHbeta'], self.obj_data['S2_abund'])
            self.obj_data['S3_Flux']    = self.metal_theo_flux('S3', self.obj_data['Te_0'], self.obj_data['n_e'], self.obj_data['cHbeta'], self.obj_data['S3_abund'])
            
            self.obj_data['H_error']    = self.obj_data['H_Flux'] * 0.01
            self.obj_data['He_error']   = self.obj_data['He_Flux'] * 0.02
            self.obj_data['S2_error']   = self.obj_data['S2_Flux'] * 0.02
            self.obj_data['S3_error']   = self.obj_data['S3_Flux'] * 0.02

            print '\nInput Parameters'
            print '-y_plus',   self.obj_data['y_plus']
            print '-n_e',      self.obj_data['n_e']
            print '-Te_0',     self.obj_data['Te_0']
            print '-cHbeta',   self.obj_data['cHbeta']
            print '-tau',      self.obj_data['tau']
            print '-xi',       self.obj_data['xi']
            print '-TOIII',    self.obj_data['TOIII'],'+/-',self.obj_data['TOIII_error']
            print '-nSII',     self.obj_data['nSII'],'+/-',self.obj_data['nSII']
            print '-S2_abund', self.obj_data['S2_abund']
            print '-S3_abund', self.obj_data['S3_abund']
            
        print '\nHydrogen emissivities'
        for i in self.nHydrogen_range:
            print '-{} {}'.format(self.obj_data['H_labels'][i], self.obj_data['H_Flux'][i])

        print '\nHelium emissivities'
        for i in self.nHelium_range:
            print '-{} {}'.format(self.obj_data['He_labels'][i], self.obj_data['He_Flux'][i])

        print '\nS2 emissivities'
        for i in self.obj_data['nS2_range']:
            print '-{} {}'.format(self.obj_data['S2_labels'][i], self.obj_data['S2_Flux'][i])

        print '\nS3 emissivities'
        for i in self.obj_data['nS3_range']:
            print '-{} {}'.format(self.obj_data['S3_labels'][i], self.obj_data['S3_Flux'][i])

        #Fluxes combined to speed the process
        self.H_He_Obs_Flux      = concatenate([self.obj_data['H_Flux'], self.obj_data['He_Flux']])                        
        self.H_He_Obs_Error     = concatenate([self.obj_data['H_error'], self.obj_data['He_error']])  
        self.S2_S3_Obs_Flux     = concatenate([self.obj_data['S2_Flux'], self.obj_data['S3_Flux']])                        
        self.S2_S3_Obs_Error    = concatenate([self.obj_data['S2_error'], self.obj_data['S3_error']]) 
       
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
        f_module    = power(10, -1 *  self.data_dic['flambda_He_vector'] * cHbeta)
             
        #Calculate theoretical helium flux for the line
        He_Flux     = y_plus * emis_module * self.data_dic['ftau_He_vector'] * CR_Module * f_module

        return He_Flux
    
    def metal_theo_flux(self, ion, Te, ne, cHbeta, ionic_abund):
        
        #Ions emissivities by Hbeta emissivity
        emis_module = self.metal_emis(ion, Te, ne) / self.Hbeta_emis
        
        
        
        #Reddening factor
        f_module    = power(10, -1 * self.data_dic['flambda_{}_vector'.format(ion)] * cHbeta)



        #Metals flux
        metal_flux  = ionic_abund * emis_module * f_module
        
        return metal_flux
    
    def metal_emis(self, ion, Te, ne):
        
        ion_emis    = self.ionDict[ion]
        pyneb_code  = '{}_pyneb_code'.format(ion)
        vector_emis =  self.data_dic['Emissivity_{}_vector'.format(ion)]
        
        for i in self.obj_data['n' + ion + '_range']:    
            vector_emis[i] = ion_emis.getEmissivity(Te, ne, wave = self.obj_data[pyneb_code][i])
            
        return vector_emis
                      
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
        
        y_plus      =   pymc2.Uniform(           'y_plus',    0.050,                       0.15)
        S2_abund    =   pymc2.Uniform(           'S2_abund',  0.000001,                    0.15)
        S3_abund    =   pymc2.Uniform(           'S3_abund',  0.000001,                    0.15)        
        ne          =   pymc2.TruncatedNormal(   'ne',      self.obj_data['nSII'],         self.obj_data['nSII_error']**-2,    a = 0.0 ,   b = 1000.0)
        Te          =   pymc2.Normal(            'Te',      self.obj_data['TOIII'],        self.obj_data['TOIII_error']**-2)
        tau         =   pymc2.TruncatedNormal(   'tau',     0.75,                          0.5**-2,    a = 0.0,    b = 7.0)
        cHbeta      =   pymc2.TruncatedNormal(   'cHbeta',  0.15,                          0.05**-2,   a = 0.0,    b = 3.0)
        xi          =   pymc2.TruncatedNormal(   'xi',      1,                             200**-2,    a = 0.0,    b = 1000.0)

        @pymc2.deterministic #Calculate Hydrogen theoretical flux
        def det_HFlux_theo(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta):
            return self.hydrogen_theo_flux(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta)
           
        @pymc2.deterministic #Calculate Helium theoretical flux 
        def det_He_Flux_theo_nof(Te=Te, ne=ne, tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus):
            return self.helium_theo_flux(Te=Te, ne=ne,  tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus)

        @pymc2.deterministic #Calculate Helium theoretical flux 
        def det_S2_Flux_theo_nof(ion='S2', Te=Te, ne=ne, tau=tau, xi=xi, cHbeta=cHbeta, abund = S2_abund):
            return self.metal_theo_flux(ion=ion, Te=Te, ne=ne, cHbeta=cHbeta, ionic_abund=abund)       

        @pymc2.deterministic #Calculate Helium theoretical flux 
        def det_S3_Flux_theo_nof(ion='S3', Te=Te, ne=ne, tau=tau, xi=xi, cHbeta=cHbeta, abund = S3_abund):
            return self.metal_theo_flux(ion=ion, Te=Te, ne=ne, cHbeta=cHbeta, ionic_abund=abund)   
       
        @pymc2.deterministic #Combine theoretical fluxes into a single array
        def H_He_TheoFlux(HFlux_theo=det_HFlux_theo, HeFlux_theo=det_He_Flux_theo_nof):
            self.H_He_Theo_Flux = empty(self.nTotal)
            self.H_He_Theo_Flux[:self.nHydrogen] = HFlux_theo[:]
            self.H_He_Theo_Flux[self.nHydrogen:] = HeFlux_theo[:]
            return self.H_He_Theo_Flux

        @pymc2.deterministic #Combine theoretical fluxes into a single array
        def S2_S3_TheoFlux(S2Flux_theo=det_S2_Flux_theo_nof, S3Flux_theo=det_S3_Flux_theo_nof):
            S2_S3_Theo_Flux = empty(5)
            S2_S3_Theo_Flux[:2] = S2Flux_theo[:]
            S2_S3_Theo_Flux[2:] = S3Flux_theo[:]
            return S2_S3_Theo_Flux
                   
        @pymc2.stochastic(observed=True) #Likelihood
        def Likelihood_model_Recomb(value = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2

        @pymc2.stochastic(observed=True) #Likelihood
        def Likelihood_model_metals(value = self.S2_S3_Obs_Flux, S2_S3_TheoFlux = S2_S3_TheoFlux, sigmaLines = self.S2_S3_Obs_Error):
            chi_F = sum(square(S2_S3_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2

        @pymc2.deterministic() #Deterministic method to track the evolution of the chi:
        def ChiSq_Recomb(H_He_ObsFlux = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
            return - chi_F / 2

        @pymc2.deterministic() #Deterministic method to track the evolution of the chi:
        def ChiSq_Metals(S2_S3_Obs_Flux = self.S2_S3_Obs_Flux, S2_S3_TheoFlux = S2_S3_TheoFlux, sigmaLines = self.S2_S3_Obs_Error):
            chi_F = sum(square(S2_S3_TheoFlux - S2_S3_Obs_Flux) / square(sigmaLines))
            return - chi_F / 2
    
        return locals()
   
class Run_MCMC(Inference_AbundanceModel):

    def __init__(self):
        
        #Import parent classes
        Inference_AbundanceModel.__init__(self)
        self.pymc_stats_keys = ['mean','95% HPD interval','standard deviation','mc error','quantiles','n']
  
    def select_inference_model(self, model_code):
        
        #Declare inference model
        if '_abs' not in model_code:
            self.inf_dict = self.helium_abundance_model()
        else:
            self.inf_dict = self.helium_abundance_model_abs()
                    
    def run_pymc2(self, db_address, iterations = 10000, burn_phase = 1500, thining=2, variables_list = None):
        
        #Ideal cond1: iter=30000, burn=5000, thin=10
        #Ideal cond2: iter=10000, burn=1500, thin=2
        
        #Prefit:
        print '\n--Starting Fmin_powell prefit'
        self.MAP_Model  = pymc2.MAP(self.inf_dict)
        self.MAP_Model.fit(method = 'fmin_powell')
        
        #Print prefit data
        self.display_run_data(self.MAP_Model, variables_list)
                
        #Launch sample
        self.pymc2_M    = pymc2.MCMC(self.MAP_Model.variables, db = 'pickle', dbname =  db_address)
        self.pymc2_M.sample(iter=iterations, burn=burn_phase, thin=thining)
        
        #Save the output csv mean data
        if variables_list != None:
            print '--Saving results in csv'
            self.csv_address = db_address + '_Parameters'
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

    def load_pymc_database(self, db_address):
                
        #Load the pymc output textfile database
        pymc_database = pymc2.database.pickle.load(db_address)
        
        #Create a dictionaries with the traces and statistics
        traces_dic = {}
        stats_dic = OrderedDict()
        
        #This variable contains all the traces from the MCMC (stochastic and deterministic)
        traces_list = pymc_database.trace_names[0] 
    
        #Get statistics from the run    
        for trace in traces_list:
            stats_dic[trace] = OrderedDict()
            traces_dic[trace] = pymc_database.trace(trace)
            
            for stat in self.pymc_stats_keys: 
                stats_dic[trace][stat] = pymc_database.trace(trace).stats()[stat] 
    
            trace_array = pymc_database.trace(trace)[:] 
            stats_dic[trace]['16th_p'] = percentile(trace_array, 16)
            stats_dic[trace]['84th_p'] = percentile(trace_array, 84)    
    
        #Generate a MCMC object to recover all the data from the run
        dbMCMC = pymc2.MCMC(traces_dic, pymc_database)
        
        return dbMCMC, stats_dic                                                                                                  

