from collections                    import OrderedDict

from numpy                          import array, arange, loadtxt, genfromtxt, sum, exp, power, zeros, ones, dot, isnan, square, concatenate, where
from numpy.random                   import normal
import pymc
from uncertainties                  import ufloat

from CodeTools.File_Managing_Tools  import Txt_Files_Manager
import pyneb                        as pn


class HeAbundance_InferenceMethods(Txt_Files_Manager):
    
    def __init__(self):

        self.Hydrogen_CollCoeff_TableAddress    = '/home/vital/git/Dazer/Dazer/dazer/lib/Astro_Libraries/Neutral_Hydrogen_Collisional_Correction_coef.txt'
        self.Helium_CollCoeff_TableAddress      = '/home/vital/git/Dazer/Dazer/dazer/lib/Astro_Libraries/Neutral_Helium_Collisional_Correction_coef.txt'
        self.Helium_OpticalDepth_TableAddress   = '/home/vital/git/Dazer/Dazer/dazer/lib/Astro_Libraries/Helium_OpticalDepthFunction_Coefficients.txt'
        
        #Import Dazer_Files Class to import data from lines_logs
        Txt_Files_Manager.__init__(self)
        
        #Declare Hydrogen and Helium lines for the analysis
        self.posHydrogen_Lines      = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
        self.Hydrogen_Wavelengths   = [4101.742,        4340.471,       4862.683,        6562.819]
        
        self.posHelium_Lines        = ['He1_3889A',  'He1_4026A',    'He1_4387A',    'He1_4471A',    'He1_4686A',    'He1_4714A',    'He1_4922A',   'He1_5876A',    'He1_6678A',   'He1_7065A',    'He1_7281A',      'He1_10830A']     
        self.Helium_Wavelengths     = [3889.0,         4026.0,         4387.0,         4471.0,         4686.0,         4714.0,         4922.0,         5876.0,         6678.0,         7065.0,         7281.0,         10830.0]
        
        self.Cand_Hydrogen_Lines    = []
        self.Cand_Helium_Lines      = []

        self.nHydrogen              = None
        self.nHelium                = None
    
        #Define indexes and labels to speed up the code
        self.H13889A_label          = 'H1_3889A'          
        self.HBeta_label            = 'H1_4861A'
        self.He3889_label           = 'He1_3889A'
        self.He3889_Check           = None
        
        #Set up the right emissivities
        pn.atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
                
        #Declare pyneb Hydrogen and Helium atoms to calculate emissivities
        self.H1                     = pn.RecAtom('H', 1)
        self.He1                    = pn.RecAtom('He', 1) 
    
        print 'Helium emissivities:\t', self.He1.printSources()

        #Import collisional coefficients table
        self.Coef_Kalpha_dict       = self.Import_Coll_Coeff_Table(self.posHydrogen_Lines, None)

        #Import Optical depth function
        self.Coef_ftau_dict         = self.Import_OpticalDepth_Coeff_Table(self.posHelium_Lines) 
               
        #Declare dictionaries to store the data
        self.Flux_dict              = OrderedDict()
        self.Error_dict             = OrderedDict()
        self.Wave_dict              = OrderedDict()
        self.PynebCode_dict         = OrderedDict()
        self.EqW_dict               = OrderedDict()
        self.EqWerror_dict          = OrderedDict()
        self.hlambda_dict           = OrderedDict()
        self.flambda_dict           = self.get_flambda_dict(self.posHydrogen_Lines + self.posHelium_Lines, self.Hydrogen_Wavelengths + self.Helium_Wavelengths) #flambda does not form part of the inference predictions
                            
    def Import_TableData(self, Address, Columns):
        
        Imported_Array  = genfromtxt(Address, dtype=float, usecols = Columns, skip_header = 2).T
        Datarray = Imported_Array[:,~isnan(Imported_Array).all(0)]
                
        return Datarray
     
    def Import_Coll_Coeff_Table(self, HydrogenLines, HeliumLines):
        
        Data_dict = OrderedDict()
                
        for i in range(len(HydrogenLines)):
            
            Line            = HydrogenLines[i]
            Data_Columns    = [0 + 3*i, 1 + 3*i, 2 + 3*i]
            Data_dict[Line] = self.Import_TableData(self.Hydrogen_CollCoeff_TableAddress, Data_Columns)
                    
        return Data_dict
     
    def Import_OpticalDepth_Coeff_Table(self, HeliumLines):
        
        Data_dict = OrderedDict()
        
        for i in range(len(HeliumLines)):
            
            Line = HeliumLines[i]
            
            if Line != self.H13889A_label:
                Data_dict[Line] = loadtxt(self.Helium_OpticalDepth_TableAddress, dtype = float, skiprows = 2, usecols = (i,))
                
            else:
                Data_dict[self.He3889_label] = loadtxt(self.Helium_OpticalDepth_TableAddress, dtype = float, skiprows = 2, usecols = (i,))
            
        return Data_dict
    
    def Import_Synthetic_Fluxes(self, Model = 'Model2'):

        self.SyntheticData_dict = OrderedDict()

        if Model == 'Model1':
            self.SyntheticData_dict['H_Flux']       = array([ 0.247269153412,   0.455769957,     1.0,            2.96628259685])
            self.SyntheticData_dict['H_wave']       = array([4101.742,      4340.471,      4862.683,       6562.819])
            self.SyntheticData_dict['H_Ions']       = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
            self.SyntheticData_dict['H_labels']     = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
            self.SyntheticData_dict['H_hlambda']    = [1.0,             1.0,            1.0,             1.0]
            self.SyntheticData_dict['H_error']      = self.SyntheticData_dict['H_Flux'] * 0.01
            
            self.SyntheticData_dict['He_Flux']      = array([ 0.0897022788179,    0.0129308815264,     0.0320439130823, 0.0992100632188, 0.0259080438348, 0.0237335057344])   
            self.SyntheticData_dict['He_wave']      = array([ 3889.0,       4026.0,         4471.0,        5876.0,     6678.0,     7065.0])
            self.SyntheticData_dict['He_labels']    = ['He1_3889A',   'He1_4026A',   'He1_4471A',   'He1_5876A','He1_6678A', 'He1_7065A']
            self.SyntheticData_dict['He_hlambda']   = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0] 
            self.SyntheticData_dict['He_error']     = self.SyntheticData_dict['He_Flux'] * 0.02
                        
            self.SyntheticData_dict['TOIII']        = 19000
            self.SyntheticData_dict['y_plus']       = 0.08
            self.SyntheticData_dict['n_e']          = 100.0
            self.SyntheticData_dict['tau']          = 0.2
            self.SyntheticData_dict['Te_0']         = 18000.0
            self.SyntheticData_dict['cHbeta']       = 0.1
            self.SyntheticData_dict['xi']           = 1.0

        if Model == 'Model2':
            self.SyntheticData_dict['H_Flux']       = array([0.24666367, 0.45494969,    1,              2.97990939])
            self.SyntheticData_dict['H_wave']       = array([4101.742,   4340.471,      4862.683,       6562.819])
            self.SyntheticData_dict['H_Ions']       = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
            self.SyntheticData_dict['H_labels']     = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
            self.SyntheticData_dict['H_hlambda']    = [1.0,             1.0,            1.0,             1.0]
            self.SyntheticData_dict['H_error']      = self.SyntheticData_dict['H_Flux'] * 0.01

            self.SyntheticData_dict['He_Flux']      = array([0.10121858,    0.01695164,     0.03928185,     0.12712208,  0.03548869,    0.0415693])
            self.SyntheticData_dict['He_wave']      = array([ 3889.0,       4026.0,         4471.0,         5876.0,      6678.0,        7065.0])
            self.SyntheticData_dict['He_labels']    = ['He1_3889A',         'He1_4026A',   'He1_4471A',   'He1_5876A',  'He1_6678A',    'He1_7065A']
            self.SyntheticData_dict['He_hlambda']   = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0] 
            self.SyntheticData_dict['He_error']     = self.SyntheticData_dict['He_Flux'] * 0.02

            self.SyntheticData_dict['TOIII']        = 17000
            self.SyntheticData_dict['y_plus']       = 0.085
            self.SyntheticData_dict['n_e']          = 500.0
            self.SyntheticData_dict['tau']          = 1.0
            self.SyntheticData_dict['Te_0']         = 16000.0
            self.SyntheticData_dict['cHbeta']       = 0.1
            self.SyntheticData_dict['xi']           = 1.0
                 
        if Model == 'Model3':
            self.SyntheticData_dict['H_Flux']       = array([0.249665005739, 0.457121205187,    1,              2.9720213093])
            self.SyntheticData_dict['H_wave']       = array([4101.742,   4340.471,      4862.683,       6562.819])
            self.SyntheticData_dict['H_Ions']       = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
            self.SyntheticData_dict['H_labels']     = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
            self.SyntheticData_dict['H_Eqw']        = [50.0,            50.0,           250,             300.0]
            self.SyntheticData_dict['H_EqwErr']     = [1.0,             1.0,            12.5,            14]
            self.SyntheticData_dict['H_hlambda']    = [1.0,             1.0,            1.0,             1.0]
            self.SyntheticData_dict['H_error']      = self.SyntheticData_dict['H_Flux'] * 0.01

            self.SyntheticData_dict['He_Flux']      = array([0.102807348487,    0.0188761305051,     0.0411173796552,     0.128607653966,  0.0373393280888,    0.04339571902,      0.761144585562])
            self.SyntheticData_dict['He_wave']      = array([ 3889.0,       4026.0,         4471.0,         5876.0,      6678.0,        7065.0,         10830.0])
            self.SyntheticData_dict['He_labels']    = ['He1_3889A',         'He1_4026A',   'He1_4471A',   'He1_5876A',  'He1_6678A',    'He1_7065A',    'He1_10830A']
            self.SyntheticData_dict['He_Eqw']       = [10.0,                10.0,           10.0,           10.0,        10.0,          10.0,           10.0]
            self.SyntheticData_dict['He_EqwErr']    = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0,            1.0]
            self.SyntheticData_dict['He_hlambda']   = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0,            1.0] 
            self.SyntheticData_dict['He_error']     = self.SyntheticData_dict['He_Flux'] * 0.02

            self.SyntheticData_dict['TOIII']        = 17000
            self.SyntheticData_dict['y_plus']       = 0.085
            self.SyntheticData_dict['n_e']          = 500.0
            self.SyntheticData_dict['tau']          = 1.0
            self.SyntheticData_dict['Te_0']         = 16000.0
            self.SyntheticData_dict['cHbeta']       = 0.1
            self.SyntheticData_dict['xi']           = 1.0  
  
        if Model == 'Model4':
            self.SyntheticData_dict['H_Flux']       = array([ 0.24685935,   0.45526236,     1.0,            2.96988502])
            self.SyntheticData_dict['H_wave']       = array([4101.742,      4340.471,      4862.683,       6562.819])
            self.SyntheticData_dict['H_Ions']       = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
            self.SyntheticData_dict['H_labels']     = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
            self.SyntheticData_dict['H_Eqw']        = [50.0,            50.0,           250,             300.0]
            self.SyntheticData_dict['H_EqwErr']     = [1.0,             1.0,            12.5,            14]
            self.SyntheticData_dict['H_hlambda']    = [1.0,             1.0,            1.0,             1.0]
            self.SyntheticData_dict['H_error']      = self.SyntheticData_dict['H_Flux'] * 0.01
            
            self.SyntheticData_dict['He_Flux']      = array([ 0.08951096,    0.01290066,     0.03201484, 0.09931435,   0.02594518, 0.02377085,      0.271656766229])   
            self.SyntheticData_dict['He_wave']      = array([ 3889.0,       4026.0,         4471.0,         5876.0,      6678.0,        7065.0,         10830.0])
            self.SyntheticData_dict['He_labels']    = ['He1_3889A',         'He1_4026A',   'He1_4471A',   'He1_5876A',  'He1_6678A',    'He1_7065A',    'He1_10830A']
            self.SyntheticData_dict['He_Eqw']       = [10.0,                10.0,           10.0,           10.0,        10.0,          10.0,   10.0]
            self.SyntheticData_dict['He_EqwErr']    = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0,    10.0]
            self.SyntheticData_dict['He_hlambda']   = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0,    10.0] 
            self.SyntheticData_dict['He_error']     = self.SyntheticData_dict['He_Flux'] * 0.02
                        
            self.SyntheticData_dict['TOIII']        = 19000
            self.SyntheticData_dict['y_plus']       = 0.08
            self.SyntheticData_dict['n_e']          = 100.0
            self.SyntheticData_dict['tau']          = 0.2
            self.SyntheticData_dict['Te_0']         = 18000.0
            self.SyntheticData_dict['cHbeta']       = 0.1
            self.SyntheticData_dict['xi']           = 1.0
                
        if Model == 'Model3_abs':
            self.SyntheticData_dict['H_Flux']       = array([0.246663665762, 0.454949690008,    1,              2.97990939454])
            self.SyntheticData_dict['H_wave']       = array([4101.742,   4340.471,      4862.683,       6562.819])
            self.SyntheticData_dict['H_Ions']       = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
            self.SyntheticData_dict['H_labels']     = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
            self.SyntheticData_dict['H_Eqw']        = [50.0,            50.0,           250,             300.0]
            self.SyntheticData_dict['H_EqwErr']     = [1.0,             1.0,            12.5,            14]
            self.SyntheticData_dict['H_hlambda']    = [1.0,             1.0,            1.0,             1.0]
            self.SyntheticData_dict['H_error']      = self.SyntheticData_dict['H_Flux'] * 0.01

            self.SyntheticData_dict['He_Flux']      = array([0.101218577881,    0.0169516350271,     0.0392818491739,     0.127122084582,  0.0354886854011,    0.0415693018961,      0.762189163904])
            self.SyntheticData_dict['He_wave']      = array([ 3889.0,       4026.0,         4471.0,         5876.0,      6678.0,        7065.0,         10830.0])
            self.SyntheticData_dict['He_labels']    = ['He1_3889A',         'He1_4026A',   'He1_4471A',   'He1_5876A',  'He1_6678A',    'He1_7065A',    'He1_10830A']
            self.SyntheticData_dict['He_Eqw']       = [10.0,                10.0,           10.0,           10.0,        10.0,          10.0,           10.0]
            self.SyntheticData_dict['He_EqwErr']    = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0,            1.0]
            self.SyntheticData_dict['He_hlambda']   = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0,            1.0] 
            self.SyntheticData_dict['He_error']     = self.SyntheticData_dict['He_Flux'] * 0.02

            self.SyntheticData_dict['TOIII']        = 17000
            self.SyntheticData_dict['y_plus']       = 0.085
            self.SyntheticData_dict['n_e']          = 500.0
            self.SyntheticData_dict['tau']          = 1.0
            self.SyntheticData_dict['Te_0']         = 16000.0
            self.SyntheticData_dict['cHbeta']       = 0.1
            self.SyntheticData_dict['xi']           = 1.0          
            self.SyntheticData_dict['abs_H']        = 1.0          
            self.SyntheticData_dict['abs_He']       = 0.5          
            self.SyntheticData_dict['Hbeta_eqw']    = 250          
         
        print 'Physical parameters:'
        self.ObjectData_dict = OrderedDict()
        self.ObjectData_dict['nSII'] = 500
        self.ObjectData_dict['nSII_error'] = 200
        self.ObjectData_dict['TOIII'] = self.SyntheticData_dict['TOIII']
        self.ObjectData_dict['TOIII_error'] = self.SyntheticData_dict['TOIII'] * 0.2
                               
        return
    
    def Calculate_Synthetic_Fluxes(self, model):
        
        self.Import_Synthetic_Fluxes(Model = model)
        
        self.Check_EmissionLinesObserved(self.SyntheticData_dict)
        
        self.Preload_Data(Hbeta_Normalized=True, Deblend_Check=False)
    
        if '_abs' not in model:
            
            H_Flux          = self.H_Flux_theo(self.SyntheticData_dict['Te_0'], self.SyntheticData_dict['n_e'], self.SyntheticData_dict['xi'], self.SyntheticData_dict['cHbeta'])
            He_Flux         = self.He_Flux_theo_nof(self.SyntheticData_dict['Te_0'], self.SyntheticData_dict['n_e'], self.SyntheticData_dict['tau'], self.SyntheticData_dict['xi'], self.SyntheticData_dict['cHbeta'], self.SyntheticData_dict['y_plus'])
    
            print '\nInput Parameters'
            print 'y_plus' ,self.SyntheticData_dict['y_plus']
            print 'n_e', self.SyntheticData_dict['n_e']
            print 'Te_0', self.SyntheticData_dict['Te_0']
            print 'cHbeta', self.SyntheticData_dict['cHbeta']
            print 'tau', self.SyntheticData_dict['tau']
            print 'xi', self.SyntheticData_dict['xi']
    
        else:
    
            H_Flux          = self.H_Flux_theo_abs(self.SyntheticData_dict['Te_0'], self.SyntheticData_dict['n_e'], self.SyntheticData_dict['xi'], self.SyntheticData_dict['cHbeta'], self.SyntheticData_dict['abs_H'], self.SyntheticData_dict['Hbeta_eqw'])
            He_Flux         = self.He_Flux_theo_abs(self.SyntheticData_dict['Te_0'], self.SyntheticData_dict['n_e'], self.SyntheticData_dict['tau'], self.SyntheticData_dict['xi'], self.SyntheticData_dict['cHbeta'], self.SyntheticData_dict['y_plus'], self.SyntheticData_dict['abs_H'], self.SyntheticData_dict['abs_He'], self.SyntheticData_dict['Hbeta_eqw'])
    
            print '\nInput Parameters'
            print 'y_plus' ,self.SyntheticData_dict['y_plus']
            print 'n_e', self.SyntheticData_dict['n_e']
            print 'Te_0', self.SyntheticData_dict['Te_0']
            print 'cHbeta', self.SyntheticData_dict['cHbeta']
            print 'tau', self.SyntheticData_dict['tau']
            print 'xi', self.SyntheticData_dict['xi']
            print 'abs_H', self.SyntheticData_dict['abs_H']     
            print 'abs_He',  self.SyntheticData_dict['abs_He']            
            print 'Hbeta_eqw', self.SyntheticData_dict['Hbeta_eqw']          
        
        print '\nHydrogen emissivities'
        for i in range(len(self.Emissivity_H_vector)):
            print self.Cand_Hydrogen_Lines[i], self.Emissivity_H_vector[i]/self.Data_Dict['Hbeta_emis'] 

        print '\nHelium emissivities'
        for i in range(len(self.Emissivity_He_vector)):
            print self.Cand_Helium_Lines[i], self.Emissivity_He_vector[i]/self.Data_Dict['Hbeta_emis'] 
        
        print '\nHydrogen fluxes'
        for flux in H_Flux:
            print flux
        
        print '\nHelium fluxes'
        for flux in He_Flux:
            print flux

        return
    
    def Import_Object_Data(self, Filefolder, CodeName, Lineslog_Extension, Deblend_Check = True):
        
        self.ObjectData_dict = OrderedDict()
            
        String_Columns      = [self.Labels_ColumnHeader, self.Ion_ColumnHeader]
        Float_Columns       = [self.Wavelength_ColumnHeader, self.GaussianFlux_ColumnHeader, self.GaussianError_ColumnHeader, self.Line_Continuum_ColumnHeader, self.EqW_ColumnHeader, self.EqW_error_ColumnHeader]
        
        Lineslog_Address    = Filefolder + CodeName + Lineslog_Extension
                
        Labels, Ions                                                = self.get_ColumnData(String_Columns, Lineslog_Address, HeaderSize = 2, StringIndexes = True, datatype = str, unpack_check = True) 
        Wavelengths, Fluxes, Fluxes_error, h_lambda, Eqw, Eqw_error = self.get_ColumnData(Float_Columns, Lineslog_Address, HeaderSize = 2, StringIndexes = True, unpack_check = True)

        self.ObjectData_dict['H_labels']        = list(Labels)
        self.ObjectData_dict['H_Ions']          = list(Ions)
        self.ObjectData_dict['H_Flux']          = Fluxes
        self.ObjectData_dict['H_error']         = Fluxes_error
        self.ObjectData_dict['H_wave']          = Wavelengths
        self.ObjectData_dict['H_Eqw']           = Eqw
        self.ObjectData_dict['H_EqwErr']        = Eqw_error
        self.ObjectData_dict['H_hlambda']       = h_lambda
        
        self.ObjectData_dict['He_labels']       = list(Labels)
        self.ObjectData_dict['He_Ions']         = list(Ions)
        self.ObjectData_dict['He_Flux']         = Fluxes
        self.ObjectData_dict['He_error']        = Fluxes_error
        self.ObjectData_dict['He_wave']         = Wavelengths
        self.ObjectData_dict['He_Eqw']          = Eqw
        self.ObjectData_dict['He_EqwErr']       = Eqw_error
        self.ObjectData_dict['He_hlambda']      = h_lambda
        
        #TOIII calculated from TSIII, if not we use TOIII, if not we use 17000
        TOIII_from_TSIII = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'TOIII_approx_TSIII', Assumption='float')
        TOIII            = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'TOIII', Assumption='float')
        
        if TOIII_from_TSIII != None:
            self.ObjectData_dict['TOIII']       = TOIII_from_TSIII.nominal_value
            self.ObjectData_dict['TOIII_error'] = TOIII_from_TSIII.std_dev
        elif TOIII != None:
            self.ObjectData_dict['TOIII']       = TOIII.nominal_value
            self.ObjectData_dict['TOIII_error'] = TOIII.std_dev
        else: 
            self.ObjectData_dict['TOIII']       = 12000.0
            self.ObjectData_dict['TOIII_error'] = 2000.0
            
        #nSII density to be used in the priors definition. If it could not be calculated the output value is 100
        ne_obj = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'nSII', Assumption = 'float')
        if (ne_obj < 0) or (ne_obj == None):
            self.ObjectData_dict['nSII']        = 50.0
            self.ObjectData_dict['nSII_error']  = 50.0
        else:
            self.ObjectData_dict['nSII']        = ne_obj.nominal_value
            self.ObjectData_dict['nSII_error']  = ne_obj.std_dev          
            
        self.ObjectData_dict['cHbeta_obs']  = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'cHBeta_stellar', Assumption = 'cHbeta_min')
        return
    
    def Check_EmissionLinesObserved(self, Data_Dict):
        
        #Empty the lists of candidate lines for Hydrogen and Helium
        self.Cand_Hydrogen_Lines = []
        self.Cand_Helium_Lines = []
        self.He3889_Check = False
        
        #Loop through the hydrogen lines to store the properties of those we can use to perform the analysis
        Obs_HydrogenLabels                                  = Data_Dict['H_labels']
        for i in range(len(self.posHydrogen_Lines)):
            Line_Label                                      = self.posHydrogen_Lines[i]
            if Line_Label in Obs_HydrogenLabels:
                Line_Index                                  = Obs_HydrogenLabels.index(Line_Label)
                self.Flux_dict[Line_Label]                  = Data_Dict['H_Flux'][Line_Index]   
                self.Error_dict[Line_Label]                 = Data_Dict['H_error'][Line_Index]
                self.Wave_dict[Line_Label]                  = Data_Dict['H_wave'][Line_Index]
                self.hlambda_dict[Line_Label]               = Data_Dict['H_hlambda'][Line_Index]
                self.PynebCode_dict[Line_Label]             = Data_Dict['H_Ions'][Line_Index][Data_Dict['H_Ions'][Line_Index].find('_')+1:len(Data_Dict['H_Ions'][Line_Index])]
                self.EqW_dict[Line_Label]                   = Data_Dict['H_Eqw'][Line_Index]
                self.EqWerror_dict[Line_Label]              = Data_Dict['H_EqwErr'][Line_Index]               
                self.Cand_Hydrogen_Lines.append(Line_Label)
                                          
        #Loop through the helium lines to store the properties of those we can use to perform the analysis
        Obs_HeliumLabels                                    = Data_Dict['He_labels']
        for i in range(len(self.posHelium_Lines)):
            Line_Label                                      = self.posHelium_Lines[i]
            if Line_Label in Obs_HeliumLabels:
                Line_Index                                  = Obs_HeliumLabels.index(Line_Label)
                self.Flux_dict[Line_Label]                  = Data_Dict['He_Flux'][Line_Index]   
                self.Error_dict[Line_Label]                 = Data_Dict['He_error'][Line_Index]
                self.Wave_dict[Line_Label]                  = Data_Dict['He_wave'][Line_Index]
                self.hlambda_dict[Line_Label]               = Data_Dict['He_hlambda'][Line_Index]
                self.PynebCode_dict[Line_Label]             = Line_Label[Line_Label.find('_')+1:-1]+'.0'
                self.EqW_dict[Line_Label]                   = Data_Dict['He_Eqw'][Line_Index]
                self.EqWerror_dict[Line_Label]              = Data_Dict['He_EqwErr'][Line_Index]                        
                self.Cand_Helium_Lines.append(Line_Label)
        
        #Print the helium lines observed for the object
        print '\nObserved hydrogen lines'
        print self.Cand_Hydrogen_Lines
        print 'Observed Helium lines'
        print self.Cand_Helium_Lines
        
        #Since the Hbeta line is not used for the fitting we remove it from the candidates lines list. However, its properties are still stored in the dictionaries
        self.Cand_Hydrogen_Lines.remove(self.HBeta_label)
        
        #Since the He_3889A line is blended convert the entry from the H1_3889A measurement
        if self.H13889A_label in self.Cand_Helium_Lines:
            
            self.Error_dict[self.He3889_label]          = self.Error_dict.pop(self.H13889A_label)
            self.Wave_dict[self.He3889_label]           = self.Wave_dict.pop(self.H13889A_label)
            self.flambda_dict[self.He3889_label]        = self.flambda_dict.pop(self.H13889A_label)
            self.hlambda_dict[self.He3889_label]        = self.hlambda_dict.pop(self.H13889A_label)
            self.PynebCode_dict[self.He3889_label]      = self.PynebCode_dict.pop(self.H13889A_label)
            self.EqW_dict[self.He3889_label]            = self.EqW_dict.pop(self.H13889A_label)
            self.EqWerror_dict[self.He3889_label]       = self.EqWerror_dict.pop(self.H13889A_label)

            He3889A_Deblended_Flux                      = self.Deblend_He3889A_line(Te = self.ObjectData_dict['TOIII'], ne = self.ObjectData_dict['nSII'], cHbeta = self.ObjectData_dict['cHbeta_obs'])
            self.Flux_dict[self.He3889_label]           = He3889A_Deblended_Flux #This does not delete the H1_3889 previous entry
            
            #We replace H13889A_on the list of Candidate helium lines
            self.Cand_Helium_Lines[self.Cand_Helium_Lines.index(self.H13889A_label)] = self.He3889_label
            self.He3889_Check = True
            
            print 'Updated observed Helium'
            print self.Cand_Helium_Lines
                                
        #Define number of lines variables to speed up the calculation
        self.nHydrogen                                      = len(self.Cand_Hydrogen_Lines)
        self.nHydrogen_range                                = arange(self.nHydrogen)
        
        self.nHelium                                        = len(self.Cand_Helium_Lines)
        self.nHelium_range                                  = arange(self.nHelium)
        
        #Convert list to array to make it more efficient
        self.Cand_Hydrogen_Lines                            = array(self.Cand_Hydrogen_Lines)
        self.Cand_Helium_Lines                              = array(self.Cand_Helium_Lines)
                
    def Preload_Data(self, Hbeta_Normalized = True, Deblend_Check = True):

        #Define the vectors to store the physical parameters
        self.Flux_H_vector                          = zeros(self.nHydrogen)
        self.Error_H_vector                         = zeros(self.nHydrogen)
        self.Emissivity_H_vector                    = zeros(self.nHydrogen)
        self.Kalpha_vector                          = zeros(self.nHydrogen)
        self.flambda_H_vector                       = zeros(self.nHydrogen)
        self.EqW_H_vector                           = zeros(self.nHydrogen)
        self.EqWerror_H_vector                      = zeros(self.nHydrogen)
        self.hlambda_H_vector                       = ones(self.nHydrogen)

        self.Flux_He_vector                         = zeros(self.nHelium)
        self.Error_He_vector                        = zeros(self.nHelium)
        self.Emissivity_He_vector                   = zeros(self.nHelium)
        self.ftau_He_vector                         = zeros(self.nHelium)
        self.flambda_He_vector                      = zeros(self.nHelium)
        self.EqW_He_vector                          = zeros(self.nHelium)
        self.EqWerror_He_vector                     = zeros(self.nHelium)
        self.hlambda_He_vector                      = ones(self.nHelium)
      
        #Define HBeta values in advance
        self.Hbeta_Flux                             = self.Flux_dict[self.HBeta_label]
        self.Hbeta_error                            = self.Error_dict[self.HBeta_label]
        self.hbeta_hlambda                          = self.hlambda_dict[self.HBeta_label]
        self.hbeta_eqw                              = self.EqW_dict[self.HBeta_label]
        
        #Load physical parameters vectors
        for i in range(self.nHydrogen):
            Line_Label                              = self.Cand_Hydrogen_Lines[i]        
            self.Flux_H_vector[i]                   = self.Flux_dict[Line_Label]
            self.Error_H_vector[i]                  = self.Error_dict[Line_Label]
            self.hlambda_H_vector[i]                = self.hlambda_dict[Line_Label]
            self.flambda_H_vector[i]                = self.flambda_dict[Line_Label]
            self.EqW_H_vector                       = self.EqW_dict[Line_Label]
            self.EqWerror_H_vector                  = self.EqWerror_dict[Line_Label]

        for i in range(self.nHelium):
            Line_Label                              = self.Cand_Helium_Lines[i]
            self.Flux_He_vector[i]                  = self.Flux_dict[Line_Label]
            self.Error_He_vector[i]                 = self.Error_dict[Line_Label]
            self.flambda_He_vector[i]               = self.flambda_dict[Line_Label]
            self.hlambda_He_vector[i]               = self.hlambda_dict[Line_Label]             
            self.EqW_H_vector                       = self.EqW_dict[Line_Label]
            self.EqWerror_H_vector                  = self.EqWerror_dict[Line_Label]

        #Combine Hydrogen and Helium values into a single vector
        if Hbeta_Normalized: 
            self.H_He_Obs_Flux                      = concatenate([self.Flux_H_vector, self.Flux_He_vector])
            self.H_He_Obs_Error                     = concatenate([self.Error_H_vector, self.Error_He_vector])

        else:
            self.H_He_Obs_Flux                      = concatenate([self.Flux_H_vector, self.Flux_He_vector])    / self.Hbeta_Flux
            self.H_He_Obs_Error                     = concatenate([self.Error_H_vector, self.Error_He_vector])  / self.Hbeta_Flux
                        
        #Dictionary with the vectors
        self.Data_Dict = {'Hbeta_Kalpha'    : None, 
                     'Hbeta_emis'           : None,
                     'Emissivity_H_vector'  : None,
                     'Kalpha_vector'        : None,
                     'Emissivity_He_vector' : None,
                     'ftau_He_vector'       : None,
                     'T_4'                  : None} 
                
    def Calculate_H_Parameters_Vectors(self, Te, ne):    
        
        #Calculate in advance the T_4 parameter
        self.Data_Dict['T_4'] = Te / 10000.0
        
        #Calculate the hbeta parameters
        self.Data_Dict['Hbeta_Kalpha']              = self.Kalpha_Ratio_H(T_4 = self.Data_Dict['T_4'], H_label=self.HBeta_label)
        self.Data_Dict['Hbeta_emis']                = self.H1.getEmissivity(Te, ne, label = self.PynebCode_dict[self.HBeta_label])

        #Calculate physical parameters for Hydrogen lines
        for i in self.nHydrogen_range:
            self.Emissivity_H_vector[i]             = self.H1.getEmissivity(Te, ne, label = self.PynebCode_dict[self.Cand_Hydrogen_Lines[i]])
            self.Kalpha_vector[i]                   = self.Kalpha_Ratio_H(T_4 = self.Data_Dict['T_4'], H_label=self.Cand_Hydrogen_Lines[i])
        
        self.Data_Dict['Emissivity_H_vector']       = self.Emissivity_H_vector
        self.Data_Dict['Kalpha_vector']             = self.Kalpha_vector
                
        return 
    
    def Calculate_He_Parameters_Vectors(self, Te, ne, tau):
    
        #Calculate physical parameters for Helium lines
        for i in self.nHelium_range:
            self.Emissivity_He_vector[i]            = self.He1.getEmissivity(Te, ne, label = self.PynebCode_dict[self.Cand_Helium_Lines[i]])
            self.ftau_He_vector[i]                  = self.OpticalDepth_He(tau = tau, T_4 = self.Data_Dict['T_4'], ne = ne, He_label = self.Cand_Helium_Lines[i])    
            
        self.Data_Dict['Emissivity_He_vector']      = self.Emissivity_He_vector
        self.Data_Dict['ftau_He_vector']            = self.ftau_He_vector
        
        return

    def H_Flux_theo(self, Te, ne, xi, cHbeta):
        
        #Hydrogen calculation parameters
        self.Calculate_H_Parameters_Vectors(Te, ne)

        #Calculate the emissivities for each of the lines for the given temperature and density
        Emissivities_module = self.Data_Dict['Emissivity_H_vector'] / self.Data_Dict['Hbeta_emis'] 
                
        #Calculate the Collisional excitation fraction
        CR_Module           = (1.0 + 0.0001* xi * self.Data_Dict['Kalpha_vector']) / (1.0 + 0.0001* xi * self.Data_Dict['Hbeta_Kalpha'] )
                
        #Calculate the reddening component
        f_module            = power(10, -1 * self.flambda_H_vector * cHbeta)
          
        #Calculate theoretical Hydrogen flux for each emission line
        H_Flux              = Emissivities_module * CR_Module * f_module
                     
        return H_Flux

    def He_Flux_theo_nof(self, Te, ne, tau, xi, cHbeta, y_plus):
        
        #Helium calculation parameters
        self.Calculate_He_Parameters_Vectors(Te, ne, tau)
             
        #Calculate the emissivities for each of the lines for the given temperature and density
        Emissivities_module = self.Data_Dict['Emissivity_He_vector'] / self.Data_Dict['Hbeta_emis']
                
        #Calculate the collisional excitation fraction
        CR_Module           = 1 / (1.0 + 0.0001* xi * self.Data_Dict['Hbeta_Kalpha'])

        #Calculate the reddening component
        f_module            = power(10, -1 * self.flambda_He_vector * cHbeta)
             
        #Calculate theoretical Hydrogen flux for each emission line
        He_Flux             = y_plus * Emissivities_module * self.Data_Dict['ftau_He_vector'] * CR_Module * f_module

        return He_Flux

    def H_Flux_theo_abs(self, Te, ne, xi, cHbeta, a_H, Eqw_hbeta):
        
        #Hydrogen calculation parameters
        self.Calculate_H_Parameters_Vectors(Te, ne)

        #Calculate the emissivities for each of the lines for the given temperature and density
        Emissivities_module = self.Data_Dict['Emissivity_H_vector'] / self.Data_Dict['Hbeta_emis'] 
                
        #Calculate the Collisional excitation fraction
        CR_Module           = (1.0 + 0.0001* xi * self.Data_Dict['Kalpha_vector']) / (1.0 + 0.0001* xi * self.Data_Dict['Hbeta_Kalpha'] )
                
        #Calculate the reddening component
        f_module            = power(10, -1 * self.flambda_H_vector * cHbeta)
          
        #Calculate the Hbeta normalization module
        EW_Hbeta_module     = (Eqw_hbeta + a_H) / Eqw_hbeta          
          
        #Calculate stellar absorption module
        a_H_module          = (a_H * self.hlambda_H_vector) / (Eqw_hbeta * self.hbeta_hlambda)          
          
        #Calculate theoretical Hydrogen flux for each emission line
        H_Flux              = Emissivities_module * CR_Module * f_module * EW_Hbeta_module - a_H_module
                     
        return H_Flux

    def He_Flux_theo_abs(self, Te, ne, tau, xi, cHbeta, y_plus, a_H, a_He, Eqw_hbeta):
        
        #Helium calculation parameters
        self.Calculate_He_Parameters_Vectors(Te, ne, tau)
             
        #Calculate the emissivities for each of the lines for the given temperature and density
        Emissivities_module = self.Data_Dict['Emissivity_He_vector'] / self.Data_Dict['Hbeta_emis']
                
        #Calculate the collisional excitation fraction
        CR_Module           = 1 / (1.0 + 0.0001* xi * self.Data_Dict['Hbeta_Kalpha'])

        #Calculate the reddening component
        f_module            = power(10, -1 * self.flambda_He_vector * cHbeta)
             
        #Calculate the Hbeta normalization module
        EW_Hbeta_module     = (Eqw_hbeta + a_H) / Eqw_hbeta
        
        #Calculate stellar absorption module
        a_He_module         = (a_He * self.hlambda_He_vector) / (Eqw_hbeta * self.hbeta_hlambda)             
             
        #Calculate theoretical Hydrogen flux for each emission line
        He_Flux             = y_plus * Emissivities_module * self.Data_Dict['ftau_He_vector'] * CR_Module * f_module  * EW_Hbeta_module - a_He_module

        return He_Flux

    def H_He_Corr_ObsFlux(self, H_He_ObsFlux, Te, ne):
                
        H_He_ObsFlux.setflags(write=True)

        if self.He3889_Check:
            H_He_ObsFlux[self.nHydrogen + 1] = H_He_ObsFlux[self.nHydrogen + 1] + 0 * Te
        else:
            H_He_ObsFlux[self.nHydrogen + 1] = H_He_ObsFlux[self.nHydrogen + 1] + 0 * Te
                        
        return H_He_ObsFlux
                 
    def get_flambda_dict(self, lines_labels, lines_wavelengths, R_v=3.2):
    
        x_true          = 1.0 / (array(lines_wavelengths) / 10000.0)
        y               = x_true - 1.82
    
        y_coeffs        = array([ones(len(y)), y, power(y, 2), power(y, 3), power(y, 4), power(y, 5), power(y, 6), power(y, 7)])
        a_coeffs        = array([1, 0.17699,    -0.50447,   -0.02427,   0.72085,    0.01979,    -0.77530,   0.32999])
        b_coeffs        = array([0, 1.41338,    2.28305,   1.07233,   -5.38434,    -0.62251,    5.30260,   -2.09002])
        
        a_x             = dot(a_coeffs,y_coeffs)
        b_x             = dot(b_coeffs,y_coeffs)
        
        X_x             = a_x + b_x / R_v
        
        y_beta          = (1 / (4862.683 / 10000)) - 1.82
        y_beta_coeffs   = array([1, y_beta, power(y_beta, 2), power(y_beta, 3), power(y_beta, 4), power(y_beta, 5), power(y_beta, 6), power(y_beta, 7)])
        
        X_x_beta        = dot(a_coeffs,y_beta_coeffs) + dot(b_coeffs,y_beta_coeffs) / R_v
    
        f               = X_x / X_x_beta - 1
        
        return dict(zip(lines_labels , f))
         
    def Kalpha_Ratio_H(self, T_4, H_label):
                
        K_alpha_Ratio   = sum(self.Coef_Kalpha_dict[H_label][0] * exp(self.Coef_Kalpha_dict[H_label][1] / T_4) * power(T_4, self.Coef_Kalpha_dict[H_label][2]))
        
        return K_alpha_Ratio

    def OpticalDepth_He(self, tau, T_4, ne, He_label):
                        
        f_tau   = 1 + (tau/2) * (self.Coef_ftau_dict[He_label][0] + (self.Coef_ftau_dict[He_label][1] + self.Coef_ftau_dict[He_label][2]*ne + self.Coef_ftau_dict[He_label][3]*ne*ne) * T_4)

        return f_tau

    def Deblend_He3889A_line(self, Te, ne, cHbeta):
        
        Emis_Hydrogen3889   = 0.104 * power(Te / 10000, 0.046)        
        Flux3389_byHbeta    = self.Flux_dict[self.H13889A_label] / self.Flux_dict[self.HBeta_label] - (Emis_Hydrogen3889) * (power(10 , - self.flambda_dict[self.He3889_label] * cHbeta))
        Flux3389            = Flux3389_byHbeta * self.Flux_dict[self.HBeta_label] 
        
        return Flux3389     
        
class HeAbundance_InferenceStructure(HeAbundance_InferenceMethods):

    def __init__(self):

        HeAbundance_InferenceMethods.__init__(self)

    def Bayesian_HeliumAbundance_Analysis(self):

        print 'Physical parameters:'
        print self.ObjectData_dict['nSII'], self.ObjectData_dict['nSII_error']
        print self.ObjectData_dict['TOIII'], self.ObjectData_dict['TOIII_error']
        
        y_plus      =   pymc.Uniform(           'He_abud',  0.050,                              0.15)
        ne          =   pymc.TruncatedNormal(   'n_e',      self.ObjectData_dict['nSII'],       self.ObjectData_dict['nSII_error']**-2,    a = 0.0 ,   b = 1000.0)
        Te          =   pymc.Normal(            'T_e',      self.ObjectData_dict['TOIII'],      self.ObjectData_dict['TOIII_error']**-2)
        tau         =   pymc.TruncatedNormal(   'Tau',      0.75,                               0.5**-2,    a = 0.0,    b = 7.0)
        cHbeta      =   pymc.TruncatedNormal(   'c_Hbeta',  0.15,                               0.05**-2,   a = 0.0,    b = 3.0)
        xi          =   pymc.TruncatedNormal(   'Xi',       1,                                  200**-2,    a = 0.0,    b = 1000.0)
#       a_H         =   pymc.TruncatedNormal(   'abs_H',    1.0,                                1**-2,      a = 0.0,    b = 14.0)
#       a_He        =   pymc.TruncatedNormal(   'abs_He',   1.0,                                0.5**-2,    a = 0.0,    b = 5.0)
#       Eqw_hbeta   =   pymc.Normal(            'Eqw_hbeta',self.EqW_dict['H1_4861A'],          self.EqWerror_dict['H1_4861A']**-2)

        #Calculate Hydrogen theoretical flux
        @pymc.deterministic
        def det_HFlux_theo(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta):
            return self.H_Flux_theo(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta)
        
        #Calculate Helium theoretical flux    
        @pymc.deterministic
        def det_He_Flux_theo_nof(Te=Te, ne=ne, tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus):
            return self.He_Flux_theo_nof(Te=Te, ne=ne,  tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus)
        
        #Combine theoretical fluxes into a single array
        @pymc.deterministic
        def H_He_TheoFlux(HFlux_theo=det_HFlux_theo, HeFlux_theo=det_He_Flux_theo_nof):
            return concatenate([HFlux_theo, HeFlux_theo])
            
        #Likelihood
        @pymc.stochastic(observed=True)
        def Likelihood_model(value = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2
        
        #Deterministic method to track the evolution of the chi:
        @pymc.deterministic()
        def ChiSq(H_He_ObsFlux = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
            return - chi_F / 2
       
        return locals()

    def Bayesian_HeliumAbundance_Analysis_obsCor(self):

        print 'Physical parameters:'
        print self.ObjectData_dict['nSII'], self.ObjectData_dict['nSII_error']
        print self.ObjectData_dict['TOIII'], self.ObjectData_dict['TOIII_error']
        
        y_plus      =   pymc.Uniform(           'He_abud',  0.050,                              0.15)
        ne          =   pymc.TruncatedNormal(   'n_e',      self.ObjectData_dict['nSII'],       self.ObjectData_dict['nSII_error']**-2,    a = 0.0 ,   b = 1000.0)
        Te          =   pymc.Normal(            'T_e',      self.ObjectData_dict['TOIII'],      self.ObjectData_dict['TOIII_error']**-2)
        tau         =   pymc.TruncatedNormal(   'Tau',      0.75,                               0.5**-2,    a = 0.0,    b = 7.0)
        cHbeta      =   pymc.TruncatedNormal(   'c_Hbeta',  0.15,                               0.05**-2,   a = 0.0,    b = 3.0)
        xi          =   pymc.TruncatedNormal(   'Xi',       1,                                  200**-2,    a = 0.0,    b = 1000.0)
#       a_H         =   pymc.TruncatedNormal(   'abs_H',    1.0,                                1**-2,      a = 0.0,    b = 14.0)
#       a_He        =   pymc.TruncatedNormal(   'abs_He',   1.0,                                0.5**-2,    a = 0.0,    b = 5.0)
#       Eqw_hbeta   =   pymc.Normal(            'Eqw_hbeta',self.EqW_dict['H1_4861A'],          self.EqWerror_dict['H1_4861A']**-2)

        #Calculate Hydrogen theoretical flux
        @pymc.deterministic
        def det_HFlux_theo(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta):
            return self.H_Flux_theo(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta)
        
        #Calculate Helium theoretical flux    
        @pymc.deterministic
        def det_He_Flux_theo_nof(Te=Te, ne=ne, tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus):
            return self.He_Flux_theo_nof(Te=Te, ne=ne,  tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus)
        
        #Combine theoretical fluxes into a single array
        @pymc.deterministic
        def H_He_TheoFlux(HFlux_theo=det_HFlux_theo, HeFlux_theo=det_He_Flux_theo_nof):
            return concatenate([HFlux_theo, HeFlux_theo])
        
        @pymc.deterministic
        def det_H_He_Corr_ObsFlux(H_He_ObsFlux = self.H_He_Obs_Flux, Te=Te, ne=ne):
            return self.H_He_Corr_ObsFlux(H_He_ObsFlux=H_He_ObsFlux, Te=Te, ne=ne)
            
        #Likelihood
        @pymc.stochastic(observed=True)
        def Likelihood_model(value = det_H_He_Corr_ObsFlux.value, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2
        
        #Deterministic method to track the evolution of the chi:
        @pymc.deterministic()
        def ChiSq(H_He_ObsFlux = det_H_He_Corr_ObsFlux.value, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
            return - chi_F / 2
       
        return locals()
   
    def Bayesian_HeliumAbundance_Analysis_abs(self):

        print 'Physical parameters:'
        print self.ObjectData_dict['nSII'], self.ObjectData_dict['nSII_error']
        print self.ObjectData_dict['TOIII'], self.ObjectData_dict['TOIII_error']
        
        y_plus      =   pymc.Uniform(           'He_abud',      0.050,                              0.15)
        ne          =   pymc.TruncatedNormal(   'n_e',          self.ObjectData_dict['nSII'],       self.ObjectData_dict['nSII_error']**-2,    a = 0.0 ,   b = 1000.0)
        Te          =   pymc.Normal(            'T_e',          self.ObjectData_dict['TOIII'],      self.ObjectData_dict['TOIII_error']**-2)
        tau         =   pymc.TruncatedNormal(   'Tau',          0.75,                               0.5**-2,    a = 0.0,    b = 7.0)
        cHbeta      =   pymc.TruncatedNormal(   'c_Hbeta',      0.15,                               0.05**-2,   a = 0.0,    b = 3.0)
        xi          =   pymc.TruncatedNormal(   'Xi',           1,                                  200**-2,    a = 0.0,    b = 1000.0)
        a_H         =   pymc.TruncatedNormal(   'abs_H',        1.0,                                1**-2,      a = 0.0,    b = 14.0)
        a_He        =   pymc.TruncatedNormal(   'abs_He',       1.0,                                0.5**-2,    a = 0.0,    b = 5.0)
        Eqw_hbeta   =   pymc.Normal(            'Eqw_hbeta',    self.EqW_dict['H1_4861A'],          self.EqWerror_dict['H1_4861A']**-2)

        #Calculate Hydrogen theoretical flux
        @pymc.deterministic
        def det_HFlux_theo_abs(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta, a_H=a_H, Eqw_hbeta=Eqw_hbeta):
            return self.H_Flux_theo_abs(Te=Te, ne=ne, xi=xi, cHbeta=cHbeta, a_H=a_H, Eqw_hbeta=Eqw_hbeta)
        
        #Calculate Helium theoretical flux    
        @pymc.deterministic
        def det_He_Flux_theo_abs(Te=Te, ne=ne, tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus, a_H=a_H, a_He=a_He, Eqw_hbeta=Eqw_hbeta):
            return self.He_Flux_theo_abs(Te=Te, ne=ne,  tau=tau, xi=xi, cHbeta=cHbeta, y_plus = y_plus, a_H=a_H, a_He=a_He, Eqw_hbeta=Eqw_hbeta)
        
        #Combine theoretical fluxes into a single array
        @pymc.deterministic
        def H_He_TheoFlux(HFlux_theo=det_HFlux_theo_abs, HeFlux_theo=det_He_Flux_theo_abs):
            return concatenate([HFlux_theo, HeFlux_theo])
        
#         @pymc.deterministic
#         def det_H_He_Corr_ObsFlux(H_He_ObsFlux = self.H_He_Obs_Flux, Te=Te, ne=ne):
#             return self.H_He_Corr_ObsFlux(H_He_ObsFlux=H_He_ObsFlux, Te=Te, ne=ne)
            
        #Likelihood
        @pymc.stochastic(observed=True)
        def Likelihood_model(value = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2
        
        #Deterministic method to track the evolution of the chi:
        @pymc.deterministic()
        def ChiSq(H_He_ObsFlux = self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
            return - chi_F / 2
       
        return locals()
    