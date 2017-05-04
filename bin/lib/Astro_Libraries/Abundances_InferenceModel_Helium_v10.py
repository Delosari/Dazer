from collections                    import OrderedDict

from numpy                          import array, loadtxt, genfromtxt, sum, exp, power, zeros, ones, dot, isnan, square, concatenate, where
from numpy.random                   import normal
import pymc
from uncertainties                  import ufloat

from CodeTools.File_Managing_Tools  import Txt_Files_Manager
import pyneb                        as pn


class HeAbundance_InferenceMethods(Txt_Files_Manager):
    
    def __init__(self):
        
        #Import Dazer_Files Class to import data from lines_logs
        Txt_Files_Manager.__init__(self)

        #Declare Hydrogen and Helium collisional lines for the analysis
        #self.posHydrogen_Ions      = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
        self.posHydrogen_Lines      = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
        self.Hydrogen_Wavelengths   = [4101.742,        4340.471,       4862.683,        6562.819]
        
        self.posHelium_Lines        = ['H1_3889A',  'He1_4026A',    'He1_4387A',    'He1_4471A',    'He1_4686A',    'He1_4714A',    'He1_4922A',   'He1_5876A',    'He1_6678A',   'He1_7065A',    'He1_7281A',      'He1_10380A']     
        self.Helium_Wavelengths     = [3889.0,         4026.0,         4387.0,         4471.0,         4686.0,         4714.0,         4922.0,         5876.0,         6678.0,         7065.0,         7281.0,         10830.0]
        
        self.Cand_Hydrogen_Lines    = []
        self.Cand_Helium_Lines      = []

        self.nHydrogen              = None
        self.nHelium                = None
    
        #Define indexes and labels to speed up the code
        self.H13889A_label          = 'H1_3889A'          
        self.HBeta_label            = 'H1_4861A'
        self.He3889_label           = 'He1_3889A'
        self.He3889blended_label    = 'He1_3889A_blended'
        self.He3889_check           = None
        
        #Set up the right emissivities
        pn.atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
        
        #Declare pyneb Hydrogen and Helium atoms to calculate emissivities
        self.H1                     = pn.RecAtom('H', 1)
        self.He1                    = pn.RecAtom('He', 1) 
        
        print 'I am using this atomic data'
        print self.He1.printSources()
        
        #Import collisional coefficients table
        self.Coef_Kalpha_dict       = self.Import_Coll_Coeff_Table(self.posHydrogen_Lines, None)

        #Import Optical depth function
        self.Coef_ftau_dict         = self.Import_OpticalDepth_Coeff_Table(self.posHelium_Lines) 
               
        #Declare dictionaries to store the data
        #WARNING: In this analysis flambda is mantained constant
        self.Flux_dict              = OrderedDict()
        self.Error_dict             = OrderedDict()
        self.Wave_dict              = OrderedDict()
        self.PynebCode_dict         = OrderedDict()
        self.EqW_dict               = OrderedDict()
        self.EqWerror_dict          = OrderedDict()
        self.hlambda_dict           = OrderedDict()
        self.flambda_dict           = self.get_flambda_dict(self.posHydrogen_Lines + self.posHelium_Lines, self.Hydrogen_Wavelengths + self.Helium_Wavelengths)
        
        #Extra
        self.EmptyRowFormat         = 'nan'            
            
    def Import_TableData(self, Address, Columns):
        
        Imported_Array  = genfromtxt(Address, dtype=float, usecols = Columns, skip_header = 2).T
        Datarray = Imported_Array[:,~isnan(Imported_Array).all(0)]
                
        return Datarray
     
    def Import_Coll_Coeff_Table(self, HydrogenLines, HeliumLines):
        
        data_dict = OrderedDict()
                
        for i in range(len(HydrogenLines)):
            
            Line            = HydrogenLines[i]
            Data_Columns    = [0 + 3*i, 1 + 3*i, 2 + 3*i]
            data_dict[Line] = self.Import_TableData(self.Hydrogen_CollCoeff_TableAddress, Data_Columns)
                    
        return data_dict
     
    def Import_OpticalDepth_Coeff_Table(self, HeliumLines):
        
        data_dict = OrderedDict()
        
        for i in range(len(HeliumLines)):
            
            Line = HeliumLines[i]
            
            if Line == self.H13889A_label:
                data_dict[self.He3889_label] = loadtxt(self.Helium_OpticalDepth_TableAddress, dtype = float, skiprows = 2, usecols = (i,))
            else:
                data_dict[Line] = loadtxt(self.Helium_OpticalDepth_TableAddress, dtype = float, skiprows = 2, usecols = (i,))
            
        return data_dict
    
    def Import_Synthetic_Fluxes(self, Model = 'Model2'):

        self.SyntheticData_dict = OrderedDict()

        if Model == 'Model1':
            #self.SyntheticData_dict['He_Flux']      = array([ 0.0891384,    0.01283333,     0.03187135, 0.09890274, 0.02582588, 0.02366021])   
            #self.SyntheticData_dict['He_Flux']      = array([ 0.0891384,    0.01283333,     0.03187135, 0.09890274, 0.02582588, 0.02366021])
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
            
            #These are the ones from the version 0
            #self.SyntheticData_dict['H_Flux']       = array([ 0.24585991,  0.45343263,  2.95803687])
            #self.SyntheticData_dict['He_Flux']      = array([ 0.10059744,  0.01684244,  0.03908421,  0.12673761,  0.03538569,  0.04145339])
            #These are the ones from the version 1
            #self.SyntheticData_dict['H_Flux']       = array([0.24625485, 0.454443,      1,              2.98352834])
            #self.SyntheticData_dict['He_Flux']      = array([0.10100783,    0.01691781,     0.03924855,     0.12725256,  0.03553524,    0.04162797])
            
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
                                
        return
    
    def Calculate_Synthetic_Fluxes(self, model):
        
        self.Import_Synthetic_Fluxes(Model = model)
        
        self.Check_EmissionLinesObserved(self.SyntheticData_dict)
        
        self.Preload_Data(Hbeta_Normalized=True, Deblend_Check=False)
        
       
        Hydrogen_Dict   = self.Calculate_H_Parameters_Vectors(self.SyntheticData_dict['Te_0'], self.SyntheticData_dict['n_e'], self.SyntheticData_dict['tau'])
        Helium_Dict     = self.Calculate_He_Parameters_Vectors(self.SyntheticData_dict['Te_0'], self.SyntheticData_dict['n_e'], self.SyntheticData_dict['tau'])
    
        H_Flux          = self.H_Flux_theo(self.SyntheticData_dict['xi'], self.SyntheticData_dict['cHbeta'], Hydrogen_Dict)
        He_Flux         = self.He_Flux_theo_nof(self.SyntheticData_dict['xi'], self.SyntheticData_dict['cHbeta'], self.SyntheticData_dict['y_plus'], Helium_Dict)
        
        print '\nInput Parameters'
        print 'y_plus' ,self.SyntheticData_dict['y_plus']
        print 'n_e', self.SyntheticData_dict['n_e']
        print 'Te_0', self.SyntheticData_dict['Te_0']
        print 'cHbeta', self.SyntheticData_dict['cHbeta']
        print 'tau', self.SyntheticData_dict['tau']
        print 'xi', self.SyntheticData_dict['xi']
        
        print '\nHydrogen emissivities'
        for i in range(len(self.Emissivity_H_vector)):
            print self.Cand_Hydrogen_Lines[i], self.Emissivity_H_vector[i]/Hydrogen_Dict['Hbeta_emis'] 

        print '\nHelium emissivities'
        for i in range(len(self.Emissivity_He_vector)):
            print self.Cand_Helium_Lines[i], self.Emissivity_He_vector[i]/Helium_Dict['Hbeta_emis'] 
        
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
        self.ObjectData_dict['H_hlambda']        = h_lambda
        
        self.ObjectData_dict['He_labels']       = list(Labels)
        self.ObjectData_dict['He_Ions']         = list(Ions)
        self.ObjectData_dict['He_Flux']         = Fluxes
        self.ObjectData_dict['He_error']        = Fluxes_error
        self.ObjectData_dict['He_wave']         = Wavelengths
        self.ObjectData_dict['He_Eqw']          = Eqw
        self.ObjectData_dict['He_EqwErr']       = Eqw_error
        self.ObjectData_dict['He_hlambda']      = h_lambda
        
        #TOIII calculated from TSIII, if not we use TOIII, if not we use 17000
        TOIII_from_TSIII                        = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'TOIII_approx_TSIII', Assumption='float').nominal_value
        TOIII_from_TSIII_error                  = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'TOIII_approx_TSIII', Assumption='float').std_dev
        TOIII                                   = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'TOIII', Assumption='float').nominal_value
        TOIII_error                             = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'TOIII', Assumption='float').std_dev
        
        if TOIII_from_TSIII != None:
            self.ObjectData_dict['TOIII']       = TOIII_from_TSIII
            self.ObjectData_dict['TOIII_error'] = TOIII_from_TSIII_error
        elif TOIII != None:
            self.ObjectData_dict['TOIII']       = TOIII
            self.ObjectData_dict['TOIII_error'] = TOIII_error
        else: 
            self.ObjectData_dict['TOIII']       = 12000.0
            self.ObjectData_dict['TOIII_error'] = 1000.0
            
        #nSII density to be used in the priors definition. If it could not be calculated the output value is 100
        ne_obj = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'nSII', Assumption = 'float')
        if (ne_obj < 0) or (ne_obj == None):
            self.ObjectData_dict['nSII']        = 50.0
            self.ObjectData_dict['nSII_error']  = 50.0
        else:
            self.ObjectData_dict['nSII']        = ne_obj.nominal_value
            self.ObjectData_dict['nSII_error']  = ne_obj.std_dev          
            
        self.ObjectData_dict['cHbeta_obs']  = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'cHBeta', Assumption = 'cHbeta_min')
        return
    
    def Check_EmissionLinesObserved(self, Data_Dict):
        
        #Empty the lists of candidate lines for Hydrogen and Helium
        del self.Cand_Hydrogen_Lines[:]
        del self.Cand_Helium_Lines[:]
                
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
                self.Cand_Helium_Lines.append(Line_Label)
        
        print 'Hydrogen lines'
        print self.Cand_Hydrogen_Lines
        print 'Helium lines'
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

            He_3889A_DeblendFlux                        = self.Deblend_He3889A_line(Te = self.ObjectData_dict['TOIII'], ne = self.ObjectData_dict['nSII'], cHbeta = self.ObjectData_dict['cHbeta_obs'])
            self.Flux_dict[self.He3889_label]           = He_3889A_DeblendFlux #This does not delete the H1_3889 previous entry
            self.Cand_Helium_Lines[self.Cand_Helium_Lines.index(self.H13889A_label)] = self.He3889_label
    
        print 'Helium lines apanadas'
        print self.Cand_Helium_Lines
        
        print '\nMy input hydrogen values:'
        for hydrogen_line in self.Cand_Hydrogen_Lines:
            print hydrogen_line, self.Flux_dict[hydrogen_line], self.Error_dict[hydrogen_line], self.PynebCode_dict[hydrogen_line]

        print '\nMy input helium values:'
        for helium_line in self.Cand_Helium_Lines:
            print helium_line, self.Flux_dict[helium_line], self.Error_dict[helium_line], self.PynebCode_dict[helium_line]
                        
        #Define number of lines variables to speed up the calculation
        self.nHydrogen                                      = len(self.Cand_Hydrogen_Lines)
        self.nHelium                                        = len(self.Cand_Helium_Lines)
        
        self.nHydrogen_range                                = range(self.nHydrogen)
        self.nHelium_range                                  = range(self.nHelium)
                
        #Oxygen temperature prior
        self.TOIIIprior                                     = Data_Dict['TOIII']
        self.TOIIIprior_sigma                               = self.TOIIIprior * 0.2
        
    def Preload_Data(self, Hbeta_Normalized = True, Deblend_Check = True):

        #Define the vectors where we store the physical parameters
        self.Flux_H_vector                          = zeros(self.nHydrogen)
        self.Error_H_vector                         = zeros(self.nHydrogen)
        self.Emissivity_H_vector                    = zeros(self.nHydrogen)
        self.Kalpha_vector                          = zeros(self.nHydrogen)
        self.flambda_H_vector                       = zeros(self.nHydrogen)
        self.hlambda_H_vector                       = ones(self.nHydrogen)

        self.Flux_He_vector                         = zeros(self.nHelium)
        self.Error_He_vector                        = zeros(self.nHelium)
        self.Emissivity_He_vector                   = zeros(self.nHelium)
        self.flambda_He_vector                      = zeros(self.nHelium)
        self.hlambda_He_vector                      = ones(self.nHelium)
      
        #Define HBeta values in advance
        self.Hbeta_Flux                             = self.Flux_dict[self.HBeta_label]
        self.Hbeta_error                            = self.Error_dict[self.HBeta_label]
        self.hbeta_hlambda                          = self.hlambda_dict[self.HBeta_label]
      
        #Calculate in advance observed physical parameters vectors for Hydrogen lines
        for i in range(self.nHydrogen):
            Line_Label                              = self.Cand_Hydrogen_Lines[i]        
            self.Flux_H_vector[i]                   = self.Flux_dict[Line_Label]
            self.Error_H_vector[i]                  = self.Error_dict[Line_Label]
            self.hlambda_H_vector[i]                = self.hlambda_dict[Line_Label]
            self.flambda_H_vector[i]                = self.flambda_dict[Line_Label]

        for i in range(self.nHelium):
            Line_Label                              = self.Cand_Helium_Lines[i]
            self.Flux_He_vector[i]                  = self.Flux_dict[Line_Label]
            self.Error_He_vector[i]                 = self.Error_dict[Line_Label]
            self.flambda_He_vector[i]               = self.flambda_dict[Line_Label]
            self.hlambda_He_vector[i]               = self.hlambda_dict[Line_Label]             

        #Combine Hydrogen and Helium values into a single vector
        if Hbeta_Normalized: 
            self.H_He_Obs_Flux                      = concatenate([self.Flux_H_vector, self.Flux_He_vector])
            self.H_He_Obs_Error                     = concatenate([self.Error_H_vector, self.Error_He_vector])

        else:
            self.H_He_Obs_Flux                      = concatenate([self.Flux_H_vector, self.Flux_He_vector])    / self.Hbeta_Flux
            self.H_He_Obs_Error                     = concatenate([self.Error_H_vector, self.Error_He_vector])  / self.Hbeta_Flux
        
        print '\nMy input 2 hydrogen values:'
        for j in range(self.nHydrogen):
            Line_Label                              = self.Cand_Hydrogen_Lines[j]                
            print Line_Label, self.H_He_Obs_Flux[j], self.H_He_Obs_Error[j]
        
        print '\nMy input 2 helium values:'
        for j in range(self.nHelium):
            Line_Label                              = self.Cand_Helium_Lines[j]
            print Line_Label, self.H_He_Obs_Flux[j + self.nHydrogen], self.H_He_Obs_Error[j]    
                
        #Dictionary with the vectors
        self.Data_Dict = {'Hbeta_Kalpha'    : None, 
                     'Hbeta_emis'           : None,
                     'Emissivity_H_vector'  : None,
                     'Kalpha_vector'        : None,
                     'Emissivity_He_vector' : None,
                     'T_4'                  : None} 
        
    def Calculate_H_Parameters_Vectors(self, T_e, n_e):    
        
        #Calculate in advance the T_4 parameter
        self.Data_Dict['T_4'] = T_e / 10000.0
        
        #Calculate the hbeta parameters
        self.Data_Dict['Hbeta_Kalpha']              = self.Kalpha_Ratio_H(T_4 = self.Data_Dict['T_4'], H_label=self.HBeta_label)
        self.Data_Dict['Hbeta_emis']                = self.H1.getEmissivity(T_e, n_e, label = self.PynebCode_dict[self.HBeta_label])

        #Calculate physical parameters for Hydrogen lines
        for i in self.nHydrogen_range:
            self.Emissivity_H_vector[i]             = self.H1.getEmissivity(T_e, n_e, label = self.PynebCode_dict[self.Cand_Hydrogen_Lines[i]])
            self.Kalpha_vector[i]                   = self.Kalpha_Ratio_H(T_4 = self.Data_Dict['T_4'], H_label=self.Cand_Hydrogen_Lines[i])
        
        self.Data_Dict['Emissivity_H_vector']       = self.Emissivity_H_vector
        self.Data_Dict['Kalpha_vector']             = self.Kalpha_vector
                
        return self.Data_Dict    
    
    def Calculate_He_Parameters_Vectors(self, T_e, n_e):
    
        #Calculate physical parameters for Helium lines
        for i in self.nHelium_range:
            self.Emissivity_He_vector[i]            = self.He1.getEmissivity(T_e, n_e, label = self.PynebCode_dict[self.Cand_Helium_Lines[i]])
        
        self.Data_Dict['Emissivity_He_vector']      = self.Emissivity_He_vector
        
        return self.Data_Dict

    def H_Flux_theo(self, xi, cHbeta, data_dict):
                
        #Calculate the emissivities for each of the lines for the given temperature and density
        Emissivities_module = data_dict['Emissivity_H_vector'] / data_dict['Hbeta_emis'] 
                
        #Calculate the Collisional excitation fraction
        CR_Module           = (1.0 + 0.0001* xi *data_dict['Kalpha_vector']) / (1.0 + 0.0001* xi * data_dict['Hbeta_Kalpha'] )
                
        #Calculate the reddening component
        f_module            = power(10, -1 * self.flambda_H_vector * cHbeta)
          
        #Calculate theoretical Hydrogen flux for each emission line
        H_Flux              = Emissivities_module * CR_Module * f_module
                     
        return H_Flux

    def He_Flux_theo_nof(self, xi, cHbeta, y_plus, data_dict):
                      
        #Calculate the emissivities for each of the lines for the given temperature and density
        Emissivities_module = data_dict['Emissivity_He_vector'] / data_dict['Hbeta_emis']
                
        #Calculate the collisional excitation fraction
        CR_Module           = 1 / (1.0 + 0.0001* xi * data_dict['Hbeta_Kalpha'])

        #Calculate the reddening component
        f_module            = power(10, -1 * self.flambda_He_vector * cHbeta)
             
        #Calculate theoretical Hydrogen flux for each emission line
        He_Flux             = y_plus * Emissivities_module * CR_Module * f_module

        return He_Flux
               
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
    
    def Deblend_He3889A_line(self, Te, ne, cHbeta):
        
        Emis_Hydrogen3889   = 0.104 * power(Te / 10000, 0.046)

#         Emis_Hbeta          = self.H1.getEmissivity(Te, ne, label = self.PynebCode_dict[self.HBeta_label]) 
#         print 'self.Flux_dict[self.H13889A_label]', self.Flux_dict[self.H13889A_label]
#         print 'self.Flux_dict[self.HBeta_label]', self.Flux_dict[self.HBeta_label]
#         print 'Flux Ratio', self.Flux_dict[self.H13889A_label] / self.Flux_dict[self.HBeta_label]
#         print 'Te', Te
#         print 'ne', ne
#         print 'Emis_Hydrogen3889', Emis_Hydrogen3889, 'comparada con', self.H1.getEmissivity(Te, ne, label = 'H1_3889A') / self.H1.getEmissivity(Te, ne, label = self.PynebCode_dict[self.HBeta_label])
#         print 'Emis_Hbeta', Emis_Hbeta        
#         print 'cHbeta', cHbeta
#         print '(power(10 , - self.flambda_dict[self.He3889_label] * cHbeta))', (power(10 , - self.flambda_dict[self.He3889_label] * cHbeta))        
        
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
        cHbeta      =   pymc.TruncatedNormal(   'c_Hbeta',  0.15,                               0.05**-2,   a = 0.0,    b = 3.0)
        xi          =   pymc.TruncatedNormal(   'Xi',       1,                                  200**-2,    a = 0.0,    b = 1000.0)
#       a_H         =   pymc.TruncatedNormal(   'abs_H',    1.0,                                1**-2,      a = 0.0,    b = 14.0)
#       a_He        =   pymc.TruncatedNormal(   'abs_He',   1.0,                                0.5**-2,    a = 0.0,    b = 5.0)
#       Eqw_hbeta   =   pymc.Normal(            'Eqw_hbeta',self.EqW_dict['H1_4861A'],          self.EqWerror_dict['H1_4861A']**-2)

        #Divide in groups the data vectors which can be calculated simultaneously
        @pymc.deterministic
        def Calculate_H_LineParameters(Te=Te, ne=ne): 
            return self.Calculate_H_Parameters_Vectors(T_e=Te, n_e=ne)            

        @pymc.deterministic
        def Calculate_He_LineParameters(Te=Te, ne=ne): 
            return self.Calculate_He_Parameters_Vectors(T_e=Te, n_e=ne)  

        #Calculate Hydrogen theoretical flux
        @pymc.deterministic
        def det_HFlux_theo(xi=xi, cHbeta=cHbeta, data_dict = Calculate_H_LineParameters):
            return self.H_Flux_theo(xi=xi, cHbeta=cHbeta, data_dict = data_dict)
        
        #Calculate Helium theoretical flux    
        @pymc.deterministic
        def det_He_Flux_theo_nof(xi=xi, cHbeta=cHbeta, y_plus = y_plus, data_dict = Calculate_He_LineParameters):
            return self.He_Flux_theo_nof(xi=xi, cHbeta=cHbeta, y_plus = y_plus, data_dict = data_dict)
        
        #Combine theoretical fluxes into a single array
        @pymc.deterministic
        def H_He_TheoFlux(HFlux_theo=det_HFlux_theo, HeFlux_theo=det_He_Flux_theo_nof):
            return concatenate([HFlux_theo, HeFlux_theo])
                
        #Likelihood
        @pymc.stochastic(observed=True)
        def Likelihood_model(value=self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2
        
        #Deterministic method to track the evolution of the chi:
        @pymc.deterministic()
        def ChiSq(H_He_ObsFlux=self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error):
            chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
            return - chi_F / 2
       
        return locals()