from collections                    import OrderedDict

from numpy                          import array, loadtxt, genfromtxt, sum, exp, power, zeros, ones, dot, isnan, square, concatenate, where
from numpy.random                   import normal
import pymc
from pyneb                          import RecAtom

from CodeTools.File_Managing_Tools  import Txt_Files_Manager


class HeAbundance_InferenceMethods(Txt_Files_Manager):
    
    def __init__(self):
        
        #Import Dazer_Files Class to import data from lines_logs
        Txt_Files_Manager.__init__(self)
        
        #Declare Hydrogen and Helium collisional lines for the analysis
        #self.posHydrogen_Ions      = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
        self.posHydrogen_Lines      = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
        self.Hydrogen_Wavelengths   = [4101.742,        4340.471,       4862.683,        6562.819]
        
        self.posHelium_Lines        = ['He1_3889A',  'He1_4026A',    'He1_4387A',    'He1_4471A',    'He1_4686A',    'He1_4714A',    'He1_4922A',   'He1_5876A',    'He1_6678A',   'He1_7065A',    'He1_7281A',      'He1_10380A']     
        self.Helium_Wavelengths     = [3889.0,        4026.0,         4387.0,         4471.0,         4686.0,         4714.0,         4922.0,         5876.0,         6678.0,         7065.0,         7281.0,         10380.0]
        
        self.Cand_Hydrogen_Lines    = []
        self.Cand_Helium_Lines      = []
        
        self.nHydrogen              = None
        self.nHelium                = None
    
        #Define indexes and labels to speed up the code               
        self.HBeta_label            = 'H1_4861A'
        self.He3889_label           = 'He1_3889A'
        self.He3889blended_label    = 'He1_3889A_blended'
        
        #Declare pyneb Hydrogen and Helium atoms to calculate emissivities
        self.H1                     = RecAtom('H', 1)
        self.He1                    = RecAtom('He', 1) 
        
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
        
        Imported_Array  = genfromtxt(Address, dtype=float, usecols = Columns, skiprows = 2).T
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
            
            Line            = HeliumLines[i]
            data_dict[Line] = loadtxt(self.Helium_OpticalDepth_TableAddress, dtype = float, skiprows = 2, usecols = (i,))
            
        return data_dict
    
    def Import_Synthetic_Fluxes(self, Model = 'Model2'):

        self.SyntheticData_dict = OrderedDict()

        if Model == 'Model1':
            #self.SyntheticData_dict['He_Flux']      = array([ 0.0891384,    0.01283333,     0.03187135, 0.09890274, 0.02582588, 0.02366021])   
            #self.SyntheticData_dict['He_Flux']      = array([ 0.0891384,    0.01283333,     0.03187135, 0.09890274, 0.02582588, 0.02366021])
            self.SyntheticData_dict['H_Flux']       = array([ 0.24685935,   0.45526236,     1.0,            2.96988502])
            self.SyntheticData_dict['H_wave']       = array([4101.742,      4340.471,      4862.683,       6562.819])
            self.SyntheticData_dict['H_Ions']       = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
            self.SyntheticData_dict['H_labels']     = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
            self.SyntheticData_dict['H_Eqw']        = [50.0,            50.0,           250,             300.0]
            self.SyntheticData_dict['H_EqwErr']     = [1.0,             1.0,            12.5,            14]
            self.SyntheticData_dict['H_hlambda']    = [1.0,             1.0,            1.0,             1.0]
            self.SyntheticData_dict['H_error']      = self.SyntheticData_dict['H_Flux'] * 0.01
            
            self.SyntheticData_dict['He_Flux']      = array([ 0.08951096,    0.01290066,     0.03201484, 0.09931435, 0.02594518, 0.02377085])   
            self.SyntheticData_dict['He_wave']      = array([ 3889.0,       4026.0,         4471.0,        5876.0,     6678.0,     7065.0])
            self.SyntheticData_dict['He_labels']    =       ['He1_3889A',   'He1_4026A',   'He1_4471A',   'He1_5876A','He1_6678A','He1_7065A']
            self.SyntheticData_dict['He_Eqw']       = [10.0,                10.0,           10.0,           10.0,        10.0,          10.0]
            self.SyntheticData_dict['He_EqwErr']    = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0]
            self.SyntheticData_dict['He_hlambda']   = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0] 
            self.SyntheticData_dict['He_error']     = self.SyntheticData_dict['He_Flux'] * 0.02
                        
            self.SyntheticData_dict['TOIII']        = 19000
            self.SyntheticData_dict['y_plus']       = 0.08
            self.SyntheticData_dict['n_e']          = 100.0
            self.SyntheticData_dict['a_He']         = 1.0
            self.SyntheticData_dict['tau']          = 0.2
            self.SyntheticData_dict['Te_0']         = 18000.0
            self.SyntheticData_dict['cHbeta']       = 0.1
            self.SyntheticData_dict['a_H']          = 1.0
            self.SyntheticData_dict['xi']           = 1.0

        if Model == 'Model2':
            
            #These are the ones from the version 0
            #self.SyntheticData_dict['H_Flux']       = array([ 0.24585991,  0.45343263,  2.95803687])
            #self.SyntheticData_dict['He_Flux']      = array([ 0.10059744,  0.01684244,  0.03908421,  0.12673761,  0.03538569,  0.04145339])
            #These are the ones from the version 1
            self.SyntheticData_dict['H_Flux']       = array([0.24625485, 0.454443,      1,              2.98352834])
            self.SyntheticData_dict['H_wave']       = array([4101.742,   4340.471,      4862.683,       6562.819])
            self.SyntheticData_dict['H_Ions']       = ['Hdelta_6_2',    'Hgamma_5_2',   'Hbeta_4_2',    'Halpha_3_2']
            self.SyntheticData_dict['H_labels']     = ['H1_4102A',      'H1_4340A',     'H1_4861A',     'H1_6563A']
            self.SyntheticData_dict['H_Eqw']        = [50.0,            50.0,           250,             300.0]
            self.SyntheticData_dict['H_EqwErr']     = [1.0,             1.0,            12.5,            14]
            self.SyntheticData_dict['H_hlambda']    = [1.0,             1.0,            1.0,             1.0]
            self.SyntheticData_dict['H_error']      = self.SyntheticData_dict['H_Flux'] * 0.01
            
            self.SyntheticData_dict['He_Flux']      = array([0.10100783,    0.01691781,     0.03924855,     0.12725256,  0.03553524,    0.04162797])
            self.SyntheticData_dict['He_wave']      = array([ 3889.0,       4026.0,         4471.0,         5876.0,      6678.0,        7065.0])
            self.SyntheticData_dict['He_labels']    = ['He1_3889A',         'He1_4026A',   'He1_4471A',   'He1_5876A',  'He1_6678A',    'He1_7065A']
            self.SyntheticData_dict['He_Eqw']       = [10.0,                10.0,           10.0,           10.0,        10.0,          10.0]
            self.SyntheticData_dict['He_EqwErr']    = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0]
            self.SyntheticData_dict['He_hlambda']   = [1.0,                 1.0,            1.0,            1.0,         1.0,           1.0] 
            self.SyntheticData_dict['He_error']     = self.SyntheticData_dict['He_Flux'] * 0.02

            self.SyntheticData_dict['TOIII']        = 17000
            self.SyntheticData_dict['y_plus']       = 0.085
            self.SyntheticData_dict['n_e']          = 500.0
            self.SyntheticData_dict['a_He']         = 0.5
            self.SyntheticData_dict['tau']          = 1.0
            self.SyntheticData_dict['Te_0']         = 16000.0
            self.SyntheticData_dict['cHbeta']       = 0.1
            self.SyntheticData_dict['a_H']          = 1.0
            self.SyntheticData_dict['xi']           = 1.0
                                   
        return
    
    def Import_Object_Data(self, Filefolder, CodeName, Lineslog_Extension):
        
        self.ObjectData_dict = OrderedDict()
            
        String_Columns      = [self.Labels_ColumnHeader, self.Ion_ColumnHeader]
        Float_Columns       = [self.Wavelength_ColumnHeader, self.GaussianFlux_ColumnHeader, self.GaussianError_ColumnHeader, self.Line_Continuum_ColumnHeader, self.EqW_ColumnHeader, self.EqW_error_ColumnHeader]
        
        Lineslog_Address    = Filefolder + CodeName + Lineslog_Extension
        
        Labels, Ions                                                = self.get_ColumnData(String_Columns, Lineslog_Address, HeaderSize = 2, StringIndexes = True, datatype = str, unpack_check = True) 
        Wavelengths, Fluxes, Fluxes_error, h_lambda, Eqw, Eqw_error = self.get_ColumnData(Float_Columns, Lineslog_Address, HeaderSize = 2, StringIndexes = True, unpack_check = True)

        self.ObjectData_dict['H_labels']     = list(Labels)
        self.ObjectData_dict['H_Ions']       = list(Ions)
        self.ObjectData_dict['H_Flux']       = Fluxes
        self.ObjectData_dict['H_error']      = Fluxes_error
        self.ObjectData_dict['H_wave']       = Wavelengths
        self.ObjectData_dict['H_Eqw']        = Eqw
        self.ObjectData_dict['H_EqwErr']     = Eqw_error
        self.ObjectData_dict['H_hlambda']    = h_lambda
        
        self.ObjectData_dict['He_labels']    = list(Labels)
        self.ObjectData_dict['He_Ions']      = list(Ions)
        self.ObjectData_dict['He_Flux']      = Fluxes
        self.ObjectData_dict['He_error']     = Fluxes_error
        self.ObjectData_dict['He_wave']      = Wavelengths
        self.ObjectData_dict['He_Eqw']       = Eqw
        self.ObjectData_dict['He_EqwErr']    = Eqw_error
        self.ObjectData_dict['He_hlambda']   = h_lambda
        
        #TOIII calculated from TSIII, if not we use TOIII, if not we use 17000
        TOIII_from_TSIII    = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'TOIII_approx_TSIII', Assumption = 'float')
        TOIII               = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'TOIII', Assumption = 'float')
        
        if TOIII_from_TSIII != None:
            self.ObjectData_dict['TOIII']   = TOIII_from_TSIII
        elif TOIII != None:
            self.ObjectData_dict['TOIII']   = TOIII
        else: 
            self.ObjectData_dict['TOIII']   = 15000.0
        
        #nSII density to be used in the priors definition. If it could not be calculated the output value is 100
        nSII               = self.GetParameter_ObjLog(CodeName, Filefolder, Parameter = 'nSII', Assumption = 'Min_Den')
        
        self.ObjectData_dict['nSII'] =  nSII
        
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
                self.EqW_dict[Line_Label]                   = Data_Dict['H_Eqw'][Line_Index]
                self.EqWerror_dict[Line_Label]              = Data_Dict['H_EqwErr'][Line_Index]
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
                self.EqW_dict[Line_Label]                   = Data_Dict['He_Eqw'][Line_Index]
                self.EqWerror_dict[Line_Label]              = Data_Dict['He_EqwErr'][Line_Index]
                self.hlambda_dict[Line_Label]               = Data_Dict['He_hlambda'][Line_Index]
                self.PynebCode_dict[Line_Label]             = Line_Label[Line_Label.find('_')+1:-1]+'.0'
                self.Cand_Helium_Lines.append(Line_Label)
        
        print 'Hydrogen lines'
        print self.Cand_Hydrogen_Lines
        print 'Helium lines'
        print self.Cand_Helium_Lines
        
        #Since the Hbeta line is not used for the fitting we remove it from the candidates lines list. However, its properties are still stored in the dictionaries
        self.Cand_Hydrogen_Lines.remove(self.HBeta_label)
        
        #Since the He_3889A line is blended we add new key He_3889A_blended which will be store the observed, blended, value
        if self.He3889_label in self.Cand_Helium_Lines:
            self.Flux_dict[self.He3889blended_label]        = self.Flux_dict[self.He3889_label]
        
        #Define number of lines variables to speed up the calculation
        self.nHydrogen                                      = len(self.Cand_Hydrogen_Lines)
        self.nHelium                                        = len(self.Cand_Helium_Lines)
        
        self.nHydrogen_range                                = range(self.nHydrogen)
        self.nHelium_range                                  = range(self.nHelium)
                
        #Oxygen temperature prior
        self.TOIIIprior                                     = Data_Dict['TOIII']
        self.TOIIIprior_sigma                               = self.TOIIIprior * 0.2
        
    def Preload_Data(self, Hbeta_Normalized = True):

        #Define the vectors where we store the physical parameters
        self.Flux_H_vector                          = zeros(self.nHydrogen)
        self.Error_H_vector                         = zeros(self.nHydrogen)
        self.Emissivity_H_vector                    = zeros(self.nHydrogen)
        self.Kalpha_vector                          = zeros(self.nHydrogen)
        self.flambda_H_vector                       = zeros(self.nHydrogen)
        #self.Eqw_H_vector                          = zeros(self.nHydrogen)
        #self.EqwError_H_vector                     = zeros(self.nHydrogen)
        self.hlambda_H_vector                       = ones(self.nHydrogen)

        self.Flux_He_vector                         = zeros(self.nHelium)
        self.Error_He_vector                        = zeros(self.nHelium)
        self.Emissivity_He_vector                   = zeros(self.nHelium)
        self.ftau_He_vector                         = zeros(self.nHelium)
        self.flambda_He_vector                      = zeros(self.nHelium)
        #self.Eqw_He_vector                         = zeros(self.nHelium)
        #self.EqwError_He_vector                    = zeros(self.nHelium)
        self.hlambda_He_vector                      = ones(self.nHelium)
      
        #Define HBeta values in advance
        self.Hbeta_Flux                             = self.Flux_dict[self.HBeta_label]
        self.Hbeta_error                            = self.Error_dict[self.HBeta_label]
        self.Hbeta_EW                               = self.EqW_dict[self.HBeta_label]
        self.Hbeta_EWerror                          = self.EqWerror_dict[self.HBeta_label]      
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
            self.H_He_Obs_Flux                      = concatenate([self.Flux_H_vector, self.Flux_He_vector]) / self.Hbeta_Flux
            self.H_He_Obs_Error                     = concatenate([self.Error_H_vector, self.Error_He_vector])  
        
        #Getting the index for the He1_3889A line in advance
        if self.He3889_label in self.Cand_Helium_Lines:
            self.He3889_index = where(self.H_He_Obs_Flux == self.Flux_dict[self.He3889blended_label] / self.Hbeta_Flux)
#             print 'flujo de la linea', self.Flux_dict[self.He3889blended_label]/ self.Hbeta_Flux
#             print 'Veamos como queda. Mi indice de He3889 es ', self.He3889_index, '. Hence:'
#             for i in range(len(self.H_He_Obs_Flux)):
#                 print 'Line', self.H_He_Obs_Flux[i]
        
        
        #Dictionary with the vectors
        self.Data_Dict = {'Hbeta_Kalpha'    : None, 
                     'Hbeta_emis'           : None,
                     'Emissivity_H_vector'  : None,
                     'Kalpha_vector'        : None,
                     'Emissivity_He_vector' : None,
                     'ftau_He_vector'       : None,
                     'T_4'                  : None} 
        
    def Calculate_Parameters_Vectors(self, T_e, n_e, tau):    
        
        #Calculate in advance the T_4 parameter
        T_4 = T_e / 10000.0
        
        #Calculate the hbeta parameters
        self.Data_Dict['Hbeta_Kalpha']              = self.Kalpha_Ratio_H(T_4 = T_4, H_label=self.HBeta_label)
        self.Data_Dict['Hbeta_emis']                = self.H1.getEmissivity(T_e, n_e, label = self.PynebCode_dict[self.HBeta_label])

        #Calculate physical parameters for Hydrogen lines
        for i in self.nHydrogen_range:
            self.Emissivity_H_vector[i]             = self.H1.getEmissivity(T_e, n_e, label = self.PynebCode_dict[self.Cand_Hydrogen_Lines[i]])
            self.Kalpha_vector[i]                   = self.Kalpha_Ratio_H(T_4 = T_4, H_label=self.Cand_Hydrogen_Lines[i])
            
        #Calculate physical parameters for Helium lines
        for i in self.nHelium_range:
            self.Emissivity_He_vector[i]            = self.He1.getEmissivity(T_e, n_e, label = self.PynebCode_dict[self.Cand_Helium_Lines[i]])
            self.ftau_He_vector[i]                  = self.OpticalDepth_He(tau = tau, T_4 = T_4, n_e = n_e, He_label = self.Cand_Helium_Lines[i])
        
        return self.Data_Dict
    
    def Calculate_H_Parameters_Vectors(self, T_e, n_e, tau):    
        
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
    
    def Calculate_He_Parameters_Vectors(self, T_e, n_e, tau):
    
        #Calculate physical parameters for Helium lines
        for i in self.nHelium_range:
            self.Emissivity_He_vector[i]            = self.He1.getEmissivity(T_e, n_e, label = self.PynebCode_dict[self.Cand_Helium_Lines[i]])
            self.ftau_He_vector[i]                  = self.OpticalDepth_He(tau = tau, T_4 = self.Data_Dict['T_4'], n_e = n_e, He_label = self.Cand_Helium_Lines[i])    
        
        self.Data_Dict['Emissivity_He_vector']      = self.Emissivity_He_vector
        self.Data_Dict['ftau_He_vector']            = self.ftau_He_vector
        
        return self.Data_Dict

    def H_Flux_theo(self, xi, cHbeta, a_H, data_dict):
                
        #Calculate the emissivities for each of the lines for the given temperature and density
        Emissivities_module = data_dict['Emissivity_H_vector'] / data_dict['Hbeta_emis'] 
        
        #Calculate the Collisional excitation fraction
        CR_Module           = (1.0 + 0.0001* xi *data_dict['Kalpha_vector']) / (1.0 + 0.0001* xi * data_dict['Hbeta_Kalpha'] )
        
        #Calculate the reddening component
        f_module            = power(10, -1 * self.flambda_H_vector * cHbeta)
        
        #Calculate the Hbeta normalization module
        EW_Hbeta_module     = (self.Hbeta_EW + a_H) / self.Hbeta_EW
        
        #Calculate stellar absorption module
        a_H_module          = (a_H * self.hlambda_H_vector) / (self.Hbeta_EW * self.hbeta_hlambda)
        
        #Calculate theoretical Hydrogen flux for each emission line
        H_Flux              = Emissivities_module * CR_Module * f_module * EW_Hbeta_module - a_H_module
                     
        return H_Flux

    def He_Flux_theo_nof(self, xi, cHbeta, a_H, a_He, y_plus, data_dict):
                      
        #Calculate the emissivities for each of the lines for the given temperature and density
        Emissivities_module = data_dict['Emissivity_He_vector'] / data_dict['Hbeta_emis']
        
        #Calculate the collisional excitation fraction
        CR_Module           = 1 / (1.0 + 0.0001* xi * data_dict['Hbeta_Kalpha'])

        #Calculate the reddening component
        f_module            = power(10, -1 * self.flambda_He_vector * cHbeta)

        #Calculate the Hbeta normalization module
        EW_Hbeta_module     = (self.Hbeta_EW + a_H) / self.Hbeta_EW
        
        #Calculate stellar absorption module
        a_He_module         = (a_He * self.hlambda_He_vector) / (self.Hbeta_EW * self.hbeta_hlambda)
             
        #Calculate theoretical Hydrogen flux for each emission line
        He_Flux             = y_plus * Emissivities_module * data_dict['ftau_He_vector'] * CR_Module * f_module * EW_Hbeta_module - a_He_module

        #Deblend the He1_3889A line 
#         self.Deblend_He3889A_line(a_H, data_dict['T_4'], cHbeta)

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
    
    def OpticalDepth_He(self, tau, T_4, n_e, He_label):
                        
        f_tau   = 1 + (tau/2) * (self.Coef_ftau_dict[He_label][0] + (self.Coef_ftau_dict[He_label][1] + self.Coef_ftau_dict[He_label][2]*n_e + self.Coef_ftau_dict[He_label][3]*n_e*n_e) * T_4)

        return f_tau
    
    def Deblend_He3889A_line(self, a_H, T_4, cHbeta):
        
        self.H_He_Obs_Flux[self.He3889_index] = (self.Flux_dict[self.He3889blended_label] / self.Flux_dict[self.HBeta_label]) * ((self.EqW_dict[self.He3889_label] + a_H) / self.EqW_dict[self.He3889_label]) - (0.104 * power(T_4, 0.046)) * (power(10 ,-self.flambda_dict[self.He3889_label] * cHbeta))  
        
        return      
        
class HeAbundance_InferenceStructure(HeAbundance_InferenceMethods):

    def __init__(self):

        HeAbundance_InferenceMethods.__init__(self)

    def Bayesian_HeliumAbundance_Analysis(self):
        
        y_plus  =   pymc.Uniform(           'He_abud',  0.060,  0.125)
        ne      =   pymc.TruncatedNormal(   'n_e',      500,    250**-2,    a = 0.0 ,   b = 1000.0, value=500)
        a_He    =   pymc.TruncatedNormal(   'abs_He',   1.0,    0.5**-2,    a = 0.0,    b = 5.0,    value=0.5)
        tau     =   pymc.TruncatedNormal(   'Tau',      0.75,   0.5**-2,    a = 0.0,    b = 7.0)
        Te      =   pymc.Normal(            'T_e',      16000.0,  2000**-2, value = 16000.0)
        cHbeta  =   pymc.TruncatedNormal(   'c_Hbeta',  0.15,   0.05**-2,   a = 0.0,    b = 3.0,    value=0.1)
        a_H     =   pymc.TruncatedNormal(   'abs_H',    1.0,    1**-2,      a = 0.0,    b = 14.0)
        xi      =   pymc.TruncatedNormal(   'Xi',       1,     200**-2,    a = 0.0,    b = 1000.0, value = 1.0)
        #xi     =   pymc.DiscreteUniform(   'Xi',       0,      1000)
                    
        #Divide in groups the data vectors which can be calculated simultaneously
        @pymc.deterministic
        def Calculate_H_LineParameters(Te=Te, ne=ne, tau=tau): 
            return self.Calculate_H_Parameters_Vectors(T_e=Te, n_e=ne, tau=tau)            

        @pymc.deterministic
        def Calculate_He_LineParameters(Te=Te, ne=ne, tau=tau): 
            return self.Calculate_He_Parameters_Vectors(T_e=Te, n_e=ne, tau=tau)  

        #Calculate Hydrogen theoretical flux
        @pymc.deterministic
        def det_HFlux_theo(xi=xi, cHbeta=cHbeta, a_H=a_H, data_dict = Calculate_H_LineParameters):
            return self.H_Flux_theo(xi=xi, cHbeta=cHbeta, a_H=a_H, data_dict = data_dict)
        
        #Calculate Helium theoretical flux    
        @pymc.deterministic
        def det_He_Flux_theo_nof(xi=xi, cHbeta=cHbeta, a_H=a_H, a_He=a_He, y_plus = y_plus, data_dict = Calculate_He_LineParameters):
            return self.He_Flux_theo_nof(xi=xi, cHbeta=cHbeta, a_H=a_H, a_He = a_He, y_plus = y_plus, data_dict = data_dict)
        
        #Combine theoretical fluxes into a single array
        @pymc.deterministic
        def H_He_TheoFlux(HFlux_theo=det_HFlux_theo, HeFlux_theo=det_He_Flux_theo_nof):
            return concatenate([HFlux_theo, HeFlux_theo])
        
        #Chi_EquivalentWidth
        @pymc.deterministic
        def ChiSq_EW(EW_Hbeta_obs=self.Hbeta_EW, EW_Hbeta_sig = self.Hbeta_EWerror):
            EW_HbetaSynth   = normal(EW_Hbeta_obs, EW_Hbeta_sig, 1)
            chi_ew          = square(EW_HbetaSynth - EW_Hbeta_obs) / square(EW_Hbeta_sig)
            return - chi_ew / 2
        
        #Chi_Temperature
        @pymc.deterministic
        def ChiSq_T(T_Synth=Te, T_Meas = self.TOIIIprior, sigma = self.TOIIIprior_sigma):
            chi_tem          = square(T_Synth - T_Meas) / square(sigma)
            return - chi_tem / 2
    
        #Likelihood
        @pymc.stochastic(observed=True)
        def Likelihood_model(value=self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error, ChiSq_EW=ChiSq_EW, ChiSq_T=ChiSq_T):
            chi_F           = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return (- chi_F / 2) + ChiSq_EW + ChiSq_T
        
        #Deterministic method to track the evolution of the chi:
        @pymc.deterministic()
        def ChiSq(H_He_ObsFlux=self.H_He_Obs_Flux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = self.H_He_Obs_Error, ChiSq_EW=ChiSq_EW, ChiSq_T=ChiSq_T):
            chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
            return (- chi_F / 2) + ChiSq_EW + ChiSq_T
       
        return locals()

