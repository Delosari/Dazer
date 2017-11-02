from os                                     import path, name
from sys                                    import argv
from pyneb                                  import atomicData, RecAtom, Atom
from collections                            import OrderedDict
from uncertainties.unumpy                   import nominal_values, std_devs
from Reddening_Corrections                  import ReddeningLaws
from lib.Astro_Libraries.Nebular_Continuum  import NebularContinuumCalculator
from lib.ssp_functions.ssp_synthesis_tools  import ssp_fitter
from DZ_LineMesurer                         import LineMesurer_v2
from numpy                                  import array, loadtxt, genfromtxt, copy, isnan, arange, insert, concatenate, mean, std, power, exp, zeros, square, empty, percentile, random, median, ones, isnan, sum as np_sum, argsort, vstack, hstack, delete, where
import pymc as pymc2
from timeit                                 import default_timer as timer
from uncertainties import ufloat
from pandas import read_excel, read_csv

class Import_model_data(ReddeningLaws):

    def __init__(self):

        #Import sub classes
        ReddeningLaws.__init__(self)

        self.paths_dict = {}
        
        #Paths for linux:
        if name == 'posix':
            self.paths_dict['inference_folder']     = '/home/vital/Astrodata/Inference_output/'
            self.paths_dict['nebular_data_folder']  = '/home/vital/Dropbox/Astrophysics/Lore/NebularContinuum/'
            self.paths_dict['Hydrogen_CollCoeff']   = '/home/vital/workspace/dazer/bin/lib/Astro_Libraries/Neutral_Hydrogen_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_CollCoeff']     = '/home/vital/workspace/dazer/bin/lib/Astro_Libraries/Neutral_Helium_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_OpticalDepth']  = '/home/vital/workspace/dazer/bin/lib/Astro_Libraries/Helium_OpticalDepthFunction_Coefficients.txt'
            self.paths_dict['lines_data_file']      = '/home/vital/workspace/dazer/bin/lib/Astro_Libraries/lines_data.xlsx'
            self.paths_dict['dazer_path']           = '/home/vital/workspace/dazer/'
            self.paths_dict['stellar_data_folder']  = '/home/vital/Starlight/'    
            self.paths_dict['lick_indexes_folder']  = '/home/vital/workspace/dazer/working_examples/'
        
        #Paths for windows:
        elif name == 'nt':
            self.paths_dict['inference_folder']     = 'D:/Inference_data/'
            self.paths_dict['nebular_data_folder']  = 'E:/Cloud Storage/Dropbox/Astrophysics/Lore/NebularContinuum/'
            self.paths_dict['Hydrogen_CollCoeff']   = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Neutral_Hydrogen_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_CollCoeff']     = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Neutral_Helium_Collisional_Correction_coef.txt'
            self.paths_dict['Helium_OpticalDepth']  = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/Helium_OpticalDepthFunction_Coefficients.txt'
            self.paths_dict['stellar_data_folder']  = 'E:/Cloud Storage/Dropbox/Astrophysics/Tools/Starlight/'  
            self.paths_dict['dazer_path']           = 'C:/Users/lativ/git/dazer/'
            self.paths_dict['lick_indexes_folder']  = 'C:/Users/lativ/git/dazer/working_examples/'
            self.paths_dict['lines_data_file']      = 'C:/Users/lativ/git/dazer/bin/lib/Astro_Libraries/lines_data.xlsx'
            
        #Lines labels and wavelength to use in pyneb
        self.lines_df = read_excel(self.paths_dict['lines_data_file'], sheetname=0, header=0, index_col=0)
        self.lines_df['pyneb_code'] = self.lines_df['pyneb_code'].astype(str)
        idx_numeric_pynebCode = ~self.lines_df['pyneb_code'].str.contains('A+')
        self.lines_df.loc[idx_numeric_pynebCode,'pyneb_code'] = self.lines_df.loc[idx_numeric_pynebCode,'pyneb_code'].astype(int)
     
    def import_table_data(self, Address, Columns):
        
        Imported_Array  = genfromtxt(Address, dtype=float, usecols = Columns, skip_header = 2).T
        Datarray = Imported_Array[:,~isnan(Imported_Array).all(0)]
                
        return Datarray
     
    def import_coll_coeff_table(self, HydrogenLines, HeliumLines):
        
        Data_dict = OrderedDict()
                
        for i in range(len(HydrogenLines)):
            
            Line            = HydrogenLines[i]
            Data_Columns    = [0 + 3*i, 1 + 3*i, 2 + 3*i]
            Data_dict[Line] = self.import_table_data(self.paths_dict['Hydrogen_CollCoeff'], Data_Columns)
                    
        return Data_dict
     
    def import_optical_depth_coeff_table(self, HeliumLines):
        
        Data_dict = OrderedDict()
        
        for i in range(len(HeliumLines)):
            
            Line = HeliumLines[i]
            
            if Line != self.H13889A_label:
                Data_dict[Line] = loadtxt(self.paths_dict['Helium_OpticalDepth'], dtype = float, skiprows = 2, usecols = (i,))
                
            else:
                Data_dict[self.He3889_label] = loadtxt(self.paths_dict['Helium_OpticalDepth'], dtype = float, skiprows = 2, usecols = (i,))
            
        return Data_dict
    
    def ready_lines_data(self, ion, lines_labels):
        
        idx_lines = (self.lines_df.index.isin(lines_labels))
        
        self.obj_data[ion + '_labels']     = array(lines_labels)
        self.obj_data[ion + '_wave']       = self.lines_df[idx_lines].wavelength.values
        self.obj_data[ion + '_pyneb_code'] = self.lines_df[idx_lines].pyneb_code.values  
        
        return
    
    def load_synthetic_data(self, model):
        
        #Dictionary with data from the object
        self.obj_data = dict()
        
        #Physical parameters from the object
        self.obj_data['z_star']         = 0.00
        self.obj_data['sigma_star']     = 1.254
        self.obj_data['Av_star']        = 0.754            

        self.obj_data['sigma_gas']      = 2.255
        
        self.obj_data['Hbeta_Flux']     = 1e4 
        
        self.obj_data['n_e']            = 120.0
        self.obj_data['tau']            = 1.0
        self.obj_data['T_low']          = 15500.0
        self.obj_data['T_high']         = (1.0807 * self.obj_data['T_low']/10000.0 - 0.0846) * 10000.0
        self.obj_data['cHbeta']         = 0.1
        self.obj_data['xi']             = 1.0  
        self.obj_data['TSIII']          = 15600.0
        self.obj_data['TSIII_error']    = 700
        self.obj_data['nSII']           = 130.0
        self.obj_data['nSII_error']     = 50.0

        self.obj_data['He1_abund']      = 0.085
        self.obj_data['He2_abund']      = 0.0

        self.obj_data['O2_abund']       = 0.00025
        self.obj_data['O3_abund']       = 0.00075

        self.obj_data['N2_abund']       = 0.00035

        self.obj_data['S2_abund']       = 0.000015
        self.obj_data['S3_abund']       = 0.000055
        
        self.obj_data['Ar3_abund']      = 0.00065
        self.obj_data['Ar4_abund']      = 0.00012
        
        print 'La direccion', self.paths_dict['lick_indexes_folder']
        self.obj_data['lick_idcs_df']   = read_csv(self.paths_dict['lick_indexes_folder'] + 'synth_lick_indeces.txt', delim_whitespace = True,\
                                          header = 0, index_col = 0, comment='L') #Dirty trick to avoid the Line_label row
        
        H1_labels = ['H1_4102A', 'H1_4341A', 'H1_6563A']
        self.ready_lines_data('H1', H1_labels)  
        
        #He1_labels = ['He1_3889A',   'He1_4026A',  'He1_4471A',  'He1_5876A', 'He1_6678A',   'He1_7065A',    'He1_10830A']
        He1_labels = ['He1_4026A',  'He1_4471A',  'He1_5876A', 'He1_6678A']                                                      
        self.ready_lines_data('He1', He1_labels)  
        
        S2_labels = ['S2_6716A', 'S2_6731A']                        
        self.ready_lines_data('S2', S2_labels)  

        S3_labels = ['S3_6312A', 'S3_9069A', 'S3_9531A']                       
        self.ready_lines_data('S3', S3_labels)
        
        O2_labels = ['O2_3726A', 'O2_3729A']
        self.ready_lines_data('O2', O2_labels)
        
        O3_labels = ['O3_4363A', 'O3_4959A', 'O3_5007A']
        self.ready_lines_data('O3', O3_labels) 
 
        N2_labels = ['N2_6548A', 'N2_6584A']
        self.ready_lines_data('N2', N2_labels)  

        Ar3_labels = ['Ar3_7136A']
        self.ready_lines_data('Ar3', Ar3_labels)               
        
        Ar4_labels = ['Ar4_4740A']
        self.ready_lines_data('Ar4', Ar4_labels)    
              
#         self.obj_data['He1_labels']     = ['He1_3889A',         'He1_4026A',    'He1_4471A',    'He1_5876A', 'He1_6678A',   'He1_7065A',    'He1_10830A']
#         self.obj_data['He1_wave']       = array([ 3889.0,       4026.0,         4471.0,         5876.0,      6678.0,        7065.0,         10830.0])
#         self.obj_data['He1_pyneb_code'] = array(['3889.0',      '4026.0',       '4471.0',       '5876.0',    '6678.0',      '7065.0',       '10830.0'])
# 
#         self.obj_data['S2_labels']      = ['S2_6716A', 'S2_6731A']
#         self.obj_data['S2_wave']        = array([6716.44, 6730.81])
#         self.obj_data['S2_pyneb_code']  = array([6716, 6730])
#         
#         self.obj_data['S3_labels']      = ['S3_6312A', 'S3_9069A', 'S3_9531A']
#         self.obj_data['S3_wave']        = array([6312.06, 9068.6, 9531.1])
#         self.obj_data['S3_pyneb_code']  = array([6312, 9069, 9531])
# 
#         self.obj_data['O2_labels']      = ['O2_3726A', 'O2_3729A']
#         self.obj_data['O2_wave']        = array([3726.032, 3728.815])
#         self.obj_data['O2_pyneb_code']  = array([3726, 3729])
#         
#         self.obj_data['O3_labels']      = ['O3_4363A', 'O3_4959A', 'O3_5007A']
#         self.obj_data['O3_wave']        = array([4363.21, 4958.911, 5006.843])
#         self.obj_data['O3_pyneb_code']  = array([4363, 4959, 5007])
# 
#         self.obj_data['N2_labels']      = ['N2_6548A', 'N2_6584A']
#         self.obj_data['N2_wave']        = array([6548.05, 6583.46])
#         self.obj_data['N2_pyneb_code']  = array([6548, 6584])
# 
#         self.obj_data['Ar3_labels']     = ['Ar3_7136A']
#         self.obj_data['Ar3_wave']       = array([7135.79])
#         self.obj_data['Ar3_pyneb_code'] = array([7136])
#         
#         self.obj_data['Ar4_labels']     = ['Ar4_4740A']
#         self.obj_data['Ar4_wave']       = array([4740])
#         self.obj_data['Ar4_pyneb_code'] = array([4740])

        #Generate synthetic observation using default values        
        self.obj_data['mask_stellar'] = OrderedDict()
        self.obj_data['mask_stellar']['He1_4026A']  = (4019, 4033)
        self.obj_data['mask_stellar']['He1_4471A']  = (4463, 4480)
        self.obj_data['mask_stellar']['He1_5876A']  = (5867, 5885)
        self.obj_data['mask_stellar']['He1_6678A']  = (6667, 6687)
        self.obj_data['mask_stellar']['H1_delta']   = (4090,4114)
        self.obj_data['mask_stellar']['H1_gamma']   = (4329,4353)
        self.obj_data['mask_stellar']['H1_beta']    = (4840,4880)
        self.obj_data['mask_stellar']['H1_alpha']   = (6551,6575)
                
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
            
class Continua_FluxCalculation(ssp_fitter):

    def __init__(self):
        
        ssp_fitter.__init__(self)
        
class Nebular_FluxCalculation(NebularContinuumCalculator):
    
    def __init__(self):
        
        NebularContinuumCalculator.__init__(self)

    def nebular_Cont(self, wave_obs, z, cHbeta, Te, He1_abund, He2_abund, Halpha_Flux):
        
        wave_obs_rest   = wave_obs / (1.0 + z)
         
        neb_gCont       = self.calculate_neb_gCont(wave_obs_rest, Te, He1_abund, He2_abund)
         
        neb_int_norm    = self.gCont_calibration(wave_obs_rest, Te, Halpha_Flux, neb_gCont)
                     
        neb_xX          = self.reddening_Xx(wave_obs_rest, 'G03_average', 3.4)
        flambda_neb     = neb_xX/self.Hbeta_xX - 1.0
         
        neb_flux_norm   = neb_int_norm * power(10, -1 * flambda_neb * cHbeta)        
        
        return neb_flux_norm

    def calculate_nebular_SED(self, wave_obs, z, cHbeta, Te, He1_abund, He2_abund, Halpha_Flux):
        
        wave_obs_rest   = wave_obs / (1.0 + z)
         
        neb_gCont       = self.calculate_neb_gCont(wave_obs_rest, Te, He1_abund, He2_abund)
         
        neb_int_norm    = self.gCont_calibration(wave_obs_rest, Te, Halpha_Flux, neb_gCont)
                     
        neb_xX          = self.reddening_Xx(wave_obs_rest, 'G03_average', 3.4)
        flambda_neb     = neb_xX/self.Hbeta_xX - 1.0
         
        neb_flux_norm   = neb_int_norm * power(10, -1 * flambda_neb * cHbeta)   
                
        nebular_SED = {}
        nebular_SED['neb_gCont']        = neb_gCont        
        nebular_SED['neb_int_norm']     = neb_int_norm                    
        nebular_SED['neb_xX']           = neb_xX
        nebular_SED['flambda_neb']      = flambda_neb        
        nebular_SED['neb_flux_norm']    = neb_flux_norm        
        
        return nebular_SED
                   
class Recombination_FluxCalibration(LineMesurer_v2):
    
    def __init__(self):
    
        #Define indexes and labels to speed up the code
        LineMesurer_v2.__init__(self,  self.paths_dict['dazer_path'] + 'format/', 'DZT_LineLog_Headers.dz')
          
        #Define indexes and labels to speed up the code
        self.Hbeta_label            = 'H1_4861A'
        self.Hbeta_wave             = 4862.683
        self.Hbeta_pynebCode        = '4_2'
        self.H13889A_label          = 'H1_3889A'          
        self.He3889_label           = 'He1_3889A'
        self.He3889_Check           = None
        
        #Set up the right emissivities
        #atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
        
        #Atoms to fit the data
        self.recombDict        = {}
        self.recombDict['H1']  = RecAtom('H', 1)
        self.recombDict['He1'] = RecAtom('He', 1)
        self.recombDict['He2'] = RecAtom('He', 2)

        #Make sure we are using the right helium emissivities
        print '--Helium emissivities: '
        self.recombDict['He1'].printSources()

        #Import collisional coefficients table #We add the Hbeta to get its coefficients
        posHydrogen_Lines       = ['H1_4102A',  'H1_4341A', 'H1_6563A']
        self.Coef_Kalpha_dict   = self.import_coll_coeff_table(posHydrogen_Lines + [self.Hbeta_label], None)

        #Import Optical depth function
        posHelium_Lines         = ['He1_3889A',  'He1_4026A',  'He1_4387A',    'He1_4471A',    'He1_4686A',    'He1_4714A',    'He1_4922A',   'He1_5876A',    'He1_6678A',   'He1_7065A',    'He1_7281A',    'He1_10830A']
        self.Coef_ftau_dict     = self.import_optical_depth_coeff_table(posHelium_Lines) 

    def H_alphaCalc(self, Te, ne, xi, cHbeta, idx_Halpha = 2):
        
        t4 = Te / 10000.0
        
        Hbeta_Kalpha        = self.Kalpha_Ratio_H(T_4 = t4, H_label = self.Hbeta_label)
        Hbeta_emis          = self.H1.getEmissivity(Te, ne, label = self.Hbeta_pynebCode)

        Halpha_Kalpha       = self.Kalpha_Ratio_H(T_4 = t4, H_label = self.obj_data['H_labels'][idx_Halpha])
        Emissivity_Halpha   = self.H1.getEmissivity(Te, ne, label = self.obj_data['H_pyneb_code'][idx_Halpha])
        
        emis_module         = Emissivity_Halpha / Hbeta_emis
        
        CR_Module           = (1.0 + 0.0001* xi * Halpha_Kalpha) / (1.0 + 0.0001* xi * Hbeta_Kalpha)
        
        f_module            = power(10, -1 * self.data_dic['flambda_H_vector'][idx_Halpha] * cHbeta)
        
        Halpha_Flux         = emis_module * CR_Module * f_module * self.Hbeta_Flux

        return Halpha_Flux

    def H_theoFlux(self, Te, ne, xi, cHbeta):
        
        #Hydrogen calculation parameters
        self.calculate_H_params(Te, ne)
        
        #Hbeta parameters
        Hbeta_emis      = self.H1.getEmissivity(Te, ne, label = self.Hbeta_pynebCode)
        Hbeta_Kalpha    = self.Kalpha_Ratio_H(T_4 = Te/10000.0, H_label = self.Hbeta_label)
        
        #Calculate the emissivities for each of the lines for the given temperature and density
        emis_module     = self.data_dic['Emissivity_H_vector'] / Hbeta_emis
                
        #Calculate the Collisional excitation fraction
        CR_Module       = (1.0 + 0.0001* xi * self.data_dic['Kalpha_vector']) / (1.0 + 0.0001* xi * Hbeta_Kalpha)
                
        #Calculate the reddening component
        f_module        = power(10, -1 * self.data_dic['flambda_H_vector'] * cHbeta)

        #Calculate theoretical Hydrogen flux for each emission line
        H_Flux          = emis_module * CR_Module * f_module
                 
        return H_Flux

    def He_theoFlux(self, Te, ne, tau, xi, cHbeta, He1_abund):
        
        #Helium calculation parameters
        self.calculate_He_params(Te, ne, tau)
        
        #Hbea kalpha
        Hbeta_emis      = self.H1.getEmissivity(Te, ne, label = self.Hbeta_pynebCode)
        Hbeta_Kalpha    = self.Kalpha_Ratio_H(T_4 = Te/10000.0, H_label = self.Hbeta_label)
        
        #Calculate the emissivities for each of the lines for the given temperature and density
        emis_module = self.data_dic['Emissivity_He_vector'] / Hbeta_emis
                
        #Calculate the collisional excitation fraction
        CR_Module   = 1 / (1.0 + 0.0001* xi * Hbeta_Kalpha)

        #Calculate the reddening component
        f_module    = power(10, -1 *  self.data_dic['flambda_He_vector'] * cHbeta)
             
        #Calculate theoretical helium flux for the line
        He_Flux     = He1_abund * emis_module * self.data_dic['ftau_He_vector'] * CR_Module * f_module

        return He_Flux
                                                                    
    def calculate_H_params(self, Te, ne):    
        
        #Calculate in advance the T_4 parameter
        T_4 = Te / 10000.0
        
        #Calculate the hbeta parameters
        #self.Hbeta_Kalpha   = self.Kalpha_Ratio_H(T_4 = T_4, H_label = self.Hbeta_label)
        #self.Hbeta_emis     = self.H1.getEmissivity(Te, ne, label = self.Hbeta_pynebCode)

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
            self.data_dic['Emissivity_He_vector'][i]    = self.He1.getEmissivity(Te, ne, label = self.obj_data['He1_pyneb_code'][i])
            self.data_dic['ftau_He_vector'][i]          = self.OpticalDepth_He(tau = tau, T_4 = T_4, ne = ne, He_label = self.obj_data['He1_labels'][i])    
                    
        return
                          
    def Kalpha_Ratio_H(self, T_4, H_label):
                      
        K_alpha_Ratio   = sum(self.Coef_Kalpha_dict[H_label][0] * exp(self.Coef_Kalpha_dict[H_label][1] / T_4) * power(T_4, self.Coef_Kalpha_dict[H_label][2]))
        
        return K_alpha_Ratio

    def OpticalDepth_He(self, tau, T_4, ne, He_label):
                        
        f_tau   = 1 + (tau/2) * (self.Coef_ftau_dict[He_label][0] + (self.Coef_ftau_dict[He_label][1] + self.Coef_ftau_dict[He_label][2]*ne + self.Coef_ftau_dict[He_label][3]*ne*ne) * T_4)

        return f_tau

    def calculate_recomb_fluxes(self, Thigh, ne, cHbeta, xi, tau, lines_abund_dict, lines_waves, lines_ions, lines_flambda):
                
        #Emissivity for the two temperature layers
        Te = Thigh
        t4 = Te / 10000.0

        #Hbeta parameters
        Emis_Hbeta                      = self.recombDict['H1'].getEmissivity(Thigh, ne, label = self.Hbeta_pynebCode) 
        Hbeta_Kalpha                    = self.Kalpha_Ratio_H(T_4 = t4, H_label = self.Hbeta_label)
        cr_Hbeta                        = (1.0 + 0.0001* xi * Hbeta_Kalpha)
        
        #Loop through the lines to calculate their emissivities
        lines_emis_vector               = empty(self.n_recombLines)
        lines_abs_vector                = empty(self.n_recombLines)
        
        for i in self.range_recombLines:
            
            ion = lines_ions[i]
            
            if ion == 'H1':
                Kalpha_i                = self.Kalpha_Ratio_H(T_4 = t4, H_label = self.Recomb_labels[i])                
                emisRatio_i             = self.recombDict[ion].getEmissivity(Te, ne, wave = self.Recomb_pynebCode[i]) / Emis_Hbeta
                cr_i                    = (1.0 + 0.0001* xi * Kalpha_i) / cr_Hbeta
                lines_emis_vector[i]    = emisRatio_i * cr_i
                
            elif ion == 'He1':
                ftau_i                  = self.OpticalDepth_He(tau = tau, T_4 = t4, ne = ne, He_label = self.Recomb_labels[i])
                emisRatio_i             = self.recombDict[ion].getEmissivity(Te, ne, wave = self.Recomb_pynebCode[i]) / Emis_Hbeta
                lines_emis_vector[i]    = lines_abund_dict[ion] * emisRatio_i * ftau_i * (1.0/cr_Hbeta)
                
            elif ion == 'He2':
                emisRatio_i             = self.recombDict['ion'].getEmissivity(Te, ne, wave = self.Recomb_pynebCode[i]) / Emis_Hbeta
                lines_emis_vector[i]    = lines_abund_dict[ion] * emisRatio_i
                        
        #Reddening factor for all the lines
        f_module = power(10, -1 * lines_flambda * cHbeta)
                
        #Metals flux
        recomb_fluxes = lines_emis_vector * f_module
                
        return recomb_fluxes
   
class Collisional_FluxCalibration(Import_model_data):
    
    def __init__(self):
                
        #Set atomic data
        atomicData.setDataFile('s_iii_coll_HRS12.dat')
        
        #Atoms to fit the data
        self.ionDict        = {}
        self.ionDict['H']   = RecAtom('H', 1)
        self.ionDict['He1'] = RecAtom('He', 1)
        self.ionDict['He2'] = RecAtom('He', 2)
        self.ionDict['S2']  = Atom('S', 2)
        self.ionDict['S3']  = Atom('S', 3) 
        self.ionDict['Ar3'] = Atom('Ar', 3)
        self.ionDict['Ar4'] = Atom('Ar', 4) 
        self.ionDict['O2']  = Atom('O', 2)
        self.ionDict['O3']  = Atom('O', 3)
        self.ionDict['N2']  = Atom('N', 2)
        
        self.highTemp_ions = array(['He', 'He1', 'He2', 'O3', 'Ar4'])
                
    def coll_Flux(self, ion, Te, ne, cHbeta, ionic_abund):
        
        #Calculate Hbeta emissivity for this ion temperatre
        Hbeta_emis      = self.H1.getEmissivity(Te, ne, label = self.Hbeta_pynebCode)

        #Ions emissivities by Hbeta emissivity
        emis_module     = self.metal_emis(ion, Te, ne) / Hbeta_emis
                
        #Reddening factor
        f_module        = power(10, -1 * self.data_dic['flambda_{}_vector'.format(ion)] * cHbeta)

        #Metals flux
        metal_flux      = ionic_abund * emis_module * f_module
        
        return metal_flux
    
    def metal_emis(self, ion, Te, ne):
                
        ion_emis    = self.ionDict[ion]
        pyneb_code  = '{}_pyneb_code'.format(ion)
        vector_emis =  self.data_dic['Emissivity_{}_vector'.format(ion)]
        
        for i in self.obj_data['n' + ion + '_range']:    
            vector_emis[i] = ion_emis.getEmissivity(Te, ne, wave = self.obj_data[pyneb_code][i])
            
        return vector_emis

    def calculate_colExcit_flux(self, Tlow, Thigh, ne, cHbeta, lines_abund_dict, lines_waves, lines_ions, lines_flambda):
        
        #Emissivity for the two temperature layers
        Emis_Hbeta_low  = self.recombDict['H1'].getEmissivity(Tlow, ne, label = self.Hbeta_pynebCode) 
        Emis_Hbeta_high = self.recombDict['H1'].getEmissivity(Thigh, ne, label = self.Hbeta_pynebCode)
        
        #Loop through the lines to calculate their emissivities
        lines_emis_vector = empty(self.n_colExcLines)
        for i in self.range_colExcLines:
            
            ion = lines_ions[i]
            
            #Set the right temperature for the ion
            if ion in self.highTemp_ions:
                Te = Thigh
                emisHbeta = Emis_Hbeta_high
            else:
                Te = Tlow
                emisHbeta = Emis_Hbeta_low
                
            #Ions emissivities by Hbeta emissivity
            ion_emis = self.ionDict[ion].getEmissivity(Te, ne, wave = self.lines_pynebCode[i])
            lines_emis_vector[i] = lines_abund_dict[ion] * ion_emis / emisHbeta

        #Reddening factor for all the lines
        f_module = power(10, -1 * lines_flambda * cHbeta)
                
        #Metals flux
        colExcit_fluxes = lines_emis_vector * f_module
                
        return colExcit_fluxes

class Inference_AbundanceModel(Import_model_data, Collisional_FluxCalibration, Recombination_FluxCalibration, Continua_FluxCalculation, Nebular_FluxCalculation):

    def __init__(self):

        #Import tools to load data
        Import_model_data.__init__(self)
        Continua_FluxCalculation.__init__(self)
        Collisional_FluxCalibration.__init__(self)
        Recombination_FluxCalibration.__init__(self)
        Nebular_FluxCalculation.__init__(self)

        #Self
        self.Rv_model = 3.4
        self.reddedning_curve_model = 'G03_average'
        self.abund_iter_dict = {}
        
        #Speed up calculation
        self.Hbeta_xX = self.reddening_Xx(array([self.Hbeta_wave]), self.reddedning_curve_model, self.Rv_model)[0]

    def measure_lines(self, line_labels, wave, flux, calibration_factor):
        
        n_lines     = line_labels.shape[0]
        range_lines = arange(n_lines)
         
        intg_flux   = empty(n_lines)
                
        for i in range_lines:
            
            self.Current_Label      = line_labels[i]
            self.Current_Ion        = self.obj_data['recombLine_ions'][i]
            self.Current_TheoLoc    = self.obj_data['recombLine_waves'][i]
            selections              = self.obj_data['lick_idcs_df'].loc[self.Current_Label][3:9].values
            line_data               = self.measure_line(wave, flux, selections, None, Measuring_Method = 'lmfit', store_data = False)
            intg_flux[i]            = line_data['flux_intg']
        
        return intg_flux  * calibration_factor   

    def calculate_simObservation(self, model, obs_lines, verbose = True):
          
        #Get physical parameters from model HII region
        self.load_synthetic_data(model = model)
        
        #Physical parameters for the nebular continua calculation
        self.load_neb_constants(self.paths_dict['nebular_data_folder'])
        
        #Recombination features                    
        self.recomb_fluxes          = self.synth_recomb_emission(obs_lines)
        self.obs_recomb_err         = self.recomb_fluxes * 0.02
                                
        #Collisinal excited features        
        self.obs_metal_fluxes       = self.synth_collisional_emission(obs_lines)
        self.obs_metal_Error        = self.obs_metal_fluxes * 0.02
                
        #-------Prepare stellar continua
        #Load stellar libraries
        default_Starlight_file      = self.paths_dict['stellar_data_folder'] + 'Dani_Bases_Extra_short.txt'
        default_Starlight_folder    = self.paths_dict['stellar_data_folder'] + 'Bases/'
        default_Starlight_coeffs    = self.paths_dict['stellar_data_folder'] + 'Bases/coeffs_sync_short.txt'  
        
        ssp_lib_dict                = self.load_stellar_bases('starlight', default_Starlight_folder, default_Starlight_file,\
                                                resample_int=1, resample_range = (3600, 6900), norm_interval = (5100,5150))
                
        self.stellar_SED            = self.calculate_synthStellarSED(self.obj_data['Av_star'], self.obj_data['z_star'], self.obj_data['sigma_star'],\
                                                default_Starlight_coeffs, ssp_lib_dict, (4000, 6900), mask_dict= self.obj_data['mask_stellar'])
        
        #Initial guess number of bases
        idx_populations             = self.stellar_SED['bases_coeff'] > 0.01
        self.population_limitss     = vstack((self.stellar_SED['bases_coeff'][where(idx_populations)]*0.95, self.stellar_SED['bases_coeff'][where(idx_populations)]*1.05)).T
        self.nBases                 = len(where(idx_populations)[0]) 
        self.range_bases            = arange(self.nBases)       
        columns_to_delete           = where(~idx_populations)
        
        #-------Prepare nebular continua
        self.idx_Halpha             = (self.obj_data['recombLine_labes'] == 'H1_6563A')
        
        self.Halpha_norm            = self.recomb_fluxes[self.idx_Halpha][0] * self.obj_data['Hbeta_Flux'] / self.stellar_SED['normFlux_stellar']
        
        self.nebular_SED            = self.calculate_nebular_SED(self.stellar_SED['stellar_wave_resam'], self.obj_data['z_star'], self.obj_data['cHbeta'],\
                                                       self.obj_data['T_low'], self.obj_data['He1_abund'], self.obj_data['He2_abund'], self.Halpha_norm)
                
        #Redshift limits for the object
        z_max_ssp                   = (self.stellar_SED['stellar_wave_resam'][0] / ssp_lib_dict['basesWave_resam'][0]) - 1.0
        z_min_ssp                   = (self.stellar_SED['stellar_wave_resam'][-1] / ssp_lib_dict['basesWave_resam'][-1]) - 1.0
        self.z_max_ssp_limit        = round(z_max_ssp - 0.001, 3)
        self.z_min_ssp_limit        = z_min_ssp

        #-------Emission sed from the object
        H_He_fluxes                 = self.recomb_fluxes * self.obj_data['Hbeta_Flux']  / self.stellar_SED['normFlux_stellar']
        self.emission_Spec          = self.calc_emis_spectrum(self.stellar_SED['stellar_wave_resam'], self.obj_data['recombLine_waves'], H_He_fluxes,\
                                                       self.obj_data['recombLine_waves'], self.obj_data['sigma_gas'], self.obj_data['z_star'])

        #Add the components to be considered in the analysis
        #obs_flux_norm              = self.stellar_SED['stellar_flux_norm']
        obs_flux_norm               = self.stellar_SED['stellar_flux_norm'] + self.nebular_SED['neb_flux_norm'] + self.emission_Spec
        
        

        
        #Store the data in the object dictionary 
        self.obj_data['normFlux_obs']           = self.stellar_SED['normFlux_stellar']
        self.obj_data['obs_wave_resam']         = self.stellar_SED['stellar_wave_resam'] 
        self.obj_data['obs_flux_norm']          = obs_flux_norm
        self.obj_data['obs_flux_norm_masked']   = self.obj_data['obs_flux_norm'] * self.stellar_SED['int_mask']
        self.obj_data['basesWave_resam']        = ssp_lib_dict['basesWave_resam'] 
        self.obj_data['bases_flux_norm']        = delete(ssp_lib_dict['bases_flux_norm'], columns_to_delete, axis=0)
        self.obj_data['int_mask']               = self.stellar_SED['int_mask']
        self.obj_data['obs_fluxEr_norm']        = self.stellar_SED['stellar_fluxEr_norm']
    
    
        self.obs_recomb_fluxes = self.measure_lines(self.Recomb_labels, self.obj_data['obs_wave_resam'], self.obj_data['obs_flux_norm'], calibration_factor = (self.obj_data['normFlux_obs'] / self.obj_data['Hbeta_Flux']))
        print 'Troleo mayor', self.obs_recomb_fluxes, type(self.obs_recomb_fluxes)
    
        #Print the model values
        if verbose:
            print '\nInput Parameters:'
            print '-He1_abund', self.obj_data['He1_abund']
            print '-n_e',       self.obj_data['n_e']
            print '-T_low',     self.obj_data['T_low']
            print '-T_high',    self.obj_data['T_high']
            print '-cHbeta',    self.obj_data['cHbeta']
            print '-tau',       self.obj_data['tau']
            print '-xi',        self.obj_data['xi']
            print '-TOIII',     self.obj_data['TSIII'],'+/-',self.obj_data['TSIII_error']
            print '-nSII',      self.obj_data['nSII'],'+/-',self.obj_data['nSII_error']
            print '\n-S2_abund',self.obj_data['S2_abund']
            print '-S3_abund',  self.obj_data['S3_abund']
            print '-O2_abund',  self.obj_data['O2_abund']
            print '-O3_abund',  self.obj_data['O3_abund']
            
            print '-N2_abund',  self.obj_data['N2_abund']
            print '-Ar3_abund',  self.obj_data['Ar3_abund']
            print '-Ar4_abund',  self.obj_data['Ar4_abund']
            
            print '\n-z_star',  self.obj_data['z_star']
            print '-sigma_star',self.obj_data['sigma_star']
            print '-Av_star',   self.obj_data['Av_star']
            
            print '\n-Wavelength ranges:'
            print '--Observation:', self.obj_data['obs_wave_resam'][0], self.obj_data['obs_wave_resam'][-1]
            print '--Bases:', ssp_lib_dict['basesWave_resam'][0], ssp_lib_dict['basesWave_resam'][-1]
            print '--z min:', self.z_min_ssp_limit
            print '--z max:', self.z_max_ssp_limit
            print '--z true:', self.obj_data['z_star']
              
        return

    def synth_collisional_emission(self, obs_lines):
        
        #Loop through the ions and read the lines (this is due to the type I format the synth data)
        lines_labes, lines_waves, lines_ions, lines_abund, lines_pynebCode = [], [], [], [], []
        abund_dict = {}
        for ion in obs_lines:
            if ion not in ['H1', 'He1', 'He2']:
                lines_ions                  += [ion] * len(self.obj_data[ion + '_wave'])
                lines_labes                 += list(self.obj_data[ion + '_labels'])
                lines_waves                 += list(self.obj_data[ion + '_wave'])
                lines_abund                 += [self.obj_data[ion + '_abund']] * len(self.obj_data[ion + '_wave']) 
                lines_pynebCode             += list(self.obj_data[ion + '_pyneb_code'])
                abund_dict[ion]             = self.obj_data[ion + '_abund']
                self.abund_iter_dict[ion]   = 0.0
                
        #Convert lists to numpy arrays
        lines_labes, lines_waves, lines_ions = array(lines_labes), array(lines_waves), array(lines_ions)
        lines_abund, lines_pynebCode = array(lines_abund), array(lines_pynebCode)
        
        #Sorting by increasing wavelength
        idx_sort        = argsort(lines_waves)
        lines_labes     = lines_labes[idx_sort]
        lines_waves     = lines_waves[idx_sort]
        lines_ions      = lines_ions[idx_sort]
        lines_abund     = lines_abund[idx_sort]
        self.lines_pynebCode = lines_pynebCode[idx_sort]

        #Calculate xX array
        lines_labes_xX  = self.reddening_Xx(lines_waves, self.reddedning_curve_model, self.Rv_model)
        lines_flambda   = lines_labes_xX/self.Hbeta_xX - 1.0
        
        #Parameters to speed the code
        self.n_colExcLines = len(lines_waves)
        self.range_colExcLines = arange(self.n_colExcLines)
        
        colExcit_fluxes = self.calculate_colExcit_flux(self.obj_data['T_low'], self.obj_data['T_high'], self.obj_data['n_e'], self.obj_data['cHbeta'], abund_dict, lines_waves, lines_ions, lines_flambda)
        
        self.obj_data['colLine_ions'] = lines_ions
        self.obj_data['colLine_labes'] = lines_labes
        self.obj_data['colLine_waves'] = lines_waves
        self.obj_data['colLine_pynebCode'] = lines_pynebCode  
        self.obj_data['colLine_flambda'] = lines_flambda
        
        return colExcit_fluxes

    def synth_recomb_emission(self, obs_lines):
        
        #Loop through the ions and read the lines (this is due to the type I format the synth data)
        lines_labes, lines_waves, lines_ions, lines_abund, lines_pynebCode = [], [], [], [], []
        abund_dict = {}
        
        for ion in obs_lines:
            if ion in ['H1', 'He1', 'He2']:
                lines_ions          += [ion] * len(self.obj_data[ion + '_wave'])
                lines_labes         += list(self.obj_data[ion + '_labels'])
                lines_waves         += list(self.obj_data[ion + '_wave'])
                lines_pynebCode     += list(self.obj_data[ion + '_pyneb_code'])
                
                if ion == 'H1':
                    lines_abund     += [1.0] * len(self.obj_data[ion + '_wave'])
                    abund_dict[ion] = 1.0                    
                else:
                    lines_abund     += [self.obj_data[ion + '_abund']] * len(self.obj_data[ion + '_wave'])
                    abund_dict[ion] = self.obj_data[ion + '_abund']
                    
                self.abund_iter_dict[ion] = 0.0    

        #Convert lists to numpy arrays
        lines_labes, lines_waves, lines_ions = array(lines_labes), array(lines_waves), array(lines_ions)
        lines_abund, lines_pynebCode = array(lines_abund), array(lines_pynebCode)

        #Sorting by increasing wavelength
        idx_sort                = argsort(lines_waves)
        lines_labes             = lines_labes[idx_sort]
        lines_waves             = lines_waves[idx_sort]
        lines_ions              = lines_ions[idx_sort]
        lines_abund             = lines_abund[idx_sort]
        self.Recomb_pynebCode   = lines_pynebCode[idx_sort]
        self.Recomb_labels      = lines_labes
        
        #Calculate xX array
        lines_labes_xX  = self.reddening_Xx(lines_waves, self.reddedning_curve_model, self.Rv_model)
        lines_flambda   = lines_labes_xX/self.Hbeta_xX - 1.0

        #Parameters to speed the code
        self.n_recombLines = len(lines_waves)
        self.range_recombLines = arange(self.n_recombLines)

        recomb_fluxes = self.calculate_recomb_fluxes(self.obj_data['T_high'], self.obj_data['n_e'],\
                                                      self.obj_data['cHbeta'], self.obj_data['xi'], self.obj_data['tau'],\
                                                      abund_dict, lines_waves, lines_ions, lines_flambda)
        
        self.obj_data['recombLine_ions'] = lines_ions
        self.obj_data['recombLine_labes'] = lines_labes
        self.obj_data['recombLine_waves'] = lines_waves
        self.obj_data['recombLine_pynebCode'] = lines_pynebCode  
        self.obj_data['recombLine_flambda'] = lines_flambda
        
        return recomb_fluxes

    def calc_emis_spectrum(self, wavelength_range, lines_waves, lines_fluxes, lines_mu, lines_sigma, redshift):
                
        lines_mu            = lines_mu * (1 + redshift)
        
        A_lines             = lines_fluxes / (lines_sigma * self.sqrt2pi)
        
        wave_matrix         = wavelength_range * ones((lines_waves.shape[0], wavelength_range.shape[0]))
        
        emis_spec_matrix    = A_lines[:,None] * exp(-(wave_matrix-lines_mu[:,None])*(wave_matrix-lines_mu[:,None])/(2 * lines_sigma * lines_sigma))
        
        combined_spectrum   = emis_spec_matrix.sum(axis=0)
        
        return combined_spectrum
     
    def He_O_S_nebStellar_model(self):
        
        ne          =   pymc2.TruncatedNormal(  'ne',           self.obj_data['nSII'],      self.obj_data['nSII_error']**-2,    a = 50.0 ,   b = 1000.0)
        cHbeta      =   pymc2.TruncatedNormal(  'cHbeta',       0.15,                       0.05**-2,   a = 0.0,    b = 3.0)
        T_low       =   pymc2.TruncatedNormal(  'T_low',        self.obj_data['TSIII'],     self.obj_data['TSIII_error']**-2,    a = 7000.0 ,   b = 20000.0)
        
        S2_abund    =   pymc2.Uniform(          'S2_abund',     0.000001,                   0.001)
        S3_abund    =   pymc2.Uniform(          'S3_abund',     0.000001,                   0.001)
        O2_abund    =   pymc2.Uniform(          'O2_abund',     0.000001,                   0.001)
        O3_abund    =   pymc2.Uniform(          'O3_abund',     0.000001,                   0.001)                        
        N2_abund    =   pymc2.Uniform(          'N2_abund',     0.000001,                   0.001)
        Ar3_abund   =   pymc2.Uniform(          'Ar3_abund',    0.000001,                   0.001)                        
        Ar4_abund   =   pymc2.Uniform(          'Ar4_abund',    0.000001,                   0.001)                        
                
        He1_abund   =   pymc2.Uniform(          'He1_abund',    0.050,                      0.15)
        tau         =   pymc2.TruncatedNormal(  'tau',          0.75,                       0.5**-2,    a = 0.0,    b = 7.0)
        cHbeta      =   pymc2.TruncatedNormal(  'cHbeta',       0.15,                       0.05**-2,   a = 0.0,    b = 3.0)
        xi          =   pymc2.TruncatedNormal(  'xi',           1,                          200**-2,    a = 0.0,    b = 1000.0)
        T_He        =   pymc2.TruncatedNormal(  'T_He',         self.obj_data['TSIII'],     self.obj_data['TSIII_error']**-2,    a = 7000.0 ,   b = 20000.0, value=14500.0)
        
        #z_star     =   pymc2.Uniform(          'z_star',       self.z_min_ssp_limit,       self.z_max_ssp_limit)
        Av_star     =   pymc2.Uniform(          'Av_star',      0.0,                        2.00)
        sigma_star  =   pymc2.Uniform(          'sigma_star',   0.0,                        2.00)
        ssp_coefs   =  [pymc2.Uniform(          'ssp_coefs_%i' % i,   self.population_limitss[i][0],   self.population_limitss[i][1])   for i in self.range_bases]

        @pymc2.deterministic()
        def calc_Thigh(Te = T_low):
            return (1.0807 * Te/10000.0 - 0.0846) * 10000.0

        @pymc2.deterministic()
        def calc_abund_dict(He1_abund=He1_abund, S2_abund=S2_abund, S3_abund=S3_abund, O2_abund=O2_abund, O3_abund=O3_abund, N2_abund=N2_abund, Ar3_abund=Ar3_abund, Ar4_abund=Ar4_abund):
            
            self.abund_iter_dict['He1'] = He1_abund
            self.abund_iter_dict['S2'] = S2_abund
            self.abund_iter_dict['S3'] = S3_abund
            self.abund_iter_dict['O2'] = O2_abund
            self.abund_iter_dict['O3'] = O3_abund
            self.abund_iter_dict['N2'] = N2_abund
            self.abund_iter_dict['Ar3'] = Ar3_abund
            self.abund_iter_dict['Ar4'] = Ar4_abund
            
            return self.abund_iter_dict

        @pymc2.deterministic()
        def calc_colExcit_fluxes(abund_dict=calc_abund_dict, T_low=T_low, T_High=calc_Thigh, ne=ne, cHbeta=cHbeta):
              
            colExcit_fluxes = self.calculate_colExcit_flux(T_low, T_High, ne, cHbeta, abund_dict, self.obj_data['colLine_waves'], self.obj_data['colLine_ions'], self.obj_data['colLine_flambda'])
            
            return colExcit_fluxes
                       
        @pymc2.deterministic
        def calc_nebular_cont(z_star=self.obj_data['z_star'], cHbeta=self.obj_data['cHbeta'], Te=self.obj_data['T_low'], He1_abund=He1_abund, He2_abund=0.0, Halpha_Flux=self.Halpha_norm):
                   
            neb_flux_norm = self.nebular_Cont(self.obj_data['obs_wave_resam'], z_star, cHbeta, Te, He1_abund, He2_abund, Halpha_Flux)

            return neb_flux_norm
        
        @pymc2.deterministic
        def calc_continuum(z_star=self.obj_data['z_star'], Av_star=Av_star, sigma_star=sigma_star, ssp_coefs=ssp_coefs, nebular_flux=calc_nebular_cont):
                                    
            ssp_grid_i      = self.physical_SED_model(self.obj_data['basesWave_resam'], self.obj_data['obs_wave_resam'], self.obj_data['bases_flux_norm'],\
                                                       Av_star, z_star, sigma_star, 3.4)
            
            fit_continuum   = ssp_grid_i.dot(ssp_coefs) + nebular_flux
            
            return fit_continuum        
        
        @pymc2.deterministic
        def calc_recomb_fluxes(abund_dict=calc_abund_dict, T_He=T_He, ne=ne, cHbeta=cHbeta, xi=xi, tau=tau):
              
            recomb_fluxes = self.calculate_recomb_fluxes(T_He, ne, cHbeta, xi, tau, abund_dict,\
                                                          self.obj_data['recombLine_waves'], self.obj_data['recombLine_ions'], self.obj_data['recombLine_flambda'])
            
            return recomb_fluxes
        
        @pymc2.deterministic
        def calc_obs_recomb_fluxes(fit_continuum=calc_continuum, emis_recomb_fluses=calc_recomb_fluxes, z_star = self.obj_data['z_star']):
            
            H_He_fluxes         = emis_recomb_fluses * self.obj_data['Hbeta_Flux']  / self.stellar_SED['normFlux_stellar']

            emiss_spectrum      = self.calc_emis_spectrum(self.obj_data['obs_wave_resam'], self.obj_data['recombLine_waves'], H_He_fluxes,\
                                                       self.obj_data['recombLine_waves'], self.obj_data['sigma_gas'], z_star)                            
            
            sym_spectrum        = fit_continuum + emiss_spectrum
            
            lines_fluxes_obs    = self.measure_lines(self.Recomb_labels, self.obj_data['obs_wave_resam'], sym_spectrum, calibration_factor = (self.obj_data['normFlux_obs'] / self.obj_data['Hbeta_Flux']))
            
            
            return lines_fluxes_obs
            
        @pymc2.stochastic(observed=True) #Likelihood
        def likelihood_ssp(value = self.obj_data['obs_flux_norm_masked'], fit_continuum=calc_continuum, sigmaContinuum=self.obj_data['obs_fluxEr_norm']):
            calc_continuum_masked = fit_continuum * self.obj_data['int_mask']
            chi_F = sum(square(calc_continuum_masked - value) / square(sigmaContinuum))
            return - chi_F / 2
 
        @pymc2.stochastic(observed=True) #Likelihood
        def likelihood_recomb(value = self.obs_recomb_fluxes, H_He_TheoFlux = calc_obs_recomb_fluxes, sigmaLines = self.obs_recomb_err):
            chi_F = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
            return - chi_F / 2

        @pymc2.stochastic(observed=True) #Likelihood
        def likelihood_colExcited(value = self.obs_metal_fluxes, theo_metal_fluzes = calc_colExcit_fluxes, sigmaLines = self.obs_metal_Error):
            chi_F = sum(square(theo_metal_fluzes - value) / square(sigmaLines))
            return - chi_F / 2
 
        return locals()
    
    def stellar_bayes(self):
                
        #z_star     =   pymc2.Uniform(    'z_star',             self.z_min_ssp_limit,       self.z_max_ssp_limit)
        Av_star     =   pymc2.Uniform(    'Av_star',            0.0,                        5.00,                   value=self.obj_data['Av_star'])
        sigma_star  =   pymc2.Uniform(    'sigma_star',         0.0,                        5.00,                   value=self.obj_data['sigma_star'])
        ssp_coefs   =  [pymc2.Uniform(    'ssp_coefs_%i' % i,   self.population_limitss[i][0],   self.population_limitss[i][1]) for i in self.range_bases]
                
        @pymc2.deterministic
        def ssp_spectrum(z_star=self.obj_data['z_star'], Av_star=Av_star, sigma_star=sigma_star, ssp_coefs = ssp_coefs):#, nebular_flux=nebular_continua):
                                    
            ssp_grid_i          = self.physical_SED_model(self.obj_data['basesWave_resam'], self.obj_data['obs_wave_resam'], self.obj_data['bases_flux_norm'],\
                                    Av_star, z_star, sigma_star, 3.4)
            
            ssp_grid_i_masked   = (self.obj_data['int_mask'] * ssp_grid_i.T).T
            
            
            fitted_spectrum     = ssp_grid_i_masked.dot(ssp_coefs)
                        
            return fitted_spectrum        
        
        @pymc2.stochastic(observed=True) #Likelihood
        def likelihood_ssp(value = self.obj_data['obs_flux_norm_masked'], StellarCont_TheoFlux=ssp_spectrum, sigmaContinuum=self.obj_data['obs_fluxEr_norm']):
            chi_F = sum(square(StellarCont_TheoFlux - value) / square(sigmaContinuum))
            return - chi_F / 2

        return locals()
 
class Run_MCMC(Inference_AbundanceModel, ssp_fitter):

    def __init__(self):
        
        #Import parent classes
        Inference_AbundanceModel.__init__(self)
        
        self.pymc_stats_keys = ['mean','95% HPD interval','standard deviation','mc error','quantiles','n']
                    
    def select_inference_model(self, model):
        
        if model == '_stellar_tests':
            self.inf_dict = self.stellar_bayes()
        else:
            self.inf_dict = self.He_O_S_nebStellar_model()
                         
    def run_pymc2(self, db_address, iterations = 10000, variables_list = None, prefit = True):
                
        #Define MCMC model        
        self.MAP_Model = pymc2.MAP(self.inf_dict)

        #Prefit:
        if prefit is not False:
            fit_method = prefit if prefit is str else 'fmin_powell'
            print '\n--Starting {} prefit'.format(fit_method)
            start = timer()
            self.MAP_Model.fit(method = fit_method)
            end = timer()
            print 'prefit interval', (end - start) 

        #Print prefit data
        print 'Initial conditions:'
        self.display_run_data(self.MAP_Model, variables_list)
                 
        #Launch sample
        print '\nInitiating fit:'
        self.pymc2_M = pymc2.MCMC(self.MAP_Model.variables, db = 'pickle', dbname =  db_address)
        self.pymc2_M.sample(iter=iterations)
         
        #Save the output csv mean data
        if variables_list != None:
            print '--Saving results in csv'
            self.csv_address = db_address + '_Parameters'
            self.pymc2_M.write_csv(self.csv_address, variables=variables_list)
             
        #Print again the output prediction for the entire trace           
        self.display_run_data(self.MAP_Model, variables_list)
 
        #Close the database
        self.pymc2_M.db.close()       

    def load_pymc_database_manual(self, db_address, burning = 0, params_list = None):
                
        #Load the pymc output textfile database
        pymc_database = pymc2.database.pickle.load(db_address)
        
        #Create a dictionaries with the traces and statistics
        traces_dic = {}
        stats_dic = OrderedDict()
        stats_dic['true_values'] = empty(len(params_list))   
        
        #This variable contains all the traces from the MCMC (stochastic and deterministic)
        traces_list = pymc_database.trace_names[0] 
    
        #Get statistics from the run
        for i in range(len(traces_list)):
            
            print i, traces_list[i]
            
            trace               = traces_list[i]
            stats_dic[trace]    = OrderedDict()
            trace_array         = pymc_database.trace(trace)[burning:]
            traces_dic[trace]   = trace_array
             
            if 'dict' not in trace:
                stats_dic[trace]['mean']                    = mean(trace_array)
                stats_dic[trace]['median']                  = median(trace_array)
                stats_dic[trace]['standard deviation']      = std(trace_array)
                stats_dic[trace]['n']                       = trace_array.shape[0]
                stats_dic[trace]['16th_p']                  = percentile(trace_array, 16)
                stats_dic[trace]['84th_p']                  = percentile(trace_array, 84)
                stats_dic[trace]['95% HPD interval']        = (stats_dic[trace]['16th_p'], stats_dic[trace]['84th_p'])
                stats_dic[trace]['trace']                   = trace_array
                 
                if trace in params_list:
                    #Get the right true value key
                    key_true = trace if trace != 'T_He' else 'T_low'
                     
                    #Special cases
                    if key_true == 'ne':
                        key_true = 'n_e'
                     
                    stats_dic[trace]['true_value'] = self.obj_data[key_true]
                 
                if params_list is not None:
                    if trace in params_list:
                        print trace, stats_dic[trace]['mean']
             
        #Generate a MCMC object to recover all the data from the run
        dbMCMC = pymc2.MCMC(traces_dic, pymc_database)
        
        return dbMCMC, stats_dic  
                            
    def display_run_data(self, database, variables_list):
        
        for param in variables_list:
            param_entry = getattr(database, param, None)
            if param_entry is not None:
                try:
                    print '-{} {}'.format(param, param_entry.value)
                except:
                    print 'I could not do it ', param

# #Observational RS3 ratio
# self.S3_6312A_idx = (self.obj_data['colLine_labes'] == 'S3_6312A')
# self.S3_9069A_idx = (self.obj_data['colLine_labes'] == 'S3_9069A')
# self.S3_9531A_idx = (self.obj_data['colLine_labes'] == 'S3_9531A')
# S3_6312A_flux, S3_9069A_flux, S3_9531A_flux = self.obs_metal_fluxes[self.S3_6312A_idx], self.obs_metal_fluxes[self.S3_9069A_idx], self.obs_metal_fluxes[self.S3_9531A_idx]
# S3_6312A_err, S3_9069A_err, S3_9531A_err =  self.obs_metal_Error[self.S3_6312A_idx], self.obs_metal_Error[self.S3_9069A_idx], self.obs_metal_Error[self.S3_9531A_idx] 
# RS3 = ufloat(S3_6312A_flux,S3_6312A_err) / (ufloat(S3_9069A_flux,S3_9069A_err) + ufloat(S3_9531A_flux,S3_9531A_err))
# self.RS3_obs = RS3.nominal_value
# self.RS3_err = RS3.std_dev