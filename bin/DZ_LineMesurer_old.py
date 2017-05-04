from bin.DZ_Logs        import LineMesurer_Log
from collections    import OrderedDict
from numpy          import array, concatenate, ones, power, searchsorted, sqrt, std, where, loadtxt, linspace, median, max as npmax, empty, full
from uncertainties  import ufloat 
from lib.Math_Libraries.FittingTools import Fitting_Gaussians, Python_linfit

class LineMesurer(Fitting_Gaussians, LineMesurer_Log):
    
    def __init__(self, Conf_Folder, LinesLogHeader_Name):

        #Import the library with the fitting algorithms with gaussians
        Fitting_Gaussians.__init__(self)
        LineMesurer_Log.__init__(self, Conf_Folder, LinesLogHeader_Name)
        
        #Line definition
        self.Current_Label          = None
        self.Current_Ion            = None
        self.Current_TheoLoc        = None

        #Textfile where we store our data
        self.Current_LinesLog       = None
                              
        #Vector we use to store the 6 wavelengths which define the emision line location and the two continuums
        self.Selections             = []
        
        #Indexes for the "Selections vector" wavelengths              
        self.ind1                   = None      #Wave1
        self.ind2                   = None      #Wave2
        self.ind3                   = None      #Wave3
        self.ind4                   = None      #Wave4
        self.ind5                   = None      #Wave5
        self.ind6                   = None      #Wave6
        
        #New design
        self.BoxSize                = 70
        
        #Dictionary to store the data
        self.Parameter_dict = OrderedDict.fromkeys(self.ColumnHeaderVector)
        
        #Extra mechanisms for testing the wide component on Halpha      
        self.force_WD = True               
                
    def Dazer_LineMesuring(self, SubWave, SubInt, Methodology = 'lmfit', store_data = True):
        
        #Clear the dictionaries with the line data
        self.Fitting_dict   = self.Fitting_dict.fromkeys(self.Fitting_dict)
        
        #Calculation and plotting  parameters of adjacent continuum     
        WPoint1, Wpoint2, FPoint1, FPoint2 = self.ContinuumRegions(SubInt, SubWave)
        
        #Check if emission or absorption line as well as line mixture
        self.Check_GaussianMixture(SubWave, Methodology = 'lmfit', force_check=True)   

        #Convert emission lines to gaussian equivalents
        GaussianWave, GaussianInt = self.Measure_LineIntensity(SubWave, SubInt, store_data = store_data)  

        #Store line data in the log                   
        if store_data and self.Fitting_dict['start treatment']:
            self.SaveLM_to_file(str(self.Current_TheoLoc), self.Fitting_dict['blended number'], SavingFile = self.Current_LinesLog)
              
        #These are the parameters which Dazer needs to make the plots
        return GaussianWave, GaussianInt, WPoint1, Wpoint2, FPoint1, FPoint2
    
    def Batch_LineMesuring(self, SubWave, SubInt, wavelengths_list, Measuring_Method = 'lmfit', store_data = False):
        
        #Clear the dictionaries with the line data
        self.Fitting_dict   = self.Fitting_dict.fromkeys(self.Fitting_dict)
        
        #Six wavelengths which mark the blue continuum, the line, and the red continuum
        self.Selections = wavelengths_list
                
        #Get indeces of line_regions in spectrum wavelength
        self.ind1, self.ind2, self.ind3, self.ind4, self.ind5, self.ind6 = searchsorted(SubWave, [self.Selections[0], self.Selections[1], self.Selections[2], self.Selections[3], self.Selections[4], self.Selections[5]])
                
        #Calculation and plotting  parameters of adjacent continuum     
        WPoint1, Wpoint2, FPoint1, FPoint2 = self.ContinuumRegions(SubInt, SubWave)
                
        #Check if emission or absorption line as well as line mixture #WE SHOULD ONLY FEED HERE THE LINE REGION NOT ALL
        self.Check_GaussianMixture(SubWave, Methodology = Measuring_Method, force_check=True)   

        #Convert emission lines to gaussian equivalents
        GaussianWave, GaussianInt = self.Measure_LineIntensity(SubWave, SubInt)  

        #Store line data in the log                   
        if store_data and self.Fitting_dict['start treatment']:
            self.SaveLM_to_file(str(self.Current_TheoLoc), self.Fitting_dict['blended number'], SavingFile = self.Current_LinesLog)
                      
        return 

    def command_LineMesuring(self, SubWave, SubInt, wavelengths_list, line_wavelength = None, line_label = None, Measuring_Method = 'lmfit', store_data = False):
        
        if self.Current_TheoLoc == None:
            self.Current_TheoLoc = line_wavelength
        
        if self.Current_Label == None:
            self.Current_Label = line_label
        
        #Clear the dictionaries with the line data
        self.Fitting_dict   = self.Fitting_dict.fromkeys(self.Fitting_dict)
        
        #Six wavelengths which mark the blue continuum, the line, and the red continuum
        self.Selections = wavelengths_list
                
        #Get indeces of line_regions in spectrum wavelength
        self.ind1, self.ind2, self.ind3, self.ind4, self.ind5, self.ind6 = searchsorted(SubWave, [self.Selections[0], self.Selections[1], self.Selections[2], self.Selections[3], self.Selections[4], self.Selections[5]])
                
        #Calculation and plotting  parameters of adjacent continuum     
        WPoint1, Wpoint2, FPoint1, FPoint2 = self.ContinuumRegions(SubInt, SubWave)
                        
        #Check if emission or absorption line as well as line mixture #WE SHOULD ONLY FEED HERE THE LINE REGION NOT ALL
        self.Check_GaussianMixture(SubWave, Methodology = Measuring_Method, force_check=True)   

        #Convert emission lines to gaussian equivalents
        GaussianWave, GaussianInt = self.Measure_LineIntensity(SubWave, SubInt)  
                
        #Store line data in the log                   
        if store_data:
            self.SaveLM_to_file(str(self.Current_TheoLoc), self.Fitting_dict['blended number'], SavingFile = self.Current_LinesLog)
                     
        return self.Fitting_dict

    def ContinuumRegions(self,SubInt, SubWave, adjacent_continuum = True):
        
        #In this case we use adjacent regions to compute the continuum level
        if adjacent_continuum == True:
        
            #We generate arrays containing the wavelength and flux values from blue and red continuum             
            FluxCont_BothSides                  = concatenate((SubInt[self.ind1:self.ind2],   SubInt[self.ind5:self.ind6]))
            WaveCont_BothSides                  = concatenate((SubWave[self.ind1:self.ind2],  SubWave[self.ind5:self.ind6]))
            
            Res = (SubWave[-1] - SubWave[0]) / len(SubWave)
            
            self.Fitting_dict['ContinuumWidth'] = Res * len(WaveCont_BothSides)
            
            #We generate and array with the standard deviation of the points on the left and the right                        
            FluxError_BothSides                 = concatenate((std(SubInt[self.ind1:self.ind2]) * ones(len(SubInt[self.ind1:self.ind2])) , std(SubInt[self.ind5:self.ind6]) * ones(len(SubInt[self.ind5:self.ind6]))))
            
            #We perform a linear regresion taking into consideration the error in the flux
            m_Continuum, n_continuum            = Python_linfit(WaveCont_BothSides, FluxCont_BothSides, FluxError_BothSides, errors_output = False)
            
            #We calculate the wavelength and flux at the middle region of each continuums for graphical purposes
            WPoint1                             = (SubWave[self.ind2] + SubWave[self.ind1]) / 2
            Wpoint2                             = (SubWave[self.ind6] + SubWave[self.ind5]) / 2
            FPoint1                             = m_Continuum * WPoint1 + n_continuum
            FPoint2                             = m_Continuum * Wpoint2 + n_continuum
                
            #We calculate the linear continuum as a vector for the whole central region and the middle point
            Linear_Flux                         = m_Continuum * SubWave + n_continuum
            self.Fitting_dict['zerolev_median'] = m_Continuum * ((SubWave[self.ind3] + SubWave[self.ind4]) / 2) + n_continuum
    
            #Finally we calculate the dispersion in the continuum using the blue and red regions only
            Dif_FluxCont                        = FluxCont_BothSides - concatenate((Linear_Flux[self.ind1:self.ind2], Linear_Flux[self.ind5:self.ind6]))
            self.Fitting_dict['zerolev_sigma']  = std(Dif_FluxCont)
            
            
            #Storing data into the fitting dictionary
            self.Fitting_dict['m_Continuum']    =  m_Continuum     #Assuming a line this is the gradient (m)
            self.Fitting_dict['n_Continuum']    =  n_continuum     #Assuming a line this is the y axis interception point (n)
            self.Fitting_dict['ContinuumFlux']  =  Linear_Flux
                                
            return WPoint1, Wpoint2, FPoint1, FPoint2  
        
        else:
            
            print 'Middle region continua'
            
            return
              
    def Check_GaussianMixture(self, SubWave, Methodology = 'lmfit', force_check = True, wide_component_present = False):
        
        #Store the procedure to use
        self.Fitting_dict['Fitting method']     = Methodology
        
        #Define dictionary initial conditions
        self.Fitting_dict['Deblend check']      = False         #This logic is true when we have a blended line
        self.Fitting_dict['start treatment']    = False         #This logic is true when a blended group is observed
        self.Fitting_dict['blended number']     = 1             #This integer describes the number of blended lines. Currently it must be equal to the number in the blended list elements.
        self.Fitting_dict['Blended list']       = self.BlendedLines_Table()
        self.Fitting_dict['Add_wideComponent']  = False
                                  
        #Select fitting type (Simple or bootstrap) 
        if 'MC' not in self.Fitting_dict['Fitting method']:
            self.Fitting_dict['MC_iteration'] = 1
        else:
            self.Fitting_dict['MC_iteration'] = 500
            
        #Check for observed wide component
        if self.force_WD:
            self.Fitting_dict['Wide component'] = True
        
        #Automatic analysis for blended lines
        if force_check: 
            #Check we are actually including the emission line intended
            if  SubWave[self.ind3] < self.Current_TheoLoc < SubWave[self.ind4]:                                   #This is a security check that should be at the top
                self.Fitting_dict['line number'] = 1
                self.Fitting_dict['start treatment'] = True 
            
                #Check we are measuring a blended line
                for Blended_group_Labels in self.Fitting_dict['Blended list'][0]:                                   
                                        
                    #Check if the line is in a blended group
                    if self.Current_Label in Blended_group_Labels:
                        self.Fitting_dict['Deblend check']  = True
                        self.Fitting_dict['line group']     = self.Fitting_dict['Blended list'][0].index(Blended_group_Labels)
                        
                        First_Group_line                    = Blended_group_Labels[0]
                                                
                        #Check whether the current line is the first in the set: Positive case proceed to deblend
                        if self.Current_Label == First_Group_line:
                                   
                                                        
                            #Security check that all the blended lines are found within the wavelength ranged marked      
                            Line_Wavelength =   array(self.Fitting_dict['Blended list'][1][self.Fitting_dict['line group']])
                                                        
                            outside_min     =   len(where(SubWave[self.ind3] > Line_Wavelength)[0])
                            outside_max     =   len(where(SubWave[self.ind4] < Line_Wavelength)[0])
                            
                            if (outside_min == 0) and (outside_max == 0):

                                self.Fitting_dict['start treatment']          = True
                                self.Fitting_dict['blended number']         = len(Blended_group_Labels)
                                self.Fitting_dict['blended label']          = "_".join(Blended_group_Labels)
                                self.Fitting_dict['blended wavelengths']    = array(self.Fitting_dict['Blended list'][1][self.Fitting_dict['line group']])
                                self.Fitting_dict['blended index']          = 0
                            else:
                                self.Fitting_dict['Deblend check'] = False  
                                
                        #Negative case use data in the log for the blended line            
                        else:                
                                        
                            Index_Labels        = where(self.ColumnHeaderVector == ('Line_Label'))[0]  #No me acaba de convencer esto....
                            Index_Blended_Set   = where(self.ColumnHeaderVector == ('Blended_Set'))[0]
                            
                            MeasuredLine_Labels = loadtxt(self.Current_LinesLog, dtype=str, skiprows = self.TableHeaderSize, usecols = Index_Labels)
                            MeasuredLine_Sets   = loadtxt(self.Current_LinesLog, dtype=str, skiprows = self.TableHeaderSize, usecols = [Index_Blended_Set])                                        
                                                
                            self.Fitting_dict['blended index']      = where(MeasuredLine_Labels == self.Current_Label)[0]
                            
                            #Check whether the lines have been measured as part of a blended set 
                            if (First_Group_line in MeasuredLine_Labels) and (self.Current_Label in MeasuredLine_Labels):                  
                                
                                if MeasuredLine_Sets[self.Fitting_dict['blended index']] !=  'None':
                                    self.Fitting_dict['Deblend check']      = True
                                    self.Fitting_dict['start treatment']    = False
#                                     self.Fitting_dict['blended label']      = "_".join(Blended_group_Labels)
                                    self.Fitting_dict['blended label']      = MeasuredLine_Sets[self.Fitting_dict['blended index']]
                                    self.Fitting_dict['blended number']     = len(Blended_group_Labels)

                                else:
                                    self.Fitting_dict['Deblend check']      = False
                                    self.Fitting_dict['start treatment']    = True
                            else:
                                self.Fitting_dict['Deblend check']          = False
                                self.Fitting_dict['start treatment']        = True
                
        #Quick check for the wide component (It must be still be set to zero at each line)
        if (self.Current_TheoLoc == 6548.05) and self.Fitting_dict['Wide component'] and self.Fitting_dict['Deblend check'] and self.Fitting_dict['start treatment']:
            self.Fitting_dict['WC_theowavelength']  = 6548.0
            self.Fitting_dict['Add_wideComponent']  = True          
#             self.Fitting_dict['blended number'] += 1
                                                                         
        return
    
    def Measure_LineIntensity(self, SubWave, SubInt, stellar_continua_removed = False, Current_Wavelength = None, TableAddress = None, TableHeaderSize = None, store_data = True):
                            
        # Check we are actually including the line intended #GET A GOOD ERROR MESSAGE FOR THE CASE THIS DOES HAPPEN                
        if self.Fitting_dict['start treatment']:
            
            #Load lmfit parameters:
            self.Load_lmfit_parameters(SubWave[self.ind3:self.ind4], SubInt[self.ind3:self.ind4], self.Fitting_dict['ContinuumFlux'][self.ind3:self.ind4], self.Fitting_dict['zerolev_sigma'], self.Fitting_dict['blended number'], wide_component = False)
 
            #Mesure the line flux and perform the fix
            self.fit_line(self.Fitting_dict['x_norm'], self.Fitting_dict['y_norm'], self.Fitting_dict['zerolev_norm'], self.Fitting_dict['sig_zerolev_norm'], self.Fitting_dict['blended number'], self.Fitting_dict['lmfit_params'], self.Fitting_dict['lmfit_params_wide'], self.Fitting_dict['MC_iteration'])
                        
            #Scale line parameters to physical units
            self.scale_lmfit_params(SubWave[self.ind3:self.ind4], SubInt[self.ind3:self.ind4], self.Fitting_dict['lmfit_output'], self.Fitting_dict['blended number'], self.Fitting_dict['x_scaler'], self.Fitting_dict['y_scaler'], self.Fitting_dict['Fitting method'])
            
            #Determine the physical parameters from the line
            self.Line_physical_data(SubWave, SubInt, self.Fitting_dict['blended number'], self.Fitting_dict['zerolev_sigma'], Current_Wavelength, TableAddress, self.Fitting_dict['Fitting method'], stellar_continua_removed)
                                 
        #In this case we get the data from the line data from the log
        elif (self.Fitting_dict['Deblend check'] == True) and (self.Fitting_dict['start treatment'] == False):
                        
            #Load the gaussian fit parameters of all the components 
            A_vector        = self.GetDataInTable('Blended_Set',"A",            self.Fitting_dict['blended label'], self.Current_LinesLog, self.TableHeaderSize, dataype='float')
            mu_vector       = self.GetDataInTable("Blended_Set","mu",           self.Fitting_dict['blended label'], self.Current_LinesLog, self.TableHeaderSize, dataype='float')
            sigma_vector    = self.GetDataInTable('Blended_Set',"sigma",        self.Fitting_dict['blended label'], self.Current_LinesLog, self.TableHeaderSize, dataype='float')
            FluxG_vector    = self.GetDataInTable('Blended_Set',"Flux_Gauss",   self.Fitting_dict['blended label'], self.Current_LinesLog, self.TableHeaderSize, dataype='float')
            Eqw_vector      = self.GetDataInTable('Blended_Set',"Eqw",          self.Fitting_dict['blended label'], self.Current_LinesLog, self.TableHeaderSize, dataype='float')
                        
            for i in range(self.Fitting_dict['blended number']):
                index = str(i)
                self.Fitting_dict['A' + index]      = ufloat(A_vector[i],         0.0)
                self.Fitting_dict['mu' + index]     = ufloat(mu_vector[i],        0.0)
                self.Fitting_dict['sigma' + index]  = ufloat(sigma_vector[i],     0.0)
                self.Fitting_dict['Eqw' + index]    = ufloat(Eqw_vector[i],       0.0)
                self.Fitting_dict['FluxG' + index]  = ufloat(FluxG_vector[i],     0.0)
 
            #Calculate the gaussian curve for plotting
            self.Fitting_dict['x_resample']             = linspace(SubWave[self.ind3], SubWave[self.ind4], 50)
            self.Fitting_dict['zerolev_resample']       = self.Fitting_dict['m_Continuum'] * self.Fitting_dict['x_resample'] + self.Fitting_dict['n_Continuum']
            
            self.Fitting_dict['peak_waves']  = mu_vector
            self.Fitting_dict['peak_Maxima'] = A_vector + self.Fitting_dict['zerolev_median']

            #The y resample is a list of arrays            
            self.Fitting_dict['y_resample'] = [] 
            for j in range(self.Fitting_dict['blended number']):
                Comp = [str(j)]
                self.Fitting_dict['y_resample'].append(self.gaussian_curve_SingleMixture_from_dict(self.Fitting_dict, self.Fitting_dict['x_resample'], self.Fitting_dict['zerolev_resample'], Ncomps = Comp))

        return self.Fitting_dict['x_resample'], self.Fitting_dict['y_resample']

    def Line_physical_data(self, SubWave, SubInt, Ncomps, sigma_continuum, Current_Wavelength, TableAddress, Methodology, stellar_continua_removed = False):

        #This is necessary in case we must include the stellar continua:
        if stellar_continua_removed:
            LocalMedian = ufloat(float(self.GetDataInTable("TheoWavelength", "Continuum_Median", Current_Wavelength, TableAddress, 2)), sigma_continuum)
        else:
            LocalMedian = ufloat(self.Fitting_dict['zerolev_median'], sigma_continuum)

        #Parameters for the classical error estimation
        D                                   = float(SubWave[-1] - SubWave[0]) / float(len(SubWave))                             #Region dispertion
        Npix_line                           = float(len(SubWave[self.ind3:self.ind4]))                                          #Number of pixels occupied by the line                                
        self.Fitting_dict['ContinuumWidth'] = float(len(SubInt[self.ind1:self.ind2]) + len(SubInt[self.ind5:self.ind6])) * D    #Number of pixels at the sides 
        
        #Calculation for each of the four conditions (Single - Blended; Simple - MC)
        
        #--Simple treatement:
        if 'MC' not in Methodology:
            
            #Single
            if Ncomps == 1:
                Flux                        = self.Fitting_dict['FluxI'].nominal_value
                Eqw                         = Flux / LocalMedian.nominal_value
                Sigma_F, Sigma_Eqw          = self.line_classic_error(Flux, Eqw, sigma_continuum, D, Npix_line)
                self.Fitting_dict['Eqw0']   = ufloat(Eqw,   Sigma_Eqw)
                self.Fitting_dict['FluxI']  = ufloat(Flux,  Sigma_F)
                
            #Blended
            else:
                for i in range(Ncomps):
                    index               = str(i)
                    Flux                = self.Fitting_dict['FluxG' + index].nominal_value
                    Eqw                 = Flux / LocalMedian.nominal_value
                    Sigma_F, Sigma_Eqw  = self.line_classic_error(Flux, Eqw, sigma_continuum, D, Npix_line)
                    self.Fitting_dict['FluxG' + index]  = ufloat(Flux, Sigma_F)
                    self.Fitting_dict['Eqw' + index]    = ufloat(Eqw, Sigma_Eqw)
                
        #--MC treatment
        else:
            #Single
            if Ncomps == 1:
                self.Fitting_dict['Eqw0']   = self.Fitting_dict['FluxI'] / LocalMedian
            
            #Blended line
            else:
                for i in range(Ncomps):
                    index = str(i)
                    self.Fitting_dict['Eqw'+index] = self.Fitting_dict['FluxG' + index] / LocalMedian        
        
        return

    def line_classic_error(self, Flux, Eqw, sigma_continuum, D, Npix_line):
        
        Sigma_F     = sigma_continuum * D * sqrt(2 * Npix_line + Eqw/D)
        Sigma_Eqw   = Eqw / Flux * sigma_continuum * D * sqrt(Eqw/D + 2 * Npix_line + power(Eqw / D, 2) / Npix_line)
               
        return Sigma_F, abs(Sigma_Eqw)
                               
    def Emission_Threshold(self, LineLoc, TotalWavelen, TotalInten):
        
        #Use this method to determine the box and location of the emission lines
        Bot                 = LineLoc - self.BoxSize
        Top                 = LineLoc + self.BoxSize
        
        indmin, indmax      = searchsorted(TotalWavelen, (Bot, Top))
        if indmax > (len(TotalWavelen)-1):
            indmax = len(TotalWavelen)-1
        
        PartialWavelength   = TotalWavelen[indmin:indmax]
        PartialIntensity    = TotalInten[indmin:indmax]
        
        Bot                 = LineLoc - 2
        Top                 = LineLoc + 2
        
        indmin, indmax      = searchsorted(PartialWavelength, (Bot, Top))
        
        LineHeight          = npmax(PartialIntensity[indmin:indmax])
        LineExpLoc          = median(PartialWavelength[where(PartialIntensity == LineHeight)])
              
        return PartialWavelength, PartialIntensity, LineHeight, LineExpLoc   
       
    def BlendedLines_Table(self):
        
        #DeblendingLines = (    (Labels), (Wavelengths), (Ions)    )
        
        DeblendingLines =  (
                            (
                            ['O2_3726A',    'O2_3729A'],
                            ['Ne3_3968A',   'H1_3970A'],
                            ['O3_5007A',    'He1_5016A'],
                            ['N2_6548A',    'H1_6563A','N2_6584A'],   
                            ['O2_7319A',    'O2_7330A'],
                            ),(
                            [3726.032,      3728.815],
                            [3968.0,        3970.0],
                            [5006.843,      5016],
                            [6548.05,       6562.819, 6583.46],
                            [7319.0,        7330.0],
                            ),( 
                            ['[OII]',       '[OII]'],
                            ['[NeIII]',     'HEps_7_2'],
                            ['[OIII]',      'HeI_5016'],
                            ['[NII]',       'Halpha_3_2','[NII]'],   
                            ['[OII]',       '[OII]'],                     
                            )
                           )
        
        return DeblendingLines
                   
    def SaveLM_to_file(self, Current_Em_Wave, Ncomps, SavingFile, ForceRemake = False):
  
        #Clear the dictionary an fill it with arrays with the same length as Ncomps        
        self.Parameter_dict = {key: ['empty'] * Ncomps for key in self.Parameter_dict.keys()}
                
        #--Confirm the current line for blended components:
        if self.Fitting_dict['Deblend check'] == False:
            Linelabel       = [self.Current_Label]
            Ion             = [self.Current_Ion]
            TheoWavelength  = [self.Current_TheoLoc]
            Blended_Set     = None   
        else:
            Linelabel       = self.Fitting_dict['Blended list'][0][self.Fitting_dict['line group']]
            Ion             = self.Fitting_dict['Blended list'][2][self.Fitting_dict['line group']]
            TheoWavelength  = self.Fitting_dict['Blended list'][1][self.Fitting_dict['line group']]
            Blended_Set     = self.Fitting_dict['blended label']

            if self.Fitting_dict['Add_wideComponent']:
                Linelabel       = Linelabel + [Linelabel[1] + '_w']
                Ion             = Ion + [Ion[1] + '_w']
                TheoWavelength  = TheoWavelength + [self.Fitting_dict['peak_waves'][1]] #Dirty trick so it is not saved in exactly the same place
                Blended_Set     = Blended_Set + '_wide'
                
        #--Confirm wide components
        Wide_components = array(['False'] * Ncomps) 
        if self.Fitting_dict['Add_wideComponent']:
            Wide_components[1] = 'True'
        
        #Case with wide component
        for i in range(Ncomps):
            index = str(i)

            self.Parameter_dict['Line_Label'][i]        = Linelabel[i]
            self.Parameter_dict['Ion'][i]               = Ion[i]
            self.Parameter_dict['TheoWavelength'][i]    = TheoWavelength[i]

            self.Parameter_dict['Wave_Peak'][i]         = self.Fitting_dict['peak_waves'][i]
            self.Parameter_dict['Flux_Int'][i]          = self.Fitting_dict['FluxI'].nominal_value
            self.Parameter_dict['Flux_Gauss'][i]        = self.Fitting_dict['FluxG' + index].nominal_value
            self.Parameter_dict['Eqw'][i]               = self.Fitting_dict['Eqw'   + index].nominal_value
            
            self.Parameter_dict['Error_FluxI'][i]       = self.Fitting_dict['FluxI'].std_dev
            self.Parameter_dict['Error_FluxG'][i]       = self.Fitting_dict['FluxG' + index].std_dev
            self.Parameter_dict['Error_Eqw'][i]         = self.Fitting_dict['Eqw'   + index].std_dev                  
                             
            self.Parameter_dict['A'][i]                 = self.Fitting_dict['A'     + index].nominal_value
            self.Parameter_dict['mu'][i]                = self.Fitting_dict['mu'    + index].nominal_value
            self.Parameter_dict['sigma'][i]             = self.Fitting_dict['sigma' + index].nominal_value
            self.Parameter_dict['h3_skewness'][i]       = 0.0
            self.Parameter_dict['h4_kurtosis'][i]       = 0.0 
                     
            self.Parameter_dict['A_err'][i]             = self.Fitting_dict['A'     + index].std_dev
            self.Parameter_dict['mu_err'][i]            = self.Fitting_dict['mu'    + index].std_dev
            self.Parameter_dict['sigma_err'][i]         = self.Fitting_dict['sigma' + index].std_dev
            self.Parameter_dict['h3_skewness_err'][i]   = 0.0
            self.Parameter_dict['h4_kurtosis_err'][i]   = 0.0 
            
            self.Parameter_dict['Continuum_Median'][i]  = self.Fitting_dict['zerolev_median']
            self.Parameter_dict['Continuum_m'][i]       = self.Fitting_dict['m_Continuum']
            self.Parameter_dict['Continuum_n'][i]       = self.Fitting_dict['n_Continuum']
            self.Parameter_dict['Continuum_Width'][i]   = self.Fitting_dict['ContinuumWidth']
            self.Parameter_dict['Continuum_sigma'][i]   = self.Fitting_dict['zerolev_sigma']
            
            self.Parameter_dict['Blended_Set'][i]       = Blended_Set
            self.Parameter_dict['Wide_component'][i]    = Wide_components[i]

            self.Parameter_dict['Wave1'][i]             = self.Selections[0] 
            self.Parameter_dict['Wave2'][i]             = self.Selections[1]
            self.Parameter_dict['Wave3'][i]             = self.Selections[2]
            self.Parameter_dict['Wave4'][i]             = self.Selections[3]
            self.Parameter_dict['Wave5'][i]             = self.Selections[4]
            self.Parameter_dict['Wave6'][i]             = self.Selections[5]
             
                        
            self.Parameter_dict['Comments'][i]          = 'None'

        #Check if the file exists, in the negative case make one
        self.CleanTableMaker(SavingFile, ForceRemake)

        #Store the data 
        for i in range(self.Fitting_dict['blended number']):
            
            Elemental_Line_String = self.CreateEmissionLineReccord(self.Parameter_dict['TheoWavelength'][i], self.Parameter_dict['Ion'][i], self.Parameter_dict['Line_Label'][i])
          
            self.InsertNewLine(Elemental_Line_String,  self.Parameter_dict['TheoWavelength'][i], self.Parameter_dict['Line_Label'][i], SavingFile)
#             self.InsertNewLine(Elemental_Line_String, self.Parameter_dict['TheoWavelength'][i], SavingFile)
            self.Replace_Row_byLabel(i, self.Parameter_dict['Line_Label'][i], self.Parameter_dict, SavingFile)
#             self.Replace_Row(i, self.Parameter_dict['TheoWavelength'][i], self.Parameter_dict, SavingFile)
             
        return
