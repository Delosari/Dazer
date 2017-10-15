from collections        import OrderedDict
from pandas             import read_csv, Series
from numpy              import concatenate, ones, power, searchsorted, sqrt, std, where, loadtxt, linspace, median, max as npmax, True_, False_, bool_, empty, isnan
from uncertainties      import ufloat 
from DZ_Logs            import LineMesurer_Log
from pandas.core.frame  import DataFrame
from lib.Math_Libraries.FittingTools import Fitting_Gaussians, Fitting_Gaussians_v2, Python_linfit, gaussian_components, gaussian_curve, gaussian_mixture
from timeit import default_timer as timer

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
                                      
        #Without a graphical interface we need a dataframe with the possible blended lines
        self.blendedLines_Table()
                
        #New design
        self.BoxSize = 70
        
        #Dictionary to store the data
        self.Parameter_dict = OrderedDict.fromkeys(self.ColumnHeaderVector)
        
        #Extra mechanisms for testing the wide component on Halpha      
        self.force_WD = True         
        
        #Parameters to save in log #WARNING: Missing emission line type
        self.saving_parameters_list = ['Ion',
                                  'lambda_theo',
                                  'lambda_obs',
                                  'flux_intg',
                                  'flux_intg_er',
                                  'flux_gauss', 
                                  'flux_gauss_er', 
                                  'eqw',
                                  'eqw_er',
                                  'A',
                                  'A_er',
                                  'mu',
                                  'mu_er',
                                  'sigma',
                                  'sigma_er',
                                  'zerolev_mean', #continuum linear value at line region center
                                  'zerolev_std',
                                  'zerolev_width',                                  
                                  'm_zerolev',
                                  'n_zerolev',                                  
                                  'Wave1',
                                  'Wave2',
                                  'Wave3',
                                  'Wave4',
                                  'Wave5',
                                  'Wave6',
                                  'blended_check',
                                  'line_number',
                                  'group_label',
                                  'add_wide_component',
                                  'fit_routine']
                             
    def measure_line(self, subWave, subFlux, wavelengths_list, lines_dataframe = None, line_label = None, Measuring_Method = 'lmfit', store_data = False):
        
        #Clear the dictionaries with the line data
        self.fit_dict = Series(index = self.fitting_parameters, dtype=object)
         
        #Get indeces of line_regions in spectrum wavelength
        self.fit_dict.loc['idx0':'idx5'] = searchsorted(subWave, wavelengths_list)
        self.fit_dict.loc['Wave1':'Wave6'] = wavelengths_list
        
        #Calculation and plotting  parameters of adjacent continuum 
        self.continuum_Regions(subWave, subFlux)

        #Check if emission or absorption line as well as line mixture #WE SHOULD ONLY FEED HERE THE LINE REGION NOT ALL
        self.check_GaussianMixture(subWave, Measuring_Method, lines_dataframe, force_check=True)

        #Convert emission lines to gaussian equivalents
        start = timer()
        self.calculate_Intensity(subWave, subFlux, lines_dataframe)
        end = timer()
        print 'calculate_Intensity', (end - start) 

        #Store line data in the log                   
        if store_data and  self.fit_dict.start_treatment:
            if lines_dataframe is None:
                self.current_df = DataFrame(columns = self.saving_parameters_list)
            else:
                self.current_df = lines_dataframe
            
            self.line_to_Logdf()
        
        #Otherwise return the data from the fit
        else:
            return self.fit_dict

    def dazer_lineMeasuring(self, subWave, subFlux, wavelengths_list, lines_dataframe, Measuring_Method = 'lmfit'):
        
        #Clear the dictionaries with the line data
        self.fit_dict = Series(index = self.fitting_parameters, dtype=object)
         
        #Get indeces of line_regions in spectrum wavelength
        self.fit_dict.loc['idx0':'idx5'] = searchsorted(subWave, wavelengths_list)
        self.fit_dict.loc['Wave1':'Wave6'] = wavelengths_list
                
        #Calculation and plotting  parameters of adjacent continuum     
        self.continuum_Regions(subWave, subFlux)
                
        #Check if emission or absorption line as well as line mixture #WE SHOULD ONLY FEED HERE THE LINE REGION NOT ALL
        self.check_GaussianMixture(subWave, Measuring_Method, lines_dataframe, force_check=True)   

        #Convert emission lines to gaussian equivalents
        self.calculate_Intensity(subWave, subFlux, lines_dataframe)  
        
        #Save the fit:              
        if self.fit_dict.start_treatment:

            #Update the data frame with the new data
            self.line_to_Logdf()
            
            #Sort the data frame
            self.current_df.sort_values(['lambda_theo'], ascending=[True], inplace=True)
            
            #Store data to dataframe
            self.save_lineslog_dataframe(self.current_df, self.lineslog_df_address)
            
        return

    def continuum_Regions(self, subWave, subFlux, adjacent_continuum = True):
        
        #Indeces from the continuum regions
        idx1, idx2, idx3, idx4, idx5, idx6 = self.fit_dict[['idx0', 'idx1', 'idx2', 'idx3', 'idx4', 'idx5']]
        
        #In this case we use adjacent regions to compute the continuum level
        if adjacent_continuum == True:
        
            #Region resolution
            region_resolution = (subWave[-1] - subWave[0]) / len(subWave) #Should we save it?

            #We generate arrays containing the wavelength and flux values from blue and red continuum             
            FluxCont_BothSides = concatenate((subFlux[idx1:idx2],   subFlux[idx5:idx6]))
            WaveCont_BothSides = concatenate((subWave[idx1:idx2],  subWave[idx5:idx6]))

            self.fit_dict['zerolev_width'] = region_resolution * len(WaveCont_BothSides)
            
            #We generate and array with the standard deviation of the points on the left and the right  #WARNING no se si esto esta bien                     
            FluxError_BothSides = concatenate((std(subFlux[idx1:idx2]) * ones(len(subFlux[idx1:idx2])), std(subFlux[idx5:idx6]) * ones(len(subFlux[idx5:idx6]))))
            
            #We perform a linear regresion taking into consideration the error in the flux #WARNING no se si esto esta bien
            m_Continuum, n_continuum = Python_linfit(WaveCont_BothSides, FluxCont_BothSides, FluxError_BothSides, errors_output = False)
            
            #We calculate the wavelength and flux at the middle region of each continuums for graphical purposes
            self.fit_dict['Blue_wave_zerolev']  = (subWave[idx2] + subWave[idx1]) / 2
            self.fit_dict['Red_wave_zerolev']   = (subWave[idx6] + subWave[idx5]) / 2
            self.fit_dict['Blue_flux_zerolev']  = m_Continuum * self.fit_dict['Blue_wave_zerolev'] + n_continuum
            self.fit_dict['Red_flux_zerolev']   = m_Continuum * self.fit_dict['Red_wave_zerolev'] + n_continuum
               
            #We calculate the linear continuum as a vector for the whole central region and the middle point
            Linear_Flux = m_Continuum * subWave + n_continuum
    
            #Finally we calculate the dispersion in the continuum using the blue and red regions only
            Dif_FluxCont = FluxCont_BothSides - concatenate((Linear_Flux[idx1:idx2], Linear_Flux[idx5:idx6]))
            
            #Storing data into the fitting dictionary
            start = timer()
            self.fit_dict['zerolev_mean']   = m_Continuum * ((subWave[idx3] + subWave[idx4]) / 2) + n_continuum
            self.fit_dict['zerolev_std']    = std(Dif_FluxCont)
            self.fit_dict['m_zerolev']      = m_Continuum #Assuming a line this is the gradient (m)
            self.fit_dict['n_zerolev']      = n_continuum #Assuming a line this is the y axis interception point (n)
            self.fit_dict['zerolev_linear'] = Linear_Flux

            end = timer()
            print 'continua', end-start
        else:
            print 'Middle region continua'
            
        return
              
    def check_GaussianMixture(self, SubWave, Methodology, lines_dataframe, force_check = True):
        
        #Store the procedure to use
        self.fit_dict['fit_routine'] = Methodology
        
        #Define dictionary initial conditions
        line_in_region                          = False
        lines_in_region                         = False
        first_blended                           = False
        in_blended_group                        = False
        
        self.fit_dict['MC_iterations']          = 1
        self.fit_dict['blended_check']          = False         #This logic is true when we have a blended line
        self.fit_dict['start_treatment']        = False         #This logic is true when a blended group is observed
        self.fit_dict['line_number']            = 0             #This integer describes the number of blended lines. Currently it must be equal to the number in the blended list elements.
        self.fit_dict['add_wide_component']     = False
        self.fit_dict['wide_component']         = False
                
        self.fit_dict['group_label']            = None
        self.fit_dict['blended_lambdas']        = None
        self.fit_dict['blended_labels']         = None
        self.fit_dict['blended_ions']           = None
        self.fit_dict['blended_index']          = 0
                
        #Select fitting type (Simple or bootstrap) 
        if 'MC' not in self.fit_dict['fit_routine']:
            self.fit_dict['MC_iterations'] = 1
        else:
            self.fit_dict['MC_iterations'] = 500
         
        #Check for observed wide component #WARNING this should come from the screen
        if self.force_WD:
            self.fit_dict['wide_component'] = True
        
        #If label is from a wide component do not measure it
        if self.Current_Label is not None:
            if '_w' in self.Current_Label:
                force_check = False
                        
        #Automatic analysis for blended lines
        if force_check:
            
            #Check we are actually including the emission line intended
            if  SubWave[self.fit_dict.idx2] < self.Current_TheoLoc < SubWave[self.fit_dict.idx3]:  #This is a security check that should be at the top
                line_in_region = True
                self.fit_dict['line_number'] = 1
                
                #Check if this is the first line in a blended group
                if self.Current_Label in self.blended_lines_df.index:
                    first_blended = True
                    
                #Check if this line is in a blended group previously measured
                elif lines_dataframe is not None:
                    if self.Current_Label in lines_dataframe.index:
                        if (lines_dataframe.loc[self.Current_Label, 'group_label'] != 'None') and (lines_dataframe.loc[self.Current_Label, 'group_label'] is not None):
                            in_blended_group = True
            
                #Check all lines are located in the selected region
                if first_blended:
                    group_wavelengths = self.blended_lines_df.loc[self.Current_Label, 'Theowavelength']
                    
                    if (SubWave[self.fit_dict.idx2] < group_wavelengths[0]) and (group_wavelengths[-1] < SubWave[self.fit_dict.idx3]):
                        lines_in_region = True
                                    
            if (line_in_region and (first_blended is False) and (in_blended_group is False)) | ((line_in_region) and (lines_in_region is False) and in_blended_group is False):
                self.fit_dict['start_treatment']    = True
                self.fit_dict['line_number']        = 1
            
            #--Case with a blended line
            else:
                self.fit_dict['blended_check'] = True
                
                #First appereance of blended group
                if first_blended and lines_in_region:
                    self.fit_dict['start_treatment'] = True
                    
                    #WARNING in reality we should load this form the screen
                    self.fit_dict['group_label']        = self.blended_lines_df.loc[self.Current_Label, 'group_label']
                    self.fit_dict['blended_lambdas']    = self.blended_lines_df.loc[self.Current_Label, 'Theowavelength']
                    self.fit_dict['blended_labels']     = self.blended_lines_df.loc[self.Current_Label, 'line_labels']
                    self.fit_dict['blended_ions']       = self.blended_lines_df.loc[self.Current_Label, 'Ion']
                    self.fit_dict['line_number']        = len(self.fit_dict['blended_lambdas'])
                    #self.fit_dict['blended_index']      = 0

                #Case with a blended line previously measured
                elif in_blended_group:
                    self.fit_dict['start_treatment'] = False 

                    self.fit_dict['group_label'] = lines_dataframe.loc[self.Current_Label, 'group_label']
                    idces_group = (lines_dataframe.group_label == self.fit_dict['group_label'])
                                   
                    self.fit_dict['blended_labels']     = lines_dataframe.loc[idces_group].index.values
                    self.fit_dict['blended_lambdas']    = lines_dataframe.loc[idces_group, 'lambda_theo'].values
                    self.fit_dict['blended_ions']       = lines_dataframe.loc[idces_group, 'Ion'].values
                    self.fit_dict['line_number']        = len(self.fit_dict['blended_lambdas'])
                    #self.fit_dict['blended_index']      = self.fit_dict['blended_labels'].index(self.Current_Label)
                    
            #Dirty trick to always add a wide components for Halpha
            if (self.Current_TheoLoc == 6548.05)  and self.fit_dict['blended_check'] and self.fit_dict['start_treatment'] and self.fit_dict['wide_component']:
                self.fit_dict['add_wide_component'] = True          
        
#         print 'Bingo', self.fit_dict['blended_check'] 
           
        return
    
    def calculate_Intensity(self, subWave, subFlux, lineslog_df):
        
        #Load indeces for he regions
        idx3, idx4 = self.fit_dict[['idx2', 'idx3']]
                
        # Check we are actually including the line intended #GET A GOOD ERROR MESSAGE FOR THE CASE THIS DOES HAPPEN                
        if self.fit_dict['start_treatment']:
          
            #Load lmfit parameters:
            start = timer()
            self.load_lmfit_parameters(subWave[idx3:idx4], subFlux[idx3:idx4], self.fit_dict.zerolev_linear[idx3:idx4], self.fit_dict.zerolev_std, self.fit_dict.line_number, wide_component = False)
            end = timer()
            print 'load_lmfit_parameters', (end - start) 
            
            
            #Mesure the line flux and perform the fix
            start = timer()
            if self.fit_dict.blended_check is False:
                self.fit_single_line_BS(self.fit_dict.x_n, self.fit_dict.y_n, self.fit_dict.zerolev_n, self.fit_dict.sigZerolev_n, self.fit_dict.params_lmfit, 
                                        bootstrap_iterations=self.fit_dict['MC_iterations'])            
            else:
                self.fit_blended_line_BS_lmfit(self.fit_dict.x_n, self.fit_dict.y_n, self.fit_dict.zerolev_n, self.fit_dict.sigZerolev_n,  self.fit_dict.line_number, self.fit_dict.params_lmfit, 
                                      self.fit_dict.add_wide_component, self.fit_dict.params_lmfit_wide, bootstrap_iterations=self.fit_dict['MC_iterations'])
            end = timer()
            print 'fit_single_line_BS', (end - start) 
            
            
            #Scale line parameters to physical units
            start = timer()
            self.rescale_lmfit_params(subWave[idx3:idx4], subFlux[idx3:idx4], self.fit_dict.line_number, self.fit_dict.x_scaler, self.fit_dict.y_scaler, self.fit_dict.fit_routine)
            end = timer()
            print 'rescale_lmfit_params', (end - start)             
            
            #Determine the physical parameters from the line
            start = timer()
            self.linePhysicalData(subWave, subFlux, self.fit_dict.line_number, self.fit_dict.zerolev_std, self.fit_dict.fit_routine)
            end = timer()
            print 'linePhysicalData', (end - start)   
                                                             
        #In this case we get the data from the line data from the log
        elif (self.fit_dict['blended_check'] == True) and (self.fit_dict['start_treatment'] == False):
            
            #Labels for the group
            idcs_group = (lineslog_df.group_label == self.fit_dict.group_label)

            self.fit_dict['flux_intg'], self.fit_dict['flux_intg_er'] = lineslog_df.loc[idcs_group, 'flux_intg'].values, lineslog_df.loc[idcs_group, 'flux_intg_er'].values
                       
            #Load the gaussian fit parameters of all the components 
            A_vector, A_er_vector           = lineslog_df.loc[idcs_group, 'A'].values, lineslog_df.loc[idcs_group, 'A_er'].values
            mu_vector, mu_er_vector         = lineslog_df.loc[idcs_group, 'mu'].values, lineslog_df.loc[idcs_group, 'mu_er'].values
            sigma_vector, sigma_er_vector   = lineslog_df.loc[idcs_group, 'sigma'].values, lineslog_df.loc[idcs_group, 'sigma_er'].values
            fluxG_vector, fluxG_er_vector   = lineslog_df.loc[idcs_group, 'flux_gauss'].values, lineslog_df.loc[idcs_group, 'flux_gauss_er'].values
            eqw_vector, eqw_er_vector       = lineslog_df.loc[idcs_group, 'eqw'].values, lineslog_df.loc[idcs_group, 'eqw_er'].values

            #Load the data into line dictionary #WARNING: Add the error to the line
            for i in range(len(mu_vector)):
                idx = str(i)
                self.fit_dict['A' + idx],           self.fit_dict['A' + idx + '_er']            = A_vector[i],      A_er_vector[i]        
                self.fit_dict['mu' + idx],          self.fit_dict['mu' + idx + '_er']           = mu_vector[i],     mu_er_vector[i]      
                self.fit_dict['sigma' + idx],       self.fit_dict['sigma' + idx + '_er']        = sigma_vector[i],  sigma_er_vector[i]    
                self.fit_dict['flux_gauss' + idx],  self.fit_dict['flux_gauss' + idx + '_er']   = fluxG_vector[i],  fluxG_er_vector[i]
                self.fit_dict['eqw' + idx],         self.fit_dict['eqw' + idx + '_er']          = eqw_vector[i],    eqw_er_vector[i]       
 
            #Calculate the gaussian curve for plotting
            self.fit_dict['x_resample'] = linspace(subWave[idx3], subWave[idx4], 50)
            self.fit_dict['zerolev_resample'] = self.fit_dict['m_zerolev'] * self.fit_dict['x_resample'] + self.fit_dict['n_zerolev']
            
            self.fit_dict['maxLambdas']  = mu_vector
            self.fit_dict['maxPeaks'] = A_vector + self.fit_dict['zerolev_mean']

            if self.fit_dict.blended_check == False:
                self.fit_dict['y_resample'] = gaussian_curve(self.fit_dict['A0'], self.fit_dict['mu0'], self.fit_dict['sigma0'], self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'])
            else:
                self.fit_dict['y_resample'] = gaussian_mixture(self.fit_dict, self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'], map(str, range(self.fit_dict.line_number)))
                self.fit_dict['y_comps']    = gaussian_components(self.fit_dict, self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'], self.fit_dict.line_number)
        return
    
    def linePhysicalData(self, subWave, subInt, Ncomps, sigma_continuum, Methodology):

        #Load indeces for he regions
        idx1, idx2, idx3, idx4, idx5, idx6 = self.fit_dict[['idx0', 'idx1', 'idx2', 'idx3', 'idx4', 'idx5']]

        #This is necessary in case we must include the stellar continua:
        local_median = self.fit_dict['zerolev_mean']

        #Parameters for the classical error estimation
        Npix_line   = len(subWave[idx3:idx4])                   #Number of pixels occupied by the lin                                
        D           = (subWave[-1] - subWave[0]) / len(subWave) #Region dispertion
        
        #Continuum width calculation
        self.fit_dict['continuum_width'] = float(len(subInt[idx1:idx2]) + len(subInt[idx5:idx6])) * D    #Number of pixels at the sides 

#         #--Simple treatement:
#         if 'MC' not in Methodology:
#             
#             #Single
#             if Ncomps == 1:
#                 Flux                    = self.fit_dict['flux_intg']
#                 Eqw                     = Flux / local_median
#                 Sigma_F, Sigma_Eqw      = self.line_classic_error(Flux, Eqw, sigma_continuum, D, Npix_line)
#                    
#             #Blended
#             else:
#                 for i in range(Ncomps):
#                     index               = str(i)
#                     Flux                = self.fit_dict['flux_gauss' + index]
#                     Eqw                 = Flux / local_median
#                     Sigma_F, Sigma_Eqw  = self.line_classic_error(Flux, Eqw, sigma_continuum, D, Npix_line)
# 
#             self.fit_dict['eqw0'] = Eqw
#             self.fit_dict['flux_intg_er'],  self.fit_dict['eqw0_er'] = Sigma_F, Sigma_Eqw
            
        #--MC treatment
        #Single
        if Ncomps == 1:
            Eqw0 = ufloat(self.fit_dict['flux_intg'], self.fit_dict['flux_intg_er']) / ufloat(local_median, sigma_continuum)
            self.fit_dict['eqw0'], self.fit_dict['eqw0_er'] = Eqw0.nominal_value, Eqw0.std_dev
            
        #Blended line
        else:
            for i in range(Ncomps):
                index = str(i)
                Eqw_idx = ufloat(self.fit_dict['flux_gauss'+index], self.fit_dict['flux_gauss'+index+'_er']) / ufloat(local_median, sigma_continuum)
                self.fit_dict['eqw'+index], self.fit_dict['eqw'+index+'_er'] = Eqw_idx.nominal_value, Eqw_idx.std_dev
       
#         #--Simple treatement:
#         if 'MC' not in Methodology:
#             
#             #Single
#             if Ncomps == 1:
#                 Flux                    = self.fit_dict['flux_intg']
#                 Eqw                     = Flux / local_median
#                 Sigma_F, Sigma_Eqw      = self.line_classic_error(Flux, Eqw, sigma_continuum, D, Npix_line)
#                    
#             #Blended
#             else:
#                 for i in range(Ncomps):
#                     index               = str(i)
#                     Flux                = self.fit_dict['flux_gauss' + index]
#                     Eqw                 = Flux / local_median
#                     Sigma_F, Sigma_Eqw  = self.line_classic_error(Flux, Eqw, sigma_continuum, D, Npix_line)
# 
#             self.fit_dict['eqw0'] = Eqw
#             self.fit_dict['flux_intg_er'],  self.fit_dict['eqw0_er'] = Sigma_F, Sigma_Eqw
#             
#         #--MC treatment
#         else:
#             #Single
#             if Ncomps == 1:
#                 Eqw0 = ufloat(self.fit_dict['flux_intg'], self.fit_dict['flux_intg_er']) / ufloat(local_median, sigma_continuum)
#                 self.fit_dict['eqw0'], self.fit_dict['eqw0_er'] = Eqw0.nominal_value, Eqw0.std_dev
#                 
#             #Blended line
#             else:
#                 for i in range(Ncomps):
#                     index = str(i)
#                     Eqw_idx = ufloat(self.fit_dict['flux_gauss'+index], self.fit_dict['flux_gauss'+index+'_er']) / ufloat(local_median, sigma_continuum)
#                     self.fit_dict['eqw'+index], self.fit_dict['eqw'+index] = Eqw_idx.nominal_value, Eqw_idx.std_dev

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
       
    def blendedLines_Table(self):
                        
        first_lines_group   = ['O2_3726A', 'Ne3_3968A', 'O3_5007A', 'N2_6548A', 'O2_7319A']
        
        self.blended_lines_df = DataFrame(index=first_lines_group)
        self.blended_lines_df['line_labels']    = [['O2_3726A', 'O2_3729A'], ['Ne3_3968A', 'H1_3970A'], ['O3_5007A', 'He1_5016A'], ['N2_6548A', 'H1_6563A', 'N2_6584A'], ['O2_7319A', 'O2_7330A']]
        self.blended_lines_df['Theowavelength'] = [[3726.032, 3728.815], [3968.0, 3970.0], [5006.843, 5016], [6548.05, 6562.819, 6583.46], [7319.0, 7330.0]]
        self.blended_lines_df['Ion']            = [['[OII]','[OII]'], ['[NeIII]','HEps_7_2'],  ['[OIII]','HeI_5016'], ['[NII]','Halpha_3_2','[NII]'],  ['[OII]', '[OII]']]
        self.blended_lines_df['group_label']    = ['O2_3726A-O2_3729A', 'Ne3_3968A-H1_3970A', 'O3_5007A-He1_5016A', 'N2_6548A-H1_6563A-N2_6584A', 'O2_7319A-O2_7330A']
        return 

    def line_to_Logdf(self):
        
        #Single line
        if self.fit_dict.blended_check is False:
            line_label      = [self.Current_Label]
            Ion             = [self.Current_Ion]
            TheoWavelength  = [self.Current_TheoLoc]
        
        #Blended line
        else: 
            line_label      = self.fit_dict.blended_labels
            Ion             = self.fit_dict.blended_ions
            TheoWavelength  = self.fit_dict.blended_lambdas
        
        #Case we have a wide component
        if self.fit_dict.add_wide_component:
            line_label      = line_label + [line_label[1] + '_w']
            Ion             = Ion + [Ion[1] + '_w']
            TheoWavelength  = TheoWavelength + [self.fit_dict['maxLambdas'][1]] #Dirty trick so it is not saved in exactly the same place
            
        #Load the data to the dataframe:
        for i in range(self.fit_dict.line_number):
            idx_str = str(i)
            
            #Insert by sections
            parameters = ['Ion', 'lambda_theo', 'lambda_obs']
            self.current_df.loc[line_label[i], parameters] = Ion[i], TheoWavelength[i], self.fit_dict['maxLambdas'][i]
            
            parameters = ['flux_intg', 'flux_intg_er']
            self.current_df.loc[line_label[i], parameters] = self.fit_dict[parameters]  
                        
            parameters = ['flux_gauss{}', 'flux_gauss{}_er', 'eqw{}', 'eqw{}_er', 'A{}', 'A{}_er', 'mu{}', 'mu{}_er', 'sigma{}', 'sigma{}_er']
            self.current_df.loc[line_label[i], [param.format('') for param in parameters]] = self.fit_dict[[param.format(idx_str) for param in parameters]].values
                        
            parameters = ['zerolev_mean', 'zerolev_std', 'zerolev_width', 'm_zerolev', 'n_zerolev']
            self.current_df.loc[line_label[i], parameters] = self.fit_dict[parameters]  
    
            parameters = ['Wave1', 'Wave2', 'Wave3', 'Wave4', 'Wave5', 'Wave6']            
            self.current_df.loc[line_label[i], parameters] = self.fit_dict[parameters].values

            parameters = ['blended_check', 'line_number', 'group_label', 'add_wide_component', 'fit_routine']
            self.current_df.loc[line_label[i], parameters] = self.fit_dict[parameters]  
            
        return
    
    def load_lineslog_dataframe(self, lineslog_address):
        
        lines_df = read_csv(lineslog_address, delim_whitespace = True, header = 0, index_col = 0, comment='L') #Dirty trick to avoid the Line_label row
        
        return lines_df
    
    def save_lineslog_dataframe(self, linelog_df, lineslog_address, columns_to_store = None):
                        
        #In case we dont want to save all the columns
        if columns_to_store != None:
            linelog_df_to_save = linelog_df[columns_to_store]
        else:
            linelog_df_to_save = linelog_df
            
        string_frame = linelog_df_to_save.to_string() 
           
        #Convert to string format and save to text file
        with open(lineslog_address, 'wb') as output_file:
            output_file.write(string_frame.encode('UTF-8'))

        return
    
    def remove_lines_df(self, lineslog_df, line_label):

        #Check if blended line:
        blended_status = lineslog_df.loc[line_label, 'blended_check']        
        
        if bool_(blended_status) is False_:
            lineslog_df.drop(line_label, inplace = True)
        else:
            label_group         = lineslog_df.loc[line_label, 'group_label']
            idcs_group          = (lineslog_df.group_label == label_group)
            labels_to_delete    = lineslog_df.loc[idcs_group].index.values
            print 'Borrame estos'
            print labels_to_delete
            lineslog_df.drop(labels_to_delete, inplace = True)
                      
class LineMesurer_v2(Fitting_Gaussians_v2, LineMesurer_Log):
    
    def __init__(self, Conf_Folder, LinesLogHeader_Name):
                
        #Import the library with the fitting algorithms with gaussians
        Fitting_Gaussians_v2.__init__(self)
        LineMesurer_Log.__init__(self, Conf_Folder, LinesLogHeader_Name)
        
        #Line definition
        self.Current_Label          = None
        self.Current_Ion            = None
        self.Current_TheoLoc        = None

        #Textfile where we store our data
        self.Current_LinesLog       = None
                                      
        #Without a graphical interface we need a dataframe with the possible blended lines
        self.blendedLines_Table()
                
        #New design
        self.BoxSize = 70
        
        #Dictionary to store the data
        self.Parameter_dict = OrderedDict.fromkeys(self.ColumnHeaderVector)
        
        #Extra mechanisms for testing the wide component on Halpha      
        self.force_WD = True         
        
        #Parameters to save in log #WARNING: Missing emission line type
        self.saving_parameters_list = ['Ion',
                                  'lambda_theo',
                                  'lambda_obs',
                                  'flux_intg',
                                  'flux_intg_er',
                                  'flux_gauss', 
                                  'flux_gauss_er', 
                                  'eqw',
                                  'eqw_er',
                                  'A',
                                  'A_er',
                                  'mu',
                                  'mu_er',
                                  'sigma',
                                  'sigma_er',
                                  'zerolev_mean', #continuum linear value at line region center
                                  'zerolev_std',
                                  'zerolev_width',                                  
                                  'm_zerolev',
                                  'n_zerolev',                                  
                                  'Wave1',
                                  'Wave2',
                                  'Wave3',
                                  'Wave4',
                                  'Wave5',
                                  'Wave6',
                                  'blended_check',
                                  'line_number',
                                  'group_label',
                                  'add_wide_component',
                                  'fit_routine']
                
        #Clear the objects to store the fitting parameter for each line
        self.fit_dict = Series(index = self.fitting_parameters, dtype=object)
                             
    def measure_line(self, subWave, subFlux, wavelengths_list, lines_dataframe = None, line_label = None, Measuring_Method = 'lmfit', store_data = False):
        
        #Reset the dictionary to None
        self.fit_dict[self.fitting_parameters]  = None
               
        #Get indeces of line_regions in spectrum wavelength
        self.fit_dict.loc['idx0':'idx5']        = searchsorted(subWave, wavelengths_list)
        self.fit_dict.loc['Wave1':'Wave6']      = wavelengths_list
                
        #Calculation and plotting  parameters of adjacent continuum
        self.continuum_Regions(subWave, subFlux, adjacent_continuum=True)

        #Check if emission or absorption line as well as line mixture #WE SHOULD ONLY FEED HERE THE LINE REGION NOT ALL
        self.check_GaussianMixture(subWave, Measuring_Method, lines_dataframe, force_check=True)

        #Convert emission lines to gaussian equivalents
        self.calculate_Intensity(subWave, subFlux, lines_dataframe)
        
        #Store line data in the log                   
        if store_data and  self.fit_dict.start_treatment:
            if lines_dataframe is None:
                self.current_df = DataFrame(columns = self.saving_parameters_list)
            else:
                self.current_df = lines_dataframe
            
            self.line_to_Logdf()
        
        #Otherwise return the data from the fit
        else:
            return self.fit_dict

    def dazer_lineMeasuring(self, subWave, subFlux, wavelengths_list, lines_dataframe, Measuring_Method = 'lmfit'):
        
        #Clear the dictionaries with the line data
        self.fit_dict = Series(index = self.fitting_parameters, dtype=object)
         
        #Get indeces of line_regions in spectrum wavelength
        self.fit_dict.loc['idx0':'idx5'] = searchsorted(subWave, wavelengths_list)
        self.fit_dict.loc['Wave1':'Wave6'] = wavelengths_list
                
        #Calculation and plotting  parameters of adjacent continuum     
        self.continuum_Regions(subWave, subFlux)
                
        #Check if emission or absorption line as well as line mixture #WE SHOULD ONLY FEED HERE THE LINE REGION NOT ALL
        self.check_GaussianMixture(subWave, Measuring_Method, lines_dataframe, force_check=True)   

        #Convert emission lines to gaussian equivalents
        self.calculate_Intensity(subWave, subFlux, lines_dataframe)  
        
        #Save the fit:              
        if self.fit_dict.start_treatment:

            #Update the data frame with the new data
            self.line_to_Logdf()
            
            #Sort the data frame
            self.current_df.sort_values(['lambda_theo'], ascending=[True], inplace=True)
            
            #Store data to dataframe
            self.save_lineslog_dataframe(self.current_df, self.lineslog_df_address)
            
        return

    def continuum_Regions(self, subWave, subFlux, adjacent_continuum = True):
                
        #Indeces from the continuum regions
        idx1, idx2, idx3, idx4, idx5, idx6 = self.fit_dict.loc['idx0':'idx5'].values

        #In this case we use adjacent regions to compute the continuum level
        if adjacent_continuum == True:
            
            #Number of continua pixels
            pixBlue             = subFlux[idx1:idx2].shape[0]
            pixRed              = subFlux[idx5:idx6].shape[0]
            
            #Empty containers to store the continua regions
            FluxCont_BothSides  = empty(pixBlue + pixRed)
            WaveCont_BothSides  = empty(pixBlue + pixRed)
            FluxError_BothSides = empty(pixBlue + pixRed)
            FluxLine_BothSides  = empty(pixBlue + pixRed)
            
            #Combine blue and red continua into one array
            FluxCont_BothSides[0:pixBlue]  = subFlux[idx1:idx2]
            FluxCont_BothSides[pixBlue:]   = subFlux[idx5:idx6]
            WaveCont_BothSides[0:pixBlue]  = subWave[idx1:idx2]
            WaveCont_BothSides[pixBlue:]   = subWave[idx5:idx6]
            FluxError_BothSides[0:pixBlue] = std(subFlux[idx1:idx2]) * ones(pixBlue)
            FluxError_BothSides[pixBlue:]  = std(subFlux[idx5:idx6]) * ones(pixRed)

            #Region resolution
            region_resolution = (subWave[-1] - subWave[0]) / len(subWave)

            #Continua width in angstroms
            zerolev_width = region_resolution * (pixBlue + pixRed)
        
            #We perform a linear regresion taking into consideration the error in the flux #WARNING no se si esto esta bien
            m_Continuum, n_continuum = Python_linfit(WaveCont_BothSides, FluxCont_BothSides, FluxError_BothSides, errors_output = False)
            
            #Security check for very low density
            if isnan(m_Continuum):
                m_Continuum = 0.0
            if isnan(n_continuum):
                n_continuum = 0.0
            
            #Features for plotting
            Blue_wave_zerolev  = (subWave[idx2] + subWave[idx1]) / 2
            Red_wave_zerolev   = (subWave[idx6] + subWave[idx5]) / 2
            Blue_flux_zerolev  = m_Continuum * Blue_wave_zerolev + n_continuum
            Red_flux_zerolev   = m_Continuum * Red_wave_zerolev + n_continuum
   
            #We calculate the linear continuum as a vector for the whole central region and the middle point
            Linear_Flux = m_Continuum * subWave + n_continuum
    
            #Finally we calculate the dispersion in the continuum using the blue and red regions only
            FluxLine_BothSides[0:pixBlue] = Linear_Flux[idx1:idx2]
            FluxLine_BothSides[pixBlue:] = Linear_Flux[idx5:idx6]            
            Dif_FluxCont = FluxCont_BothSides - FluxLine_BothSides
            
            #Storing data into the fitting dictionary
            self.fit_dict['zerolev_mean']       = m_Continuum * ((subWave[idx3] + subWave[idx4]) / 2) + n_continuum
            self.fit_dict['zerolev_std']        = std(Dif_FluxCont) if std(Dif_FluxCont) != 0.0 else 1.0e-5
            self.fit_dict['m_zerolev']          = m_Continuum #Assuming a line this is the gradient (m)
            self.fit_dict['n_zerolev']          = n_continuum #Assuming a line this is the y axis interception point (n)
            self.fit_dict['zerolev_linear']     = Linear_Flux
            self.fit_dict['Blue_wave_zerolev']  = Blue_wave_zerolev
            self.fit_dict['Red_wave_zerolev']   = Red_wave_zerolev
            self.fit_dict['Blue_flux_zerolev']  = Blue_flux_zerolev
            self.fit_dict['Red_flux_zerolev']   = Red_flux_zerolev
            self.fit_dict['zerolev_width']      = zerolev_width
            
#             print 'zerolev_mean', self.fit_dict['zerolev_mean'], m_Continuum, n_continuum
#             print 'zerolev_std', self.fit_dict['zerolev_std']
#             print 'm_zerolev', self.fit_dict['m_zerolev']  
#             print 'n_zerolev', self.fit_dict['n_zerolev']
#             print 'zerolev_linear', self.fit_dict['zerolev_linear'] 
#             print 'Blue_wave_zerolev', self.fit_dict['Blue_wave_zerolev']
#             print 'Red_wave_zerolev', self.fit_dict['Red_wave_zerolev']
#             print 'Blue_flux_zerolev', self.fit_dict['Blue_flux_zerolev']
#             print 'Red_flux_zerolev', self.fit_dict['Red_flux_zerolev']                            
#             print 'zerolev_width', self.fit_dict['zerolev_width']                            
                                    
        else:

            #Storing data into the fitting dictionary
            self.fit_dict['zerolev_mean']       = (subFlux[idx3] + subFlux[idx4]) / 2 
            self.fit_dict['zerolev_std']        = 0.0
            self.fit_dict['m_zerolev']          = (subFlux[idx4] - subFlux[idx3]) / (subWave[idx4] - subWave[idx3])            
            self.fit_dict['n_zerolev']          = subFlux[idx3] - self.fit_dict['m_zerolev'] * subWave[idx3]
            self.fit_dict['zerolev_linear']     = self.fit_dict['m_zerolev'] * subWave + self.fit_dict['n_zerolev']
            self.fit_dict['Blue_wave_zerolev']  = 0.0
            self.fit_dict['Red_wave_zerolev']   = 0.0
            self.fit_dict['Blue_flux_zerolev']  = 0.0
            self.fit_dict['Red_flux_zerolev']   = 0.0
            self.fit_dict['zerolev_width']      = self.fit_dict['zerolev_linear'].shape[0]
                       
        return
              
    def check_GaussianMixture(self, SubWave, Methodology, lines_dataframe, force_check = True):
        
        #Store the procedure to use
        self.fit_dict['fit_routine'] = Methodology
        
        #Define dictionary initial conditions
        line_in_region                          = False
        lines_in_region                         = False
        first_blended                           = False
        in_blended_group                        = False
        
        self.fit_dict['MC_iterations']          = 1
        self.fit_dict['blended_check']          = False         #This logic is true when we have a blended line
        self.fit_dict['start_treatment']        = False         #This logic is true when a blended group is observed
        self.fit_dict['line_number']            = 0             #This integer describes the number of blended lines. Currently it must be equal to the number in the blended list elements.
        self.fit_dict['add_wide_component']     = False
        self.fit_dict['wide_component']         = False
                
        self.fit_dict['group_label']            = None
        self.fit_dict['blended_lambdas']        = None
        self.fit_dict['blended_labels']         = None
        self.fit_dict['blended_ions']           = None
        self.fit_dict['blended_index']          = 0
                
        #Select fitting type (Simple or bootstrap) 
        if 'MC' not in self.fit_dict['fit_routine']:
            self.fit_dict['MC_iterations'] = 1
        else:
            self.fit_dict['MC_iterations'] = 500
         
        #Check for observed wide component #WARNING this should come from the screen
        if self.force_WD:
            self.fit_dict['wide_component'] = True
        
        #If label is from a wide component do not measure it
        if self.Current_Label is not None:
            if '_w' in self.Current_Label:
                force_check = False
                
        #Automatic analysis for blended lines
        if force_check:

            #Check we are actually including the emission line intended
            if  SubWave[self.fit_dict.idx2] < self.Current_TheoLoc < SubWave[self.fit_dict.idx3]:  #This is a security check that should be at the top
                line_in_region = True
                self.fit_dict['line_number'] = 1
                
                #Check if this is the first line in a blended group
                if self.Current_Label in self.blended_lines_df.index:
                    first_blended = True
                    
                #Check if this line is in a blended group previously measured
                elif lines_dataframe is not None:
                    if self.Current_Label in lines_dataframe.index:
                        if (lines_dataframe.loc[self.Current_Label, 'group_label'] != 'None') and (lines_dataframe.loc[self.Current_Label, 'group_label'] is not None):
                            in_blended_group = True
            
                #Check all lines are located in the selected region
                if first_blended:
                    group_wavelengths = self.blended_lines_df.loc[self.Current_Label, 'Theowavelength']
                    
                    if (SubWave[self.fit_dict.idx2] < group_wavelengths[0]) and (group_wavelengths[-1] < SubWave[self.fit_dict.idx3]):
                        lines_in_region = True
                                    
            if (line_in_region and (first_blended is False) and (in_blended_group is False)) | ((line_in_region) and (lines_in_region is False) and in_blended_group is False):
                self.fit_dict['start_treatment']    = True
                self.fit_dict['line_number']        = 1
            
            #--Case with a blended line
            else:
                self.fit_dict['blended_check'] = True
                
                #First appereance of blended group
                if first_blended and lines_in_region:
                    self.fit_dict['start_treatment'] = True
                    
                    #WARNING in reality we should load this form the screen
                    self.fit_dict['group_label']        = self.blended_lines_df.loc[self.Current_Label, 'group_label']
                    self.fit_dict['blended_lambdas']    = self.blended_lines_df.loc[self.Current_Label, 'Theowavelength']
                    self.fit_dict['blended_labels']     = self.blended_lines_df.loc[self.Current_Label, 'line_labels']
                    self.fit_dict['blended_ions']       = self.blended_lines_df.loc[self.Current_Label, 'Ion']
                    self.fit_dict['line_number']        = len(self.fit_dict['blended_lambdas'])
                    #self.fit_dict['blended_index']      = 0

                #Case with a blended line previously measured
                elif in_blended_group:
                    self.fit_dict['start_treatment'] = False 

                    self.fit_dict['group_label'] = lines_dataframe.loc[self.Current_Label, 'group_label']
                    idces_group = (lines_dataframe.group_label == self.fit_dict['group_label'])
                                   
                    self.fit_dict['blended_labels']     = lines_dataframe.loc[idces_group].index.values
                    self.fit_dict['blended_lambdas']    = lines_dataframe.loc[idces_group, 'lambda_theo'].values
                    self.fit_dict['blended_ions']       = lines_dataframe.loc[idces_group, 'Ion'].values
                    self.fit_dict['line_number']        = len(self.fit_dict['blended_lambdas'])
                    #self.fit_dict['blended_index']      = self.fit_dict['blended_labels'].index(self.Current_Label)
                    
            #Dirty trick to always add a wide components for Halpha
            if (self.Current_TheoLoc == 6548.05)  and self.fit_dict['blended_check'] and self.fit_dict['start_treatment'] and self.fit_dict['wide_component']:
                self.fit_dict['add_wide_component'] = True          
        
#         print 'Bingo', self.fit_dict['blended_check'] 
           
        return
    
    def calculate_Intensity(self, subWave, subFlux, lineslog_df):
        
        #Load indeces for he regions
        idx3, idx4 = self.fit_dict.idx2, self.fit_dict.idx3
                
        # Check we are actually including the line intended #GET A GOOD ERROR MESSAGE FOR THE CASE THIS DOES HAPPEN                
        if self.fit_dict['start_treatment']:
            
            #Load lmfit parameters:
            self.load_lmfit_parameters(subWave[idx3:idx4], subFlux[idx3:idx4], self.fit_dict.zerolev_linear[idx3:idx4], self.fit_dict.zerolev_std, self.fit_dict.line_number, wide_component = False)
                   
            #Mesure single line 
            if self.fit_dict.blended_check is False:
                
                #One iteration
                if self.fit_dict.MC_iterations == 1:
                    self.fit_single_line(self.fit_dict.x_n, self.fit_dict.y_n, self.fit_dict.zerolev_n, self.fit_dict.sigZerolev_n, self.fit_dict.params_lmfit)
                
                #Several interations
                else:
                    self.fit_single_line_BS(self.fit_dict.x_n, self.fit_dict.y_n, self.fit_dict.zerolev_n, self.fit_dict.sigZerolev_n, self.fit_dict.params_lmfit, self.fit_dict['MC_iterations'])            
            
            #Mesure blended line
            else:
                self.fit_blended_line_BS_lmfit(self.fit_dict.x_n, self.fit_dict.y_n, self.fit_dict.zerolev_n, self.fit_dict.sigZerolev_n,  self.fit_dict.line_number, self.fit_dict.params_lmfit, 
                                      self.fit_dict.add_wide_component, self.fit_dict.params_lmfit_wide, bootstrap_iterations=self.fit_dict['MC_iterations'])

            #Scale line parameters to physical units
            self.rescale_lmfit_params(subWave[idx3:idx4], subFlux[idx3:idx4], self.fit_dict.line_number, self.fit_dict.x_scaler, self.fit_dict.y_scaler, self.fit_dict.fit_routine)
                  
            #Determine the physical parameters from the line
            self.linePhysicalData(subWave, subFlux, self.fit_dict.line_number, self.fit_dict.zerolev_std, self.fit_dict.fit_routine)
                                                             
        #In this case we get the data from the line data from the log
        elif (self.fit_dict['blended_check'] == True) and (self.fit_dict['start_treatment'] == False):
                        
            #Labels for the group
            idcs_group = (lineslog_df.group_label == self.fit_dict.group_label)

            self.fit_dict['flux_intg'], self.fit_dict['flux_intg_er'] = lineslog_df.loc[idcs_group, 'flux_intg'].values, lineslog_df.loc[idcs_group, 'flux_intg_er'].values
                       
            #Load the gaussian fit parameters of all the components 
            A_vector, A_er_vector           = lineslog_df.loc[idcs_group, 'A'].values, lineslog_df.loc[idcs_group, 'A_er'].values
            mu_vector, mu_er_vector         = lineslog_df.loc[idcs_group, 'mu'].values, lineslog_df.loc[idcs_group, 'mu_er'].values
            sigma_vector, sigma_er_vector   = lineslog_df.loc[idcs_group, 'sigma'].values, lineslog_df.loc[idcs_group, 'sigma_er'].values
            fluxG_vector, fluxG_er_vector   = lineslog_df.loc[idcs_group, 'flux_gauss'].values, lineslog_df.loc[idcs_group, 'flux_gauss_er'].values
            eqw_vector, eqw_er_vector       = lineslog_df.loc[idcs_group, 'eqw'].values, lineslog_df.loc[idcs_group, 'eqw_er'].values

            #Load the data into line dictionary #WARNING: Add the error to the line
            for i in range(len(mu_vector)):
                idx = str(i)
                self.fit_dict['A' + idx],           self.fit_dict['A' + idx + '_er']            = A_vector[i],      A_er_vector[i]        
                self.fit_dict['mu' + idx],          self.fit_dict['mu' + idx + '_er']           = mu_vector[i],     mu_er_vector[i]      
                self.fit_dict['sigma' + idx],       self.fit_dict['sigma' + idx + '_er']        = sigma_vector[i],  sigma_er_vector[i]    
                self.fit_dict['flux_gauss' + idx],  self.fit_dict['flux_gauss' + idx + '_er']   = fluxG_vector[i],  fluxG_er_vector[i]
                self.fit_dict['eqw' + idx],         self.fit_dict['eqw' + idx + '_er']          = eqw_vector[i],    eqw_er_vector[i]       
 
            #Calculate the gaussian curve for plotting
            self.fit_dict['x_resample'] = linspace(subWave[idx3], subWave[idx4], 50)
            self.fit_dict['zerolev_resample'] = self.fit_dict['m_zerolev'] * self.fit_dict['x_resample'] + self.fit_dict['n_zerolev']
            
            self.fit_dict['maxLambdas']  = mu_vector
            self.fit_dict['maxPeaks'] = A_vector + self.fit_dict['zerolev_mean']

            if self.fit_dict.blended_check == False:
                self.fit_dict['y_resample'] = gaussian_curve(self.fit_dict['A0'], self.fit_dict['mu0'], self.fit_dict['sigma0'], self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'])
            else:
                self.fit_dict['y_resample'] = gaussian_mixture(self.fit_dict, self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'], map(str, range(self.fit_dict.line_number)))
                self.fit_dict['y_comps']    = gaussian_components(self.fit_dict, self.fit_dict['x_resample'], self.fit_dict['zerolev_resample'], self.fit_dict.line_number)
        return
    
    def linePhysicalData(self, subWave, subInt, Ncomps, sigma_continuum, Methodology):
                
        #Load indeces for he regions
        idx1, idx2, idx3, idx4, idx5, idx6 = self.fit_dict.loc['idx0':'idx5'].values

        #This is necessary in case we must include the stellar continua:
        local_median = self.fit_dict['zerolev_mean']

        #Parameters for the classical error estimation
        Npix_line   = subWave[idx3:idx4].shape[0]                       #Number of pixels occupied by the lin                                
        D           = (subWave[-1] - subWave[0]) / subWave.shape[0] #Region dispertion
        
        #Continuum width calculation
        self.fit_dict['continuum_width'] = float(len(subInt[idx1:idx2]) + len(subInt[idx5:idx6])) * D    #Number of pixels at the sides 
            
        #--MC treatment
        #Single
        if Ncomps == 1:
            if local_median > 1e-50:
                Eqw0 = ufloat(self.fit_dict['flux_intg'], self.fit_dict['flux_intg_er']) / ufloat(local_median, sigma_continuum)
                self.fit_dict['eqw0'], self.fit_dict['eqw0_er'] = Eqw0.nominal_value, Eqw0.std_dev
            else:
                self.fit_dict['eqw0'], self.fit_dict['eqw0_er'] = 0.0, 0.0
            
        #Blended line
        else:
            for i in range(Ncomps):
                index = str(i)
                Eqw_idx = ufloat(self.fit_dict['flux_gauss'+index], self.fit_dict['flux_gauss'+index+'_er']) / ufloat(local_median, sigma_continuum)
                self.fit_dict['eqw'+index], self.fit_dict['eqw'+index+'_er'] = Eqw_idx.nominal_value, Eqw_idx.std_dev
       
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
       
    def blendedLines_Table(self):
                        
        first_lines_group   = ['O2_3726A', 'Ne3_3968A', 'O3_5007A', 'N2_6548A', 'O2_7319A']
        
        self.blended_lines_df = DataFrame(index=first_lines_group)
        self.blended_lines_df['line_labels']    = [['O2_3726A', 'O2_3729A'], ['Ne3_3968A', 'H1_3970A'], ['O3_5007A', 'He1_5016A'], ['N2_6548A', 'H1_6563A', 'N2_6584A'], ['O2_7319A', 'O2_7330A']]
        self.blended_lines_df['Theowavelength'] = [[3726.032, 3728.815], [3968.0, 3970.0], [5006.843, 5016], [6548.05, 6562.819, 6583.46], [7319.0, 7330.0]]
        self.blended_lines_df['Ion']            = [['[OII]','[OII]'], ['[NeIII]','HEps_7_2'],  ['[OIII]','HeI_5016'], ['[NII]','Halpha_3_2','[NII]'],  ['[OII]', '[OII]']]
        self.blended_lines_df['group_label']    = ['O2_3726A-O2_3729A', 'Ne3_3968A-H1_3970A', 'O3_5007A-He1_5016A', 'N2_6548A-H1_6563A-N2_6584A', 'O2_7319A-O2_7330A']
        return 

    def line_to_Logdf(self):
        
        #Single line
        if self.fit_dict.blended_check is False:
            line_label      = [self.Current_Label]
            Ion             = [self.Current_Ion]
            TheoWavelength  = [self.Current_TheoLoc]
        
        #Blended line
        else: 
            line_label      = self.fit_dict.blended_labels
            Ion             = self.fit_dict.blended_ions
            TheoWavelength  = self.fit_dict.blended_lambdas
        
        #Case we have a wide component
        if self.fit_dict.add_wide_component:
            line_label      = line_label + [line_label[1] + '_w']
            Ion             = Ion + [Ion[1] + '_w']
            TheoWavelength  = TheoWavelength + [self.fit_dict['maxLambdas'][1]] #Dirty trick so it is not saved in exactly the same place
            
        #Load the data to the dataframe:
        for i in range(self.fit_dict.line_number):
            idx_str = str(i)
            
            #Insert by sections
            parameters = ['Ion', 'lambda_theo', 'lambda_obs']
            self.current_df.loc[line_label[i], parameters] = Ion[i], TheoWavelength[i], self.fit_dict['maxLambdas'][i]
            
            parameters = ['flux_intg', 'flux_intg_er']
            self.current_df.loc[line_label[i], parameters] = self.fit_dict[parameters]  
                        
            parameters = ['flux_gauss{}', 'flux_gauss{}_er', 'eqw{}', 'eqw{}_er', 'A{}', 'A{}_er', 'mu{}', 'mu{}_er', 'sigma{}', 'sigma{}_er']
            self.current_df.loc[line_label[i], [param.format('') for param in parameters]] = self.fit_dict[[param.format(idx_str) for param in parameters]].values
                        
            parameters = ['zerolev_mean', 'zerolev_std', 'zerolev_width', 'm_zerolev', 'n_zerolev']
            self.current_df.loc[line_label[i], parameters] = self.fit_dict[parameters]  
    
            parameters = ['Wave1', 'Wave2', 'Wave3', 'Wave4', 'Wave5', 'Wave6']            
            self.current_df.loc[line_label[i], parameters] = self.fit_dict[parameters].values

            parameters = ['blended_check', 'line_number', 'group_label', 'add_wide_component', 'fit_routine']
            self.current_df.loc[line_label[i], parameters] = self.fit_dict[parameters]  
            
        return
    
    def load_lineslog_dataframe(self, lineslog_address):
        
        lines_df = read_csv(lineslog_address, delim_whitespace = True, header = 0, index_col = 0, comment='L') #Dirty trick to avoid the Line_label row
        
        return lines_df
    
    def save_lineslog_dataframe(self, linelog_df, lineslog_address, columns_to_store = None):
                        
        #In case we dont want to save all the columns
        if columns_to_store != None:
            linelog_df_to_save = linelog_df[columns_to_store]
        else:
            linelog_df_to_save = linelog_df
            
        string_frame = linelog_df_to_save.to_string() 
           
        #Convert to string format and save to text file
        with open(lineslog_address, 'wb') as output_file:
            output_file.write(string_frame.encode('UTF-8'))

        return
    
    def remove_lines_df(self, lineslog_df, line_label):

        #Check if blended line:
        blended_status = lineslog_df.loc[line_label, 'blended_check']        
        
        if bool_(blended_status) is False_:
            lineslog_df.drop(line_label, inplace = True)
        else:
            label_group         = lineslog_df.loc[line_label, 'group_label']
            idcs_group          = (lineslog_df.group_label == label_group)
            labels_to_delete    = lineslog_df.loc[idcs_group].index.values
            print 'Borrame estos'
            print labels_to_delete
            lineslog_df.drop(labels_to_delete, inplace = True)
                      
       