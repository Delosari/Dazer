from numpy                                      import max as npmax, min as npmin, median, ones, searchsorted, array_equal
from pandas                                     import DataFrame
from os.path                                    import isfile
from bin.DZ_Configuration                           import ImportConfiguration
from bin.DZ_LineMesurer                             import LineMesurer
from lib.Plotting_Libraries.dazer_plotter import Plot_Conf
from lib.CodeTools.File_Managing_Tools    import File_Manager

class Plots_Manager(ImportConfiguration, File_Manager, Plot_Conf, LineMesurer):
    
    def __init__(self):
        
        #Import the data catalogue and GUI conf
        ImportConfiguration.__init__(self)
        
        #Getting the data
        File_Manager.__init__(self)
        
        #For generating graphs
        Plot_Conf.__init__(self)
        
        #Tools to measure lines
        LineMesurer.__init__(self, self.Conf_Folder, self.LinesLogHeader_Name)
        
        #Spectra properties        
        self.Wave                       = None
        self.Int                        = None
        self.ExtraData                  = None
        self.Selections                 = []

        #This is a list with the lines indexes of the emission lines on the plot.
        self.EM_in_plot                 = []

        #This list contains the number of screens available
        self.ScreenCodes                = ['OBJECTS',  'TREATMENT',  'LINEMESURER']
        
        #This variable details the current screen (its value should be one from those in the self.ScreenCodes variable)
        self.CurrentScreen              = None
        
        #These indexes detail the current screen we should display
        self.Folder_i                   = None
        self.Treatment_i                = None
        self.EmLine_i                   = None
        
        #These are the names of the files which are obtained from the previous indexes and which are feed into the line mesurer
        self.Current_Folder             = None
        self.Current_Treatment          = None
        self.Current_Spec               = None
        self.Current_Code               = None
                                                
    def ManagePlots(self):
        
#-------Cases with introduction screens------------------------------------------------
        if (self.CurrentScreen == None):         
            Text = 'Initiating Visualizer'
            self.Axis.text(0.95, 0.05, Text, verticalalignment = 'bottom', horizontalalignment ='right', transform = self.Axis.transAxes, fontsize=15)
            
        elif (self.Folder_i == None) and (self.CurrentScreen == 'OBJECTS'):           
            Text = 'Analyzed objects'
            self.Axis.text(0.95, 0.05, Text, verticalalignment = 'bottom', horizontalalignment ='right', transform = self.Axis.transAxes, fontsize=15)
            
        elif (self.Treatment_i == None) and (self.CurrentScreen == 'TREATMENT'):           
            Text = 'Treatment steps'
            self.Axis.text(0.95, 0.05, Text, verticalalignment = 'bottom', horizontalalignment='right', transform = self.Axis.transAxes, fontsize=15)
            
        elif (self.EmLine_i == None) and (self.CurrentScreen == 'LINEMESURER'):
            Text = 'Emission line measurement'
            self.Axis.text(0.95, 0.05, Text, verticalalignment='bottom', horizontalalignment='right', transform=self.Axis.transAxes, fontsize=15)
            
            #At this point we load the data from the spectra (COULD WE MAKE IT LOAD DIRECTLY?)
            self.load_SpectraData()

#-------Cases with image screen------------------------------------------------
        elif self.CurrentScreen == 'OBJECTS':
            
            #Store the object data at each new key presssing            
            self.Current_Folder         = self.FilesList[0][self.Folder_i]
            self.Current_Code           = self.Current_Folder[self.Current_Folder[0:-1].rfind("/")+1:len(self.Current_Folder)-1]               
            self.Current_Spec           = self.SpectraPreffix + self.Current_Code + self.SpectraSuffix
            self.lineslog_df_address    = self.Current_Folder + self.Current_Code + '_WHT_linesLog_reduc.txt'
            
            #Load the lines dataframe if it exists:
            if isfile(self.lineslog_df_address):
                self.current_df = self.load_lineslog_frame(self.lineslog_df_address) #WARNING, there should be a better way to choose the file
            else:
                self.current_df = DataFrame(columns=self.saving_parameters_list)
            
            print '-- Displaying Object', self.Current_Code
            
            #Scaling the axes length to 1,1
            self.Axis.set_xlim(0,1)
            self.Axis.set_ylim(0,1)
            
            #Set Plot title
            self.Axis.set_title('Object ' + self.Current_Code, fontsize=35)   

            if isfile(self.Current_Folder + self.Current_Code + '.png'): 
                
                #Insert image of the plot
                self.insert_image(self.Current_Folder + self.Current_Code + '.png')
                                
            #Load data from the object log file
            self.Display_PhysicalData(self.Current_Folder, self.Current_Code)

#-------Case with plot screen------------------------------------------------
        elif (self.Treatment_i != None) and (self.CurrentScreen == 'TREATMENT'):
            
            #Store treatment to be displayed with each key pressing
            self.Current_Treatment = self.FilesList[1][self.Folder_i][self.Treatment_i]
            
            print '--- Displaying step', self.Current_Treatment
            
            #Load figure data
            #self.RestoreFig_v2(self.Current_Folder + self.Current_Treatment)

            self.restore_from_pickle(self.Current_Folder + self.Current_Treatment)
            #self.Fig.draw()


#-------Case for emission lines measurement screen------------------------------------------------
        elif (self.EmLine_i != None) and (self.CurrentScreen == 'LINEMESURER'):            
            
            #Determine current line
            self.Current_Label      = self.Labels_Total[self.EM_in_plot[self.EmLine_i]]
            self.Current_Ion        = self.Ions_Total[self.EM_in_plot[self.EmLine_i]]
            self.Current_TheoLoc    = self.Wavelengths_Total[self.EM_in_plot[self.EmLine_i]]
                       
            print '--- Displaying Object', self.Current_Label
                        
            #Check if line has been measure previously and in that case load data from text file
            if len(self.Selections) < 2:
                if self.Current_Label in self.current_df.index:
                    self.Selections = list(self.current_df.loc[self.Current_Label, 'Wave1':'Wave6'].values)
                           
            self.PlottingInLineMeasuring()
                                
        return

    def Display_PhysicalData(self, FileFolder, CodeName):
        
        #Display data if in screen if it is contained in the log:
        if isfile(FileFolder + CodeName + '_log.txt'):
        
            Log_Labels              = ['Obj_RA',  'Obj_DEC',  'Obj_z_blue'] + ['cHBeta_red', 'TOIII_pn', 'nSII_pn']   + ['OI_HI_pn','NI_HI_pn','SI_HI_pn','HeI_HI_from_O_pn']                                                                                                             +   ['Spectra_Meet']
            ScreenText_Parameters   = ['RA (Hours)', 'DEC(Deg)',  'z']      + [r'$c(H\beta)$', '$T_{e}$', r'$n_{e}$'] + [r'$\frac{O}{H}$', r'$\frac{N}{H}$', r'$\frac{S}{H}$',r'$\frac{He}{H}$'] + ['Arms match?']        
        
            #These coordinates assume x and y ranges going from 0 to 1
            y_min   = 0.05
            y_max   = 0.90
            y_step  = (y_min - y_max) / len(Log_Labels)
            x_label = 0.40
            x_Magni = 0.70
            
            for i in range(len(Log_Labels)):
                if ScreenText_Parameters[i] == 'RA (Hours)':
                    magnitude = self.GetParameter_ObjLog(CodeName, FileFolder, "Obj_RA")        #This needs to be updated to use pandas
                                
                elif ScreenText_Parameters[i] == 'DEC(Degrees)':
                    magnitude = self.GetParameter_ObjLog(CodeName, FileFolder, "Obj_DEC")       #This needs to be updated to use pandas                        

                else:
                    magnitude = self.GetParameter_ObjLog(CodeName, FileFolder, Log_Labels[i])   #This needs to be updated to use pandas                                        

                #Write the line
                y_text = y_max + i * y_step                    
                self.plot_text(x_label, y_text, ScreenText_Parameters[i])
                self.plot_text(x_Magni, y_text, magnitude, fontsize=16)    

        return

    def load_SpectraData(self):
        
        self.Wave, self.Int, self.ExtraData = self.get_spectra_data(self.Current_Folder + self.Current_Spec)
        
        Wmin = npmin(self.Wave)
        Wmax = npmax(self.Wave)
        
        #We erase the stored emision lines list
        del self.EM_in_plot[:]          
                                         
        for i in range(len(self.Wavelengths_Total)):
            if Wmin <= self.Wavelengths_Total[i] <= Wmax:
                    self.EM_in_plot.append(i)
                               
    def PlottingInLineMeasuring(self, store_data = True):
        
        # Plot main features-------------------------------------------------------------------
        n_selections = len(self.Selections)

        #Scale factors for the plot
        SubWave, SubInt, EM_Height, EM_ExpLoc = self.Emission_Threshold(self.Current_TheoLoc, self.Wave, self.Int)
        
        #Median flux in region    
        OverallMedian = median(SubInt)

        #Get neighboring lines
        minLim, maxLim = searchsorted(self.Wavelengths_Total, [SubWave[0], SubWave[-1]])
        neighbours_Lines = self.Wavelengths_Total[minLim:maxLim] 
        
        #Plot markers with near lines
        self.Axis.plot(neighbours_Lines, ones(len(neighbours_Lines)) * OverallMedian, 'o', color = self.colorVector['grey'])
        
        #Plotting current wavelength at different color
        self.Axis.plot(self.Current_TheoLoc, OverallMedian,'o', color = self.colorVector['olive'])            
        
        #Make a spec plot as a polyline or as a histogram               
        if self.HistogramTypeSpec:
            self.Axis.step(SubWave, SubInt, color = self.colorVector['silver'], where='mid')
        else:
            self.Axis.plot(SubWave, SubInt, color = self.colorVector['silver'])

        #Setting plot format
        self.Axis.set_xlim([SubWave[0],SubWave[-1]])
        self.Axis.set_ylim(OverallMedian/10, EM_Height*2)
        self.Axis.set_xlabel( r'Wavelength $(\AA)$')
        self.Axis.set_ylabel( r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})$')
        self.Axis.set_title('Measuring lines in ' +  self.Current_Code)   

        # Plot conditional features-------------------------------------------------------------------   
        if n_selections < 6: 
            
            #Box Text
            box_text = '{:s} emission {:0.2f} $\AA$'.format(self.Current_Ion, self.Current_TheoLoc)
            box_format = self.PlotFormat_Vector[0]
            self.Axis.text(0.65, 0.80, box_text, transform = self.Axis.transAxes, fontsize=18, verticalalignment='top', bbox = box_format)
            
            #Plot the region global median line
            self.Axis.axhline(y=OverallMedian, xmin=0, xmax=1, color=self.colorVector['grey'], linestyle = "--")           
        
        # Case with 1/3 regions    
        if n_selections == 2:
            self.ind1, self.ind2 = searchsorted(SubWave, self.Selections)
            self.Axis.fill_between(SubWave[self.ind1:self.ind2], OverallMedian, SubInt[self.ind1:self.ind2], facecolor = self.colorVector['green'], step='mid', alpha=0.5)
        
        # Case with 2/3 regions 
        if n_selections == 4:
            self.ind1, self.ind2, self.ind3, self.ind4 = searchsorted(SubWave, self.Selections)  
            self.Axis.fill_between(SubWave[self.ind1:self.ind2], OverallMedian, SubInt[self.ind1:self.ind2], facecolor =self.colorVector['green'],  step='mid', alpha=0.5)
            self.Axis.fill_between(SubWave[self.ind3:self.ind4], OverallMedian, SubInt[self.ind3:self.ind4], facecolor =self.colorVector['green'],  step='mid', alpha=0.5)
        
        # Case with 3/3 regions    
        if n_selections == 6: 
                        
            self.ind1, self.ind2,self.ind3, self.ind4, self.ind5, self.ind6 = searchsorted(SubWave, self.Selections)

            self.dazer_lineMeasuring(SubWave, SubInt, self.Selections, self.current_df, Measuring_Method = 'lmfit')
            
            self.Axis.scatter(self.fit_dict['maxLambdas'], self.fit_dict['maxPeaks'], color = self.colorVector['skin'])

            #------------Plotting line continuum elements----------------------
            self.Axis.plot(self.fit_dict.Blue_wave_zerolev, self.fit_dict.Blue_flux_zerolev, 'x', color=self.colorVector['orangish']) 
            self.Axis.plot(self.fit_dict.Red_wave_zerolev, self.fit_dict.Red_flux_zerolev, 'x', color=self.colorVector['orangish']) 
            self.Axis.plot(SubWave, self.fit_dict['zerolev_linear'] ,'--', color=self.colorVector['olive']) 
            self.Axis.axhline(y = self.fit_dict['zerolev_mean'], xmin=0, xmax=1, color =  self.colorVector['grey'], linestyle = "--")
            
            self.Axis.fill_between(SubWave[self.ind1:self.ind2], self.fit_dict['zerolev_linear'][self.ind1:self.ind2], SubInt[self.ind1:self.ind2],  step='mid', color=self.colorVector['green'], alpha=0.20)
            self.Axis.fill_between(SubWave[self.ind3:self.ind4], self.fit_dict['zerolev_linear'][self.ind3:self.ind4], SubInt[self.ind3:self.ind4],  step='mid', color=self.colorVector['cyan'],  alpha=0.40)
            self.Axis.fill_between(SubWave[self.ind5:self.ind6], self.fit_dict['zerolev_linear'][self.ind5:self.ind6], SubInt[self.ind5:self.ind6],  step='mid', color=self.colorVector['green'], alpha=0.20)

            #------------Plotting Gaussian components----------------------
            #Case of a single line
            if (self.fit_dict['line_number'] == 1) and (self.fit_dict['blended_check'] == False):
                
                self.Axis.plot(self.fit_dict['x_resample'], self.fit_dict['y_resample'], '--', color = self.colorVector['yellow'], linewidth=1.5)
                
                #Inserting text box          
                Text        =   "%s %0.2f $\AA$ (%0.2f$\AA$) \nFlux %0.5e \nEqW %0.2f $\AA$ \nContinuum %0.1f pixels"
                Text_Format =   (self.Current_Ion, self.Current_TheoLoc, max(SubInt[self.ind3:self.ind4]) - self.Current_TheoLoc, self.fit_dict['flux_intg'], self.fit_dict['eqw0'] * -1, self.fit_dict['zerolev_width'])
                self.Axis.text(0.65, 0.80, Text % Text_Format, transform=self.Axis.transAxes, fontsize=18, verticalalignment='top',bbox=self.PlotFormat_Vector[0])
            
            #Case with a blended line 
            elif (self.fit_dict['line_number'] == 1) and (self.fit_dict['blended_check'] == True):

                    #Constant text format for both cases
                    Text = "%s %0.2f $\AA$ (%0.2f$\AA$) \nFlux %0.5e \nEqW %0.2f $\AA$ \nContinuum %0.1f pixels"
    
                    #First component: We plot all global component
                    if (self.fit_dict['blended_check'] == True) and (self.fit_dict['start_treatment'] == True):
                        
                        #Plot Global component
                        self.Axis.plot(self.fit_dict['x_resample'], self.fit_dict['y_resample'], '-', color='yellow', linewidth=2.0)
                        
                        #Plot individual components
                        for i in range(1, self.fit_dict['line_number'] + 1):
                            self.Axis.plot(self.fit_dict['x_resample'], self.fit_dict['y_comps'], '--', color=self.colorVector['orangish'], linewidth=1.0)

                        Text_Format =   (self.Current_Ion, self.Current_TheoLoc, max(SubInt[self.ind3:self.ind4]) - self.Current_TheoLoc, self.fit_dict['flux_intg'], self.fit_dict['eqw0'], self.fit_dict['zerolev_width'])
    
                    #Subsequent compontents: We Plot all individual gaussians
                    else:
                        for i in range(self.fit_dict['line_number']):
                            index = str(i)
                            self.Axis.plot(self.fit_dict['x_resample'], self.fit_dict['y_comps'], '-', color=self.colorVector['orangish'], linewidth=2, linestyle='--')
                            Text_Format =   (self.Current_Ion, self.Current_TheoLoc, max(SubInt[self.ind3:self.ind4]) - self.Current_TheoLoc, self.Fitting_dict['flux_gauss' + index], self.fit_dict['eqw' + index] * -1, self.fit_dict['zerolev_width'])
                                        
                    self.Axis.text(0.65, 0.80, Text % Text_Format, transform=self.Axis.transAxes, fontsize=18, verticalalignment='top',bbox=self.PlotFormat_Vector[0])
                    
            else:
                Text = "%s %0.2f $\AA$ (%0.2f$\AA$) \nFlux %0.5e \nEqW %0.2f $\AA$ \nContinuum %0.1f pixels"
                for i in range(self.fit_dict['line_number']):
                    index = str(i)
                    self.Axis.plot(self.fit_dict['x_resample'], self.fit_dict['y_comps'][i], '-', color='red', linewidth=2, linestyle='--')
                    Text_Format = (self.Current_Ion, self.Current_TheoLoc, max(SubInt[self.ind3:self.ind4]) - self.Current_TheoLoc, self.fit_dict['flux_intg'], self.fit_dict['eqw' + index] * -1, self.fit_dict['zerolev_width'])
                    
            #------------Visual signals for special cases----------------------

            #Case with rare very small lines
            if len(SubInt[self.ind3:self.ind4])<4:
                self.Axis.plot(EM_ExpLoc, EM_Height, 'o', color='orange')
       
        return
    
    def Manage_DepthIndexes(self, Command):
            
        if Command == 'NextObject':
            self.CurrentScreen = self.Screen_keySwitcher(self.CurrentScreen, self.ScreenCodes, Command, Type = 'Word_Search')
        
        if Command == 'PreviousObject':
            self.CurrentScreen = self.Screen_keySwitcher(self.CurrentScreen, self.ScreenCodes, Command, Type = 'Word_Search')
            
            #WARNING: This loop is necessary to clear the position index of the screens below it. It can bring troubles we have more screens
            if self.CurrentScreen == 'OBJECTS':
                self.Treatment_i = None
                self.EmLine_i = None
                
            if self.CurrentScreen == 'TREATMENT':
                self.EmLine_i = None
                
        return
    
    def Manage_SubSectionIndexes(self, Command):
                
        if self.CurrentScreen == 'OBJECTS':  
                                                        #INDEX                 #LIST                                    #DIRECTION
            self.Folder_i = self.Screen_keySwitcher(self.Folder_i,             self.FilesList[0],                       Command)
            
        elif self.CurrentScreen == 'TREATMENT':
                                                        #INDEX                #LIST                                     #DIRECTION
            self.Treatment_i = self.Screen_keySwitcher(self.Treatment_i,      self.FilesList[1][self.Folder_i],         Command)
            
        elif self.CurrentScreen == 'LINEMESURER':
                                                        #INDEX                #LIST                                     #DIRECTION             
            self.EmLine_i = self.Screen_keySwitcher(self.EmLine_i,            self.EM_in_plot,                          Command)
            del self.Selections[:] #WE WANT TO RESET THIS VECTOR AFTER EACH KEY PRESSED... NO SURE IF BEST PLACE
            
        return 

    def Screen_keySwitcher(self, ScreenName, List, InputCommand, Next_Command = "NextObject", Previous_Command = "PreviousObject", Type = 'Num_Search'):
        
        #NOT THE CLEANEST ARRAY... MAYBE WE SHOULD MAKE TWO...
        if Type == 'Word_Search':
            if ScreenName == None:
                Current_Position = None
            else:
                Current_Position = List.index(ScreenName)
        
        else:
            Current_Position = ScreenName
                    
        #Warning the first value of Current_Position must be None
        if Current_Position == None:         # This works as a flag for the first position on the list 
            Index = 0
            Current_Position = 0
            
        else:
            
            if InputCommand == Next_Command:
                Index = 1
            
            if InputCommand == Previous_Command:
                Index = -1
                
            if (InputCommand != Previous_Command) and (InputCommand != Next_Command) :        
                Index = 0
                print "Key has not action assigned"
        
        NewPosition = Current_Position + Index
        List_Length = len(List)
    
        if NewPosition >= List_Length: 
            NewPosition = List_Length -1
            print "End of the list reached"
        
        if NewPosition < 0:
            NewPosition = 0
            print "First point of the list"

        if Type == 'Word_Search':
            Output =List[NewPosition]
        
        else:
            Output = NewPosition

        return Output
