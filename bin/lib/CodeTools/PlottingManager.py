from collections                    import OrderedDict, Sequence
from operator                       import itemgetter
from os                             import walk
from traceback                      import format_exc

from matplotlib._png                import read_png
from matplotlib.offsetbox           import OffsetImage, AnnotationBbox
from mpl_toolkits.mplot3d           import Axes3D
from numpy                          import array, sum, delete, char, genfromtxt, isnan, ndarray, arange, zeros, log10, where, loadtxt

from Astro_Libraries.AstroClass     import SlOutput_2_StarData
from Astro_Libraries.AstroMethods   import Fits2Data, StarlightFileManager
from CodeTools.File_Managing_Tools  import File_Manager
from CodeTools.various              import ColorsSet, TextFileChecker, FamilyOfItemsInArray
import matplotlib.pyplot            as plt

class Error_manager():
    
    def __init__(self):
        
        self.objec_error_dict   = OrderedDict()
        self.Lines_With_Errors  = []
        
    def log_error(self, object_name, verbose = True):
        
        self.objec_error_dict[object_name] = format_exc()
        
        if verbose:
            print '\n-- Error logged'
    
    def log_emLine_error(self, CodeName, Label):
        
        #Case for the first error log for a given object 
        if Label not in self.objec_error_dict:
            LinesList   = [[],[]]
        
        #Case the objects has more lines with an error
        else:
            LinesList   = self.objec_error_dict[CodeName] 

        LinesList[0].append([Label, format_exc()])
                
        self.objec_error_dict[CodeName] = LinesList
        
    def display_errors(self, type = 'general', extra_verbose = False):
        
        if type == 'general':
            
            if len(self.objec_error_dict) != 0:
                print '\n-These objects treatment produced a script error:\n'
                print self.objec_error_dict.keys(),'\n'
                for i in range(len(self.objec_error_dict)):
                    print '--Object', self.objec_error_dict.keys()[i], 
                    print self.objec_error_dict.values()[i], '\n'
                print 'hola'
                #We empty the dictionary as we don not need it anymore
                self.objec_error_dict.clear()
        
        elif type == 'emLine measurement':
            
            if len(self.objec_error_dict) != 0:
                print '\n-These objects treatment produced a script error:\n'
                print self.objec_error_dict.keys(),'\n'
                
                for i in range(len(self.objec_error_dict)):
                    print '--Object', self.objec_error_dict.keys()[i] 
                    print self.objec_error_dict.values()[i][0], '\n'
                    
                if extra_verbose == True:
                    for i in range(len(self.objec_error_dict)):
                        print '--Object',   self.objec_error_dict.keys()[i] 
                        print '--Label',    self.objec_error_dict.values()[i][0], '\n'
                        print '--Errors',   self.objec_error_dict.values()[i][1], '\n'  
                        
        return ' '
     
class myPickle(File_Manager, Error_manager):
    ''' Class employed to generate plots while saving the graph data and format into a text file '''
        
    def __init__(self):
        
        #Launch File_Manager
        File_Manager.__init__(self)
        Error_manager.__init__(self)
        
        #Files for managing the graph-to-text functionality
        self.Tags               = []                    #Terminal, file type, folder,
        self.FigConf            = None                  #ColorConf, FigWidth, FigHeight, AxisFormat, AxisFormat, TickLength, TickWidth, myLabelSize
        self.FileAddress        = []                    #FileAddress1, FileAddres2, FileAddress,
        self.DataAxis           = 1                     #AxisNumber1, AxisNumber2, AxisNumber3
        self.DataLineFormat     = []                    #[[LineStyle1, MarkerStyle1, LineColor1]; [LineStyle2, MarkerStyle2, LineColor2]; [LineStyle3, MarkerStyle3, LineColor3]]
        self.DataLineLabel      = []                    #[LineLabel1; LineLabel2; LineLabel3]
        self.DataXValues        = []                    #[(X1); (X2); (X3)]       
        self.DataYValues        = []                    #[(Y1); (Y2); (Y3)]
        self.DataXError         = []                    #[(ErrX1), (ErrX2), ErrX3)]
        self.DataYError         = []                    #[(ErrY1), (ErrY2), ErrY3)] 
        
        #Declaring the canvas
        self.FigCanvas          = None
   
        self.HistogramFormat    = []                    #[Label, Color, Transparency, Width, Fill, Edgecolor, Linestyle, Linewidth, Log])        
        self.HistogramData_X    = []
        self.HistogramData_Y    = []
   
        self.DataFilledFormat   = []                    #[DataLabel, FillColor, transparency, Fill_Style, BorderColor, Type]
        self.DataFilled_Initial = []
        self.DataFilled_Final   = []
             
        self.DataFilledBetweem_Format   = []            #[[Label, Where, FillColor, transparency, Fill_Style, BorderColor, Type]]
        self.DataFilledBetweem_Xrange   = []            #[[XRange1],[XRange2],[XRange3]]
        self.DataFilledBetweem_YLow     = []            #[[YLow1],[YLow2],[YLow3]]
        self.DataFilledBetweem_YHigh   = []             #[[Y_High1],[Y_High2],[Y_High]]
        
        self.PlotNotation       = []                    #[Plot_Title, Plot_xlabel, Plot_ylabel, XlabelSize, YLabelSize, TitleLabelSize, LegendLocation, Expand, XLabelPad, YLabelPad, Y_TitlePad]
        
        self.Text_X             = []
        self.Text_Y             = []
        self.Text_Values        = []
        self.TextFormat         = []                    #[[x_coords1, y_coords1, Text1, fontsize1, x_pad1, y_pad1, VerticalAligment1, PlottingType], [...,...,...]]
        
        self.ImageVector        = []
        
        self.FileSeparator      = '| '
        self.ParameterSeparator = '; '
        self.EmptyRowFormat     = 'nan'
        
        self.HeadersData        = []
        self.MatrixData         = []
        
        self.Fig                = None
        self.Color_Vector       = None
        self.Axis1              = None
        self.AxHor1             = None
        self.AxHor2             = None        
        
        self.MissingFiles       = []
             
    def FigFormat_One(self, ColorConf = 'Night', FigWidth = 12, FigHeight = 8, AxisFormat = 111, TickLength = 7, TickWidth = 2, myLabelSize = 14, StoreParameters = True):
        
        if AxisFormat != None:  

            self.Fig    = plt.figure(figsize = (FigWidth, FigHeight))  
            self.Axis1  = self.Fig.add_subplot(AxisFormat)
        
        self.Color_Vector       =  ColorsSet(ColorConf)
    
        self.Fig.set_facecolor(self.Color_Vector[0])
        self.Fig.set_edgecolor(self.Color_Vector[0])
        
        self.Axis1.set_facecolor(self.Color_Vector[0])
                
        self.Axis1.spines['bottom'].set_color(self.Color_Vector[1])
        self.Axis1.spines['top'].set_color(self.Color_Vector[1]) 
        self.Axis1.spines['right'].set_color(self.Color_Vector[1])
        self.Axis1.spines['left'].set_color(self.Color_Vector[1])
        
        self.Axis1.yaxis.label.set_color(self.Color_Vector[1])
        self.Axis1.xaxis.label.set_color(self.Color_Vector[1])    
        self.Axis1.yaxis.get_offset_text().set_color(self.Color_Vector[1])
        self.Axis1.xaxis.get_offset_text().set_color(self.Color_Vector[1])
    
        self.Axis1.tick_params(axis='x', length=TickLength, width=TickWidth, labelsize=myLabelSize, colors=self.Color_Vector[1])
        self.Axis1.tick_params(axis='y', length=TickLength, width=TickWidth, labelsize=myLabelSize, colors=self.Color_Vector[1])
        
        if StoreParameters == True:
            #                0             1          2            3            4        5            6     
            self.FigConf = [ColorConf, FigWidth, FigHeight, AxisFormat , TickLength, TickWidth, myLabelSize]
            self.Tags.append('SingleAxis')

    def GenerateFigure(self):
        
        self.Fig = plt.figure(figsize = (8, 12))  

    def DataPloter_One(self, X, Y, LineLabel, LineColor = None, LineStyle='-', LineWidth = 1, MarkerColor=None, MarkerStyle = None, MarkerSize  = None, XError = None, YError = None, ErrorBarsColor = None, ErrorBarsWidth = 1, ExtraData = None, StoreParameters = True):
        if LineColor == None:
            LineColor = self.Color_Vector[2][0]             #Clever but not clean improve
            
        if MarkerColor == None:
            MarkerColor = LineColor     
        
        if ErrorBarsColor == None:
            ErrorBarsColor = self.Color_Vector[1]
        
        if (YError == None) and (XError == None):           #Case without error bars
            #Plotting lines
            if LineStyle != None:
                Line = self.Axis1.plot(X, Y, linestyle=LineStyle, marker=MarkerStyle, color=LineColor, label=LineLabel, linewidth = LineWidth)
            
            #Plotting data points    
            else:
                if MarkerStyle == None:
                    MarkerStyle = 'o'
                if MarkerSize == None:
                    MarkerSize = 50
                Line = self.Axis1.scatter(X, Y, label=LineLabel, marker=MarkerStyle, color=MarkerColor, s=MarkerSize)

        else:                                               #Case with errorbars            
            if MarkerSize == None:
                MarkerSize = 10
            if MarkerStyle == None:
                MarkerStyle = 'o'
            Line = self.Axis1.errorbar(X, Y, label=LineLabel, xerr = XError, yerr = YError, ecolor=ErrorBarsColor, elinewidth = ErrorBarsWidth, fmt=MarkerStyle, ms=MarkerSize, c=MarkerColor) 
        
        if StoreParameters == True:
        
            self.DataLineLabel.append(LineLabel)
    #                                     1            2            3        4            5            6            7                8
            self.DataLineFormat.append([LineColor, LineStyle, LineWidth, MarkerColor, MarkerStyle, MarkerSize, ErrorBarsColor, ErrorBarsWidth])
            self.DataXValues.append(X)
            self.DataYValues.append(Y)
            self.DataXError.append(XError)
            self.DataYError.append(YError)
                          
        return

    def DataArea_One(self, Points_O, Points_F, DataLabel, FillColor = None, transparency = 1, Fill_Style = None, BorderColor = None, Type = 'infinite', StoreParameters = True):
        #IT CANNOT STORE LABELS WITH SPACES. We must move to a " | and ; " structure
        if FillColor == None:
            FillColor = self.Color_Vector[1]
        
        if Type == 'infinite':
            if isinstance(Points_O, Sequence) or isinstance(Points_O, ndarray):
                for k in range(len(Points_O)):
                    Filled = self.Axis1.axvspan(Points_O[k], Points_F[k], facecolor=FillColor, alpha=transparency, label=DataLabel)
            else:
                Filled = self.Axis1.axvspan(Points_O, Points_F, facecolor=FillColor, alpha=transparency, label=DataLabel)
 
        if StoreParameters == True:
            
            #                                0         1               2           3           4         5
            self.DataFilledFormat.append([DataLabel, FillColor, transparency, Fill_Style, BorderColor, Type])        
            self.DataFilled_Initial.append(Points_O)
            self.DataFilled_Final.append(Points_F)
                
        return
    
    def DataAreaFilled_One(self, X_Range, Y_Low, Y_High, Label, Where = None, FillColor = None, transparency = 1, Fill_Style = None, BorderColor = None, Type = 'between',  StoreParameters = True):
        
        #WHY IS THIS A PROBLEM ONLY FOR THIS ONE ??????
        if Where == 'None':
            Where = None
        
        #IT CANNOT STORE LABELS WITH SPACES. We must move to a " | and ; " structure
        if FillColor == None:
            FillColor = self.Color_Vector[1]
        
        if Type == 'between':
            if isinstance(X_Range, Sequence):
                for k in range(len(X_Range)):
                    Filled_between  = self.Axis1.fill_between(X_Range[k], Y_Low[k], Y_High[k], label=Label, facecolor=FillColor, alpha=transparency, where=Where)
            else:
                Filled_between      = self.Axis1.fill_between(X_Range, Y_Low, Y_High, label=Label, facecolor=FillColor, alpha=transparency, where=Where)
 
            if StoreParameters == True:
                #                                        0     1       2            3           4            5          6
                self.DataFilledBetweem_Format.append([Label, Where, FillColor, transparency, Fill_Style, BorderColor, Type])        
                self.DataFilledBetweem_Xrange.append(X_Range)
                self.DataFilledBetweem_YLow.append(Y_Low)
                self.DataFilledBetweem_YHigh.append(Y_High)
                               
        return    
        
    def TextPlotter(self, x_coords, y_coords, Text,  fontsize = 10, Color = None, x_pad = 1, y_pad = 1, VerticalAligment = 'top', PlottingType = 'Standard', StoreParameters = True):
    
        if Color == None:
            Color = self.Color_Vector[1]
        
        if PlottingType == 'Standard':
        
            if isinstance(Text, list):
                for i in range(len(Text)):
                    Label_point = Text[i]
                    self.Axis1.text(x_coords[i] * x_pad, y_coords[i] * y_pad, Label_point, fontsize=fontsize, color = Color)
                    
                #Some Tuppling
                x_coords    = x_coords           
                y_coords    = y_coords
                Text        = Text
                
            else:
                Label_point = Text                
                self.Axis1.text(x_coords * x_pad, y_coords * y_pad, Label_point, fontsize=fontsize, color = Color)
        
                #Some  more Tuppling
                x_coords    = [x_coords]           
                y_coords    = [y_coords]
                Text        = [Text]

        if PlottingType == 'fraction':        
            if isinstance(Text, list):          #NOT SURE IF IT IS ALWAYS LIST...
                for i in range(len(Text)):
                    Label_point = Text[i]
                    self.Axis1.text(x_coords[i] * x_pad, y_coords[i] * y_pad, Label_point,  fontsize=fontsize, color = Color, horizontalalignment='center', verticalalignment='center',transform=self.Axis1.transAxes)
                    
                #Some Tuppling
                x_coords    = x_coords           
                y_coords    = y_coords
                Text        = Text
                
            else:
                Label_point = Text                
                self.Axis1.text(x_coords * x_pad, y_coords * y_pad, Label_point,  fontsize=fontsize, color = Color, horizontalalignment='center', verticalalignment='center',transform=self.Axis1.transAxes)
        
                #Some  more Tuppling
                x_coords    = [x_coords]           
                y_coords    = [y_coords]
                Text        = [Text]
        
        if StoreParameters:
            self.TextFormat.append([fontsize, Color, x_pad, y_pad, VerticalAligment, PlottingType])
            self.Text_X.append(x_coords)
            self.Text_Y.append(y_coords)
            self.Text_Values.append(Text)

        return
        
    def CumulativeHistogram_One(self, Folder, File, PlotType, Color_Vector):
    
        #Vector where we save the histograms data
        OutPutData = []
        
        #Extract the data from the starlight ouuput
        index, x_j, Mini_j, Mcor_j, age_j, Z_j, LbyM, Mstars    = SlOutput_2_StarData(Folder, File)
        
        #Case light fraction plot
        if PlotType == 'Light_fraction':    
            ParameterVector = x_j
        
        #Case mass fraction plot    
        if PlotType == 'Mass_fraction':
            ParameterVector = Mcor_j
 
        #Extract metallicity values
        zValues = FamilyOfItemsInArray(Z_j)   
        
        #Define age bins width
#         AgeBins_HW = 0.20
        AgeBins_HW = 0.10
        AgeBins = arange(4.80, 10.60, AgeBins_HW)
        
        #Array to store intermediate data
        DataVector = []
                        
        for i in range(len(zValues)):
            z               = zValues[i]
            MassSumVector   = zeros(len(AgeBins))
            LightMassVector = zeros(len(AgeBins))
            for j in range(len(Z_j)):
                if Z_j[j] == z:
                    for k in range(len(AgeBins)):
                        Bin = AgeBins[k]
                        if (Bin - AgeBins_HW/2) <= log10(age_j[j]) < (Bin + AgeBins_HW/2):
                            MassSumVector[k] = MassSumVector[k] + ParameterVector[j]
    
            DataVector.append([z, Color_Vector[2][i], MassSumVector, LightMassVector])
            
        for i in range(len(AgeBins)):
            Total = 0  
            ListMasses = []
            for j in range(len(zValues)):
                ListMasses.append(DataVector[j][2][i])
                Total = Total +DataVector[j][2][i] 
                          
            Indexes_Small2Big = [n[0] for n in sorted(enumerate(ListMasses), key=lambda x:x[1])]
            Indexes_Small2Big.reverse()
     
            for Index in Indexes_Small2Big:
                if DataVector[Index][2][i] > 0:
#                         self.Histogram_One([AgeBins[i]- AgeBins_HW/4], [DataVector[Index][2][i]], str(DataVector[Index][0]), DataVector[Index][1], 1.0, Width = AgeBins_HW/2, Edgecolor = DataVector[Index][1], Linewidth = 1)

                    X_Coord         = [AgeBins[i]- AgeBins_HW/4]
                    Y_Coord         = [DataVector[Index][2][i]]
                    Label           = str(DataVector[Index][0])
                    Color           = DataVector[Index][1]
                    Transparency    = 1.0
                    Width           = AgeBins_HW/2,
                    Fill            = True
                    Edgecolor       = DataVector[Index][1]
                    Linestyle       = 'solid'
                    Linewidth       = 1
                    Log             = False
                    
                    OutPutData.append([X_Coord, Y_Coord, Label, Color, Transparency, Width, Fill, Edgecolor, Linestyle, Linewidth, Log])
                    
#                         self.Histogram_One(OutPutData[-1][0], OutPutData[-1][1], OutPutData[-1][2], OutPutData[-1][3], OutPutData[-1][4], OutPutData[-1][5], OutPutData[-1][6], OutPutData[-1][7], OutPutData[-1][8], OutPutData[-1][9], OutPutData[-1][10])

                    
            if Total > 0:
                self.Histogram_One([AgeBins[i]- AgeBins_HW/2], [Total], 'Bin total', Color_Vector[0], Width = AgeBins_HW, Fill = False, Edgecolor = self.Color_Vector[1], Linestyle='dotted')
    
                X_Coord         = [AgeBins[i]- AgeBins_HW/2]
                Y_Coord         = [Total]
                Label           =  'Bin total'
                Color           = Color_Vector[0]
                Transparency    = 1
                Width           = AgeBins_HW
                Fill            = False
                Edgecolor       = Color_Vector[1]
                Linestyle       = 'dotted'
                Linewidth       = 1
                Log             = False
                
                OutPutData.append([X_Coord, Y_Coord, Label, Color, Transparency, Width, Fill, Edgecolor, Linestyle, Linewidth, Log])

#                     self.Histogram_One(OutPutData[-1][0], OutPutData[-1][1], OutPutData[-1][2], OutPutData[-1][3], OutPutData[-1][4], OutPutData[-1][5], OutPutData[-1][6], OutPutData[-1][7], OutPutData[-1][8], OutPutData[-1][9], OutPutData[-1][10])        
#                     self.Histogram_One([AgeBins[i]- AgeBins_HW/2], [Total], 'Bin_total', self.Color_Vector[0], Width = AgeBins_HW, Fill = False, Edgecolor = self.Color_Vector[1], Linestyle='dotted')        
        
        return OutPutData
    
    def Histogram_One(self, X_Coord, Y_Coord, Label, Color = None, Transparency = 1.0, Width=0.8, Fill = True, Edgecolor = None, Linestyle='solid', Linewidth=1, Log=False, StoreParameters = True):
        
#         print 'Los valores', X_Coord, Y_Coord, Label, Color, Transparency, Width, Fill, Edgecolor, Linestyle, Linewidth, Log, StoreParameters
        Bar = self.Axis1.bar(X_Coord, Y_Coord, label=Label, color=Color, width=Width, alpha=Transparency, fill = Fill, edgecolor=Edgecolor, linestyle = Linestyle, linewidth = Linewidth, log = Log)
        
        if StoreParameters == True:

            self.HistogramFormat.append([Label, Color, Transparency, Width, Fill, Edgecolor, Linestyle, Linewidth, Log])       
            self.HistogramData_X.append(X_Coord)
            self.HistogramData_Y.append(Y_Coord)
        
        return
    
    def File2Data(self, FileFolder, FileName):

        #Fits files     
        if (".fit"  in FileName) or ('.fits' in FileName):
            Int, Wave, ExtraData = Fits2Data(FileFolder,FileName)
            
            # Clumsy need to clear this
            if type(Int[0]) == type(Wave):                                      
                Y = Int[0]
                X = Wave
            else:
                X, Y = Wave, Int
            
            return X, Y, ExtraData
          
        #Textfile         
        elif TextFileChecker(FileFolder,FileName):
            if ".slOutput" in FileName:
                Input_Wavelength, Input_Flux, Output_Flux, MaskPixels, ClippedPixels, FlagPixels, Parameters = StarlightFileManager(FileFolder, FileName)
                return Input_Wavelength, Input_Flux, Output_Flux, MaskPixels, ClippedPixels, FlagPixels, Parameters
            
            elif '.txt' in FileName:
                x, y        = self.get_ColumnData([0, 1], FileFolder + FileName, HeaderSize = 0, StringIndexes = False, unpack_check = True)
                ExtraData   = None
            
            else:
                
                #This method is obsolete... try changing this with loadtxt
                #x, y, ExtraData = TextManager(FileFolder, FileName)
                
                x, y, Extradata = None, None, None
            
            return x, y, ExtraData
        
    def Labels_Legends_One(self, Plot_Title=None, Plot_xlabel=None, Plot_ylabel=None, XlabelSize=20, YLabelSize=20, TitleLabelSize= 25, LegendLocation = 1, Expand= False,  XLabelPad = 0, YLabelPad = 0, Y_TitlePad = 1.02,  StoreParameters = True):
        
        if self.Axis1.get_legend_handles_labels()[1] != None:
            
            if  (self.AxHor1 == None) and (self.AxHor2 == None):

                Old_Handles, Old_Labels = self.Axis1.get_legend_handles_labels()
                Handles_by_Label        = OrderedDict(zip(Old_Labels, Old_Handles))
                
                Hl = sorted(zip(Handles_by_Label.values(), Handles_by_Label.keys()), key=itemgetter(1))
                        
                New_Handles, New_labels = zip(*Hl)        
                
                if Expand:
                    myLegend = self.Axis1.legend(New_Handles, New_labels, loc=LegendLocation, prop={'size':20}, ncol = len(New_labels)/2+1, scatterpoints=1, numpoints=1)   
                else:
                    myLegend = self.Axis1.legend(New_Handles, New_labels, loc=LegendLocation, prop={'size':20}, scatterpoints=1, numpoints=1) 
    
    
            else: 
                
                Old_Handles1, Old_Labels1 = self.Axis1.get_legend_handles_labels()
                Old_Handles2, Old_Labels2 = self.AxHor1.get_legend_handles_labels()
                Old_Handles3, Old_Labels3 = self.AxHor2.get_legend_handles_labels()
    
                Total_Handles   = Old_Handles1 + Old_Handles2 + Old_Handles3
                Total_labels    = Old_Labels1 + Old_Labels2 + Old_Labels3
    
                Handles_by_Label            = OrderedDict(zip(Total_labels, Total_Handles))
    
                Hl = zip(Handles_by_Label.values(), Handles_by_Label.keys())

                New_Handles, New_labels = zip(*Hl)        
                
                if Expand:
                    myLegend = self.Axis1.legend(New_Handles, New_labels, loc=LegendLocation, prop={'size':20}, ncol = len(New_labels)/2+1, scatterpoints=1, numpoints=1)   
                else:
                    myLegend = self.Axis1.legend(New_Handles, New_labels, loc=LegendLocation, prop={'size':20}, scatterpoints=1, numpoints=1) 
    
            Leg_Frame = myLegend.get_frame()
            Leg_Frame.set_facecolor(self.Color_Vector[0])
            Leg_Frame.set_edgecolor(self.Color_Vector[1])
        
            for label in myLegend.get_texts():
                label.set_fontsize('large')
         
            for label in myLegend.get_lines():
                label.set_linewidth(1.5)
             
            for text in myLegend.get_texts():
                text.set_color(self.Color_Vector[1])
            
            #WARNING!= THIS DOES NOT WORK ANYMORE
            if (Plot_xlabel == None) and (Plot_ylabel == None): 
                for i in range(len(self.FileAddress)):
                    if (".fits" or '.fit') in self.FileAddress[i]:
                        Plot_xlabel, Plot_ylabel =  r'Wavelength $(\AA)$', 'Flux' + r'$(erg\,cm^{-2} s^{-1} \AA^{-1})$' 

        self.Axis1.set_xlabel(Plot_xlabel,  fontsize = XlabelSize,      color = self.Color_Vector[1])
        self.Axis1.set_ylabel(Plot_ylabel,  fontsize = YLabelSize,      color = self.Color_Vector[1])
        self.Axis1.set_title(Plot_Title,    fontsize = TitleLabelSize,  color = self.Color_Vector[1], y = Y_TitlePad)    

        self.Axis1.xaxis.labelpad = XLabelPad
        self.Axis1.yaxis.labelpad = YLabelPad
        
        if  StoreParameters == True:
            #                        0            1            2            3            4        5                6            7            8        9            10    
            self.PlotNotation = [Plot_Title, Plot_xlabel, Plot_ylabel, XlabelSize, YLabelSize, TitleLabelSize, LegendLocation, Expand, XLabelPad, YLabelPad, Y_TitlePad]
        
        
        # LEGEND String locations
        # upper right        1
        # upper left         2
        # lower left         3
        # lower right        4
        # center right       5
        # center left        6
        # center right       7
        # lower center       8    
        # upper center       9
        # center            10
            
        return
    
    def ClearFigure(self):
        plt.cla()
        
    def SaveManager(self, SavingName = None, SavingFolder="/home/vital/Desktop/", SavingExtension='.png', ForceSave=False, ForceDisplay=None, savevectorfile = True):  
        
        if (ForceSave == True) or ('SingleFile' not in self.Tags):
           
            if SavingName == None:
                SavingName = " ".join(self.FileAddress)
            
            
            Output_Address = SavingFolder + SavingName + SavingExtension
        
               
            plt.savefig(Output_Address, dpi=300, facecolor=self.Color_Vector[0], bbox_inches='tight', pad_inches=0.2)
        
            print '--'+SavingName + SavingExtension, 'saved at', Output_Address
            
            if savevectorfile:
            
                self.SaveVectorFile(Output_Address.replace(SavingExtension,'.plot'))
                
        if ('Terminal' in self.Tags) or (ForceDisplay==True):
            plt.show()  
        
        elif (('SingleFile' in self.Tags) and ForceDisplay != False):
            plt.show()  
                        
        return
      
    def Conf_2_Line(self, single_check, line_name, line_parameters, plot_textfile):
        
        #THIS ONE SHOULD BE INCLUDED AS PART OF MY BASIC METHODS LIBRARY        
        if single_check:
            plot_textfile.write(line_name + self.ParameterSeparator.join(map(str, line_parameters)) + '\n')
        
        else:
            for i in range(len(line_parameters)):
                Component = self.ParameterSeparator.join(map(str, line_parameters[i])).replace(', ', ',') 
                line_parameters[i] = Component
            plot_textfile.write(line_name + self.FileSeparator.join(map(str, line_parameters)) + '\n')
            
        return
    
    def SaveVectorFile(self, SavingAddress, StringFormat = '{:>20}', ScientificNotation = '{:>20.5e}', ScientificNotationOld = '%20.5e'):
        
        TextFile = open(SavingAddress,"w")
        
        #---------------------Saving configuration of plot--------------------------------------------------
                
        self.Conf_2_Line(True, '-Tags:\n',          self.Tags,          TextFile) #0 Tag line            1 Index
        self.Conf_2_Line(True, '-FigConf:\n',       self.FigConf,       TextFile) #1 FigConf line        3 Index
        self.Conf_2_Line(True, '-Files:\n',         self.FileAddress,   TextFile) #2 FileAddress line    5 Index
        self.Conf_2_Line(True, '-Plot Notation:\n', self.PlotNotation,  TextFile) #3 PlotNotation line   7 Index

        #---------------------Saving data from the plot--------------------------------------------------
        LenMax          = 0
        LengthVector    = []    #It would be more efficient for this one to be a numpy array
        Data_Format     = []    
                
        #LINE DATA        
        if len(self.DataLineFormat) > 0:
            self.Conf_2_Line(True,  '-Data Labels:\n',   self.DataLineLabel, TextFile) #4 Labels line         9 Index
            self.Conf_2_Line(False, '-Data Plot Conf:\n', self.DataLineFormat, TextFile)
 
            for i in range(len(self.DataLineLabel)):
                ArrayLength = len(self.DataXValues[i])
                if ArrayLength > LenMax:
                    LenMax = ArrayLength
                    
                XColumn     = "X_" + str(i)
                self.HeadersData.append(XColumn)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.DataXValues[i])
                LengthVector.append(len(self.DataXValues[i]))
                
                YColumn = 'Y_' + str(i)
                self.HeadersData.append(YColumn)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.DataYValues[i])
                LengthVector.append(len(self.DataYValues[i]))
                
                if self.DataXError[i] != None:
                    XError_Column =  'XE_' + str(i)
                    self.HeadersData.append(XError_Column)
                    Data_Format.append(ScientificNotation)
                    self.MatrixData.append(self.DataXError[i])
                    LengthVector.append(len(self.DataXError[i]))
                    
                if self.DataYError[i] != None:
                    YError_Column =  'YE_' + str(i)
                    self.HeadersData.append(YError_Column)
                    Data_Format.append(ScientificNotation)            
                    self.MatrixData.append(self.DataYError[i])
                    LengthVector.append(len(self.DataYError[i]))
        
        #IMAGES DATA
        if len(self.ImageVector) > 0:            
            self.Conf_2_Line(False, '-Inserted images format:\n', self.ImageVector, TextFile)

        #FILLED REGIONS DATA
        if len(self.DataFilledFormat) > 0:
            self.Conf_2_Line(False, '-Filled data format:\n', self.DataFilledFormat, TextFile)
            
            for i in range(len(self.DataFilledFormat)):
                ArrayLength = len(self.DataFilled_Initial[i])
                if ArrayLength > LenMax:                    
                    LenMax = ArrayLength
                        
                XColumn     = "FO_" + str(i)
                self.HeadersData.append(XColumn)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.DataFilled_Initial[i])
                LengthVector.append(len(self.DataFilled_Initial[i]))
                
                YColumn = 'FF_' + str(i)
                self.HeadersData.append(YColumn)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.DataFilled_Final[i])
                LengthVector.append(len(self.DataFilled_Final[i]))
            
        #FILLED BETWEEN REGIONS DATA
        if len(self.DataFilledBetweem_Format) > 0:
            self.Conf_2_Line(False, '-Filled between data format:\n', self.DataFilledBetweem_Format, TextFile)
            
            for i in range(len(self.DataFilledBetweem_Format)):
                ArrayLength = len(self.DataFilledBetweem_Xrange[i])       
                if ArrayLength > LenMax:                    
                    LenMax = ArrayLength       
 
                XColumn     = "FBX_" + str(i)
                self.HeadersData.append(XColumn)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.DataFilledBetweem_Xrange[i])
                LengthVector.append(len(self.DataFilledBetweem_Xrange[i]))        
      
                YColumn_1     = "FBO_" + str(i)
                self.HeadersData.append(YColumn_1)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.DataFilledBetweem_YLow[i])
                LengthVector.append(len(self.DataFilledBetweem_YLow[i]))             

                YColumn_2     = "FBF_" + str(i)
                self.HeadersData.append(YColumn_2)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.DataFilledBetweem_YHigh[i])
                LengthVector.append(len(self.DataFilledBetweem_YHigh[i]))       
       
        #HISTOGRAMS DATA
        if len(self.HistogramFormat) > 0:
            self.Conf_2_Line(False, '-Histogram data format:\n', self.HistogramFormat, TextFile)
            
            for i in range(len(self.HistogramFormat)):
                ArrayLength = len(self.HistogramData_X[i])       
                if ArrayLength > LenMax:                    
                    LenMax = ArrayLength    

                XColumn     = "HX_" + str(i)
                self.HeadersData.append(XColumn)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.HistogramData_X[i])
                LengthVector.append(len(self.HistogramData_X[i]))        
      
                YColumn_1     = "HY_" + str(i)
                self.HeadersData.append(YColumn_1)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.HistogramData_Y[i])
                LengthVector.append(len(self.HistogramData_Y[i]))              
        
        # PLOTED TEXT DATA (The coordinates and text go into columns with the scientific data)         
        if len(self.TextFormat) > 0:
            self.Conf_2_Line(False, '-Text in plot:\n', self.TextFormat, TextFile)
            
            for i in range(len(self.TextFormat)):
                ArrayLength = len(self.TextFormat[i][0])       
                if ArrayLength > LenMax:                    
                    LenMax = ArrayLength       
 
                XColumn         = "TX_" + str(i)
                self.HeadersData.append(XColumn)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.Text_X[i])
                LengthVector.append(len(self.Text_X[i]))        
      
                YColumn_1       = "TY_" + str(i)
                self.HeadersData.append(YColumn_1)
                Data_Format.append(ScientificNotation)
                self.MatrixData.append(self.Text_Y[i])
                LengthVector.append(len(self.Text_Y[i]))             

                YColumn_2       = "TV_" + str(i)
                self.HeadersData.append(YColumn_2)
                Data_Format.append(StringFormat)
                self.MatrixData.append(self.Text_Values[i])
                LengthVector.append(len(self.Text_Values[i]))  
           
        #STORE SCIENTIFIC DATA
        TitleFormat = []
        for Header in self.HeadersData:
            TitleFormat.append(StringFormat)
        
        TextFile.write('-Data Matrix:\n')       #Headers cannot contain spaced strings
        TextFile.write(' '.join(TitleFormat).format(*tuple(self.HeadersData)) + '\n')

        #Columns with data
        for i in range(LenMax): 
            LineFormat  = []
            LineValue   = []  
            for k in range(len(self.MatrixData)):
                if i < LengthVector[k]:
                    LineFormat.append(Data_Format[k])
                    LineValue.append(self.MatrixData[k][i])
                else:
                    LineFormat.append(StringFormat)
                    LineValue.append(self.EmptyRowFormat)            

            TextFile.write(' '.join(LineFormat).format(*tuple(LineValue)) + '\n')

        TextFile.close()
             
        return
            
    def ResetPlot(self): 
        plt.cla()
        
        del self.HeadersData[:]
        del self.MatrixData[:]         
        
        del self.FileAddress[:]
        del self.PlotNotation[:]
        
        del self.DataLineLabel[:]
        del self.DataLineFormat[:]
        del self.DataXValues[:]
        del self.DataYValues[:]
        del self.DataXError[:]
        del self.DataYError[:]
        
        del self.DataFilledFormat[:]        
        del self.DataFilled_Initial[:] 
        del self.DataFilled_Final[:]
        
        del self.DataFilledBetweem_Format[:]  
        del self.DataFilledBetweem_Xrange[:] 
        del self.DataFilledBetweem_YLow[:]    
        del self.DataFilledBetweem_YHigh[:]         
        
        del self.HistogramFormat[:]
        del self.HistogramData_X[:]
        del self.HistogramData_Y[:]
        
        del self.ImageVector[:]
        
        del self.Text_X[:]
        del self.Text_Y[:]
        del self.Text_Values[:]
        del self.TextFormat[:]
        
    def DisplayFigure(self):
        plt.show()
        
    def DeleteAxis(self):
        self.Axis1.clear()
    
    def InsertFigure(self, ImageFolder, ImageName, Image_Coordinates = [0.25,0.5], Zoom = 1, Image_Box=(20,-20), Image_xyCoords = 'data', Image_BoxCoords = "offset points", StoreParameters = True):
        
        #THIS MODULE IS NOT SAVED INTO THE pyRestore
        
        arr_hand = read_png(ImageFolder + ImageName)
        Image_Frame = OffsetImage(arr_hand, zoom=Zoom)
        ab = AnnotationBbox(Image_Frame, Image_Coordinates,
            xybox=Image_Box,
            xycoords=Image_xyCoords,
            boxcoords=Image_BoxCoords)
        self.Axis1.add_artist(ab)
        
        if StoreParameters == True:
            #WARNING: WE ARE NOT SAVING THE IMAGE_BOXCOORDS..., WE MUST IMPLEMENT THE (; |) TEXT FILE DESIGN TO INSERT STRINGS WITH SPACES 
            self.ImageVector.append([ImageFolder, ImageName, Image_Coordinates, Zoom, Image_Box, Image_xyCoords])
        
        return    

    def Triple_Axis(self):
        
        self.AxHor1 = self.Axis1.twiny()
        self.AxHor2 = self.AxHor1.twiny()
        
        self.Axis1.ticklabel_format(style='sci', axis='x', scilimits=(5, 5))
        self.Axis1.tick_params(axis='both', labelsize=20)
        
        self.AxHor1.tick_params(axis='x', length=7, width=2, labelsize=20, colors=self.Color_Vector[1])
        self.AxHor1.set_frame_on(True)
        self.AxHor1.patch.set_visible(False)
        self.AxHor1.xaxis.set_ticks_position('bottom')
        #AxHor1.xaxis.set_ticks_color(Pv.Color_Vector[1])
        
        self.AxHor1.xaxis.set_label_position('bottom')
        self.AxHor1.spines['bottom'].set_position(('outward', 45))
        self.AxHor1.spines['bottom'].set_color(self.Color_Vector[1])
        self.AxHor1.ticklabel_format(style='sci')
        self.AxHor1.ticklabel_format(style='sci', axis='x', scilimits=(5, 5))
        self.AxHor1.xaxis.get_offset_text().set_color(self.Color_Vector[1])
        
        self.AxHor2.tick_params(axis='x', length=7, width=2, labelsize=20, colors=self.Color_Vector[1])
        self.AxHor2.set_frame_on(True)
        self.AxHor2.patch.set_visible(False)
        self.AxHor2.xaxis.set_ticks_position('bottom')
        self.AxHor2.xaxis.set_label_position('bottom')
        self.AxHor2.spines['bottom'].set_position(('outward', 100))
        self.AxHor2.spines['bottom'].set_color(self.Color_Vector[1])
        self.AxHor2.xaxis.get_offset_text().set_color(self.Color_Vector[1])
        self.AxHor2.ticklabel_format(style='sci', axis='x', scilimits=(6, 6))

        return
              
class LoadPickle(myPickle):                             #INTEGRATE LOAD PICKLE INTO PICKLE
    
    def __init__(self, CheckComputer = True):
        
        myPickle.__init__(self)
              
        return

    def RestoreFig(self, PlotFileAddress, Figure, Axis):
        
        self.Fig = Figure
        self.Axis1 = Axis
    
        FilesFolders = self.ComputerRoot + PlotFileAddress           
                       
        TextFile = open(FilesFolders,"r")
        PlotLogLines = TextFile.readlines()
        TextFile.close()
        
        self.ImportData(PlotLogLines, FilesFolders)
        
        self.ExportDataToFigure()
    
        return
    
    def RestoreFig_v2(self, PlotFileAddress):
            
        FilesFolders = PlotFileAddress           
                       
        TextFile = open(FilesFolders,"r")
        PlotLogLines = TextFile.readlines()
        TextFile.close()
        
        self.ImportData(PlotLogLines, FilesFolders)
        
        self.ExportDataToFigure()
    
        return    

    def ImportData(self, PlotLogLines, FilesFolders):
        
        #First we recover the name of the vectors stored on the text file
        PlotParameters = []
                
        Count = 0
        while PlotLogLines[Count] != '-Data Matrix:\n':
            if Count % 2 == 0:
                PlotParameters.append(PlotLogLines[Count][1:-2])
            Count += 1
                
        
        #Second we load the text data on the text vectors (this will only load them if they exist
        self.Tags                       = self.ImportLineData(PlotParameters, 'Tags', PlotLogLines) 
 
        self.FigConf                    = self.ImportLineData(PlotParameters, 'FigConf', PlotLogLines)
        
        self.FileAddress                = self.ImportLineData(PlotParameters, 'Files', PlotLogLines)  
        
        self.PlotNotation               = self.ImportLineData(PlotParameters, 'Plot Notation', PlotLogLines)  
        
        self.DataLineLabel              = self.ImportLineData(PlotParameters, 'Data Labels', PlotLogLines)   
        
        self.DataLineFormat             = self.ImportLineData(PlotParameters, 'Data Plot Conf', PlotLogLines)
        
        self.HistogramFormat            = self.ImportLineData(PlotParameters, 'Histogram data format', PlotLogLines)
                
        self.DataFilledFormat           = self.ImportLineData(PlotParameters, 'Filled data format', PlotLogLines)
        
        self.DataFilledBetweem_Format   = self.ImportLineData(PlotParameters, 'Filled between data format', PlotLogLines)

        self.TextFormat                 = self.ImportLineData(PlotParameters, 'Text in plot', PlotLogLines)
        
        self.ImageVector                = self.ImportLineData(PlotParameters, 'Inserted images format', PlotLogLines)
        
        self.HeadersData                = PlotLogLines[Count + 1].strip('\n').split()
        
        if self.DataFilledFormat            != None:          #A NEW DESIGN WOULD BE NICE IN WHICH ALL THE VECTORS ARE UNIFORMILY CHECKED AND NOT JUST THE LESS COMMON
            self.DataFilled_Initial         = [None] * len(self.DataFilledFormat)
            self.DataFilled_Final           = [None] * len(self.DataFilledFormat)
        
        if self.DataFilledBetweem_Format    != None:
            self.DataFilledBetweem_Xrange   = [None] * len(self.DataFilledBetweem_Format)
            self.DataFilledBetweem_YLow     = [None] * len(self.DataFilledBetweem_Format)
            self.DataFilledBetweem_YHigh    = [None] * len(self.DataFilledBetweem_Format)          
        
        if self.HistogramFormat != None:
            self.HistogramData_X            = [None] * len(self.HistogramFormat)
            self.HistogramData_Y            = [None] * len(self.HistogramFormat)
        
        if self.TextFormat != None:
            self.Text_X                              = [None] * len(self.TextFormat)
            self.Text_Y                              = [None] * len(self.TextFormat)
            self.Text_Values                         = [None] * len(self.TextFormat)
        
        #This is not nice... the second condition should not be necessary
        if self.DataLineLabel != None:
            DataColumns = len(self.DataLineLabel)
        elif self.HistogramFormat != None:
            DataColumns = len(self.HistogramFormat)
        
        if DataColumns > 0:
            #We create empty lists to store the matrix data 'SHOULD WE CHANGE THIS BY NUMPY arrays????????
            self.DataXValues                    = [None] * DataColumns
            self.DataYValues                    = [None] * DataColumns
            self.DataXError                     = [None] * DataColumns
            self.DataYError                     = [None] * DataColumns
            
            #We recover the scientific data and we put it on the corresponding data vectors
            for j in range(len(self.HeadersData)):                      # No bad but it would be more elegant to actually use load text only once for each Data Type
                Header      =   self.HeadersData[j]
                DataType    =   Header[0:Header.find('_')]
                DataIndex   =   int(Header[Header.find('_')+1:len(Header)])

                if DataType != 'TV':
                    
                    DataArray = genfromtxt(FilesFolders, dtype=float, usecols = [j], skiprows = Count + 2)
                    
                    if DataType == 'X':
                        self.DataXValues[DataIndex]                     = DataArray[~isnan(DataArray)]
                    elif DataType == 'Y':
                        self.DataYValues[DataIndex]                     = DataArray[~isnan(DataArray)]
                    elif DataType == 'XE':
                        self.DataXError[DataIndex]                      = DataArray[~isnan(DataArray)]
                    elif DataType == 'YE':
                        self.DataYError[DataIndex]                      = DataArray[~isnan(DataArray)]
                    elif DataType == 'FO':
                        self.DataFilled_Initial[DataIndex]              = DataArray[~isnan(DataArray)]
                    elif DataType == 'FF':
                        self.DataFilled_Final[DataIndex]                = DataArray[~isnan(DataArray)]
                    elif DataType == 'FBX':
                        self.DataFilledBetweem_Xrange[DataIndex]        = DataArray[~isnan(DataArray)]
                    elif DataType == 'FBO':                   
                        self.DataFilledBetweem_YLow[DataIndex]          = DataArray[~isnan(DataArray)]
                    elif DataType == 'FBF':
                        self.DataFilledBetweem_YHigh[DataIndex]         = DataArray[~isnan(DataArray)]
                    elif DataType == 'HX':
                        #THIS ONLY WORKS FOR 1 POINT bars... (BUT WE CAN PUT MANY)
                        self.HistogramData_X[DataIndex]                 = DataArray
                        
                    elif DataType == 'HY':
                        self.HistogramData_Y[DataIndex]                 = DataArray
                        
                    elif DataType == 'TX':
                        self.Text_X[DataIndex]          = DataArray[~isnan(DataArray)]
                    elif DataType == 'TY':                   
                        self.Text_Y[DataIndex]          = DataArray[~isnan(DataArray)]
                        
                else:
    
                    DataArray = genfromtxt(FilesFolders, dtype=str, usecols = [j], skiprows = Count + 2)
                    Index_nan = where(char.find(DataArray, self.EmptyRowFormat) > -1)
                    self.Text_Values[DataIndex] = delete(DataArray, Index_nan).tolist()
        
        #Finally we load the text vector data
        return   
       
    def ExportDataToFigure(self):
                
        self.FigFormat_One(ColorConf = self.FigConf[0], 
                           FigWidth = None,
                           FigHeight = None,
                           AxisFormat = None,                       #Doing this the Fig and Axis1 will not be generated again 
                           TickLength = int(self.FigConf[4]), 
                           TickWidth = int(self.FigConf[5]), 
                           myLabelSize = int(self.FigConf[6]),  
                           StoreParameters = False)
       
        if self.DataLineLabel != None: #WE SHOULD CHANGE THE DEFAULT EMPTY VALUE TO NONE   
            for i in range(len(self.DataLineLabel)):
                self.DataPloter_One(X               = self.DataXValues[i], 
                                    Y               = self.DataYValues[i], 
                                    LineLabel       = self.DataLineLabel[i],
                                    LineColor       = self.Tupling(self.DataLineFormat[i][0]), 
                                    LineStyle       = self.DataLineFormat[i][1], 
                                    LineWidth       = int(self.DataLineFormat[i][2]),
                                    MarkerColor     = self.Tupling(self.DataLineFormat[i][3]),
                                    MarkerStyle     = self.DataLineFormat[i][4],
                                    MarkerSize      = int(self.DataLineFormat[i][5]),
                                    XError          = self.DataXError[i],
                                    YError          = self.DataYError[i],
                                    ErrorBarsColor  = self.Tupling(self.DataLineFormat[i][6]),
                                    ErrorBarsWidth  = int(self.DataLineFormat[i][7]),
                                    ExtraData       = None,
                                    StoreParameters = False)
            
        if self.TextFormat != None:
            for i in range(len(self.TextFormat)):
                self.TextPlotter(x_coords           = self.Text_X[i],
                                 y_coords           = self.Text_Y[i],
                                 Text               = self.Text_Values[i], 
                                 fontsize           = int(self.TextFormat[i][0]),
                                 Color              = self.Tupling(self.TextFormat[i][1]), 
                                 x_pad              = float(self.TextFormat[i][2]),
                                 y_pad              = float(self.TextFormat[i][3]),
                                 VerticalAligment   = self.TextFormat[i][4],
                                 PlottingType       = self.TextFormat[i][5], 
                                 StoreParameters    = False)
                
                
                #[DataLabel, FillColor, transparency, Fill_Style, BorderColor, Type]
        
        if self.ImageVector != None:
            for i in range(len(self.ImageVector)):
                self.InsertFigure(ImageFolder = self.ImageVector[i][0],
                                  ImageName = self.ImageVector[i][1],
                                  Image_Coordinates = self.Tupling(self.ImageVector[i][2]),
                                  Zoom = float(self.ImageVector[i][3]),
                                  Image_Box = self.Tupling(self.ImageVector[i][4]),
                                  Image_xyCoords = self.ImageVector[i][5],
                                  StoreParameters = False)
              
        if self.DataFilledFormat != None:
            for i in range(len(self.DataFilledFormat)):
                self.DataArea_One(Points_O          =self.DataFilled_Initial[i], 
                                  Points_F          =self.DataFilled_Final[i], 
                                  DataLabel         =self.DataFilledFormat[i][0],
                                  FillColor         =self.Tupling(self.DataFilledFormat[i][1]), 
                                  transparency      =float(self.DataFilledFormat[i][2]), 
                                  Fill_Style        =self.DataFilledFormat[i][3], 
                                  BorderColor       =self.DataFilledFormat[i][4], 
                                  Type              =self.DataFilledFormat[i][5], 
                                  StoreParameters   =False)

        if self.DataFilledBetweem_Format != None:
            for i in range(len(self.DataFilledBetweem_Format)):
                self.DataAreaFilled_One(X_Range         = self.DataFilledBetweem_Xrange[i],
                                        Y_Low           = self.DataFilledBetweem_YLow[i],
                                        Y_High          = self.DataFilledBetweem_YHigh[i], 
                                        Label           = self.DataFilledBetweem_Format[i][0], 
                                        Where           = self.DataFilledBetweem_Format[i][1], 
                                        FillColor       = self.Tupling(self.DataFilledBetweem_Format[i][2]), 
                                        transparency    = float(self.DataFilledBetweem_Format[i][3]), 
                                        Fill_Style      = self.DataFilledBetweem_Format[i][4], 
                                        BorderColor     = self.Tupling(self.DataFilledBetweem_Format[i][5]), 
                                        Type            = self.DataFilledBetweem_Format[i][6], 
                                        StoreParameters = False)
                
        if self.HistogramFormat != None:
            for i in range(len(self.HistogramFormat)):
                self.Histogram_One(X_Coord              = self.HistogramData_X[i],
                                   Y_Coord              = self.HistogramData_Y[i], 
                                   Label                = self.HistogramFormat[i][0], 
                                   Color                = self.Tupling(self.HistogramFormat[i][1]), 
                                   Transparency         = float(self.HistogramFormat[i][2]), 
                                   Width                = float(self.HistogramFormat[i][3]), 
                                   Fill                 = (self.HistogramFormat[i][4] == 'True'), 
                                   Edgecolor            = self.Tupling(self.HistogramFormat[i][5]), 
                                   Linestyle            = self.HistogramFormat[i][6], 
                                   Linewidth            = self.HistogramFormat[i][7], 
                                   Log                  = (self.HistogramFormat[i][8] == 'True'), 
                                   StoreParameters      = False)

        if self.PlotNotation != [''] :  #NOT SURE WHY THIS IS WORKING... IT SHOULD BE CHANGE TO A NONE CHEK
            self.Labels_Legends_One(Plot_Title  = self.PlotNotation[0], 
                           Plot_xlabel          = self.PlotNotation[1], 
                           Plot_ylabel          = self.PlotNotation[2], 
                           XlabelSize           = int(self.PlotNotation[3]), 
                           YLabelSize           = int(self.PlotNotation[4]), 
                           TitleLabelSize       = int(self.PlotNotation[5]), 
                           LegendLocation       = int(self.PlotNotation[6]), 
                           Expand               = self.PlotNotation[7] == True, 
                           XLabelPad            = int(self.PlotNotation[8]), 
                           YLabelPad            = int(self.PlotNotation[9]), 
                           Y_TitlePad           = float(self.PlotNotation[10]),  
                           StoreParameters      = False)
        
        return
        
    def DisplayFigure(self):
        
        plt.show()
        
        return
        
    def ImportLineData(self, PlotParameters, Parameter, PlotLogLines):
        
        if Parameter in PlotParameters:
            Index = PlotParameters.index(Parameter) * 2 + 1
            ListElements = [i.split(self.ParameterSeparator)        for i in      PlotLogLines[Index].strip('\n').split(self.FileSeparator)]
        
            if Parameter in ['Tags', 'FigConf', 'Files', 'Plot Notation', 'Data Labels']:
                ListElements = ListElements[0]
        
        else:
            ListElements = None
        
        return ListElements
    
    def Tupling(self, Parameter):
        Output = None
    
        if Parameter != 'None':
            Output = tuple(map(float, Parameter.translate(None,'()').split(',')))
            
        return Output
                
    def FindAndOrganize(self, FilePattern, FilesFolders, unpack = False, CheckComputer=False, Sort_Output = 'Alpha'):
 
        if CheckComputer == True:
            FilesFolders = self.ComputerRoot + FilesFolders
         
        #Checking for arguments in terminal
        ArgumentsCheck, Arguments = Arguments_Checker(argv) 
        print "Initiating " + Arguments[0][Arguments[0].rfind("/")+1:len(Arguments[0])] + " task\n"
            
        print "Files found meeting the pattern: " + str(FilePattern), "@", FilesFolders, ':\n' 
            
        #Command from terminal
        if ArgumentsCheck == True:
            self.Tags.append('Terminal')
            FilesList, FlagsList = Arguments_Handler(ArgumentsCheck, Arguments)
             
            if len(FilesList) == 1:
                self.Tags.append('SingleFile')
             
            self.ListFiles = (FilesList,)
            return self.ListFiles
     
        #Command from eclipse
        else:
            self.Tags.append('Editor')
            if unpack == False:     
                FilesList = FilesFinder(FilesFolders,FilePattern)
 
                Single_List = []
                Combined_List = []
                 
                LeftPattern = '_@'
                RightPattern = '@_'
 
                for FileAddress in FilesList:
                    FileName = FileAddress[FileAddress.rfind("/")+1:len(FileAddress)]
                    if (LeftPattern in FileName) and (RightPattern in FileName):
                        FileLabel = FileName[FileName.find(LeftPattern)+2:FileName.find(RightPattern)]
                        Matching_List = []
                        Already_Seen = False
                         
                        for item in Combined_List:
                            for family in item:
                                if family == FileAddress:
                                    Already_Seen = True
                         
                        if Already_Seen == False:
                            for FileLocation in FilesList:
                                FilNam = FileLocation[FileLocation.rfind("/")+1:len(FileLocation)]
                                if FileLabel in FilNam:
                                    Matching_List.append(FileLocation)
                                    
                            Combined_List.append(tuple(Matching_List))
                         
                    else:
                        Single_List.append((FileAddress,))
                 
                self.ListFiles = tuple(Single_List + Combined_List)    
                 
                if len(self.ListFiles) == 1:
                    self.Tags.append('SingleFile')
            else:
                FoldersList = []
                ArchivesList = []
                 
                if type(FilePattern) is not list:
                    myPatternList = [FilePattern]
                      
                else:
                    myPatternList = FilePattern
                   
                for Root, Dirs, Archives in walk(FilesFolders):
                    ValidArchives = []
                    for Archive in Archives:
                        Meets_one_Pattern = False
                        for i in range(len(myPatternList)):
                            if (myPatternList[i] in Archive):
                                if "~" not in Archive:                   
                                    Meets_one_Pattern = True
                                    print '--- File', Archive, '@', Dirs, Root               
                 
                        if Meets_one_Pattern:        
                            if Root.endswith("/"):
                                FinalName = Root
                            else:
                                FinalName = Root + "/"
                 
                            if FinalName in FoldersList:
                                ValidArchives.append(Archive)
                            else:
                                FoldersList.append(FinalName)
                                ValidArchives.append(Archive)
                      
                    if len(ValidArchives) > 0:
                        ValidArchives.sort()
                        ArchivesList.append(ValidArchives)
             
                if Sort_Output == 'Alpha':
                    FoldersList, ArchivesList = zip(*sorted(zip(FoldersList, ArchivesList), key=itemgetter(0), reverse=False))
                     
                self.ListFiles = (FoldersList, ArchivesList)    
 
        return self.ListFiles

