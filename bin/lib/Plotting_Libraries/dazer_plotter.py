import pickle
import corner
from uncertainties                      import UFloat
from numpy                              import array, arange, zeros, log10, ndarray, percentile, median, string_, where, argsort, sort, unique, max, sum as np_sum
from collections                        import OrderedDict, Sequence
from matplotlib                         import image, colors, cm, rcParams, pyplot as plt
from matplotlib._png                    import read_png
from matplotlib.offsetbox               import OffsetImage, AnnotationBbox
from matplotlib import gridspec
from lib.Math_Libraries.sigfig    import round_sig

class Fig_Conf():
    
    def __init__(self):
        
        self.Fig                = None
        self.Axis               = None
        self.ColorVector        = None
        
    def Import_FigConf(self, Fig = None, Axis = None):
        
        if Fig != None:
            self.Fig = Fig
            
        if Axis != None:
            self.Axis = Axis

    def define_ColorVector(self, PlotStyle, n_colors, color_map):
        
        self.ColorVector = ['white', 'black', 'blue']
    
    def gen_colorList(self, vmin=0.0, vmax=1.0, color_palette = None):
        
        self.colorNorm = colors.Normalize(vmin, vmax)
        self.cmap = cm.get_cmap(name = color_palette)

    def get_color(self, idx):
        
        return self.cmap(self.colorNorm(idx))
          
    def define_format(self, plotStyle, plotSize):

        #Default sizes for computer
        sizing_dict = {}
        sizing_dict['figure.figsize'] = (14, 8)
        sizing_dict['legend.fontsize'] = 15
        sizing_dict['axes.labelsize'] = 20
        sizing_dict['axes.titlesize'] = 24
        sizing_dict['xtick.labelsize'] = 14
        sizing_dict['ytick.labelsize'] = 14
        
        self.colorVector = {
        'iron':'#4c4c4c',
        'silver':'#cccccc',                  
        'dark blue':'#0072B2',
        'green':'#009E73', 
        'orangish':'#D55E00',
        'pink':'#CC79A7',
        'yellow':'#F0E442',
        'cyan':'#56B4E9',
        'olive':'#bcbd22',
        'grey':'#7f7f7f',
        'skin':'#FFB5B8'}

        #sizing_dict['text.usetex'] = True
        
        #--Update the colors/format
        if plotStyle == None:
            self.ColorVector = [None, None, None]
        
        elif plotStyle == 'dark':
            plt.style.use('dark_background')

        elif plotStyle == 'night':
            
            plt.style.use('seaborn-colorblind')
            
            iron_color = '#4c4c4c' #Iron: (76 76 76)
            silver_color = '#cccccc' #Silver: (204 204 204) 
            sizing_dict['axes.facecolor']   = iron_color
            sizing_dict['figure.facecolor'] = iron_color
            sizing_dict['axes.edgecolor']   = silver_color
            sizing_dict['text.color']       = silver_color
            sizing_dict['axes.labelcolor']  = silver_color
            sizing_dict['xtick.color']      = silver_color
            sizing_dict['ytick.color']      = silver_color
            sizing_dict['axes.edgecolor']   = silver_color

            
            #'plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) + cycler('linestyle', ['-', '--', ':', '-.'])))'
            #This should be the set up for the cycler we just need to add the colors
            #axes.prop_cycle : cycler('color', 'bgrcmyk')

        else:
            plt.style.use(plotStyle)
        
        #--Load particular configuration for this plot
        if plotSize == 'medium':            
            rcParams.update(sizing_dict)
        
        elif type(plotSize) is dict:
            sizing_dict.update(plotSize)
            rcParams.update(sizing_dict)

        '''
        Seaborn color blind
        #0072B2 dark blue
        #009E73 green 
        #D55E00 orangish
        #CC79A7 pink
        #F0E442 yellow
        #56B4E9 cyan
        #bcbd22 olive #adicional
        #7f7f7f grey
        #FFB5B8 skin
        '''
 
    
        '''
        Matplotlib default palete
        #17becf dark blue
        #bcbd22 orange
        #2ca02c green
        #e377c2 red
        #8c564b purple
        #9467bd brown
        #d62728 pink
        #7f7f7f grey
        #ff7f0e olive
        #1f77b4 cyan
        '''


  
        '''
        --These are matplotlib styles
        seaborn-darkgrid
        seaborn-notebook
        classic
        seaborn-ticks
        grayscale
        bmh
        seaborn-talk
        dark_background
        ggplot
        fivethirtyeight
        seaborn-colorblind
        seaborn-deep
        seaborn-whitegrid
        seaborn-bright
        seaborn-poster
        seaborn-muted
        seaborn-paper
        seaborn-white
        seaborn-pastel
        seaborn-dark
        seaborn
        seaborn-dark-palette
        '''
                  
    def FigConf(self, plotStyle = None, plotSize = 'medium', Figtype = 'Single', AxisFormat = 111, n_columns = None, n_rows = None, n_colors = None, color_map = None, axis_not = None):
                
        #Set the figure format before creating it
        self.define_format(plotStyle, plotSize)
        
        if Figtype == 'Single':
            
            if AxisFormat == 111:
                self.Fig                = plt.figure()  
                self.Axis               = self.Fig.add_subplot(AxisFormat)
            else:
                self.Fig, self.Axis    = plt.subplots(n_rows, n_columns)  
                self.Axis              = self.Axis.ravel()
        
        elif Figtype == 'Posteriors':
            self.Fig                = plt.figure()
            AxisFormat              = int(str(n_colors)  + '11')
            self.Axis              = self.Fig.add_subplot(AxisFormat)
        
        elif Figtype == 'grid':
            self.Fig, self.Axis = plt.subplots(n_rows, n_columns) 
            if (n_rows * n_columns) != 1:
                self.Axis = self.Axis.ravel()
            else:
                self.Axis = [self.Axis]  
                  
        elif Figtype == 'Grid':
            frame1 = plt.gca()
            frame1.axes.xaxis.set_visible(False)
            frame1.axes.yaxis.set_visible(False)
            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
            self.Fig, self.Axis = plt.subplots(n_rows, n_columns)  
            self.Axis = self.Axis.ravel()

        elif Figtype == 'Grid_size':
            self.Fig = plt.figure()
            gs = gridspec.GridSpec(n_rows, n_columns, height_ratios=[2.5, 1])
            self.ax1 = self.Fig.add_subplot(gs[0,:])
            self.ax2 = self.Fig.add_subplot(gs[1,:])

        if axis_not == 'sci': 
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
        #self.define_ColorVector(PlotStyle, n_colors, color_map)
   
        return
   
    def FigWording(self, xlabel, ylabel, title, loc = 'best', Expand = False,  XLabelPad = 0.0, YLabelPad = 0.0, Y_TitlePad = 1.02, cb_title = None, sort_legend = False, ncols_leg = 1, graph_axis = None):

        if graph_axis == None:
            Axis = self.Axis
        else:
            Axis = graph_axis

        Axis.set_xlabel(xlabel)
        Axis.set_ylabel(ylabel)
        Axis.set_title(title, y = Y_TitlePad)    
        
        if (XLabelPad != 0) or (YLabelPad !=0):
            Axis.xaxis.labelpad = XLabelPad
            Axis.yaxis.labelpad = YLabelPad
           
        if cb_title != None:
            self.cb.set_label(cb_title, fontsize=18)
        
        self.legend_conf(Axis, loc, sort_legend=sort_legend, ncols=ncols_leg)
      
    def legend_conf(self, Axis = None, loc = 'best', sort_legend = False, ncols = 1):
        
        if Axis == None:
            Axis = self.Axis
            
        Axis.legend(loc = loc, ncol = ncols)
 
        #Security checks to avoid empty legends
        if Axis.get_legend_handles_labels()[1] != None:
             
            if len(Axis.get_legend_handles_labels()[1]) != 0:
                Old_Handles, Old_Labels = Axis.get_legend_handles_labels()
                 
                if sort_legend:
                    labels, handles = zip(*sorted(zip(Old_Labels, Old_Handles), key=lambda t: t[0]))
                    Handles_by_Label = OrderedDict(zip(labels, handles))
                    Axis.legend(Handles_by_Label.values(), Handles_by_Label.keys(), loc=loc, ncol = ncols)
                else:
                    Handles_by_Label = OrderedDict(zip(Old_Labels, Old_Handles))
                    Axis.legend(Handles_by_Label.values(), Handles_by_Label.keys(), loc=loc, ncol = ncols)
 
        
#         #Security checks to avoid empty legends
#         if self.Axis.get_legend_handles_labels()[1] != None:
#             
#             if len(self.Axis.get_legend_handles_labels()[1]) != 0:
#                 Old_Handles, Old_Labels = Axis.get_legend_handles_labels()
#                 Handles_by_Label = OrderedDict(zip(Old_Labels, Old_Handles))
#                 
#                 Hl          = zip(Handles_by_Label.values(), Handles_by_Label.keys())
#             
#                 New_Handles, New_labels = zip(*Hl)        
#     
#                 myLegend    = Axis.legend(New_Handles, New_labels, loc=loc, prop={'size':12}, scatterpoints=1, numpoints=1)
#     
#                 for i in range(len(New_labels)):
#                     myLegend.legendHandles[i]._sizes = [30]
#                     
#                 if fontize != None:
#     
#                     Leg_Frame = myLegend.get_frame()
# #                     Leg_Frame.set_facecolor(self.ColorVector[0])
# #                     Leg_Frame.set_edgecolor(self.ColorVector[1])
#         
#                     for label in myLegend.get_texts():
#                         label.set_fontsize(fontize)
#                  
#                     for label in myLegend.get_lines():
#                         label.set_linewidth(1.5)
#                      
# #                     for text in myLegend.get_texts():
# #                         text.set_color(self.ColorVector[1])   
                    
        return

    def bayesian_legend_conf(self, Axis, loc = 'best', fontize = None, edgelabel = False):
        #WARNING: THIS DOES NOT WORK WITH LEGEND RAVELIN
   
        if Axis.get_legend_handles_labels()[1] != None:
            Old_Handles, Old_Labels = Axis.get_legend_handles_labels()
            Handles_by_Label = OrderedDict(zip(Old_Labels, Old_Handles))
            
            Hl          = zip(Handles_by_Label.values(), Handles_by_Label.keys())

            New_Handles, New_labels = zip(*Hl)        

            myLegend    = Axis.legend(New_Handles, New_labels, loc=loc, prop={'size':12}, scatterpoints=1, numpoints=1)

            if fontize != None:

                Leg_Frame = myLegend.get_frame()
                Leg_Frame.set_facecolor(self.Color_Vector[0])
                Leg_Frame.set_edgecolor(self.Color_Vector[1])
    
                for label in myLegend.get_texts():
                    label.set_fontsize('large')
             
                for label in myLegend.get_lines():
                    label.set_linewidth(1.5)
                 
                for text in myLegend.get_texts():
                    text.set_color(self.Color_Vector[1])   
            
            if edgelabel:
                Leg_Frame = myLegend.get_frame()
                Leg_Frame.set_edgecolor('black')
                                
    def savefig(self, output_address, extension = '.png', reset_fig = True, pad_inches=0.2):

        plt.savefig(output_address + extension, dpi=350.0, bbox_inches='tight', pad_inches=pad_inches)
        
        if reset_fig:
            self.reset_fig()
    
        return
    
    def reset_fig(self):
        
        plt.cla()
        
        return
    
    def display_fig(self):
        
        plt.show()

    def save_manager(self, output_address, save = None, display = None, reset_fig = True, save_pickle = None, pad_inches = 0.2, extension = '.png'):
  
        if (save == None) and (display == None):

            #If only one file
            if 'SingleFile' in self.Tags:    
                display = True
            
            #For multiple files
            else:
                save = True
                    
        if save:
            plt.savefig(output_address + extension, dpi=150, bbox_inches='tight', pad_inches=pad_inches, facecolor=self.Fig.get_facecolor())
            
        if save_pickle:            
            with open(output_address + '.dz_pickle','wb') as fid:
                pickle.dump(self.Fig, fid)
        
        if display:
            self.Fig.tight_layout()
            plt.show()

        if reset_fig:
            plt.cla()

        return

    def restore_from_pickle(self, pickle_address, Fig=None):

        if Fig == None:
            with open(pickle_address, 'rb') as fid:
                self.Axis = pickle.load(fid)
                
        return

class Plot_Conf(Fig_Conf):
    
    def __init__(self):
        
        Fig_Conf.__init__(self)
        
    def data_plot(self, x, y, label = '', color = None, linestyle = None, markerstyle = None, linewidth = None, markersize = None, x_error=None, y_error=None, cmap=None, graph_axis = None):
        
        if graph_axis == None:
            Axis = self.Axis
        else:
            Axis = graph_axis
        
        #If not symbol style defined the default will be a line
        if (linestyle == None) and (markerstyle == None):
            linestyle = '-'
        
        #Define data style:
        if linestyle != None:
            
            if linestyle == 'step':
                Axis.step(x, y, label=label, color=color, where = 'mid')
            
            elif (x_error == None) and (x_error == None): 
                Axis.plot(x, y, label=label, color=color, linestyle=linestyle, marker=markerstyle, linewidth = linewidth)
#                 print line[-1].get_color()
            else:
                Axis.errorbar(x, y, label=label, xerr = x_error, yerr = y_error, elinewidth = 1, fmt=markerstyle, ms=markersize, c=color)
        
        elif markerstyle != None:
            if (x_error is None) and (y_error is None): 
                if cmap == None:
                    
                    Axis.scatter(x, y, label=label, color=color, marker=markerstyle, s=markersize)
                else:
                    scatter_data = self.Axis.scatter(x, y, label=label, c=color, marker=markerstyle, s=markersize, cmap=cmap)
            else:
                Axis.errorbar(x, y, label=label, xerr = x_error, yerr = y_error, elinewidth = 2, fmt=markerstyle, ms=10, c=color)

        if cmap != None:
            
            self.cb = plt.colorbar(scatter_data, cmap=cmap)

    def surface_plot(self, x_range, y_range, data_matrix, plotAxis = None, interpolation_type = 'bilinear'):
        
        #Decide axis for the image to be plotted
        if plotAxis == None:
            plotAxis = self.Axis
        
        im = image.NonUniformImage(ax=plotAxis, interpolation=interpolation_type, cmap=cm.gist_earth)
        
        #Input data
        im.set_data(x_range, y_range, data_matrix)
        
        #Plot the image
        plotAxis.images.append(im)
        
        #Set the axis limit
        plotAxis.set_xlim(x_range[0],x_range[-1])
        plotAxis.set_ylim(y_range[0],y_range[-1])

    def insert_image(self, image_address, Image_Coordinates = [0.20,0.5], Zoom = 1, Image_Box=(20,-20), Image_xyCoords = 'data', Image_BoxCoords = "offset points"):
                
        arr_hand = read_png(image_address)
        Image_Frame = OffsetImage(arr_hand, zoom=Zoom)
        ab = AnnotationBbox(Image_Frame, Image_Coordinates,
            xybox=Image_Box,
            xycoords=Image_xyCoords,
            boxcoords=Image_BoxCoords)
        self.Axis.add_artist(ab)

    def plot_text(self, x_coords, y_coords, text, x_pad = 1, y_pad = 1, color = 'black', fontsize = 12, axis_plot=None):
        
        if axis_plot == None:
            axis_plot = self.Axis
        
        #Single string go ahead
        if isinstance(text, (str, string_)):
            axis_plot.text(x_coords * x_pad, y_coords * y_pad, text, fontsize=fontsize)
        
        #Uncertainties entry
        elif isinstance(text, UFloat):
            entry = r'${:L}$'.format(text)
            axis_plot.text(x_coords * x_pad, y_coords * y_pad, entry, fontsize=fontsize)
        
        #else it is an array loop through
        else:
            for i in range(len(text)):
                if isinstance(text[i], UFloat):
                    entry = r'{:L}'.format(text[i])  
                else:
                    entry = text[i]
            
                axis_plot.text(x_coords[i] * x_pad, y_coords[i] * y_pad, entry, fontsize=fontsize, color = color)
            
    def Histogram_One(self, X_Coord, Y_Coord, Label, Color = None, Transparency = 1.0, Width=0.8, Fill = True, Edgecolor = None, Linestyle='solid', Linewidth=1, Log=False, StoreParameters = True):
        
        #print 'Los valores', X_Coord, Y_Coord, Label, Color, Transparency, Width, Fill, Edgecolor, Linestyle, Linewidth, Log, StoreParameters
        Bar = self.Axis.bar(X_Coord, Y_Coord, label=Label, color=Color, width=Width, alpha=Transparency, fill = Fill, edgecolor=Edgecolor, linestyle = Linestyle, linewidth = Linewidth, log = Log)
                
        return

    def populations_histogram(self, folder_starlight, file_starlight, plot_type):
                
        #Extract the data from the starlight ouuput
        index, x_j, Mini_j, Mcor_j, age_j, Z_j, LbyM, Mstars = self.SlOutput_2_StarData(folder_starlight, file_starlight)
        
        x_j     = array(x_j)
        age_j   = array(age_j)
        Mcor_j  = array(Mcor_j)
        Z_j     = array(Z_j)
        
        #Set fractioin type
        if plot_type == 'Light_fraction':
            fraction_param = x_j
        elif plot_type == 'Mass_fraction':
            fraction_param = Mcor_j
        
        #Get metallicities from calculation
        zValues = sort(array(unique(Z_j)))
        
        #Get color list for metallicities
        self.gen_colorList(0, len(zValues))
        
        #Define age bins width
        log_agej    = log10(age_j)
        ageBins_HW  = 0.10      
        log_ageBins = arange(4.80, 10.60, ageBins_HW)
            
        #Plot the age bins per metallicity:
        for log_t in log_ageBins:
            
            idx_age     = (log_agej >= log_t - ageBins_HW/2) & (log_agej <= log_t + ageBins_HW/2)
            idx_param   = (fraction_param > 0.0)
            idx_total   = idx_age & idx_param
               
            #Check populatioins were found
            if np_sum(idx_total) > 1:
                                                                             
                #Load x and y data for the plot
                z_array     = Z_j[idx_total]
                param_array = fraction_param[idx_total]
                
                #Sort by decreasing param value
                idx_sort    = argsort(param_array)[::-1]
                z_sort      = z_array[idx_sort]
                param_sort  = param_array[idx_sort]
                
                #Plot the individual bars            
                for i in range(len(param_sort)):                    
                    x, y        = log_t, param_sort[i]
                    label       = 'z = {}'.format(z_sort[i])
                    idx_color   = where(zValues == z_sort[i])[0][0] 
                    color       = self.get_color(idx_color)
                    self.Axis.bar(x, y, label=label, color=color, width=ageBins_HW/2, fill = True, edgecolor=color, log = True)
                    
                #Plot age bin total
                param_total = sum(param_array) 
                self.Axis.bar(log_t, param_total, label='Age bin total',  width=ageBins_HW, fill = False, edgecolor='black', linestyle= '--', log = True)
                
                
                
#         #Generate metallicity matrix
#         z_AgeBins_Matrix = zeros((len(zValues), len(ageBins))) 
#         print z_AgeBins_Matrix.shape
#         for i in range(len(zValues)):
#             
#             #Get all populations matching current metallicity 
#             idcs_Z      = (Z_j == zValues[i])
#             
#             #For z_i get time bin idcs
#             age_array = log10(age_j[idcs_Z])
#             idcs_hist = digitize(age_array, ageBins)
#             print idcs_hist.shape
#             
#             #Get parameter for given metallicty
#             param_array = fraction_param[idcs_Z]
#             
#             #Organize paramb by increasing age bin
#             print 'este', param_array[idcs_hist]
#             z_AgeBins_Matrix[:,i] = param_array[idcs_hist]
#             
#         #Loop throught the matrix columns to plot lower values in front
#         for i_t in range(len(ageBins_HW)):
#             age_column  = z_AgeBins_Matrix[:,i_t]
#             idx_size    = age_column.argsort()[::-1]
#             for idx in idx_size:
#                 x, y = ageBins_HW[i_t]
#                 label = age_column[idx], 'z = {}'.format()
#                 self.Axis.bar(x, y, label=label, color=self.get_color(idx), width=ageBins_HW/2, fill = True, edgecolor=self.get_color(idx), log = True)
                
            
        return

    def CumulativeHistogram_One(self, Folder, File, PlotType, Color_Vector):
    
        #Vector where we save the histograms data
        OutPutData = []
        
        #Extract the data from the starlight ouuput
        index, x_j, Mini_j, Mcor_j, age_j, Z_j, LbyM, Mstars    = self.SlOutput_2_StarData(Folder, File)
        
        #Case light fraction plot
        if PlotType == 'Light_fraction':    
            ParameterVector = x_j
        
        #Case mass fraction plot    
        if PlotType == 'Mass_fraction':
            ParameterVector = Mcor_j
 
        #Extract metallicity values
        zValues = self.FamilyOfItemsInArray(Z_j)   
        
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
    
            DataVector.append([z, Color_Vector[i], MassSumVector, LightMassVector])
            
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
                self.Histogram_One([AgeBins[i]- AgeBins_HW/2], [Total], 'Bin total', self.ColorVector[0], Width = AgeBins_HW, Fill = False, Edgecolor = self.ColorVector[1], Linestyle='dotted')
    
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
    
    def area_fill(self, Points_O, Points_F, DataLabel, color = None, alpha = 1, Fill_Style = None, BorderColor = None, graph_axis = None):

        if graph_axis == None:
            Axis = self.Axis
        else:
            Axis = graph_axis
        
        if isinstance(Points_O, Sequence) or isinstance(Points_O, ndarray):
            for k in range(len(Points_O)):
                Axis.axvspan(Points_O[k], Points_F[k], facecolor = color, alpha=alpha, label=DataLabel)
        else:
            Axis.axvspan(Points_O, Points_F, facecolor = color, alpha=alpha, label=DataLabel)
                
        return
    
    def plot_Ownposteriors_histagram(self, Traces, labels, sigfig_n = 4, xlim = None, ylim = None):
        
        plt.rcParams['ytick.labelsize'] = 20
        plt.rcParams['xtick.labelsize'] = 20

        n_traces = len(Traces)
                
        for i in range(len(Traces)):
            Trace           = Traces[i]
            Axis_Location   = int(str(n_traces) + '1' +  str(i + 1))
            Axis            = plt.subplot(Axis_Location) 
            HDP_coords      = [percentile(Trace, 16), percentile(Trace, 84)]
            mean_value      = median(Trace)
            
            for HDP in HDP_coords:
                Axis.axvline(x = HDP, label = 'Percentiles: 16th-84th: ' + round_sig(HDP_coords[0], n=sigfig_n, scien_notation=False) + ' - ' + round_sig(HDP_coords[1], n=sigfig_n, scien_notation=False), color='grey', linestyle = 'dashed')
            
 
                    
            Median_text_legend = r'Median value: $' + round_sig(mean_value, n=sigfig_n, scien_notation=False) + '_{-' + round_sig(mean_value - HDP_coords[0], n=sigfig_n, scien_notation=False) + '}^{+' + round_sig(HDP_coords[1] - mean_value, n=sigfig_n, scien_notation=False) + '}$'                               
            Axis.axvline(x = mean_value, label = Median_text_legend, color='grey', linestyle = 'solid')
            Axis.hist(Trace, histtype='stepfilled', bins=35, alpha=.7, color=self.ColorVector[2][i], normed=False)
            Axis.set_ylabel(labels[i], fontsize=20)
                        
            self.bayesian_legend_conf(Axis, loc='best', edgelabel=True)

            
            if xlim != None:
                Axis.set_xlim(xlim[0], xlim[1])
                
            if ylim != None:
                Axis.set_ylim(ylim[0], ylim[1])
        return

    def plot_triangle_histContours(self, Traces, labels, pymc_database, true_values = None):
        n_traces = len(Traces)
        
        #Make plots
        List_Arrays = []
        for i in range(n_traces):
            List_Arrays.append(pymc_database.trace(Traces[i])[:])
        
        Samples = array(List_Arrays).T
                    
        #Make 
        if true_values != None:
            self.figure = corner.corner(Samples[:,:], fontsize=30, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_args={"fontsize": 200}, truths = true_values, title_fmt='0.3f')
        else:
            self.figure = corner.corner(Samples[:,:],fontsize=30, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_args={"fontsize": 200}, title_fmt='0.3f')
        
        self.Axis = self.figure.get_axes()
        
        for item in self.Axis:
            item.title.set_fontsize(14)
            
        for item in  self.Axis:
            item.yaxis.label.set_fontsize(20)

        for item in  self.Axis:
            item.xaxis.label.set_fontsize(20)

           
#         for item in  self.Axis:
#             item.get_xticklabels(15)            
# 
#         for item in  self.Axis:
#             item.get_yticklabels(15)     
           
#         for item in ([, self.Axis.xaxis.label, self.Axis.yaxis.label] + self.Axis.get_xticklabels() + self.Axis.get_yticklabels()):
        
    
        rcParams.update({'font.size': 30})        




