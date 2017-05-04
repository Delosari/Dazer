from collections            import OrderedDict
from itertools              import cycle

import corner
from matplotlib             import rc
from matplotlib.mlab        import detrend_mean
from numpy                  import log as nplog, max as npmax, min as npmin, array, argmin, histogram2d, argsort, reshape, arange, percentile, mean, median
from pymc                   import database, MCMC

from Math_Libraries.sigfig  import round_sig
import matplotlib.pyplot    as plt
import seaborn              as sb


class bayes_plotter():
    
    def __init__(self):
        
        self.Fig                = None
        self.Axis               = None
        self.Valid_Traces       = None
        
        self.pymc_database      = None
        self.dbMCMC             = None

        self.Traces_filter      = None
        self.pymc_stats_keys    = ['mean', '95% HPD interval', 'standard deviation',  'mc error', 'quantiles', 'n']

    def load_pymc_database(self, Database_address):
        
        #In case the database is open from a previous use
        if self.pymc_database != None:
            self.pymc_database.close()
        
        #Load the pymc output textfile database
        self.pymc_database  = database.pickle.load(Database_address)
        
        #Create a dictionary with the bases to 
        self.Traces_dict = {}
        self.traces_list = self.pymc_database.trace_names[0] #This variable contains all the traces from the MCMC (stochastic and deterministic)
        
        for trace in self.traces_list:
            self.Traces_dict[trace] = self.pymc_database.trace(trace)
    
        #Generate a MCMC object to recover all the data from the run
        self.dbMCMC      = MCMC(self.Traces_dict, self.pymc_database)
        
        return

    def extract_traces_statistics(self, traces_list = None):
        
        #     traces_yplus = bp.pymc_database.trace('He_abud')[:]
        #     print 'The y_plus trace evolution\n'
        #     print 'Mean inf',      statistics_dict['He_abud']['mean']
        #     print 'Median numpy',  median(traces_yplus)
        #     print 'Mean numpy',    mean(traces_yplus), '\n'
        # 
        #     print 'percentiles: 25, 50, 75' 
        #     print percentile(traces_yplus,25), percentile(traces_yplus,50), percentile(traces_yplus,75),'\n'
        #     print 'percentiles: 16, 84' 
        #     print percentile(traces_yplus,16), percentile(traces_yplus,84),'\n'
        #     print 'percentiles: 37.73, 68.27' 
        #     print percentile(traces_yplus,37.73), percentile(traces_yplus,68.27),'\n'
        #     print 'Standard deviation'
        #     print percentile(traces_yplus,4.55), percentile(traces_yplus,95.45)
        #     print 'HUD 95', statistics_dict['He_abud']['95% HPD interval'],'\n'
        #     
        #     print 'PYMC std', statistics_dict['He_abud']['standard deviation']
        #     print 'Numpy std', std(traces_yplus), '\n'
        
        
        self.statistics_dict = OrderedDict()
                
        #If no list input we extract all the traces from the analysis        
        if traces_list == None:
            traces_list = self.traces_list
                
        for trace in traces_list:
            self.statistics_dict[trace] = OrderedDict()
            
            for stat in self.pymc_stats_keys:
                self.statistics_dict[trace][stat] = self.dbMCMC.trace(trace).stats()[stat]                
            
            Trace_array = self.pymc_database.trace(trace)[:] 
            self.statistics_dict[trace]['16th_p'] = percentile(Trace_array, 16)
            self.statistics_dict[trace]['84th_p'] = percentile(Trace_array, 84)
        
        return self.statistics_dict
        
    def close_database(self):
        
        self.pymc_database.close()
        self.pymc_database = None
        
        return

    def Import_FigConf(self, Fig = None, Axis = None):
        
        if Fig != None:
            self.Fig = Fig
            
        if Axis != None:
            self.Axis = Axis
    
    def FigConf(self, Figtype = 'Single', FigWidth = 16, FigHeight = 9, AxisFormat = 111, fontsize = 8, PlotStyle = 'Night', n_columns = None, n_rows = None, n_colors = None, color_map = 'colorblind'):
        
        self.Triangle_Saver = False
        
        if Figtype == 'Single':
            self.Fig                = plt.figure(figsize = (FigWidth, FigHeight))  
            self.Axis1              = self.Fig.add_subplot(AxisFormat)
            #fig.subplots_adjust(hspace = .5, wspace=.001) 
        
        elif Figtype == 'Posteriors':
            self.Fig                = plt.figure(figsize = (FigWidth, FigHeight))
            AxisFormat              = int(str(n_colors)  + '11')
            self.Axis1              = self.Fig.add_subplot(AxisFormat)
            
        elif Figtype == 'Grid':
            self.Fig, self.Axis1    = plt.subplots(n_rows, n_columns, figsize=(FigWidth, FigHeight))  
            self.Axis1              = self.Axis1.ravel()
        
        elif Figtype == 'triangle':
            self.Triangle_Saver = True
            return
        
        rc('legend',    fontsize=fontsize + 3,      handlelength=3, frameon = True)
        rc('axes',      titlesize=fontsize+5)
        rc('axes',      labelsize=fontsize+5)
        rc('xtick',     labelsize=fontsize+2)
        rc('ytick',     labelsize=fontsize+2)
        rc('text',      usetex=True)
        rc('font',      size=fontsize, style = 'normal', variant='normal', stretch='normal', weight='normal')
   
#         params = {'backend': 'wxAgg', 'lines.markersize' : 2, 'axes.labelsize': fontlabel_size, 'text.fontsize': fontlabel_size, 'legend.fontsize': fontlabel_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'text.usetex': True, 'figure.figsize': fig_size}
#         plt.rcParams.update(params)
   
        self.define_ColorVector(PlotStyle, n_colors, color_map)
   
        return
   
    def legend_conf(self, Axis, loc = 'best', fontize = None, edgelabel = False):
        #WARNING: THIS DOES NOT WORK WITH LEGEND RAVELIN
   
        if self.Axis1.get_legend_handles_labels()[1] != None:
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

        return
   
    def define_ColorVector(self, PlotStyle, n_colors, color_map):
   
        self.ColorVector = [None, None, None]
        
        self.ColorVector[0]        = 'w'
        self.ColorVector[1]        = 'b'
        
        if n_colors == None:
            self.ColorVector[2]    = cycle(sb.color_palette(color_map))
        else:
            self.ColorVector[2]    = sb.color_palette(color_map, n_colors)
        
    def Compute_SigmaLevel(self, trace1, trace2, nbins=20):
        
        # From a set of traces, bin by number of standard deviations
        L, xbins, ybins = histogram2d(trace1, trace2, nbins)
        L[L == 0]       = 1E-16
        logL            = nplog(L)
    
        shape           = L.shape
        L               = L.ravel()
    
        #Obtain the indices to sort and unsort the flattened array
        i_sort          = argsort(L)[::-1]
        i_unsort        = argsort(i_sort)
    
        L_cumsum        = L[i_sort].cumsum()
        L_cumsum        /= L_cumsum[-1]
        
        xbins           = 0.5 * (xbins[1:] + xbins[:-1])
        ybins           = 0.5 * (ybins[1:] + ybins[:-1])
    
        return xbins, ybins, L_cumsum[i_unsort].reshape(shape)
    
    def plot_SigmaLevels_MCMC(self, xdata, ydata, trace_x, trace_y, plot_Scatter = True, **kwargs):
        
        #"""Plot traces and contours"""
        xbins, ybins, sigma = self.Compute_SigmaLevel(trace_x, trace_y)
        self.Axis.contour(xbins, ybins, sigma.T, levels=[0.683, 0.955], **kwargs)
        
        if plot_Scatter:
            self.Axis.plot(trace_x, trace_y, ',k', alpha=0.1)

    def plot_posteriors_histagram(self, Traces, labels):
        
        n_traces = len(Traces)
        
        self.FigConf(Figtype = 'Posteriors', n_colors = n_traces)
        
        for i in range(len(Traces)):
            Trace           = Traces[i]
            Axis_Location   = int(str(n_traces) + '1' +  str(i + 1))
            Axis            = plt.subplot(Axis_Location) 
            HDP_coords      = self.statistics_dict[Trace]['95% HPD interval']
            
            for HDP in HDP_coords:
                Axis.axvline(x = HDP, label = 'HPD interval: ' + round_sig(HDP_coords[0], 4) + ' - ' + round_sig(HDP_coords[1], 4), color='grey', linestyle = 'dashed')
                                                
            Axis.axvline(x = self.statistics_dict[Trace]['mean'], label = 'Mean value: ' + round_sig(self.statistics_dict[Trace]['mean'], 4), color='grey', linestyle = 'solid')
            Axis.hist(self.Traces_dict[Trace][:], histtype='stepfilled', bins=35, alpha=.7, color=self.ColorVector[2][i], normed=False)
            Axis.set_ylabel(labels[i],fontsize=20)
                        
            self.legend_conf(Axis, loc='best', edgelabel=True)

#             Axis.yaxis.set_ticks(arange(0, 250, 50))

    def plot_Ownposteriors_histagram(self, Traces, labels, sigfig_n = 4, xlim = None, ylim = None):
        
        n_traces = len(Traces)
        
        self.FigConf(Figtype = 'Posteriors', n_colors = n_traces)
        
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
            Axis.set_ylabel(labels[i],fontsize=20)
                        
            self.legend_conf(Axis, loc='best', edgelabel=True)

#             Axis.yaxis.set_ticks(arange(0, 250, 50))
            
            if xlim != None:
                Axis.set_xlim(xlim[0], xlim[1])
                
            if ylim != None:
                Axis.set_ylim(ylim[0], ylim[1])

              
        return
    
    def plot_tracers(self, Traces, labels):
        
        n_traces = len(Traces)
        
        self.FigConf(Figtype = 'Posteriors', n_colors = n_traces) 
        
        for i in range(len(Traces)):
            Trace           = Traces[i]
            Axis_Location   = int(str(n_traces) + '1' +  str(i + 1))
            Axis            = plt.subplot(Axis_Location)         
            Variable_value  = self.statistics_dict[Trace]['mean']
            Variable_error  = self.statistics_dict[Trace]['standard deviation']
            label = labels[i] + ': ' + round_sig(Variable_value, 4) + r'$\pm$' + round_sig(Variable_error,4)            
            Axis.plot(self.pymc_database.trace(Trace)[:], label = label, color=self.ColorVector[2][i])
            Axis.axhline(y = self.statistics_dict[Trace]['mean'],  color=self.ColorVector[2][i], linestyle = '--' )
            Axis.set_ylabel(labels[i],fontsize=20)
            self.legend_conf(Axis, loc=1, edgelabel=True)
            
            plt.locator_params(axis = 'y', nbins = 5)
           
    def plot_acorrelation(self, Traces, labels, n_columns=4, n_rows=2):
        
        n_traces = len(Traces)
        
        self.FigConf(Figtype = 'Grid', n_colors = n_traces, n_columns=n_columns, n_rows=n_rows, FigHeight=9, FigWidth=16)         
        
        for i in range(len(Traces)):
            Trace   = Traces[i]                
            maxlags = min(len(self.pymc_database.trace(Trace)[:]) - 1, 100)
            label   = labels[i]
            if Trace != 'ChiSq':
                self.Axis1[i].acorr(x = self.pymc_database.trace(Trace)[:], color=self.ColorVector[2][i], detrend= detrend_mean, maxlags=maxlags)
            else:
                #Apano momentaneo
                chisq_adapted = reshape(self.pymc_database.trace(Trace)[:], len(self.pymc_database.trace(Trace)[:]))
                maxlags = min(len(chisq_adapted) - 1, 100)
                self.Axis1[i].acorr(x = chisq_adapted, color=self.ColorVector[2][i], detrend= detrend_mean, maxlags=maxlags)
            self.Axis1[i].set_xlim(0, maxlags)
            self.Axis1[i].set_title(label)
                
        return

    def plot_chiSq_Behaviour(self, Traces,  labels):
    
        n_traces = len(Traces)
        
        self.FigConf(Figtype = 'Grid', n_colors = n_traces, n_columns=4, n_rows=2, FigHeight=9, FigWidth=16)
        
        chisq_adapted   = reshape(self.pymc_database.trace('ChiSq')[:], len(self.pymc_database.trace('ChiSq')[:])) * -2
        y_lim           = 30
        min_chi_index   = argmin(chisq_adapted)
        
        for i in range(len(Traces)):
            Trace   = Traces[i]
            label   = labels[i]       

            if Trace != 'ChiSq':
                self.Axis1[i].scatter(x = self.pymc_database.trace(Trace)[:], y = chisq_adapted, color=self.ColorVector[2][i])
                x_min           = npmin(self.pymc_database.trace(Trace)[:]) 
                x_max           = npmax(self.pymc_database.trace(Trace)[:])

                self.Axis1[i].axvline(x = self.statistics_dict[Trace]['mean'], label = 'Inference value: ' + round_sig(self.statistics_dict[Trace]['mean'], 4,scien_notation=False), color='grey', linestyle = 'solid')
                self.Axis1[i].scatter(self.pymc_database.trace(Trace)[:][min_chi_index], chisq_adapted[min_chi_index], color='Black', label = r'$\chi^{2}_{min}$ value: ' + round_sig(self.pymc_database.trace(Trace)[:][min_chi_index],4,scien_notation=False))
                
                self.Axis1[i].set_ylabel(r'$\chi^{2}$',fontsize=20)
                self.Axis1[i].set_ylim(0, y_lim)
                self.Axis1[i].set_xlim(x_min, x_max)
                self.Axis1[i].set_title(label,fontsize=20)
                legend_i = self.Axis1[i].legend(loc='best', fontsize='x-large')
                legend_i.get_frame().set_facecolor('white')

    def plot_triangle_histContours(self, Traces, labels, true_values = None):
        
        n_traces = len(Traces)

        self.FigConf(Figtype = 'triangle')
        
        #Make plots
        List_Arrays = []
        for i in range(n_traces):
            List_Arrays.append(self.pymc_database.trace(Traces[i])[:])
        
        Samples     = array(List_Arrays).T
                    
        #Make 
        if true_values != None:
            self.figure = corner.corner(Samples[:,:], labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_args={"fontsize": 100}, truths = true_values, title_fmt='0.3f')
        else:
            self.figure = corner.corner(Samples[:,:], labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_args={"fontsize": 100}, title_fmt='0.3f')
                     
    def savefig(self, output_address, extension = '.png', reset_fig = True):
    
        if self.Triangle_Saver:
            self.figure.savefig(output_address + extension, dpi=500, bbox_inches='tight', pad_inches=0.2)
            
        else:
            plt.tight_layout()
            plt.savefig(output_address + extension, dpi=500, facecolor=self.ColorVector[0], bbox_inches='tight', pad_inches=0.2)
        
        if reset_fig:
            self.reset_fig()
    
        return
    
    def reset_fig(self):
        
        plt.cla()
        
        return
    
    
    
    
# def GetLower_ChiSquare(MCMC_Database, variable, mean = None, min_value = None, max_value = None, nbins = 100):
#     
#     print 'Treating variable', variable
#     
#     
#     ChiArray        =  Database_MCMC.trace('ChiSq')[:] * -1
#     VariableArray   = Database_MCMC.trace(variable)[:]
#         
#     if mean == None:
#         mean = average(VariableArray)
# 
#     if min_value == None:
#         min_value = 0.8
# 
#     if max_value == None:
#         max_value = 1.2
#         
#     Bins                        = linspace(mean*min_value, mean*max_value, nbins)
#     Points_Index                = digitize(VariableArray, Bins)
#          
#     Unique_Bins                 = unique(Points_Index)
#     Min_Index, Max_Index        = min(Unique_Bins), max(Unique_Bins)
#      
#     Min_values_Chi              = zeros(len(Points_Index))
#     Min_values_Variable         = zeros(len(Points_Index))
#      
#     for i in range(Min_Index, Max_Index):
#         inBinIndex              = where(Points_Index == i)[0]
#  
#         if len(ChiArray[inBinIndex]) != 0:
#             Min_values_Chi[i]       = np_min(ChiArray[inBinIndex])
#             Min_values_Variable[i]  = Bins[i]
#  
#     return Min_values_Variable, Min_values_Chi, mean