from numpy import arange, sin, pi
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
import wx
import wx.lib.mixins.inspection as WIT
import pickle

class myCanvasFrame(wx.Frame):

    def __init__(self):
        
        wx.Frame.__init__(self, None, -1, 'myCanvasFrame', size=(550, 350))

    def example_calculation_plot(self):

        self.figure = Figure((8,6))
        self.axes = self.figure.add_subplot(111)
        t = arange(0.0, 3.0, 0.01)
        s = sin(2 * pi * t)
        self.axes.plot(t, s)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()
        self.add_toolbar()  # comment this out for no toolbar      
          
    def pickle_to_windows(self, pickle_address):
        
        self.figure = self.load_pickle_figure(pickle_address)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()

        self.add_toolbar()  # comment this out for no toolbar        
        
    def load_pickle_figure(self, pickle_address):
        
        with open(pickle_address, 'rb') as fid:
            fig = pickle.load(fid)
        
        return fig

    def add_toolbar(self):
        
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.update() # update the axes menu on the toolbar

class App_pickleWindow(wx.App): 
    
    def __init__(self, pickle_address): 
        wx.App.__init__(self) 
        frame = myCanvasFrame()
        #frame.example_calculation_plot()
        frame.pickle_to_windows(pickle_address)
        frame.Show(True)

if __name__ == '__main__':

    app = App_pickleWindow('/home/vital/Dropbox/Astrophysics/Data/WHT_observations_BackUp/data/oxygen_vs_sulfur_plot.dz_pickle')
    app.MainLoop()
