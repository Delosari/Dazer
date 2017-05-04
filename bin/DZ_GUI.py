'''
Created on May 2, 2015

@author: vital
'''

import ttk
import Tkinter as tk
from matplotlib.backends.backend_tkagg  import FigureCanvasTkAgg, NavigationToolbar2TkAgg

LARGE_FONT  = ("Ubuntu", 12)
NORM_FONT   = ("Ubuntu", 10)
SMALL_FONT  = ("Ubuntu", 8)

class Tk_GUI(tk.Tk):
        
    def __init__(self, PlottingVector,*args, **kwargs):
        
        self.Bg_Color = 'white'        
        self.Fg_Color = 'black'
                
        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.wm_title(self, "Dazer")
        
        tk.Tk.configure(self,bg =self.Bg_Color)
        
        container = tk.Frame(self)

        container.grid(row = 0)
        
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.Menubar_Definition(container)

        self.frames = {}

        #for Screen in [PageInitial, ]:
        for Screen in [Page_NebularContinuum, PageInitial]:
            frame = Screen(container, self, PlottingVector)

            self.frames[Screen] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(PageInitial)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()
        
    def Menubar_Definition(self, container):
        
        menubar = tk.Menu(container)
        
        File_Menu = tk.Menu(menubar, tearoff=0)
        File_Menu.add_command(label="Save settings", command = lambda: popupmsg("Not supported yet"))
        File_Menu.add_separator()
        File_Menu.add_command(label="Exit", command=quit)
        menubar.add_cascade(label="File", menu=File_Menu)
        
        Catalogue_Menu = tk.Menu(menubar, tearoff=0)
        #Catalogue_Menu.add_command(label="Catalogue view", command = lambda: self.show_frame(PageInitial))        
        Catalogue_Menu.add_command(label="Define Catalogue", command = lambda: popupmsg("Not supported yet"))
        Catalogue_Menu.add_command(label="Switch Catalogue", command = lambda: popupmsg("Not supported yet"))
        menubar.add_cascade(label="Catalogue", menu=Catalogue_Menu)       
 
        Nebular_Menu = tk.Menu(menubar, tearoff=0)
        Nebular_Menu.add_command(label="Calculate object nebular continuum", command = lambda: popupmsg("Not supported yet"))
        Nebular_Menu.add_command(label="Simulate nebular emission", command = lambda: self.show_frame(Page_NebularContinuum))
        menubar.add_cascade(label="Nebular Emission", menu=Nebular_Menu)       

        Starlight_Menu = tk.Menu(menubar, tearoff=0)
        Starlight_Menu.add_command(label="Resample spectrum",       command = lambda: popupmsg("Not supported yet"))
        Starlight_Menu.add_command(label="Check object masks",      command = lambda: popupmsg("Not supported yet"))
        Starlight_Menu.add_command(label="Set configuration file",  command = lambda: popupmsg("Not supported yet"))        
        Starlight_Menu.add_command(label="Launch Starlight",        command = lambda: popupmsg("Not supported yet"))
        Starlight_Menu.add_command(label="Check output file",       command = lambda: popupmsg("Not supported yet"))
        menubar.add_cascade(label="Starlight", menu=Starlight_Menu)  

        Starlight_Menu = tk.Menu(menubar, tearoff=0)
        Starlight_Menu.add_command(label="Check bases",       command = lambda: popupmsg("Not supported yet"))
        Starlight_Menu.add_command(label="Define catalogue library bases",      command = lambda: popupmsg("Not supported yet"))
        menubar.add_cascade(label="PopStar", menu=Starlight_Menu)  

        EmLineMesurer_Menu = tk.Menu(menubar, tearoff=0)
        EmLineMesurer_Menu.add_command(label="Define emission lines location",       command = lambda: popupmsg("Not supported yet"))
        EmLineMesurer_Menu.add_command(label="Import emission lines location",       command = lambda: popupmsg("Not supported yet"))
        EmLineMesurer_Menu.add_command(label="Set measuring scheme",                 command = lambda: popupmsg("Not supported yet"))
        EmLineMesurer_Menu.add_command(label="Measure emission lines for catalgue",  command = lambda: popupmsg("Not supported yet"))
        EmLineMesurer_Menu.add_command(label="Check object lines log",               command = lambda: popupmsg("Not supported yet"))
        menubar.add_cascade(label="Emission line measurer", menu=EmLineMesurer_Menu)  


        Pyneb_Menu = tk.Menu(menubar, tearoff=0)
        Pyneb_Menu.add_command(label="Define atomic parameters",                    command = lambda: popupmsg("Not supported yet"))
        Pyneb_Menu.add_command(label="Generate diagnostics plot",                   command = lambda: popupmsg("Not supported yet"))
        Pyneb_Menu.add_command(label="Define Electronic temperature/density",       command = lambda: popupmsg("Not supported yet"))
        Pyneb_Menu.add_command(label="Measure Abundances",                          command = lambda: popupmsg("Not supported yet"))
        menubar.add_cascade(label="PyNeb", menu=Pyneb_Menu)  

        Cloudy_Menu = tk.Menu(menubar, tearoff=0)
        Cloudy_Menu.add_command(label="Launch script",                               command = lambda: popupmsg("Not supported yet"))
        Cloudy_Menu.add_command(label="Plot script",                                 command = lambda: popupmsg("Not supported yet"))
        menubar.add_cascade(label="Cloudy", menu=Cloudy_Menu)  

        tk.Tk.config(self, menu=menubar)   
             
class PageInitial(tk.Frame):

    def __init__(self, parent, controller, PlottingVector):
        
        #self.Bg_Color = "#%02x%02x%02x" % (int(PlottingVector.Color_Vector[0][0]*255), int(PlottingVector.Color_Vector[0][1]*255), int(PlottingVector.Color_Vector[0][2]*255))        
        #self.Fg_Color = "#%02x%02x%02x" % (int(PlottingVector.Color_Vector[1][0]*255), int(PlottingVector.Color_Vector[1][1]*255), int(PlottingVector.Color_Vector[1][2]*255))

        self.Bg_Color = 'white'       
        self.Fg_Color = 'black'
                
        tk.Frame.__init__(self, parent,)
        #label = tk.Label(self, text="Catalogue view", bg = self.Bg_Color, fg = self.Fg_Color , font=LARGE_FONT)

        #label.grid(row=0, pady=10,padx=10)

        PlottingVector.FigCanvas = FigureCanvasTkAgg(PlottingVector.Fig, self)
             
        PlottingVector.FigCanvas.get_tk_widget().grid(row=1)
              
class Page_NebularContinuum(tk.Frame):
    
    def __init__(self, parent, controller, PlottingVector):

        #self.Bg_Color = "#%02x%02x%02x" % (int(PlottingVector.Color_Vector[0][0]*255), int(PlottingVector.Color_Vector[0][1]*255), int(PlottingVector.Color_Vector[0][2]*255))        
        #self.Fg_Color = "#%02x%02x%02x" % (int(PlottingVector.Color_Vector[1][0]*255), int(PlottingVector.Color_Vector[1][1]*255), int(PlottingVector.Color_Vector[1][2]*255))

        self.Bg_Color = 'white'       
        self.Fg_Color = 'black'
                
        tk.Frame.__init__(self, parent,)
        label = tk.Label(self, text="Nebular continuum calculation", bg = self.Bg_Color, fg = self.Fg_Color , font=LARGE_FONT)
        label.grid(row = 0, pady=10,padx=10)

        PlottingVector.FigCanvas = FigureCanvasTkAgg(PlottingVector.Fig, self)
                
        PlottingVector.FigCanvas.get_tk_widget().grid(row = 1)
        
def popupmsg(msg):
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    popup.mainloop()
