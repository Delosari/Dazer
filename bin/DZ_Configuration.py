from numpy import loadtxt
from os import path, getcwd

class ImportConfiguration():
    
    def __init__(self):
        
        Root_Folder             = '/home/vital/'
 
        self.Catalogue_Folder = Root_Folder + "Dropbox/Astrophysics/Data/WHT_observations/objects/"
        self.Flow_Configuration(Root_Folder, DefaultSpectra_Preffix = '', DefaultSpectra_Suffix='_WHT.fits')
        self.DataType = 'WHT'
        self.Pattern_PlotFiles  = ".dz_pickle" 

#         self.Catalogue_Folder   = Root_Folder + "Dropbox/Astrophysics/Data/WHT_CandiatesObjects/"
#         self.Flow_Configuration(Root_Folder, DefaultSpectra_Preffix = 'spec', DefaultSpectra_Suffix='_dr10.fits')
#         self.DataType = 'dr10'

#         self.Catalogue_Folder   = Root_Folder + "Dropbox/Astrophysics/Data/WHT_Catalogue_SulfurRegression/Objects/"
#         self.Flow_Configuration(Root_Folder, DefaultSpectra_Preffix = 'obj', DefaultSpectra_Suffix='_WHT.fits')
#         self.DataType = 'WHT'
#         self.Pattern_PlotFiles  = ".plot"

#         self.Catalogue_Folder   = Root_Folder + "Dropbox/Astrophysics/Data/Fabian_Catalogue/27/"
#         self.Flow_Configuration(Root_Folder, DefaultSpectra_Preffix = 'spec', DefaultSpectra_Suffix='_dr10.fits')
#         self.DataType = 'dr10'

#         self.Catalogue_Folder   = Root_Folder + "Dropbox/Astrophysics/Data/ToCompare/"
#         self.Flow_Configuration(Root_Folder, DefaultSpectra_Preffix = 'spec', DefaultSpectra_Suffix='_dr10.fits')
#         self.DataType = 'dr10'

#         self.Catalogue_Folder   = Root_Folder + 'Dropbox/Astrophysics/Data/Starburst_Spectra_z0.004/'
#         self.Flow_Configuration(Root_Folder, DefaultSpectra_Preffix = '', DefaultSpectra_Suffix='.txt')
#         self.DataType           = 'txt'
#         self.Pattern_PlotFiles  = ".txt"

#         self.Catalogue_Folder   = Root_Folder + 'Dropbox/Astrophysics/Data/Popstar_spectra_z0.004/'
#         self.Flow_Configuration(Root_Folder, DefaultSpectra_Preffix = '', DefaultSpectra_Suffix='.txt')
#         self.DataType           = 'txt'
#         self.Pattern_PlotFiles  = ".txt"

    def Flow_Configuration(self, Root_Folder, Configuration ='default', DefaultSpectra_Preffix = 'obj', DefaultSpectra_Suffix = '_WHT.fits'):

        self.SpectraPreffix                 = DefaultSpectra_Preffix
        self.SpectraSuffix                  = DefaultSpectra_Suffix

        #WARNING need a better way to find bin folder
        __location__ = path.realpath(path.join(getcwd(), path.dirname(__file__)))
        root_folder  = __location__[0:__location__.find('/bin')+1]

        self.Conf_Folder                    = root_folder + 'format/'
        self.LinesLogHeader_Name            = 'DZT_LineLog_Headers.dz'
        self.LinesLogExtension_Name         = '_linesLog_reduc.txt'

        self.RemakeFiles                    = False
        self.HistogramTypeSpec              = True
                                
        self.PlotFormat_Vector              = [{'boxstyle':'round', 'facecolor':'#4c4c4c'}, {'boxstyle':'round', 'facecolor':"red"}, {"alpha":0.5}]
        
        EmissionLineDatabaseFolder          = root_folder + 'format/'
        EmissionLineDatabaseFile            = 'emlines_pyneb_optical_infrared.dz'

        self.Labels_Total, self.Ions_Total, self.labels_format = loadtxt(EmissionLineDatabaseFolder + EmissionLineDatabaseFile, dtype=str, usecols = (0,1,3), unpack = True)
        self.Wavelengths_Total              = loadtxt(EmissionLineDatabaseFolder + EmissionLineDatabaseFile, dtype=float, usecols = [2])        
        
        
        
        
        
        
        
        
        
        