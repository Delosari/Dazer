from DZ_Configuration                             import ImportConfiguration
from DZ_LineMesurer                               import LineMesurer
from DZ_Logs                                      import LineMesurer_Log
from user_conf.ManageFlow                         import DataToTreat
from lib.CodeTools.various                        import vitools, Error_manager
from lib.Math_Libraries.fitting_methods           import Fitter
from lib.CodeTools.File_Managing_Tools            import File_Manager
from lib.Plotting_Libraries.dazer_plotter         import Plot_Conf
from lib.Astro_Libraries.Reddening_Corrections    import ReddeningLaws
from lib.Astro_Libraries.Abundances_Class         import Chemical_Analysis_pyneb
from DZ_observation_reduction                     import spectra_reduction


class Dazer(Plot_Conf, File_Manager, Fitter, ImportConfiguration, LineMesurer, LineMesurer_Log, vitools, Error_manager, ReddeningLaws, Chemical_Analysis_pyneb, spectra_reduction):

    def __init__(self):
                
        #----Load lib:
        
        #For treating files
        File_Manager.__init__(self)
        
        #For generating graphs
        Plot_Conf.__init__(self)
                
        #Scientific operations
        Fitter.__init__(self)

        #Astro_libraries
        ReddeningLaws.__init__(self)
        Chemical_Analysis_pyneb.__init__(self)
        #spectra_reduction.__init__(self) #WARNING: I need a new place to put this
        
        #Load user configuration
        ImportConfiguration.__init__(self)
        LineMesurer.__init__(self, self.Conf_Folder, self.LinesLogHeader_Name)        
        
        #Extra functions
        vitools.__init__(self)
        Error_manager.__init__(self)
             
    def import_catalogue(self):
        
        Catalogue = DataToTreat()
        
        return Catalogue