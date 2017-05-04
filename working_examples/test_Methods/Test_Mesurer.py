from DZ_DataExplorer                import Plots_Manager

from numpy                          import searchsorted as np_searchsorted, where


def Detect_LineData(CodeName, Folder, LineLabel):
       
    LineLogName = Folder +  CodeName + '_LickIndexes.txt'
       
    Labels = pv.get_ColumnData(['Line_Label'], LineLogName, datatype = str)
    TheoWavelengthC, Wave1C, Wave2C, Wave3C, Wave4C, Wave5C, Wave6C  = pv.get_ColumnData(['TheoWavelength', 'Wave1', 'Wave2', 'Wave3', 'Wave4', 'Wave5', 'Wave6'], LineLogName, datatype=float)
    
    Str_ind = where(Labels == LineLabel)
    
    return TheoWavelengthC[Str_ind][0], Wave1C[Str_ind][0], Wave2C[Str_ind][0], Wave3C[Str_ind][0], Wave4C[Str_ind][0], Wave5C[Str_ind][0], Wave6C[Str_ind][0]
  
# We declare the folder and log file to drop the lines data

pv = Plots_Manager()
  
# Forcing the remake of new files
pv.RemakeFiles = True
  
#Object and line to treat
Obj_Folder      = '/home/vital/Dropbox/Astrophysics/Data/WHT_Catalogue/Objects/SHOC579/'
Obj_Pattern     = 'objSHOC579_WHT.fits'
 
#Find and organize files from terminal command or .py file
FilesList = pv.Folder_Explorer(Obj_Pattern, Obj_Folder, CheckComputer=False)
  
#Set color format
pv.FigFormat_One(ColorConf='Night1')
  
#Identify file
CodeName, FileName, FileFolder = pv.Analyze_Address(FilesList[0])
  
#Load fits file data
Wave, Int, ExtraData = pv.File2Data(FileFolder, FileName)


#---------------------------Complete Mode--------------------------------------------------------------------
     
pv.RemakeFiles = False
pv.Wave = Wave
pv.Int = Int
pv.Current_Code = CodeName
Selections = []
 
Line_Label          =  'N2_6548A'
pv.Current_LinesLog = '/home/vital/Dropbox/Astrophysics/Data/WHT_Catalogue/Objects/51959-092/' +  CodeName + '_WHT_LinesLog_v3.txt'
 
Ion             = pv.GetDataInTable('Line_Label',"Ion", Line_Label, pv.Current_LinesLog, pv.TableHeaderSize)
TheoWavelength  = pv.GetDataInTable('Line_Label',"TheoWavelength", Line_Label, pv.Current_LinesLog, pv.TableHeaderSize, dataype='float')
 
pv.Current_Label, pv.Current_Ion, pv.Current_TheoLoc =  Line_Label, Ion, TheoWavelength
pv.Selections = pv.RangesR(Selections, pv.Current_TheoLoc, pv.Current_LinesLog)
   
pv.PlottingInLineMeasuring()
 
pv.DisplayFigure()
  
pv = Plots_Manager()


# #---------------------------Simple Mode--------------------------------------------------------------------

# #Load data into the class
# Line_Label      =  'O2_3726A'
# Line_Label    = 'He1_4472A'
# #["Line_Label", "Ion", "TheoWavelength", "Blended_Set","Wave1", "Wave2", "Wave3", "Wave4", "Wave5", "Wave6"]
# pv.Current_Label, pv.Current_Ion                                                = Line_Label, 'NoImporta'
# pv.Current_TheoLoc, pv.Wave1, pv.Wave2, pv.Wave3, pv.Wave4, pv.Wave5, pv.Wave6  = Detect_LineData(CodeName, FileFolder, Line_Label)
# pv.Selections                                                                   = [pv.Wave1, pv.Wave2, pv.Wave3, pv.Wave4, pv.Wave5, pv.Wave6]
# Ind_min                                                                         = np_searchsorted(Wave, pv.Wave1)
# Ind_max                                                                         = np_searchsorted(Wave, pv.Wave6)
# SubWave                                                                         = Wave[Ind_min - 1 : Ind_max + 1]             
# SubInt                                                                          = Int[Ind_min - 1 : Ind_max + 1]
# pv.ind1, pv.ind2, pv.ind3, pv.ind4, pv.ind5, pv.ind6                            = np_searchsorted(SubWave, pv.Wave1), np_searchsorted(SubWave, pv.Wave2), np_searchsorted(SubWave, pv.Wave3), np_searchsorted(SubWave, pv.Wave4), np_searchsorted(SubWave, pv.Wave5), np_searchsorted(SubWave, pv.Wave6)
# #  
# #Calculation and plotting of continuum properties     
# pv.LocalMedian, pv.ContinuumFlux, pv.SigmaContinuum, pv.Continuum_Gradient, pv.Continuum_n, WPoint1, Wpoint2, FPoint1, FPoint2 = pv.ContinuumRegions(SubInt, SubWave, adjacent_continuum = True)
#    
# #Continuum calculation properties
# pv.Check_GaussianMixture(SubWave, force_check = True)   
#    
# #Convert emission lines to gaussian equivalents
# GaussianWave, GaussianInt = pv.Measure_LineIntensity(SubWave, SubInt, Methodology = 'MCMC_kmpfit', SaveData = False)  
#    
# pv.SaveLM_to_file(pv.Current_TheoLoc, '/home/vital/Desktop/SHOC579_linelog.txt', True)
#    
# #Plot the fitted line
# pv.DataPloter_One(SubWave, SubInt, LineLabel='Obj SHOC579', LineColor = 'y')
# #   
# # if pv.Fitting_dict['blended number'] > 1:
# #     for i in range(pv.Fitting_dict['blended number']):
# #         pv.DataPloter_One(GaussianWave, pv.Fitting_dict['y_resample'][i], LineLabel='Gaussian fitting', LineColor = 'r')
# # else:
# #     pv.DataPloter_One(GaussianWave, pv.Fitting_dict['y_resample'], LineLabel='Gaussian fitting', LineColor = 'r')
#  
# pv.DataPloter_One(GaussianWave, pv.Fitting_dict['y_resample_total'], LineLabel='Gaussian fitting', LineColor = 'r')
#  
# pv.DisplayFigure()



  
print 'Se acabo'