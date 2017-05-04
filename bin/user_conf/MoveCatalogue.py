'''
Created on Jul 23, 2015

@author: vital
'''

from numpy import where, array

from CodeTools.PlottingManager              import myPickle
from PipeLineMethods.ManageFlow             import DataToTreat


#Generate dazer object
pv                                          = myPickle()

Catalogue_Dic                               = DataToTreat()

#Define home and destination folders
OldFolder                                   = '/home/vital/Dropbox/Astrophysics/Data/WHT_HII_Galaxies_BackUp/'     
NewFolder                                   = Catalogue_Dic['Obj_Folder']

#Find and organize files according to the pattern you want to move
ImagePattern                                = '.png'
ImagesList                                  = pv.Folder_Explorer(ImagePattern,      OldFolder, CheckComputer=False)

LinesLogPattern                             = '_LickIndexes.txt' 
LickIndexesList                             = pv.Folder_Explorer(LinesLogPattern,   OldFolder, CheckComputer=False)

Masks_Pattern                               = '_Mask.lineslog'
MasksList                                   = pv.Folder_Explorer(Masks_Pattern,     OldFolder, CheckComputer=False)

#Get the list of scientific objects
Scientific_ObjectList                       = pv.get_TableColumn(['Objects'], TableAddress = Catalogue_Dic['Data_Folder'] + 'WHT_Galaxies_properties.txt', datatype = str)

#Loop through scientific objects
for i in range(len(Scientific_ObjectList)):
    
    CodeName = Scientific_ObjectList[i]
    
    ImageFolder = OldFolder + CodeName + '/'
    ImageName = CodeName + ImagePattern
    
#     print 'Looking for image', ImageName, '@', ImageFolder
    
    if (ImageFolder + ImageName) in ImagesList:
    
#         print 'Found copying to', Catalogue_Dic['Folder'] + 'Objects/' + CodeName + '/'
        pv.copyFile(ImageName, ImageFolder, Catalogue_Dic['Folder'] + 'Objects/' + CodeName + '/')
        
        
    LickFolder = OldFolder + CodeName + '/'
    LickName = CodeName + LinesLogPattern
    
#     print 'Looking for lickindexes', LickName, '@', LickFolder
    
    if (LickFolder + LickName) in LickIndexesList:
    
#         print 'Found copying to', Catalogue_Dic['Folder'] + 'Objects/' + CodeName + '/'
        pv.copyFile(LickName, LickFolder, Catalogue_Dic['Folder'] + 'Objects/' + CodeName + '/')
        
    MaskFolder = OldFolder + CodeName + '/'
    MaskName = CodeName + Masks_Pattern
    
#     print 'Looking for lickindexes', LickName, '@', LickFolder
    
    if (MaskFolder + MaskName) in MasksList:
    
#         print 'Found copying to', Catalogue_Dic['Folder'] + 'Objects/' + CodeName + '/'
        pv.copyFile(MaskName, MaskFolder, Catalogue_Dic['Folder'] + 'Objects/' + CodeName + '/')        
        
        
        
        