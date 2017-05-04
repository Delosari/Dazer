'''
Created on Dec 16, 2014

@author: delosari
'''

from os import remove
from os.path import isfile

import CodeTools.PlottingManager as plotMan
from PipeLineMethods.ManageFlow import DataToTreat


Pv  = plotMan.myPickle()

LogFiles_Extension  = ".plot"
CombinedFits        = "WHT.fits"
NebularFits         = "_Neb.fits"
StellarRemovedFits  = "_WHT_Neb_Star.fits"

StellarContinuum    = "_StellarContinuum.fits"

MaskFile            = "_mask.txt"

TexFiles            = ['.tex']

OldFiles            = ['_LinesLog_v2.txt', '_LinesLog_v3.txt', '.LinesLog_v3']

RootFolder          = DataToTreat()
Pattern             = [LogFiles_Extension, NebularFits, StellarContinuum, StellarRemovedFits, CombinedFits] + OldFiles + TexFiles

LogImages_Extension = ".png"

ForceDelete = True

#Find and organize files from terminal command or .py file
FilesList           = Pv.FindAndOrganize(Pattern, RootFolder, CheckComputer=True)

#Loop through files
for m in range(len(FilesList)):
    for j in range(len(FilesList[m])):
        
        CodeName, FileName, FileFolder = Pv.FileAnalyzer(FilesList[m][j], Texting=False)

        #Case of logs
        if LogFiles_Extension in FileName:
            
            LogName     = FileName
            ImageName   = LogName.replace(LogFiles_Extension, LogImages_Extension)
            
            #Deleting log
            if isfile(FileFolder + LogName):
                print '--', LogName
                
                if ForceDelete == True:
                    remove(FileFolder + LogName)                
            #Deleting images    
            if isfile(FileFolder + ImageName):
                print '--', ImageName
                
                if ForceDelete == True:
                    remove(FileFolder + ImageName)

        #Case of fits file        
        if 'fit' in FileName:
            print '\n-Fits Found'
            FitsFile = FileName
            print '--',FitsFile,'\n'
            
            #Deleting Fits file
            if ForceDelete == True:
                remove(FileFolder + FitsFile)
            
        #Case of line logs file        
        if ('LinesLog' in FileName) or ('.LinesLog_v3' in FileName):
            print '\n-Lines Log Found2'
            LinesLog = FileName
            print '--',LinesLog,'\n'
            
            #Deleting Fits file
            if ForceDelete == True:
                remove(FileFolder + LinesLog)     
                
                
        #Case of many tex files:
#         for tex_extension in TexFiles:
#             if tex_extension in FileName:
#                 print 'Tex file found3:'
#                 Tex_File = FileName
#                 print '--', Tex_File, '\n'
#                 
#                 #Deleting Fits file
#                 if ForceDelete == True:
#                     remove(FileFolder + Tex_File)
                      
            