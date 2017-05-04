import os

import pyfits
from uncertainties import unumpy, ufloat  # Array manipulation

from CodeTools.vitools import GetParameterFromDigitalLines
import Scientific_Lib.AstroMethods as astro
import math as m
import vitools as vit


def SetLogFile(LogFile, LogFolder, LogFileFormat = "/Users/INAOE_Vital/Dropbox/Astrophysics/Data/WHT_HII_Galaxies/WHT_Log_Format.txt", ForceRemake = False):

    Structure_Log = vit.File2Table(LogFileFormat, "")

    #We import the data from previous log. If ForceRemake flag is on, we make new files
    if os.path.isfile(LogFolder + LogFile) and (ForceRemake == False):
        
        Obj_Log = vit.File2Table(LogFolder, LogFile)

        for i in range(1,len(Structure_Log)):
            for j in range(len(Obj_Log)):
                if Structure_Log[i][0] == Obj_Log[j][0]:
                    Structure_Log[i][1] = Obj_Log[j][1]
                    Structure_Log[i][2] = Obj_Log[j][2]
                    Structure_Log[i][3] = Obj_Log[j][3]
                    
    OutputFile = open(LogFolder + LogFile, 'w')
    ColumnSizes = []
    
    for i in range(len(Structure_Log[0])):
        Biggest_Length = 0
        for j in range(1,len(Structure_Log)):
            if len(Structure_Log[j][i]) > Biggest_Length:
                Biggest_Length = len(Structure_Log[j][i])
        ColumnSizes.append(Biggest_Length + 2)
        
    NewFormatLine = "%" + str(ColumnSizes[0]) + "s" + "%" + str(ColumnSizes[1]) + "s" + "%" + str(ColumnSizes[2]) + "s" + "%" + str(ColumnSizes[3]) + "s"
        
    for z in range(1, len(Structure_Log)):  
        NewLine = NewFormatLine % (Structure_Log[z][0],Structure_Log[z][1],Structure_Log[z][2],Structure_Log[z][3])
        OutputFile.write(NewLine+"\n")

    OutputFile.close()
                
    return

def Parameter_2_Log(FileFolder, ObjectCode, Parameter, Parameter_Value, Parameter_Error = 'none'):
    
    Magnitude = str(Parameter_Value)
    
    LogExtension = "_log.txt"
    
    Log_File_Address = FileFolder + ObjectCode + LogExtension
    
    Obj_Log = vit.File2Table(Log_File_Address, "")
        
    Index = "None"

    for i in range(len(Obj_Log)):
        Item = Obj_Log[i][0]
        if Parameter == Item:
            Index = i
    
    if Index == "None":
        print "WARNING: The parameter cannot be saved since its category cannot be found in object log"
        print "Parameter name ", Parameter
        print "Parameter Magnitude", Magnitude
        print "Object Log"
        for i in range(len(Obj_Log)):
            print Obj_Log[i][0]
           
    Obj_Log[Index][1] = Magnitude
         
    OutputFile = open(Log_File_Address, 'w')
    
    
    ColumnSizes = []
    for i in range(len(Obj_Log[0])):
        Biggest_Length = 0
        for j in range(1,len(Obj_Log)):
            if len(Obj_Log[j][i]) > Biggest_Length:
                Biggest_Length = len(Obj_Log[j][i])
        ColumnSizes.append(Biggest_Length + 2)
        
    NewFormatLine = "%" + str(ColumnSizes[0]) + "s" + "%" + str(ColumnSizes[1]) + "s" + "%" + str(ColumnSizes[2]) + "s" + "%" + str(ColumnSizes[3]) + "s"
    
    for z in range(len(Obj_Log)):  
        NewLine = NewFormatLine % (Obj_Log[z][0],Obj_Log[z][1],Obj_Log[z][2],Obj_Log[z][3])
        OutputFile.write(NewLine+"\n")

    OutputFile.close()
                
    return

def Log_2_Parameter(FileFolder, ObjectCode, ParameterToFind):
    
    LogExtension = "_log.txt"
     
    Log_File_Address = FileFolder + ObjectCode + LogExtension
    
    print 'Analyzing object:'
    print FileFolder,'|',  ObjectCode,'|', LogExtension
    
    Obj_Log = vit.File2Table(Log_File_Address, "")
    
    Index = "None"

    for i in range(len(Obj_Log)):
        Item = Obj_Log[i][0]
        if ParameterToFind == Item:
            Index = i
    
    if Index == "None":
        print "WARNING: The parameter cannot be found in log"
        print ParameterToFind
        print Obj_Log
        print FileFolder + ObjectCode + LogExtension
        
    Magnitude = Obj_Log[Index][1]
                
    return Magnitude

def Data_2_Fits(FileFolder, FitsName, Header, Wavelength, Intensity, NewKeyWord = None):
    
    Primary_HDU = pyfits.PrimaryHDU(header=Header)

    if NewKeyWord != None:
        Primary_HDU.header.update(NewKeyWord[0], NewKeyWord[1])

    Column1 = pyfits.Column(name='Wave', format='E',array=Wavelength)
    Column2 = pyfits.Column(name='Int', format='E', array=Intensity)
    Columns = pyfits.ColDefs([Column1, Column2])
    Table_HDU = pyfits.new_table(Columns)
    
    HDU_list = pyfits.HDUList([Primary_HDU, Table_HDU])
        
    HDU_list.writeto(FileFolder + FitsName, clobber = True)
    
    return
    
def GetLineFlux(Line_Label, FolderName, FileName, CodeName, LogFileFormat = '_WHT_Neb_Star_LinesLog_v3.txt'):
    
    R_V = 3.2
    
    Valid_Flag = True

    LogFileName = CodeName + LogFileFormat
    LogLines = vit.File2Lines(FolderName, LogFileName)
    
    Flux_Brute = GetParameterFromDigitalLines(LogLines,Line_Label,"FluxBrute",2)
    Flux_Error = GetParameterFromDigitalLines(LogLines,Line_Label,"ErrorEL",2)
    
    Wavelength = GetParameterFromDigitalLines(LogLines,Line_Label,"TheoWavelength",2)
    Wave_Error = 1.0

    if (Flux_Brute != None) and (Flux_Error != None):
        Flux_V =  ufloat(float(Flux_Brute), float(Flux_Error))
        Wave_V =  ufloat(float(Wavelength), float(Wave_Error))
         
        f_HBeta = astro.Reddening_f_Cal(4861.333, R_V)
        f_lambda = astro.Reddening_f_Cal(4861.333, R_V)
    
        cHBeta = float(astro.Log_2_Parameter(FolderName, CodeName, "cHBeta"))
        
        if cHBeta < 0:
            Flux_Deredd = Flux_V
        else:
            Flux_Deredd = Flux_V * 10**(cHBeta*(f_lambda-f_HBeta))
        
        LineData = Flux_Deredd
        
#         LineData = [Line_Label, Valid_Flag, Flux_Deredd, Wave_V]
        
    else:
        Valid_Flag = False
#         LineData = [Line_Label, Valid_Flag, (Flux_Brute,Flux_Error), (Wavelength,Wave_Error)]
        LineData = None
        
    return LineData

    
    
    
    
    
    
    
