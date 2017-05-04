'''
Created on Jan 29, 2014

@author: vital
'''

# import mimetypes
# import os
# import sys
# 
# from numpy import loadtxt, where, delete

# def Arguments_Checker(LineaCommando):
#     Num_Argum = len(LineaCommando)
#     
#     Arguments = []
#     
#     if Num_Argum == 1:
#         ArgumentsCheck = False
#         Arguments.append(LineaCommando[0])
#         return ArgumentsCheck, Arguments
#     
#     if Num_Argum > 1:
#         for i in range(len(LineaCommando)):
#             ArgumentsCheck = True
#             Arguments.append(LineaCommando[i])
#         return ArgumentsCheck, Arguments
#     
# def Arguments_Handler(ArgumentsCheck, Arguments):   
#     Files = []
#     FlagsList = []
#     Folder = os.getcwd()
#     
#     if ArgumentsCheck:
#         for i in range(1, len(Arguments)):
#             if "--" in Arguments[i]:
#                 FlagsList.append(Arguments[i])
#             else:
#                 Files.append(Folder + '/' + Arguments[i])
#         
#         Num_Files = len(Files)
#         Num_Flags = len(FlagsList)
#          
#         if Num_Files > 0:
#             print "-Files to treat:"
#             print Files
#         if Num_Flags > 0:
#             print "-FlagsList activated:"
#             print FlagsList
#         
#         return tuple(Files), FlagsList
#     
#     else:
#         return tuple(Files), FlagsList
# 
# def SaveXYTable_2_Text(X,Y,AddressAndName):
# 
#     if len(X) == len(Y):
#     
#         File = open(AddressAndName,"w")
#         
#         for i in range(len(X)):
#             x_i = X[i]
#             y_i = Y[i]
#             Element = str(x_i) + " " + str(y_i) + "\n"
#             File.write(Element)
#         
#         File.close()
#         
#     else:
#         
#         print "WARNING: File ", AddressAndName
#         print 'cannot be created: X and Y parameters do not have the same length'
# 
#     return
# 
# def AddressAnalyzer(FileAddress):    
#     FolderName = FileAddress[0:FileAddress.rfind("/")+1]
#     FileName = FileAddress[FileAddress.rfind("/")+1:len(FileAddress)]
#     CodeName = FolderName[FolderName[0:-1].rfind("/")+1:len(FolderName)-1]
# 
#     if FileName.startswith("obj") or FileName.startswith("std"):   
#         ObjName = FileName[3:FileName.find("_")]
#         return ObjName, FileName, FolderName
#     else:
#         ObjName = "None"
#         return CodeName, FileName, FolderName
#              
# def FilesFinder(Folder,myPattern):   
#     FileList = []
#         
#     if type(myPattern) is not list:
#         myPatternList = [myPattern]
#         
#     else:
#         myPatternList = myPattern
#     
#     for Root, Dirs, Archives in os.walk(Folder):
#         for Archive in Archives:
#             Meets_one_Pattern = False
#             for i in range(len(myPatternList)):
#                 if (myPatternList[i] in Archive):
#                     if "~" not in Archive:                   
#                         Meets_one_Pattern = True
#                         
#             if Meets_one_Pattern:        
#                 if Root.endswith("/"):
#                     FileList.append(Root + Archive)
#                     print "- " + Archive
#                 else:
#                     FileList.append(Root + "/" + Archive)
#                     print "- " + Archive
#                  
#     print "\n"  
#                       
#     return FileList
# 
# def FilesFinder_OR(Folder,myPattern):   
#     FileList = []
#         
#     print "Files found meeting one of the patterns: " + ", ".join(myPattern) + "\n"
#     print "in folder: " + Folder
#     
#     for Root, Dirs, Archives in os.walk(Folder):
#         for Archive in Archives:
#             Meets_one_Pattern = False
#             for i in range(len(myPattern)):
#                 if (myPattern[i] in Archive):
#                     if "~" not in Archive:                   
#                         Meets_one_Pattern = True
#                         
#             if Meets_one_Pattern:        
#                 if Root.endswith("/"):
#                     FileList.append(Root + Archive)
#                     print "- " + Archive
#                 else:
#                     FileList.append(Root + "/" + Archive)
#                     print "- " + Archive
#     
#                         
#     print "\n"  
#     return FileList
# 
# def FilesCounter(Folder,myPattern): 
#     
#     Counter = 0
# 
#     if myPattern is not list:         
#              
#         for Root, Dirs, Archives in os.walk(Folder):
#             for Archive in Archives:
#                 if (myPattern in Archive) and ("~" not in Archive):
#                     Counter = Counter + 1
#                     
#                                                  
#     print "Files meeting conditions: ", Counter  
#     return Counter
# 
# def FileRenamer(Input_FileAddress,Old_Pattern,New_Pattern):
#     
#     InObjName, InFileName, InFolderName = AddressAnalyzer(str(os.path.abspath(Input_FileAddress)))
# 
#     if Old_Pattern != "None":
#         OutFileName = InFileName.replace(Old_Pattern,New_Pattern)
#     
#     else:
#         OutFileName = New_Pattern
#             
#     Output_FileAddress = InFolderName + OutFileName
#     
# #     print "Old name", Input_FileAddress
# #     print "New name", Output_FileAddress
#     
#     os.rename(Input_FileAddress,Output_FileAddress)
#     
#     return
# 
# def File2Lines(Address,File):
#     File = Address + File
#     TextFile = open(File,"r")
#     Lines = TextFile.readlines()
#     TextFile.close()
#     return Lines
# 
# def File2Table(Address,File):
#     File = Address + File
#     TextFile = open(File,"r")
#     Lines = TextFile.readlines()
#     TextFile.close()
#     
#     Table = []
#     for line in Lines:
#         Table.append(line.split())
#             
#     return Table    
# 
# def List2Columns(FileLines,ColumnName,HeaderSize):
#     
#     HeadersLine = FileLines[HeaderSize].split()
#     HeaderIndex = "None"
#     
#     for i in range(len(HeadersLine)):
#         if ColumnName in HeadersLine[i]:
#             HeaderIndex = i
#     
#     if HeaderIndex == "None":
#         print "WARNING: La columna --" + ColumnName + "-- no aparece en el archivo"
#     
#     Column = []
#     
#     for i in range(HeaderSize,len(FileLines)):
#         Line = FileLines[i].split()
#         Column.append(Line[HeaderIndex])
#     return Column
# 
# def List2Text(FileAddress,List):
# 
#     print 'Nombre del archivo', FileAddress
#     print 'Y es:', type(FileAddress)
#     File = open(FileAddress,"w")
# 
#     for element in range(len(List)):
#         File.write(str(List[element]))
#         
#     File.close()
# 
# def ListTraveler(Current_Position, List, InputKey, Next_Key = "right", Previous_Key = "left"):
#     #Warning the first value of Current_Position must be None
#     
#     if Current_Position == None:         # This works as a flag for the first position on the list 
#         Index = 0
#         Current_Position = 0
#         
#     else:
#         
#         if InputKey == Next_Key:
#             Index = 1
#         
#         if InputKey == Previous_Key:
#             Index = -1
#             
#         if (InputKey != Previous_Key) and (InputKey != Next_Key) :        
#             Index = 0
#             print "Key", InputKey, "has not action assigned"
#     
#     NewPosition = Current_Position + Index
#     List_Length = len(List)
# 
#     if NewPosition >= List_Length: 
#         NewPosition = List_Length -1
#         print "End of the list reached"
#     
#     if NewPosition < 0:
#         NewPosition = 0
#         print "First point of the list"
#     
#     return NewPosition
# 
# def ListTraveler_2(Current_Position, List, InputCommand, Next_Command = "NextObject", Previous_Command = "PreviousObject"):
#     
#     #Warning the first value of Current_Position must be None
#     if Current_Position == None:         # This works as a flag for the first position on the list 
#         Index = 0
#         Current_Position = 0
#         
#     else:
#         
#         if InputCommand == Next_Command:
#             Index = 1
#         
#         if InputCommand == Previous_Command:
#             Index = -1
#             
#         if (InputCommand != Previous_Command) and (InputCommand != Next_Command) :        
#             Index = 0
#             print "Key has not action assigned"
#     
#     NewPosition = Current_Position + Index
#     List_Length = len(List)
# 
#     if NewPosition >= List_Length: 
#         NewPosition = List_Length -1
#         print "End of the list reached"
#     
#     if NewPosition < 0:
#         NewPosition = 0
#         print "First point of the list"
#     
# #     print NewPosition, '/', len(List)
#     return NewPosition
# 
# def Lines2NewTextFile(FileAdress,FileName,Lines):
#     TextFile = open(FileAdress + FileName,"w")
#     for element in range(len(Lines)):
#         TextFile.write(str(Lines[element]))    
#     TextFile.close()
#     return   
#               
# def CheckListForItem(List,Item):
#     ItemFound = False
#     for Value in List:
#         if Item in Value:
#             ItemFound = True
#     return ItemFound
# 
# def CheckTableForItem(ListHeader,Header,ListValues,NonDesiredValue):
#     ItemFound = False
#     RowIndex = "None"
#     for i in range(len(ListHeader)):
#         if Header in ListHeader[i]:
#             RowIndex = i
#     
#     if RowIndex == NonDesiredValue:
#         print "WARNING: La fila --" + Header + "-- no existe"
#     
#     
#     if ListValues[RowIndex] != "None":
#         ItemFound = True
#         
#     return ItemFound
# 
# def GetValueFromTable(ListHeader,Header,ListValues):
#     RowIndex = "None"
#     for i in range(len(ListHeader)):
#         if Header in ListHeader[i]:
#             RowIndex = i
#     
#     if RowIndex == "None":
#         print "WARNING: La fila --" + Header + "-- no existe"
#     
#     Item = ListValues[RowIndex]
#     
#     if Item == None or Item == "None":
#         print "WARNING: No se puede extraer el --" + Header +"-- de la tabla" 
#  
#         
#     return Item
# 
# def GetParameterFromTable(FolderName,FileName,Row_Label,Column_Label,HeaderSize,Row_Existance_Necesary=False):
#     
#     TextLines = File2Lines(FolderName,FileName)
#     
#     Parameter = None
#     HeaderIndex = None
#     
#     HeaderLine = TextLines[HeaderSize-1].split()
#     
#     for i in range(len(HeaderLine)):
#         item = HeaderLine[i]
#         if item == Column_Label:
#             HeaderIndex = i
#     
#     if HeaderIndex == None:
#         print "WARNING: Header not found"
#         print Column_Label
#         print HeaderLine
#         
#     for i in range(HeaderSize,len(TextLines)):
#         LineElements = TextLines[i].split()
#         if LineElements[0] == Row_Label:
#             Parameter = LineElements[HeaderIndex]
# 
#     return Parameter
#     
# def GetParameterFromDigitalLines(LogLines,Row_Label,Column_Label,HeaderSize,Row_Existance_Necesary=False):
#         
#     Parameter = None
#     HeaderIndex = None
#     
#     HeaderLine = LogLines[HeaderSize-1].split()
#     
#     for i in range(len(HeaderLine)):
#         item = HeaderLine[i]
#         if item == Column_Label:
#             HeaderIndex = i
#     
#     if HeaderIndex == None:
#         print "WARNING: Header not found"
#         print Column_Label
#         print HeaderLine
#         
#     for i in range(HeaderSize,len(LogLines)):
#         LineElements = LogLines[i].split()
#         if LineElements[0] == Row_Label:
#             Parameter = LineElements[HeaderIndex]
# 
#     return Parameter
# 
# # def GetLineIntensityFromTable(FileFolder,FileName,Line_Label):
# #     
# #     TextLines = File2Lines(FileFolder,FileName)
# #     
# #     HeaderSize = 2
# #     
# #     IntensityValue = None
# #     ErrorValue = None
# #     
# #     for i in range(2,len(TextLines)):
# #         LineElements = TextLines[i].split()
# #         
# #         if LineElements[0] == Line_Label:
# #             IntensityValue = LineElements[4]
# #             ErrorValue = LineElements[8]
# #             Wave = LineElements[2]
# #     
# #     return (float(IntensityValue), float(ErrorValue)), float(Wave)
# 
# def RGB_2_MatploblibColor(RGB_Colors):
#     MatColors = [] 
#     for Color in RGB_Colors:
#         NewColor = (Color[0]/255.0, Color[1]/255.0, Color[2]/255.0)
#         MatColors.append(NewColor)
#     return MatColors
#     
# def ColorsSet(Instruction):    
#     
#     if Instruction == 'Mine':
#         Iron = (76,76,76)
#         Silver = (204,204,204)      
#         ONE = (0, 128, 128)
#         TWO = (204, 147, 147)
#         THREE = (141, 203, 226)
#         FOUR = (202, 203, 130)
#         FIVE = (234, 184, 130)
#         SIX = (252, 141, 98)
#         SEVEN = (102, 194, 165)
#         EIGHT = (115, 135, 155)              
# 
#         AxisColor = RGB_2_MatploblibColor([Silver])[0]
#         FondoColor = RGB_2_MatploblibColor([Iron])[0]
#         Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT])  
#     
#     if Instruction == 'Night':
#         Iron = (76,76,76)
#         Silver = (204,204,204)
#         ONE = (255, 217, 47)
#         TWO = (102, 194, 165)
#         THREE = (252, 141, 98)
#         FOUR = (229, 196, 148)
#         FIVE = (230, 138, 195)
#         SIX = (166, 216, 84)
#         SEVEN = (179, 179, 179)
#         EIGHT = (141, 203, 226)        
#         
#         AxisColor = RGB_2_MatploblibColor([Silver])[0]
#         FondoColor = RGB_2_MatploblibColor([Iron])[0]
#         Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
#  
#     if Instruction == 'Night1':
#         Iron = (76,76,76)
#         Silver = (204,204,204)
#         ONE = (0, 128, 128)
#         TWO = (204, 147, 147)
#         THREE = (141, 203, 226)
#         FOUR = (202, 203, 130)
#         FIVE = (234, 184, 130)
#         SIX = (115, 135, 155)
#         SEVEN = (240, 130, 37)
#         EIGHT = (80, 159, 69)        
#         
#         AxisColor = RGB_2_MatploblibColor([Silver])[0]
#         FondoColor = RGB_2_MatploblibColor([Iron])[0]
#         Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT])  
#  
#    
#     if Instruction == 'Day1':
#         Black = (0,0,0)
#         White = (255,255,255)   
#         ONE = (27, 158, 119)
#         TWO = (230, 171, 32)
#         THREE = (231, 41, 138)
#         FOUR = (217, 95, 20)
#         FIVE = (102, 166, 30)
#         SIX = (117, 112, 179)
#         SEVEN = (166, 118, 29)
#         EIGHT = (102, 102, 102)
#     
#         AxisColor = RGB_2_MatploblibColor([Black])[0]
#         FondoColor = RGB_2_MatploblibColor([White])[0]
#         Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
# 
#     if Instruction == 'Day2':
#         Black = (0,0,0)
#         White = (255,255,255)   
#         ONE = (68, 170, 153)
#         TWO = (170, 68, 153)
#         THREE = (221, 204, 119)
#         FOUR = (204, 102, 119)
#         FIVE = (51, 34, 136)
#         SIX = (136, 204, 238)
#         SEVEN = (136, 34, 85)
#         EIGHT = (16, 119, 51)
#     
#         AxisColor = RGB_2_MatploblibColor([Black])[0]
#         FondoColor = RGB_2_MatploblibColor([White])[0]
#         Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
#         
#     if Instruction == 'Day_Elemental':
#         Black = (0,0,0)
#         White = (255,255,255)   
#         ONE = (240, 130, 37)            
#         TWO = (80, 159, 69)
#         THREE = (68, 87, 159)
#         FOUR = (222, 194, 35)
#         FIVE = (51, 34, 136)
#         SIX = (136, 204, 238)
#         SEVEN = (136, 34, 85)
#         EIGHT = (16, 119, 51)
#     
#         AxisColor = RGB_2_MatploblibColor([Black])[0]
#         FondoColor = RGB_2_MatploblibColor([White])[0]
#         Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
#         
#     if Instruction == 'Night_Elemental':
#         Iron = (76,76,76)
#         Silver = (204,204,204)
#         ONE = (240, 130, 37)            
#         TWO = (80, 159, 69)
#         THREE = (68, 87, 159)
#         FOUR = (222, 194, 35)
#         FIVE = (230, 138, 195)
#         SIX = (166, 216, 84)
#         SEVEN = (179, 179, 179)
#         EIGHT = (141, 203, 226)        
#         
#         AxisColor = RGB_2_MatploblibColor([Silver])[0]
#         FondoColor = RGB_2_MatploblibColor([Iron])[0]
#         Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
# 
#     #Color_Vector = [BackgroundColor, ForegroundColor, PaletteColors]      
#     #                        0                1            2
# 
#     Color_Vector = [FondoColor, AxisColor, Jet_Colors]
# 
#     return Color_Vector
# 
# 
# #     
# #     Black = (0,0,0)
# #     White = (255,255,255)
# #     
# #     Iron = (76,76,76)
# #     Silver = (204,204,204)
# #     Steel = (102,102,102)
# #     Magnesium = (179,179,179)  
# #     
# #     Teal = (0,128,128)
# #     mySalmon = (204,147,147)
# #     mySky = (141,203,226)
# #     myAsparagus = (202,203,130)  
# #     myCantaloupe = (234,184,130)  
# #     myMidnight = (115,135,155) 
# #     
# #     Helium = (240,130,37)  
# #     Oxygen = (80,159,69)  
# #     Nitrogen = (68,87,159) 
# #     Sulfur = (222,194,35) 
# #     
# #     Light Color blind
# #     ONE = (255, 217, 47)
# #     TWO = (102, 194, 165)
# #     THREE = (252, 141, 98)
# #     FOUR = (229, 196, 148)
# #     FIVE = (230, 138, 195)
# #     SIX = (166, 216, 84)
# #     SEVEN = (179, 179, 179)
# #     EIGHT = (141, 203, 226)
# 
# #     Dark Color blind 1
# #     ONE = (27, 158, 119)
# #     TWO = (230, 171, 32)
# #     THREE = (231, 41, 138)
# #     FOUR = (217, 95, 20)
# #     FIVE = (102, 166, 30)
# #     SIX = (117, 112, 179)
# #     SEVEN = (166, 118, 29)
# #     EIGHT = (102, 102, 102)
# 
# #     Dark Color blind 2
# #     ONE = (68, 170, 153)
# #     TWO = (170, 68, 153)
# #     THREE = (221, 204, 119)
# #     FOUR = (204, 102, 119)
# #     FIVE = (51, 34, 136)
# #     SIX = (136, 204, 238)
# #     SEVEN = (136, 34, 85)
# #     EIGHT = (16, 119, 51)
# #       
# #     TEAL = (60/255.0,60/255.0,60/255.0)
# #     Pearl_Axis = (206/255.0,204/255.0,247/255.0)
# #     Soft_Red = (165/255.0,42/255.0,42/255.0)
# #     Blue_Light = (0/255.0,127/255.0,174/255.0)
# #     Blue_Dark = (0/255.0,0/255.0,128/255.0)
# #     Yellow_Dark = (140/255.0,110/255.0,0/255.0)
# #     Orange_light = (149/255.0,125/255.0,71/255.0)
# #     Green_sea = (46/255.0,139/255.0,87/255.0)
#     
# def GenerateOneFrameFigure_WithColors(Title,Abcisas,Ordenadas,Figure,Axes1,BackgroundColor,Axes_Color):  
#     Figure.set_facecolor(BackgroundColor)
# 
#     Axes1.spines['bottom'].set_color(Axes_Color)
#     Axes1.spines['top'].set_color(Axes_Color) 
#     Axes1.spines['right'].set_color(Axes_Color)
#     Axes1.spines['left'].set_color(Axes_Color)
#     Axes1.yaxis.label.set_color(Axes_Color)
#     Axes1.xaxis.label.set_color(Axes_Color)    
#     
#     Axes1.tick_params(axis='x', length=7,width=2,labelsize = 18,colors=Axes_Color)
#     Axes1.tick_params(axis='y', length=7,width=2,labelsize = 18,colors=Axes_Color)
#         
#     Axes1.set_xlabel(Abcisas,fontsize=20)
#     Axes1.set_ylabel(Ordenadas,fontsize=20)
#     Axes1.set_title(Title, fontsize=25,color=Axes_Color)
#     
# def ListsDivider(Dividendo,Divisor):
#     Result = []
#     
#     if len(Dividendo) == len(Divisor):
#         for i in range(len(Dividendo)):
#             Value = Dividendo[i] / Divisor[i]
#             print str(Dividendo[i]) + "/" + str(Divisor[i]) + "=" + str(Value)
#             Result.append(Value)
#     else:
#         print "WARNING: Lists do not have the same lengths"
#     
#     return Result
# 
# def Single_Header_Check(FileLines):
#     
#     Pure_Table = True
#     
#     FullEmptyVector = []
#     
#     for Line in FileLines:
#         if not Line.strip():
#             FullEmptyVector.append(0)
#             Pure_Table = False
#         else:
#             FullEmptyVector.append(1)
#     
#     if Pure_Table:
#         return Pure_Table
# 
# def plot_FlagsChecker(Flag, FlagsList):
#     
#     Flag_Exists = False
#     
#     for Item in FlagsList:
#         
#         if Flag in Item:
#     
#             Flag_Exists = True
#     
#     return Flag_Exists
#     
# def plot_FlagsRecover(Flag, FlagsList):
#     
#     FlagValue = "None"
#         
#     for Item in FlagsList:
#         
#         if Flag in Item:
#     
#             FlagValue = int(Item[Item.find("_")+1:len(Item)])
#     
#     if FlagValue == "None":
#         
#         print "WARNING: Expected flag value not found"
#     
#     return FlagValue
# 
# def query_yes_no(question, default="yes"):
#     """Ask a yes/no question via raw_input() and return their answer.
# 
#     "question" is a string that is presented to the user.
#     "default" is the presumed answer if the user just hits <Enter>.
#         It must be "yes" (the default), "no" or None (meaning
#         an answer is required of the user).
# 
#     The "answer" return value is one of "yes" or "no".
#     """
#     valid = {"yes":True,   "y":True,  "ye":True,
#              "no":False,     "n":False}
#     if default == None:
#         prompt = " [y/n] "
#     elif default == "yes":
#         prompt = " [Y/n] "
#     elif default == "no":
#         prompt = " [y/N] "
#     else:
#         raise ValueError("invalid default answer: '%s'" % default)
# 
#     while True:
#         sys.stdout.write(question + prompt)
#         choice = raw_input().lower()
#         if default is not None and choice == '':
#             return valid[default]
#         elif choice in valid:
#             return valid[choice]
#         else:
#             sys.stdout.write("Please respond with 'yes' or 'no' "\
#                              "(or 'y' or 'n').\n")
# 
# def FilledTablesIndexes(Table,TableHeader):
#         FirstValue = 0
#         FinalValue = 0
#         for m in range(TableHeader,len(Table)):
#             Line = Table[m].split()
#             if Line[3] != "None":
#                 if FirstValue ==0 and FinalValue == 0:
#                     FirstValue = m
#                 if m > FinalValue:
#                     FinalValue = m
#                    
#         return FirstValue,FinalValue     
# 
# def LineFinder(myFile,myText):
#     for i in range(len(myFile)):
#         if myText in myFile[i]:
#             return i
#         
# def Sublist(array,inf,sup):
#     subArray = []
#     for i in range(inf,sup):
#         subArray.append(array[i])
#     return subArray
# 
# def TextFileChecker(FileFolder, FileName):
#     FileCheck = False
#         
#     if (mimetypes.guess_type(FileFolder + FileName)[0] == 'text/plain') or (mimetypes.guess_type(FileFolder + FileName)[0] == 'application/x-ns-proxy-autoconfig') or (mimetypes.guess_type(FileFolder + FileName)[0] == None):
#         FileCheck = True
# 
#     return FileCheck
# 
# def FamilyOfItemsInArray(Array):   
#     d = {}
#     for i in Array:
#         if i in d:
#             d[i] = d[i]+1
#         else:
#             d[i] = 1
#     a = d.keys()
#     a.sort()
#     a[::-1]
#     return list(a)    
#        
#     return
# 
# def replace_line(file_name, line_num, text):
#     Input_File = open(file_name, 'r')
#     lines = Input_File.readlines()
#     lines[line_num] = text
#     
#     Output_File = open(file_name, 'w')
#     Output_File.writelines(lines)
#     
#     Input_File.close()
#     Output_File.close()
# 
# def replace_line_v2(FileAddress, Index, NewText):
#     s = open(FileAddress,"r+")
#     Text_Lines = s.readlines()
#     for i in range(len(Text_Lines)):
#         if i == Index: 
#             print 'old', Text_Lines[i]
#             Text_Lines[i].replace(Text_Lines[i], NewText)
#             print 'New', NewText
#     s.close()
#     return
# 
# def Parameter_2_Log(FileFolder, ObjectCode, Parameter, Parameter_Value, Parameter_Error = 'none'):
#     
#     Magnitude = str(Parameter_Value)
#     
#     LogExtension = "_log.txt"
#     
#     Log_File_Address = FileFolder + ObjectCode + LogExtension
#     
#     Obj_Log = File2Table(Log_File_Address, "")
#         
#     Index = "None"
# 
#     for i in range(len(Obj_Log)):
#         Item = Obj_Log[i][0]
#         if Parameter == Item:
#             Index = i
#     
#     if Index == "None":
#         print "WARNING: The parameter cannot be saved since its category cannot be found in object log"
#         print "Parameter name ", Parameter
#         print "Parameter Magnitude", Magnitude
#         print "Object Log"
#         for i in range(len(Obj_Log)):
#             print Obj_Log[i][0]
#            
#     Obj_Log[Index][1] = Magnitude
#          
#     OutputFile = open(Log_File_Address, 'w')
#     
#     
#     ColumnSizes = []
#     for i in range(len(Obj_Log[0])):
#         Biggest_Length = 0
#         for j in range(1,len(Obj_Log)):
#             if len(Obj_Log[j][i]) > Biggest_Length:
#                 Biggest_Length = len(Obj_Log[j][i])
#         ColumnSizes.append(Biggest_Length + 2)
#         
#     NewFormatLine = "%" + str(ColumnSizes[0]) + "s" + "%" + str(ColumnSizes[1]) + "s" + "%" + str(ColumnSizes[2]) + "s" + "%" + str(ColumnSizes[3]) + "s"
#     
#     for z in range(len(Obj_Log)):  
#         NewLine = NewFormatLine % (Obj_Log[z][0],Obj_Log[z][1],Obj_Log[z][2],Obj_Log[z][3])
#         OutputFile.write(NewLine+"\n")
# 
#     OutputFile.close()
#                 
#     return
# 
# def Log_2_Parameter(FileFolder, ObjectCode, ParameterToFind):
#     
#     LogExtension = "_log.txt"
#      
#     Log_File_Address = FileFolder + ObjectCode + LogExtension
#     
#     Obj_Log = File2Table(Log_File_Address, "")
#     
#     Index = "None"
#     
#     for i in range(len(Obj_Log)):
#         Item = Obj_Log[i][0]
#         if ParameterToFind == Item:
#             Index = i
#     
#     if Index == "None":
#         print "WARNING: The parameter cannot be found in log"
#         print ParameterToFind
#         print Obj_Log
#         print FileFolder + ObjectCode + LogExtension
#         
#     Magnitude = Obj_Log[Index][1]
#                 
#     return Magnitude
# 
# def GetParameter_LineLog(CodeName, FileFolder, LineLabel, Parameter_Header, LinesLog_suffix, LinesLogHeader_Address = "/home/vital/Workspace/Dazer_Algorithms/DZT_LineLog_Headers.dz"):
#         
#     #The lines log includes the datatype suffix... is this right????
#     #LinesLog_suffix = _WHT_LinesLog_v3
#         
#     Textfile_Headers                = loadtxt(LinesLogHeader_Address, dtype=str, usecols = [1], unpack = True)
#     Header_Index                    = where(Textfile_Headers==Parameter_Header)[0][0]
#         
#     Labels_Column                   = loadtxt(FileFolder + CodeName + LinesLog_suffix, dtype = str, skiprows = 2, usecols = [0] , unpack = True) 
#     Parameter_Column                = loadtxt(FileFolder + CodeName + LinesLog_suffix, skiprows = 2, usecols = [Header_Index] , unpack = True) 
#     
#     Label_Index                     = where(Labels_Column==LineLabel)[0][0]
#         
#     Parameter                       = Parameter_Column[Label_Index]
#     
#     return Parameter
# 
# def GetParameter_ObjLog(CodeName, FileFolder, Parameter, Assumption = None):
#     
#     ObjLog_Address              = FileFolder + CodeName + '_log.txt'
#     ObjLog_Parameters           = loadtxt(ObjLog_Address, dtype=str, skiprows = 2, usecols = [0], unpack = True)
#     ObjLog_Magnitudes           = loadtxt(ObjLog_Address, dtype=str, skiprows = 2, usecols = [1], unpack = True)
#     
#     Parameter_Index             = where(ObjLog_Parameters == Parameter)[0][0]
#     Parameter_Magnitude         = ObjLog_Magnitudes[Parameter_Index]
#     
#     if Assumption != None:
#         
#         #Temperature needed case for HII region
#         if Assumption == 'Min_Temp':
#             if (Parameter_Magnitude != 'None') and (Parameter_Magnitude != 'None') and (Parameter_Magnitude != '-'):
#                 Converted_Parameter = float(Parameter_Magnitude[0:Parameter_Magnitude.find('+/-')])
#             else:
#                 Converted_Parameter = 10000.0
# 
#         #Density needed case for HII region
#         elif Assumption == 'Min_Den':
#             if (Parameter_Magnitude != 'None') and (Parameter_Magnitude != 'None') and (Parameter_Magnitude != '-'):
#                 Converted_Parameter = float(Parameter_Magnitude[0:Parameter_Magnitude.find('+/-')])
#                 if Converted_Parameter < 75.0:
#                     Converted_Parameter = 50.0
#                 
#             else:
#                 Converted_Parameter = 100.0                
#         
#         elif Assumption == "MatchingSpectra":
#             if Parameter_Magnitude == 'True':
#                 Wave_Index          = where(ObjLog_Parameters == 'WHT_Wmatch')[0][0]
#                 Converted_Parameter = float(ObjLog_Magnitudes[Wave_Index])
#     
#             else:
#                 Blue_lambdaF_ind    = ObjLog_Magnitudes[where(ObjLog_Parameters == 'WHT_Wmax_Blue')[0][0]]
#                 Red_lambdaF_ind     = ObjLog_Magnitudes[where(ObjLog_Parameters == 'WHT_Wmin_Red')[0][0]]
#                 
#                 Converted_Parameter = (float(Blue_lambdaF_ind) + float(Red_lambdaF_ind)) / 2
#         
#         elif Assumption == 'float':
#             if (Parameter_Magnitude != 'None') and (Parameter_Magnitude != 'None') and (Parameter_Magnitude != '-'):
#                 Converted_Parameter = float(Parameter_Magnitude)
#             else:
#                 Converted_Parameter = None
#             
#         
#         return Converted_Parameter
#     
#     else:
#         return Parameter_Magnitude
# 
# def SaveParameter_ObjLog(CodeName, FileFolder, Parameter,  Magnitude, Assumption = None):
#     
#     ObjLog_Address              = FileFolder + CodeName + '_log.txt'
#     ObjLog_Parameters           = loadtxt(ObjLog_Address, dtype=str, skiprows = 2, usecols = [0], unpack = True)
#     ObjLog_Magnitudes           = loadtxt(ObjLog_Address, dtype=str, skiprows = 2, usecols = [1], unpack = True)
#     
#     Parameter_Index             = where(ObjLog_Parameters == Parameter)[0][0]
#     Parameter_Magnitude         = ObjLog_Magnitudes[Parameter_Index]
#     
#     return Parameter_Magnitude


