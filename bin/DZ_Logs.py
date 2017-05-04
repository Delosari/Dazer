from bisect                     import bisect
from os.path                    import exists
from numpy                      import where, loadtxt, in1d, searchsorted, asscalar

class LineMesurer_Log():
    
    def __init__(self, Conf_Folder, LinesLogHeader_Name):
        
        self.ColumnHeaderVector, LineLog_Width, LineLog_Format = loadtxt(Conf_Folder + LinesLogHeader_Name, dtype=str, usecols = [1, 2, 3], skiprows=1, unpack=True)
        self.LineLog_WidthtDict = dict(zip(self.ColumnHeaderVector, LineLog_Width))
        self.LineLog_FormatDict = dict(zip(self.ColumnHeaderVector, LineLog_Format))
        
        self.TableHeaderSize    = 2
        self.ColumnWidth        = str(50 + 2)    #len(max(self.ColumnHeaderVector,key=len)) + 2
    
#     def CleanTableMaker(self, TableAddress, RemakeLog):
#         
#         #Check if file exists or if remake instruction is enforced
#         if (exists(TableAddress) == False) or (RemakeLog == True):
#             
#             #Create the table
#             myTextFile = open(TableAddress, "w")   
#             
#             #Write the context text
#             IntroductionLines = ["-File address: " + TableAddress]
#             myTextFile.write(IntroductionLines[0]+'\n')
#         
#             #Generate the headers rows with the width from configuration            
#             HeadersRow = "".join(('{'+':>' + self.LineLog_WidthtDict[self.ColumnHeaderVector[i]] + 's' + '}').format(self.ColumnHeaderVector[i]) for i in range(len(self.ColumnHeaderVector)))
#                                                 
#             #Write the headers row and close the file
#             myTextFile.write(HeadersRow + "\n")
#             myTextFile.close() 
#     
#     def GetDataInTable(self, ColumnInput, ColumnOutput, RowOutput, TableAddress, TableHeaderSize, dataype = 'str'):
#         
#         #This method could be expanded to load data with the same type.. maybe using a dictionary
#         #Add WARNING IF ROW IS NOT IN FILE
#         
#         #Confirm you are asking the right header
#         if in1d(ColumnOutput, self.ColumnHeaderVector)[0]:
#             
#             #Find Column Indeces
#             ColumnInput = where(in1d(self.ColumnHeaderVector, [ColumnInput]))[0][0]
#             ColumnOuput = where(in1d(self.ColumnHeaderVector, [ColumnOutput]))[0][0]
# 
#             #Load Columns 
#             Columns = loadtxt(TableAddress, usecols = [ColumnInput, ColumnOuput], skiprows = TableHeaderSize, dtype=str, unpack=True, ndmin=1)
#                         
#             if len(Columns[0]) > 0:
#                 InputColumn = Columns[0]
#                 OutputColumn = Columns[1]
# 
#                 #Find target line row
#                 target_row = where(InputColumn.astype(type(RowOutput)) == RowOutput)
#                                 
#                 #Get parameter
#                 Parameter = OutputColumn[target_row]
#                 
#                 #If only one entry meets condition:
#                 if len(Parameter) == 1:
#                     if dataype != 'str':
#                         Parameter = Parameter.astype(float)
#                     return Parameter[0]
#                     
#                 #If several entries in table
#                 else:
#                     if dataype != 'str':
#                         Parameter = Parameter.astype(float)
#                     return Parameter
#             
#         else:
#             print 'WARNING: Parameter', ColumnInput, 'could not be found for ', RowOutput
#             print 'in', TableAddress
#             Parameter = None
#             
#         return
#     
#     def InsertNewLine(self, Elemental_Line_String, Current_Em_El, Current_line_label, TableAddress):
# 
#         #Load wavelength of observed lines
#         LabelColumn         = loadtxt(TableAddress, usecols = [0], skiprows = self.TableHeaderSize, ndmin = 1, dtype = str)
#         WavelengthColumn    = loadtxt(TableAddress, usecols = [2], skiprows = self.TableHeaderSize, ndmin = 1, dtype = float)
#         boolean_vector      = in1d(Current_line_label, LabelColumn)
#                
#         #Check if the current line has not been observed.
#         if boolean_vector[0] == False:
#         
#             #Load the table
#             TableFile   = open(TableAddress,"r")
#             TableLines  = TableFile.readlines()
#             TableFile.close()
#             
#             #Check we are not adding the first line
#             if len(range(self.TableHeaderSize,len(TableLines))) > 0:
#                         
#                 #Find best index to insert the table
#                 LineLocation = bisect(WavelengthColumn, Current_Em_El)
#                 
#                 #Insert new row
#                 TableLines.insert(LineLocation + self.TableHeaderSize ,Elemental_Line_String+"\n")
#                
#                 #Save the table
#                 OutputFile = open(TableAddress, 'w')
#                 OutputFile.writelines(TableLines)
#                 OutputFile.close()
#             
#             else:
#                 #Case we are writting the first line
#                 TableLines.append(Elemental_Line_String + "\n")
#                 OutputFile = open(TableAddress, 'w')
#                 OutputFile.writelines(TableLines)
#                 OutputFile.close()
#                 
#         return
#             
#     def Replace_Row(self, Index, RowName, DataDict, TableAddress):
# 
#         #Load the table
#         TableFile   = open(TableAddress,"r")
#         TableLines  = TableFile.readlines()
#         TableFile.close()
# 
#         #Find row Index
# 
#         WavelengthColumn = loadtxt(TableAddress, usecols = [2], skiprows = self.TableHeaderSize, dtype = float)
#         
#         print 'por que no casa'
#         print WavelengthColumn
#         print RowName
#         print where(WavelengthColumn == RowName)
#         
#         RowIndex = where(WavelengthColumn == RowName)[0][0] + self.TableHeaderSize
#         
#         #Generate the line record        
#         NewRow = ""
#         for i in range(len(self.ColumnHeaderVector)):
#             
#             Parameter   =   DataDict[self.ColumnHeaderVector[i]][Index]
#             
#             #Generate the formated entry
#             if Parameter != None:
#                 style    =   '{'+':>' + self.LineLog_WidthtDict[self.ColumnHeaderVector[i]] + self.LineLog_FormatDict[self.ColumnHeaderVector[i]] + '}'
#                 
#             #In case it is not calculated change it to 'None'
#             else:
#                 Parameter = 'None'
#                 style =  '{'+':>' + self.LineLog_WidthtDict[self.ColumnHeaderVector[i]] + 's' + '}'
#                 
#             Entry       =   style.format(Parameter)
# 
#             NewRow      =   NewRow + Entry
#                 
#         #Replace the row   
#         TableLines[RowIndex] = NewRow + "\n"
#         
#         #Load the column
#         out = open(TableAddress, 'w')
#         out.writelines(TableLines)
#         out.close()
# 
#         return
# 
# 
#     def Replace_Row_byLabel(self, Index, Label, DataDict, TableAddress):
# 
#         #Load the table
#         TableFile   = open(TableAddress,"r")
#         TableLines  = TableFile.readlines()
#         TableFile.close()
# 
#         #Find row Index
# 
#         LabelColumn = loadtxt(TableAddress, usecols = [0], skiprows = self.TableHeaderSize, dtype = str)
#                 
#         RowIndex = where(LabelColumn == Label)[0][0] + self.TableHeaderSize
#         
#         #Generate the line record        
#         NewRow = ""
#         for i in range(len(self.ColumnHeaderVector)):
#             
#             Parameter   =   DataDict[self.ColumnHeaderVector[i]][Index]
#             
#             #Generate the formated entry
#             if Parameter != None:
#                 style    =   '{'+':>' + self.LineLog_WidthtDict[self.ColumnHeaderVector[i]] + self.LineLog_FormatDict[self.ColumnHeaderVector[i]] + '}'
#                 
#             #In case it is not calculated change it to 'None'
#             else:
#                 Parameter = 'None'
#                 style =  '{'+':>' + self.LineLog_WidthtDict[self.ColumnHeaderVector[i]] + 's' + '}'
#                 
#             Entry       =   style.format(Parameter)
# 
#             NewRow      =   NewRow + Entry
#                 
#         #Replace the row   
#         TableLines[RowIndex] = NewRow + "\n"
#         
#         #Load the column
#         out = open(TableAddress, 'w')
#         out.writelines(TableLines)
#         out.close()
# 
#         return
# 
#     def CreateEmissionLineReccord(self, Current_Em_El, Current_Ion,  Current_Label):
#         
#         #Generate empty list
#         Empty_List = ['None'] * len(self.ColumnHeaderVector)
#         
#         #Insert the identifiying entries on the empty list    
#         Empty_List[0], Empty_List[1], Empty_List[2] = Current_Label, Current_Ion, str(Current_Em_El)
#                               
#         #Generate the headers rows with the width from configuration            
#         Elemental_Line_String = "".join(('{'+':>' + self.LineLog_WidthtDict[self.ColumnHeaderVector[i]] + 's' + '}').format(Empty_List[i]) for i in range(len(self.ColumnHeaderVector)))
#            
#         return Elemental_Line_String
#     
#     def DeleteLine(self,  Current_Em_Wave, TableAddress, TableHeaderSize, Deblend_Check, List_BlendedLines, Current_BlendedGroup):
#         #WARNING CORRECT THIS BLENDED LINES
# 
#         #Open text file and import the lines
#         TableFile = open(TableAddress,"r")
#         TableLines = TableFile.readlines()
#         TableFile.close()
#         
#         #Index for the lines to delete (CHANGE THIS TO THE ENUMERATE SCHEME
#         LocationLineToDelete = None        
#         
#         #FOR THE TIME BEING THIS ONLY WORKS IF WE ARE CURRENTLY MEASURING A LINE
#         if Deblend_Check == False:
#             Measured_Wavelengths = loadtxt(TableAddress, dtype=float, skiprows= TableHeaderSize, usecols = [2])
#             LocationLineToDelete = int(where(Measured_Wavelengths==Current_Em_Wave)[0]) + TableHeaderSize
#             del TableLines[LocationLineToDelete]
#                     
#         else:
#             print 'Current_BlendedGroup', Current_BlendedGroup
#             grouped_lines_wavelengths   = List_BlendedLines[1][Current_BlendedGroup]
#             remove_list = []
# 
#             for Current_Em_Wave in grouped_lines_wavelengths:
#                 for j in range(TableHeaderSize,len(TableLines)):
#                     Wave = float(TableLines[j].split()[2])
#                     if Wave == Current_Em_Wave:
#                         remove_list.append(j)
#                         
#             TableLines = [v for i, v in enumerate(TableLines) if i not in remove_list]
#                                              
#         out = open(TableAddress, 'w')
#         out.writelines(TableLines)
#         out.close()
#         
#         return 
# 
#     def Get_DataLog_Parameter(self, Row_Label, Column_Label, TableAddress, TableHeaderSize):
#         
#         TableFile = open(TableAddress,"r")
#         LogLines = TableFile.readlines()
#         TableFile.close()
#             
#         Parameter = None
#         HeaderIndex = None
#         
#         HeaderLine = LogLines[TableHeaderSize-1].split()
#         
#         for i in range(len(HeaderLine)):
#             item = HeaderLine[i]
#             if item == Column_Label:
#                 HeaderIndex = i
#         
#         if HeaderIndex == None:
#             print "WARNING: Header not found"
#             print Column_Label
#             print HeaderLine
#             
#         for i in range(TableHeaderSize, len(LogLines)):
#             LineElements = LogLines[i].split()
#             if LineElements[0] == Row_Label:
#                 Parameter = LineElements[HeaderIndex]
#     
#         return Parameter
# 
#     def RangesR(self, Selections, RowName, TableAddress):                                                       	              
#         
#         #Same wavelengths units everywhere???
#         #Use this trick to check we are not measuring a line (and not load new data):       
#         if len(Selections) == 0:
#             
#             #Load wavelength of observed lines #THERE IS A METHOD TO DISABLE THIS WARNING ON AN EMPTY FILE
#             WavelengthColumn = loadtxt(TableAddress, usecols = [2], skiprows = self.TableHeaderSize, ndmin = 1)            
#             
#            
#             #Check if the current line has not been observed.
#             row_check = in1d(RowName, WavelengthColumn)
#                        
#             if row_check[0] == True:
#                 
#                 #Get index of observed wavelength
#                 index_row           = where(WavelengthColumn == RowName)[0][0]
#                                 
#                 #Get from indeces of the columns with the wave points               
#                 indeces_waveC       = where(in1d(self.ColumnHeaderVector, ['Wave1', 'Wave2', 'Wave3', 'Wave4', 'Wave5', 'Wave6']))[0]
#                 
#                 #Load the columns
#                 WaveIndexesColumns  = loadtxt(TableAddress, usecols = indeces_waveC, skiprows = self.TableHeaderSize, ndmin = 2)
#                                 
#                 for i in range(6):
#                     wave_value = WaveIndexesColumns[index_row][i]
#                     if wave_value != 'None':
#                         Selections.append(wave_value)
#                         
#         return Selections