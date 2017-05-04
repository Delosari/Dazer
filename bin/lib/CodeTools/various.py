from numpy          import argsort, diff, searchsorted
from mimetypes      import guess_type
from itertools      import izip
from traceback      import format_exc
from collections    import OrderedDict

# from Scientific_Lib.AstroMethods    import Fits2Data, StarlightFileManager
class Error_manager():
    
    def __init__(self):
        
        self.objec_error_dict   = OrderedDict()
        self.Lines_With_Errors  = []
        
    def log_error(self, object_name, verbose = True):
        
        self.objec_error_dict[object_name] = format_exc()
        
        if verbose:
            print '\n-- Error logged'
    
    def log_emLine_error(self, CodeName, Label):
        
        #Case for the first error log for a given object 
        if Label not in self.objec_error_dict:
            LinesList   = [[],[]]
        
        #Case the objects has more lines with an error
        else:
            LinesList   = self.objec_error_dict[CodeName] 

        LinesList[0].append([Label, format_exc()])
                
        self.objec_error_dict[CodeName] = LinesList
        
    def display_errors(self, error_type = 'general', extra_verbose = False):
        
        if error_type == 'general':
            
            if len(self.objec_error_dict) != 0:
                print '\n-These objects treatment produced a script error:\n'
                print self.objec_error_dict.keys(),'\n'
                for i in range(len(self.objec_error_dict)):
                    print '--Object', self.objec_error_dict.keys()[i],
                    print self.objec_error_dict.values()
                    
                #We empty the dictionary as we don not need it anymore
                self.objec_error_dict.clear()
        
        elif error_type == 'emLine measurement':
            
            if len(self.objec_error_dict) != 0:
                print '\n-These objects treatment produced a script error:\n'
                print self.objec_error_dict.keys(),'\n'
                
                for i in range(len(self.objec_error_dict)):
                    print '--Object', self.objec_error_dict.keys()[i] 
                    print self.objec_error_dict.values()[i][0], '\n'
                    
                if extra_verbose == True:
                    for i in range(len(self.objec_error_dict)):
                        print '--Object',   self.objec_error_dict.keys()[i] 
                        print '--Label',    self.objec_error_dict.values()[i][0], '\n'
                        print '--Errors',   self.objec_error_dict.values()[i][1], '\n'  
                        
        return ' '

class vitools():

    def __init__(self):
        '''
        Constructor
        '''
        
    def ufloatDict_nominal(self, ufloat_dict):
        'This gives us a dictionary of nominal values from a dictionary of uncertainties'
        return OrderedDict(izip(ufloat_dict.keys(), map(lambda x: x.nominal_value, ufloat_dict.values())))
    
     
    def ufloatDict_stdev(self, ufloat_dict):
        'This gives us a dictionary of nominal values from a dictionary of uncertainties'
        return OrderedDict(izip(ufloat_dict.keys(), map(lambda x: x.std_dev, ufloat_dict.values())))

    def search_targets_in_array(self, known_array, test_array):
        #This function gives the indeces of the closest values within a sorted array
        
        index_sorted        = argsort(known_array)
        known_array_sorted  = known_array[index_sorted]
        known_array_middles = known_array_sorted[1:] - diff(known_array_sorted.astype('f'))/2
        idx1                = searchsorted(known_array_middles, test_array)
        indices             = index_sorted[idx1]
        
        return indices

def ColorsSet(Instruction):    
     
    if Instruction == 'Mine':
        Iron = (76,76,76)
        Silver = (204,204,204)      
        ONE = (0, 128, 128)
        TWO = (204, 147, 147)
        THREE = (141, 203, 226)
        FOUR = (202, 203, 130)
        FIVE = (234, 184, 130)
        SIX = (252, 141, 98)
        SEVEN = (102, 194, 165)
        EIGHT = (115, 135, 155)              
 
        AxisColor = RGB_2_MatploblibColor([Silver])[0]
        FondoColor = RGB_2_MatploblibColor([Iron])[0]
        Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT])  
     
    if Instruction == 'Night':
        Iron = (76,76,76)
        Silver = (204,204,204)
        ONE = (255, 217, 47)
        TWO = (102, 194, 165)
        THREE = (252, 141, 98)
        FOUR = (229, 196, 148)
        FIVE = (230, 138, 195)
        SIX = (166, 216, 84)
        SEVEN = (179, 179, 179)
        EIGHT = (141, 203, 226)        
         
        AxisColor = RGB_2_MatploblibColor([Silver])[0]
        FondoColor = RGB_2_MatploblibColor([Iron])[0]
        Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
  
    if Instruction == 'Night1':
        Iron = (76,76,76)
        Silver = (204,204,204)
        ONE = (0, 128, 128)
        TWO = (204, 147, 147)
        THREE = (141, 203, 226)
        FOUR = (202, 203, 130)
        FIVE = (234, 184, 130)
        SIX = (115, 135, 155)
        SEVEN = (240, 130, 37)
        EIGHT = (80, 159, 69)        
         
        AxisColor = RGB_2_MatploblibColor([Silver])[0]
        FondoColor = RGB_2_MatploblibColor([Iron])[0]
        Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT])  
  
    
    if Instruction == 'Day1':
        Black = (0,0,0)
        White = (255,255,255)   
        ONE = (27, 158, 119)
        TWO = (230, 171, 32)
        THREE = (231, 41, 138)
        FOUR = (217, 95, 20)
        FIVE = (102, 166, 30)
        SIX = (117, 112, 179)
        SEVEN = (166, 118, 29)
        EIGHT = (102, 102, 102)
     
        AxisColor = RGB_2_MatploblibColor([Black])[0]
        FondoColor = RGB_2_MatploblibColor([White])[0]
        Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
 
    if Instruction == 'Day2':
        Black = (0,0,0)
        White = (255,255,255)   
        ONE = (68, 170, 153)
        TWO = (170, 68, 153)
        THREE = (221, 204, 119)
        FOUR = (204, 102, 119)
        FIVE = (51, 34, 136)
        SIX = (136, 204, 238)
        SEVEN = (136, 34, 85)
        EIGHT = (16, 119, 51)
     
        AxisColor = RGB_2_MatploblibColor([Black])[0]
        FondoColor = RGB_2_MatploblibColor([White])[0]
        Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
         
    if Instruction == 'Day_Elemental':
        Black = (0,0,0)
        White = (255,255,255)   
        ONE = (240, 130, 37)            
        TWO = (80, 159, 69)
        THREE = (68, 87, 159)
        FOUR = (222, 194, 35)
        FIVE = (51, 34, 136)
        SIX = (136, 204, 238)
        SEVEN = (136, 34, 85)
        EIGHT = (16, 119, 51)
     
        AxisColor = RGB_2_MatploblibColor([Black])[0]
        FondoColor = RGB_2_MatploblibColor([White])[0]
        Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
         
    if Instruction == 'Night_Elemental':
        Iron = (76,76,76)
        Silver = (204,204,204)
        ONE = (240, 130, 37)            
        TWO = (80, 159, 69)
        THREE = (68, 87, 159)
        FOUR = (222, 194, 35)
        FIVE = (230, 138, 195)
        SIX = (166, 216, 84)
        SEVEN = (179, 179, 179)
        EIGHT = (141, 203, 226)        
         
        AxisColor = RGB_2_MatploblibColor([Silver])[0]
        FondoColor = RGB_2_MatploblibColor([Iron])[0]
        Jet_Colors = RGB_2_MatploblibColor([ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT]) 
 
    #Color_Vector = [BackgroundColor, ForegroundColor, PaletteColors]      
    #                        0                1            2
 
    Color_Vector = [FondoColor, AxisColor, Jet_Colors]
 
    return Color_Vector

def TextFileChecker(FileFolder, FileName):
    FileCheck = False
         
    if (guess_type(FileFolder + FileName)[0] == 'text/plain') or (guess_type(FileFolder + FileName)[0] == 'application/x-ns-proxy-autoconfig') or (guess_type(FileFolder + FileName)[0] == None):
        FileCheck = True
 
    return FileCheck

def RGB_2_MatploblibColor(RGB_Colors):
    MatColors = [] 
    for Color in RGB_Colors:
        NewColor = (Color[0]/255.0, Color[1]/255.0, Color[2]/255.0)
        MatColors.append(NewColor)
    return MatColors

def FamilyOfItemsInArray(Array):   
    d = {}
    for i in Array:
        if i in d:
            d[i] = d[i]+1
        else:
            d[i] = 1
    a = d.keys()
    a.sort()
    a[::-1]
    return list(a)    
        
    return
