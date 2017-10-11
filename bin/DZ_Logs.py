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