'''
Created on Jun 3, 2016

@author: vital
'''

from bin.dazer_methods import Dazer
from collections import OrderedDict

dz = Dazer()

FileFolder  = '/home/vital/Desktop/mariluz/leecloudys/cloudys_10/'
FileName    = '00500.in.out'

FilesList               = dz.Folder_Explorer('.in.out', FileFolder, CheckComputer=False, Sort_Output='alphabetically')
size_dict, photons_dict = OrderedDict(), OrderedDict()

for file in FilesList:
    
    CodeName, FileName, FileFolder = dz.Analyze_Address(file)
    CodeName = FileName[0:FileName.find('.')]

    File            = open(file)
    File_text       = File.readlines()
    
    j = 0
    found = False
    while found == False:
        
        if '* radius=' in  File_text[j]:
            print CodeName, j, File_text[j]
            found = True
        j += 1
        
    j = 0
    found = False
    while found == False:
        
        if '* Q(H)=' in  File_text[j]:
            print CodeName, j, File_text[j]
            found = True
        j += 1        
    
# for i in range(len(File_text)):
#     if '* radius=' in File_text[i]:
#         print i, File_text[i]

print 'se acabo'
# 
# print 'Line 1494'
# print File_text[1494]
# print type(File_text[1494])
# 
#     print 'bingo'
#     print File_text[1494][File_text[1494].find('diel  1260A'):len(File_text[1494])].split()
    

    
    
