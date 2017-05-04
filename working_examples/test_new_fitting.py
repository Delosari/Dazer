'''
Created on Oct 19, 2016

@author: vital
'''

import pandas as pd

lines_log_address   = '/home/vital/Dropbox/Astrophysics/Data/WHT_Catalogue_SulfurRegression/Objects/70/70_WHT_LinesLog_v3.txt'

line_df             = pd.read_csv(lines_log_address, skiprows = [0], delim_whitespace = True, header = 0, index_col = 0)

print line_df
'H1_6563A'