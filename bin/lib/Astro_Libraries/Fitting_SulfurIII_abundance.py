'''
Created on Nov 17, 2015

@author: vital
'''

from Plotting_Libraries.dazer_plotter           import Plot_Conf
import numpy as np
import pyneb as pn


#Declare script classes
dz = Plot_Conf()

#Define figure format
dz.FigConf(FigWidth = 16, FigHeight = 9)

Te_range = np.linspace(5000, 25000, 1000)
Ne_range = np.linspace(100, 1000, 10)