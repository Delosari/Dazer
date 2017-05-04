'''
Created on Aug 25, 2015

@author: vital
'''
from numpy import array

from Astro_Libraries.Reddening_Corrections              import ReddeningLaws    
from CodeTools.PlottingManager                          import myPickle


#Generate dazer object
pv = myPickle()
rd = ReddeningLaws()

#Define plot frame and colors
pv.FigFormat_One(ColorConf='Night1')

#Import from cardelli
x_true, X_Angs, Xx, f_lambda = rd.Epm_ReddeningPoints()

# #Cardelli
# Hbeta_Angs              = array([4862.683])
# Xx_Hbeta_Cardelli       = rd.Reddening_f(Hbeta_Angs, R_v=3.1, curve_parametrization = 'Cardelli1989')
f_lambda_cardelli             = rd.Reddening_f(X_Angs, R_v=3.1, curve_parametrization = 'Cardelli1989')
# f_lambda_cardelli       = Xx_Cardelli / Xx_Hbeta_Cardelli - 1
pv.DataPloter_One(X = X_Angs, Y = f_lambda_cardelli, LineLabel = 'Cardelli - flambda',  MarkerColor=pv.Color_Vector[2][3], LineStyle=None)


Xx_Millers              = rd.Reddening_f(X_Angs, R_v=3.1, curve_parametrization = 'Miller&Mathews1972')


pv.DataPloter_One(X = X_Angs, Y = f_lambda, LineLabel = 'EPM - flambda',          MarkerColor=pv.Color_Vector[2][1], LineStyle=None)
pv.DataPloter_One(X = X_Angs, Y = Xx_Millers, LineLabel = 'Miller - flambda',     MarkerColor=pv.Color_Vector[2][2], LineStyle=None)

pv.Labels_Legends_One()

pv.DisplayFigure()