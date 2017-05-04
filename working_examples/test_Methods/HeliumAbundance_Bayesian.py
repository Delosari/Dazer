#!/usr/bin/python

from numpy import linspace, zeros, exp

from Astro_Libraries.he_bayesian         import He_Inference_Abundance
from CodeTools.PlottingManager                      import myPickle


#Generate dazer object
pv = myPickle()

#Define plot frame and colors
pv.FigFormat_One(ColorConf='Night1')

He_abun_tools = He_Inference_Abundance()


print 'the matrix'
for matrix in He_abun_tools.H_CollCoeff_Matrix:
    print He_abun_tools.H_CollCoeff_Matrix[0]

#Display tables
# LineCoefficients = He_abun_tools.He_CollCoeff_Matrix[5]
# for i in range(len(LineCoefficients[0])):
#     print LineCoefficients[0][i], LineCoefficients[1][i], LineCoefficients[2][i]
#     
# Line_Optical    = He_abun_tools.He_OpticalDepth_Matrix
# LineDepth_Coeff = Line_Optical[4]
# for i in range(len(LineDepth_Coeff)):
#     print LineDepth_Coeff[i]


#Plot Neutral Heliuum Collisional Correction
# Te_range    = linspace(10000, 25000, 50)
# y_values    = zeros(len(Te_range))
#  
# print 'The temperature range is', Te_range
# for i in range(len(Te_range)):
#     T_4         = Te_range[i] / 10000  
#     y_values[i] = 6.7937 * exp(-0.1116 / T_4) * (T_4**-0.1116)
#      
# label = 'Halpha'
# pv.DataPloter_One(X = Te_range, Y = y_values, LineLabel = label,LineColor = pv.Color_Vector[2][1])

# #Plot Neutral Hydrogen Collisional Correction
xi          = 1e-4
Te_range    = linspace(10000, 25000, 50)
LinesList   = He_abun_tools.Hydrogen_Lines
Coeff_Table = zeros(shape=(len(LinesList), len(Te_range)))
   
for i in range(len(LinesList)):
    for j in range(len(Te_range)):
        Temp                = Te_range[j]
        i_K_alpha_Ratio     = He_abun_tools.H_Kalpha_Ratio(i, Temp)
        Coeff_Table[i,j]    = i_K_alpha_Ratio
           
for i in range(len(LinesList)):
    label = LinesList[i]
    pv.DataPloter_One(X = Te_range, Y = Coeff_Table[i], LineLabel = label,LineColor = pv.Color_Vector[2][i])
   
title       = 'Neutral Hydrogen Collisonal Correction'
xtitle      = r'Te $(K)$'
ytitle      = r'$C/R(\lambda)$'
 
pv.Labels_Legends_One(Plot_Title = title, Plot_xlabel = xtitle, Plot_ylabel = ytitle, LegendLocation = 'best')
pv.Axis1.set_xlim(10000,25000)
pv.Axis1.set_ylim(-0.01,0.08)
pv.DisplayFigure()
 
#Plot the one Neutral Helium Collisional fraction
# xi          = 1e-4
# Te_range    = linspace(0, 25000, 100)
# den         = 1000
#   
# LinesList   = He_abun_tools.Helium_Lines[3]
# CR_Coef     = zeros(len(Te_range))
#   
# for j in range(len(Te_range)):
#     Temp                = Te_range[j]
#     CR_Ratio            = He_abun_tools.He_CR_ratio(3, Temp, den)
#     CR_Coef[j]          = CR_Ratio
#           
# label = LinesList 
# pv.DataPloter_One(X = Te_range, Y =CR_Coef, LineLabel = label,LineColor = pv.Color_Vector[2][3])
#   
# title       = 'Neutral Helium Collisonal Correction'
# xtitle      = r'Te $(K)$'
# ytitle      = r'$C/R(\lambda)$'
#   
# pv.Labels_Legends_One(Plot_Title = title, Plot_xlabel = xtitle, Plot_ylabel = ytitle, LegendLocation = 'best')
# pv.Axis1.set_xlim(0,27500)
# pv.Axis1.set_ylim(0,0.7)
# pv.DisplayFigure()

#Plot Neutral Helium Collisional fracion
# pv.DisplayFigure()
#  
# xi          = 1e-4
# Te_range    = linspace(10000, 25000, 50)
# den         = 100
#    
# LinesList   = He_abun_tools.Helium_Lines
# Coeff_Table = zeros(shape=(len(LinesList), len(Te_range)))
#    
# for i in range(len(LinesList)):
#     for j in range(len(Te_range)):
#         Temp                = Te_range[j]
#         CR_Ratio            = He_abun_tools.He_CR_ratio(i, Temp, den)
#         Coeff_Table[i,j]    = CR_Ratio
#            
# for i in range(len(LinesList)):
#     label = LinesList[i]
#     pv.DataPloter_One(X = Te_range, Y = Coeff_Table[i], LineLabel = label,LineColor = pv.Color_Vector[2][i])
#    
# title       = 'Neutral Helium Collisonal Correction'
# xtitle      = r'Te $(K)$'
# ytitle      = r'$C/R(\lambda)$'
#    
# pv.Labels_Legends_One(Plot_Title = title, Plot_xlabel = xtitle, Plot_ylabel = ytitle, LegendLocation = 'best')
#   
# pv.DisplayFigure()

#Plot Helium lines optical depth
 
# Tau_range   = linspace(0, 100, 50)
# Temp        = 10000
# den         = 100
#  
# LinesList   = He_abun_tools.Helium_Lines
# Coeff_Table = zeros(shape=(len(LinesList), len(Tau_range)))
#  
# for i in range(len(He_abun_tools.He_OpticalDepth_Matrix)):
#     print He_abun_tools.He_OpticalDepth_Matrix[i]
#      
# for i in range(len(LinesList)):
#     for j in range(len(Tau_range)):
#         tau                = Tau_range[j]
#         f_tau               = He_abun_tools.He_OpticalDepth(tau, Temp, den, i)
#         Coeff_Table[i,j]    = f_tau
#              
# for i in range(len(LinesList)):
#     label = LinesList[i]
#     pv.DataPloter_One(X = Tau_range, Y = Coeff_Table[i], LineLabel = label,LineColor = pv.Color_Vector[2][i])
#      
# title       = 'Optical depth correction factor for He lines'
# xtitle      = r'$(\tau_{\lambda})$'
# ytitle      = r'$(f_{\tau})$'
#      
# pv.Labels_Legends_One(Plot_Title = title, Plot_xlabel = xtitle, Plot_ylabel = ytitle, LegendLocation = 'best')
#     
# pv.DisplayFigure()

# Tau_range   = linspace(0, 100, 50)
# Temp        = 1
# den         = 0.100
#  
# y_f_tau     = zeros(len(Tau_range))
#     
# for j in range(len(Tau_range)):
#     y_f_tau[j] = 1 + (Tau_range[j]/2) * (0.359 + (-3.46e-2 - 1.84e-4*den + 3.039e-7*den**2 )*Temp)
# #                 1 + (P[1]/2)         * [0.359 + (-3.46e-2 - 1.84e-4*P[0]+ 3.039e-7*P[0]^2) *t]
# pv.DataPloter_One(X = Tau_range, Y = y_f_tau, LineLabel = '7065',LineColor = pv.Color_Vector[2][2])
#      
# title       = 'Optical depth correction factor for He lines'
# xtitle      = r'$(\tau_{\lambda})$'
# ytitle      = r'$(f_{\tau})$'
#      
# pv.Labels_Legends_One(Plot_Title = title, Plot_xlabel = xtitle, Plot_ylabel = ytitle, LegendLocation = 'best')
#     
# pv.DisplayFigure()





