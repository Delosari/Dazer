'''
Created on May 23, 2016

@author: vital
'''
from lmfit          import Parameters, minimize, report_fit
from lmfit.models   import LinearModel
from numpy          import array, loadtxt, searchsorted, where, logical_and, hstack, linspace, min as np_min, max as np_max, fliplr, trapz
from scipy.integrate import simps

from CodeTools.PlottingManager              import myPickle
from Plotting_Libraries.dazer_plotter       import Plot_Conf


def determine_linear_continua(x, y, indices_list):
    
    lineal_mod = LinearModel(prefix='lineal_')

#     x_combine = x[indices_list[0]]
#     y_combine = y[indices_list[0]]
    x_combine = array([])
    y_combine = array([])    

    for i in range(len(indices_list)):

        x_combine = hstack([x_combine, x[indices_list[i]]])
        y_combine = hstack([y_combine, y[indices_list[i]]])
        
    Lineal_parameters = lineal_mod.guess(y_combine, x=x_combine)

    x_lineal    = linspace(np_min(x_combine), np_max(x_combine), 100)
    y_lineal    = Lineal_parameters['lineal_slope'].value * x_lineal + Lineal_parameters['lineal_intercept'].value

    return x_lineal, y_lineal, Lineal_parameters

#Import the classes
pv      = myPickle()
dz      = Plot_Conf() 

#Define figure format
dz.FigConf()

#Declare the data
FilesFolder         = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Total_Grid_with_continua/' 
File_con            = 'S_Ar_Test_age5.0_zStar-2.4_zGas0.1_u-1.5.con'
File_trans_punch    = 'S_Ar_Test_age5.0_zStar-2.4_zGas0.31_u-1.5.transContinuum'
File_inci_punch    = 'S_Ar_Test_age5.0_zStar-2.4_zGas0.31_u-1.5.inciContinuum'

#Importing the .con columns
lambda_angs, incident, trans, diff_Out, net_trans, reflc, total = loadtxt(FilesFolder + File_con, skiprows = 0, usecols = (0, 1, 2, 3, 4, 5, 6), unpack = True)
lambda_angs, incident, trans, diff_Out, net_trans, reflc, total = lambda_angs[::-1], incident[::-1], trans[::-1], diff_Out[::-1], net_trans[::-1], reflc[::-1], total[::-1]

#Plotting .con continua
dz.data_plot(lambda_angs, incident, label='Incident continuum',linestyle='-')
# dz.data_plot(lambda_angs, trans, label='Transmitted continuum')
dz.data_plot(lambda_angs, net_trans, label='Net transmitted continuum', linestyle='-')
# dz.data_plot(lambda_angs, reflc, label='Reflected continuum', linestyle=':')
# dz.data_plot(lambda_angs, total, label='Total continuum', linestyle='--')

#Importing the punch emission
# lambda_angs, Flux_trans = loadtxt(FilesFolder + File_trans_punch, skiprows = 9, usecols = (0, 1), unpack = True)
# lambda_angsI, Flux_inci = loadtxt(FilesFolder + File_inci_punch, skiprows = 9, usecols = (0, 1), unpack = True)
# dz.data_plot(lambda_angs, Flux_trans, label='Transmitted punch', linestyle='--')
# dz.data_plot(lambda_angsI, Flux_inci, label='Incident punch', linestyle='--')

#Declaring Halpha, Hbeta regions and continua
ind_Hbeta   = [4480, 4600, 4830, 4888, 5080, 5160]
Blue_Hbeta, line_Hbeta, Red_Hbeta = where(logical_and(lambda_angs>=ind_Hbeta[0], lambda_angs<=ind_Hbeta[1])),where(logical_and(lambda_angs>=ind_Hbeta[2], lambda_angs<=ind_Hbeta[3])), where(logical_and(lambda_angs>=ind_Hbeta[4], lambda_angs<=ind_Hbeta[5]))
ind_Halpha  = [6390, 6528, 6527, 6634, 6788, 7020]
Blue_Halpha, line_Halpha, Red_Halpha = where(logical_and(lambda_angs>=ind_Halpha[0], lambda_angs<=ind_Halpha[1])),where(logical_and(lambda_angs>=ind_Halpha[2], lambda_angs<=ind_Halpha[3])), where(logical_and(lambda_angs>=ind_Halpha[4], lambda_angs<=ind_Halpha[5]))

#Calculating each line linear continuum
hbeta_con_wave, hbeta_con_flux, hbeta_linealParameters = determine_linear_continua(x = lambda_angs, y = net_trans, indices_list = [Blue_Hbeta, Red_Hbeta])
halpha_con_wave, halpha_con_flux, halpha_linealParameters = determine_linear_continua(x = lambda_angs, y = net_trans, indices_list = [Blue_Halpha, Red_Halpha])
hbeta_inci_wave, hbeta_inci_flux, hbeta_inci_linealParameters = determine_linear_continua(x = lambda_angs, y = incident, indices_list = [Blue_Hbeta, Red_Hbeta])


#Calculate lines flux
Hbeta_TL = simps(y = net_trans[line_Hbeta], x = lambda_angs[line_Hbeta])
Hbeta_Cont =  simps(y = hbeta_linealParameters['lineal_slope'].value * lambda_angs[line_Hbeta] + hbeta_linealParameters['lineal_intercept'].value, x = lambda_angs[line_Hbeta])


Halpha_Tl = simps(y = net_trans[line_Halpha], x = lambda_angs[line_Halpha])
Halpha_Cont = simps(y = halpha_linealParameters['lineal_slope'].value * lambda_angs[line_Halpha] + halpha_linealParameters['lineal_intercept'].value, x = lambda_angs[line_Halpha])

Halpha_flux = Halpha_Tl - Halpha_Cont
Hbeta_Flux = Hbeta_TL - Hbeta_Cont

print 'Halpha, line/continuum', Halpha_Tl,Halpha_Cont
print 'Hbeta, line/continuum', Hbeta_TL, Hbeta_Cont
print 'Hbeta/Hbeta + cont', Halpha_Tl / Hbeta_TL
print 'Halpha/Hbeta', Halpha_flux / Hbeta_Flux

Halpha_conlevel =  halpha_linealParameters['lineal_slope'].value * 6563 + halpha_linealParameters['lineal_intercept'].value
Hbeta_conlevel = hbeta_linealParameters['lineal_slope'].value * 4860 + hbeta_linealParameters['lineal_intercept'].value

print 'Halpha eqw', Halpha_flux / Halpha_conlevel
print 'Hbeta eqw',  Hbeta_Flux / Hbeta_conlevel
print 'Hbeta eqw inci',  Hbeta_Flux / (hbeta_inci_linealParameters['lineal_slope'].value * 4860 + hbeta_inci_linealParameters['lineal_intercept'].value)
print 'Halpha eqw inci',  Halpha_flux / (hbeta_inci_linealParameters['lineal_slope'].value * 4860 + hbeta_inci_linealParameters['lineal_intercept'].value)

#Plotting Eqw calculation
dz.Axis.fill_between(lambda_angs[Blue_Hbeta],   0, net_trans[Blue_Hbeta], facecolor = 'blue', alpha = 0.5)
dz.Axis.fill_between(lambda_angs[line_Hbeta],   0, net_trans[line_Hbeta], facecolor = 'green', alpha = 0.5)
dz.Axis.fill_between(lambda_angs[Red_Hbeta],   0, net_trans[Red_Hbeta], facecolor = 'red', alpha = 0.5)
dz.Axis.fill_between(lambda_angs[Blue_Halpha],   0, net_trans[Blue_Halpha], facecolor = 'blue', alpha = 0.5)
dz.Axis.fill_between(lambda_angs[line_Halpha],   0, net_trans[line_Halpha], facecolor = 'green', alpha = 0.5)
dz.Axis.fill_between(lambda_angs[Red_Halpha],    0, net_trans[Red_Halpha], facecolor = 'red', alpha = 0.5)
dz.data_plot(6563, halpha_linealParameters['lineal_slope'].value * 6563 + halpha_linealParameters['lineal_intercept'].value, markerstyle='o')
dz.data_plot(4860, hbeta_linealParameters['lineal_slope'].value * 4860 + hbeta_linealParameters['lineal_intercept'].value, markerstyle='o')
dz.data_plot(4860, hbeta_inci_linealParameters['lineal_slope'].value * 4860 + hbeta_inci_linealParameters['lineal_intercept'].value, markerstyle='o')

dz.data_plot(hbeta_con_wave, hbeta_con_flux, label='Lineal fitting', linestyle='-')
dz.data_plot(halpha_con_wave, halpha_con_flux, label='Lineal fitting', linestyle='-')

xtitle  = r'lambda $(\AA)$'
ytitle  = r'Flux'
title   = 'Transmitted spectrum: $log(U) = -1.5$, $log(age) = 5.0$, $z_{gas}=0.00132$, $z_{star} = -2.4$'
dz.FigWording(xtitle, ytitle, title, axis_Size=30, title_Size=30, legend_size=25, legend_loc='best')
# dz.Axis.set_yscale('log')
# dz.Axis.set_xscale('log')
dz.Axis.set_xlim(3000.0, 10000.0)

dz.display_fig()
