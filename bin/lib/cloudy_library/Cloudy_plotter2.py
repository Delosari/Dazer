from collections                        import OrderedDict
from lmfit                              import Parameters, minimize, report_fit
from lmfit.models                       import LinearModel
from numpy                              import log10 as nplog10, zeros, min, max, linspace, array, concatenate, isfinite, greater, hstack, min as np_min, max as np_max
from pandas                             import Series
from Math_Libraries.bces_script         import bces
from Plotting_Libraries.dazer_plotter   import Plot_Conf
from cloudy_library.cloudy_methods      import Cloudy_Tools

def Ar_S_model(Line_dict, threshold = 0.0):
    
#     Hbeta           = Line_dict['4861.36A']
#     x_axis = Line_dict['10.51m'] / Hbeta
#     y_axis = Line_dict['4740.12A'] / Hbeta
#     return x_axis, y_axis

#     S2_S3_ratio     = (Line_dict['6716A'] + Line_dict['6731A']) / (Line_dict['9068.62A'] + Line_dict['9532A'])
#      
#     x_axis = nplog10(S2_S3_ratio)
#           
#     return x_axis, None

    S2_S3_ratio     = (Line_dict['9068.62A'] + Line_dict['9532A']) / Line_dict['10.51m']
    Ar2_Ar3_ratio   = (Line_dict['7135A']) / (Line_dict['4740.12A'])

#     S2_S3_ratio     = (Line_dict['6716A'] + Line_dict['6731A']) / (Line_dict['9068.62A'] + Line_dict['9532A'])
#     Ar2_Ar3_ratio   = (Line_dict['3727A']) / (Line_dict['4959A'] + Line_dict['5007A'] )
      
    x_axis = nplog10(S2_S3_ratio)
    y_axis = nplog10(Ar2_Ar3_ratio)
      
    if (x_axis < threshold):
      
        print 'values',x_axis, y_axis
           
        return x_axis, y_axis
      
    else:
          
        return None, None

dz              = Plot_Conf()
ct              = Cloudy_Tools()

#Declare the number of colors
dz.FigConf(n_colors = 7)

#Define script name and location
ScriptFolder    = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Total_Grid_SAL2/'

#Orion nebula has a metallicity of 2/3 solar
Grid_Values = OrderedDict()

Grid_Values['age']      = ['5.0', '5.25', '5.5', '5.75', '6.0', '6.25', '6.5', '6.75', '7.0', '7.5'] 
Grid_Values['zStars']   = ['-2.1'] 
Grid_Values['zGas']     = ['0.05', '0.2', '0.4', '1'] 
Grid_Values['u']        = ['-3.5', '-3.0', '-2.5', '-2.0', '-1.5'] 

# Grid_Values['age']      = ['5.0', '5.25', '5.5', '5.75', '6.0', '6.25', '6.5', '6.75'] 
# Grid_Values['zStars']   = ['-2.1'] 
# Grid_Values['zGas']     = ['0.05', '0.2', '0.4', 1] 
# Grid_Values['u']        = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5'] 

#Dictionary of dictionaries
Grid_frame =  ({k : Series(v) for k, v in Grid_Values.iteritems()})

x_linealFitting = array([])
y_linealFitting = array([])

#Fore the each grid point run a cloudy simulation
Model_dict = OrderedDict()
for age in Grid_Values['age']:
    for zStar in Grid_Values['zStars']:
        for zGas in Grid_Values['zGas']:                                        
            for ionization_parameter in Grid_Values['u']: 
                
                Model_dict['Name']      = 'S_Ar_Test_' + 'age'+ age + '_zStar' + zStar + '_zGas' + zGas + '_u' + ionization_parameter
                Model_dict['age']       = age
                Model_dict['zStars']    = zStar
                Model_dict['zGas']      = zGas
                Model_dict['u']         = ionization_parameter
    
                print 'Current file', Model_dict['Name'] + '.in'
                ScriptName              = Model_dict['Name'] + '.in'
                Line_dict               = ct.load_predicted_lines_individual(ScriptName, ScriptFolder)
                
                #Coloring scheme
                parameter_divider = zGas
                color = Grid_Values['zGas'].index(parameter_divider)
                label = r'$Z = {logage}$'.format(logage = round(float(parameter_divider) * 0.02, 3))
                
                x_values, y_values      = Ar_S_model(Line_dict, threshold = 4)
                
                if (x_values != None) and (y_values != None):
                
                    dz.data_plot(x_values, y_values, color=dz.ColorVector[2][color], label=label, markerstyle='o')
    
                    x_linealFitting = hstack([x_linealFitting, x_values])
                    y_linealFitting = hstack([y_linealFitting, y_values])

#Lineal model
lineal_mod          = LinearModel(prefix='lineal_')
Lineal_parameters   = lineal_mod.guess(y_linealFitting, x=x_linealFitting)
x_lineal            = linspace(np_min(x_linealFitting), np_max(x_linealFitting), 100)
y_lineal            = Lineal_parameters['lineal_slope'].value * x_lineal + Lineal_parameters['lineal_intercept'].value
dz.data_plot(x_lineal, y_lineal, label='Lineal fitting', color = 'black', linestyle='-')

# #Plot fitting formula
formula = r"$log\left(Ar^{{+2}}/Ar^{{+3}}\right) = {m} \cdot log\left(S^{{+2}}/S^{{+3}}\right) + {n}$".format(m=round(Lineal_parameters['lineal_slope'].value,3), n=round(Lineal_parameters['lineal_intercept'].value, 3))
# formula2 = r"$m = {m} \pm {merror}; n = {n} \pm {nerror}$".format(m=round(m[0],3), merror=round(m_err[0],3), n=round(n[0],3), nerror=round(n_err[0],3))
dz.Axis.text(0.50, 0.15, formula, transform=dz.Axis.transAxes, fontsize=20) 
# dz.Axis.text(0.50, 0.08, formula2, transform=dz.Axis.transAxes, fontsize=20) 

#Plot wording
xtitle  = r'$log([SIII]/[SIV])$'
ytitle  =  r'$log([ArIII]/[ArIV])$'
title   = 'Argon - Sulfur ionic relation in Cloudy photoionization models'
dz.FigWording(xtitle, ytitle, title, axis_Size = 20.0, title_Size = 20.0, legend_size=20.0, legend_loc='upper left')
# dz.Axis.set_xlim(0,6)
# dz.Axis.set_ylim(0,6)

#Display figure
dz.display_fig()

print 'Data treated'