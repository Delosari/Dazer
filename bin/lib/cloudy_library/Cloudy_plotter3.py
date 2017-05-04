from collections                        import OrderedDict
from lmfit                              import Parameters, minimize, report_fit
from lmfit.models                       import LinearModel
from numpy                              import log10 as nplog10, zeros, min, max, linspace, array, concatenate, isfinite, greater, hstack, min as np_min, max as np_max, power, pi, isinf
import pandas                           as pd
from Math_Libraries.bces_script         import bces
from Plotting_Libraries.dazer_plotter   import Plot_Conf
from cloudy_library.cloudy_methods      import Cloudy_Tools
import pyneb as pn

# Tem = 10000
# ne = 100
# 
# S4 = pn.Atom('S', 4)
# 
# Emis105000um = S4.getEmissivity(Tem, ne, wave=105000)
# 
# print Emis105000um

def Ar_S_model(Line_dict, threshold = 0.0, z = None):
     
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
     
    print 'Entran', x_axis, y_axis, x_axis < threshold
     
    if (y_axis < threshold):
                 
        return x_axis, y_axis
       
    else:
           
        return None, None
 
def import_popstar_data():
     
    FilesFolder                 = '/home/vital/Dropbox/Astrophysics/Lore/PopStar/' 
    TableAddressn100            = 'mnras0403-2012-SD2_clusterTable2.txt'
     
    Frame_MetalsEmission        = pd.read_csv(FilesFolder + TableAddressn100, delim_whitespace = True)
     
    nH          = 10.0 #cm^-3
    c           = 29979245800.0 #cm/s
    pc_to_cm    = 3.0856776e18 #cm/pc
     
    Frame_MetalsEmission['logR_cm'] = nplog10(Frame_MetalsEmission['logR'] * pc_to_cm)
    Frame_MetalsEmission['Q']       = nplog10(power(10, Frame_MetalsEmission['logU']) * 4 * pi * c * nH * power(Frame_MetalsEmission['logR'] * pc_to_cm, 2))
 
    return Frame_MetalsEmission
 
def Ar_S_model_pyneb(Line_dict, diags, Ar3_atom, Ar4_atom, S3_atom, S4_atom):
 
    TSIII, NSII = diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+', diag_den = '[SII] 6731/6716',
                                            value_tem = Line_dict['6312A'] / (Line_dict['9068.62A'] + Line_dict['9532A']),
                                            value_den = Line_dict['6731A']/Line_dict['6716A']) 
 
    TOIII, NSII_2 = diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+', diag_den = '[SII] 6731/6716',
                                            value_tem = Line_dict['4363A']/(Line_dict['5007A'] + Line_dict['4959A']),
                                            value_den = Line_dict['6731A']/Line_dict['6716A']) 
 
         
    Ar3         = Ar3_atom.getIonAbundance(int_ratio = Line_dict['7751.11A'] + Line_dict['7135A'], tem=TSIII, den=NSII, to_eval = 'L(7751) + L(7136)', Hbeta = Line_dict['4861.36A'])
    Ar4         = Ar4_atom.getIonAbundance(int_ratio = Line_dict['4740.12A'] + Line_dict['4711.26A'], tem=TOIII, den=NSII, to_eval = 'L(4740) + L(4711)', Hbeta = Line_dict['4861.36A'])                             
          
    S3          = S3_atom.getIonAbundance(int_ratio = (Line_dict['9068.62A'] + Line_dict['9532A']), tem=TSIII, den=NSII, to_eval = 'L(9069)+L(9531)', Hbeta = Line_dict['4861.36A'])
    S4          = S4_atom.getIonAbundance(int_ratio = (Line_dict['10.51m']), tem=TOIII, den=NSII, wave = 105000., Hbeta = Line_dict['4861.36A'])
    x_axis      = nplog10(Ar3) - nplog10(Ar4)
    y_axis      = nplog10(S3) - nplog10(S4)
     
    indexes     = x_axis>0.0
  
    x_values, y_values = x_axis[indexes][0], y_axis[indexes][0]
    
    if (isinf(x_values) == False) & (isinf(y_values) == False):
        
        if (x_values < 3) and (y_values < 3):
            return x_values, y_values
        else:
            return None, None    
    else:
        return None, None
    
pn.atomicData.setDataFile('s_iii_coll_HRS12.dat')
 
diags  = pn.Diagnostics()
Ar3 = pn.Atom('Ar', 3)
Ar4 = pn.Atom('Ar', 4)
S3 = pn.Atom('S', 3)
S4 = pn.Atom('S', 4)
 
dz     = Plot_Conf()
ct     = Cloudy_Tools()
diags  = pn.Diagnostics()
 
#Declare the number of colors
dz.FigConf(n_colors = 7)
 
#Define script name and location
ScriptFolder    = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Total_Q_R_Grid/'
 
Grid_Values                 = OrderedDict()
# [ 5.    5.48  5.7   5.85  6.    6.1   6.18  6.24  6.3   6.35  6.4   6.44 6.48  6.51  6.54  6.57  6.6   6.63  6.65  6.68  6.7   6.72]
Grid_Values['age']          = ['5.', '5.48', '5.7', '5.85',  '6.', '6.18', '6.3', '6.4', '6.51', '6.57', '6.63', '6.68', '6.72' ] 
Grid_Values['clus_mass']    = ['12000.', '20000.', '40000.', '60000.', '100000.', '200000.'] 
Grid_Values['zGas']         = ['0.0001', '0.0004', '0.004', '0.008'] 
Grid_Values['zStars']       = ['-2.1'] 
 
#Data from popstar
Frame_MetalsEmission = import_popstar_data()
 
x_linealFitting = array([])
y_linealFitting = array([])
 
#Fore the each grid point run a cloudy simulation
Model_dict = OrderedDict()
for age in Grid_Values['age']:
    for mass in Grid_Values['clus_mass']:
        for zGas in Grid_Values['zGas']:                                        
            for zStar in Grid_Values['zStars']:
                                 
                Model_dict['Name']          = 'S_Ar_Test_' + 'age'+ age + '_zStar' + zStar + '_clusMass' + mass + '_zGas' + zGas
                Model_dict['age']           = age
                Model_dict['zGas']          = zGas
                Model_dict['Metals_frac']   = str(float(zGas) / 0.02)
                Model_dict['zGas']          = zGas
                Model_dict['zStars']        = zStar
     
                index = (Frame_MetalsEmission["Z"] == float(zGas)) & (Frame_MetalsEmission["M_Msun"] == float(mass)) & (Frame_MetalsEmission["t"] == float(age))
                Model_dict['Q']             = Frame_MetalsEmission.loc[index, 'Q'].values[0]
                Model_dict['R']             = Frame_MetalsEmission.loc[index, 'logR_cm'].values[0]    
     
                ScriptName              = Model_dict['Name'] + '.in'
                Line_dict               = ct.load_predicted_lines_individual(ScriptName, ScriptFolder)
                 
#                 for coso in Line_dict:
#                     print coso
                 
                #Coloring scheme
                parameter_divider = zGas
                color = Grid_Values['zGas'].index(parameter_divider)
                label = r'$Z = {logage}$'.format(logage = Model_dict['zGas'])
                 
#                 x_values, y_values      = Ar_S_model(Line_dict, threshold = 4, z = float(Model_dict['zGas']))
                x_values, y_values      = Ar_S_model_pyneb(Line_dict, diags, Ar3, Ar4, S3, S4)
                 
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