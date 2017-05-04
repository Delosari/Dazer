'''
Created on Apr 17, 2016

@author: vital
'''
from collections                        import OrderedDict

from numpy                              import log10 as nplog10, zeros, min, max, linspace, array, concatenate, isfinite, greater
from pandas import Series
from uncertainties.unumpy               import uarray, nominal_values, std_devs, log10 as umlog10

from CodeTools.PlottingManager          import myPickle
from ManageFlow                         import DataToTreat
from Math_Libraries.bces_script         import bces
from Plotting_Libraries.dazer_plotter   import Plot_Conf
from cloudy_library.cloudy_methods      import Cloudy_Tools
import pyneb as pn


def import_data_from_objLog_triple(FilesList, pv):
    
    Valid_objects       = []
    List_Abundances     = ['ArIII_HII', 'ArIV_HII', 'TSIII', 'SIII_HII']
    Empty_Array         = [0] * (len(List_Abundances) + 1)
    
    #Dictionary of dictionaries to store object abundances
    Abund_dict = OrderedDict()
    for abund in List_Abundances:
        Abund_dict[abund] = []

    #Loop through files
    for i in range(len(FilesList)):
        
        #Analyze file address
        CodeName, FileName, FileFolder  = pv.Analyze_Address(FilesList[i])  
      
        #Loop through abundances in the log
        All_observed    = True
        Empty_Array[0]  = CodeName
        for j in range(len(List_Abundances)):
            abundance           = List_Abundances[j]
            Empty_Array[j+1]    = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter = abundance, Assumption = 'float')
            
            #If the abundance was measure store it 
            if Empty_Array[j+1] == None:
                All_observed = False
                
        if All_observed:
            Valid_objects.append(array(Empty_Array, copy=True))
            
    return array(Valid_objects)

def Ar_S_model(Line_dict, threshold = 0.0):
    
    R3          = (Line_dict['9068.62A'] + Line_dict['9532A']) / Line_dict['6312A']
    TSIII_4     = (0.5147 + 0.0003187 * R3 + (23.6404 / R3)) 
    TOIII_4     = ((TSIII_4 + 0.0846)/1.0807)

    TSIII       = TSIII_4 * 10000
    TOIII       = TOIII_4 * 10000
        
    logS3HI     = nplog10(Line_dict['10.51m']/Line_dict['4861.36A']) + 6.3956 + 0.0416/TOIII_4 - 0.4216 * nplog10(TOIII_4) - 12
    logS2HI     = nplog10((Line_dict['9068.62A'] + Line_dict['9532A']) / Line_dict['4861.36A']) + 6.012 + 0.6309 / TSIII_4 - 0.5722 * nplog10(TSIII_4)  -12 

    logAr3HI    = nplog10(Line_dict['4740.12A'] / Line_dict['4861.36A']) +  5.705 + 1.246/TOIII_4 - 0.156 * nplog10(TOIII_4) - 12 
    logAr2HI    = nplog10(Line_dict['7135A'] / Line_dict['4861.36A']) +  6.157 + 0.808/TSIII_4 - 0.508 * nplog10(TSIII_4) - 12 

    x_axis      = logS2HI - logS3HI
    y_axis      = logAr2HI - logAr3HI
    
    indexes     = x_axis>0.0

    return x_axis[indexes], y_axis[indexes], TSIII[indexes], TOIII[indexes]


pv      = myPickle()
dz      = Plot_Conf()
ct      = Cloudy_Tools()
diags   = pn.Diagnostics()

#Define data type and location
Catalogue_Dic                   = DataToTreat()
Pattern                         = Catalogue_Dic['Datatype'] + '.fits'

FilesList                   = pv.Folder_Explorer(Pattern, Catalogue_Dic['Obj_Folder'], CheckComputer=False)
  
Abundances_Matrix           = import_data_from_objLog_triple(FilesList, pv)
  
Objects                     = Abundances_Matrix[:,0]
ArIII_HII_array             = Abundances_Matrix[:,1]
ArIV_HII_array              = Abundances_Matrix[:,2]
Temps                       = Abundances_Matrix[:,3]    
SIII_HII_array              = Abundances_Matrix[:,4]   
    
logArII_ArIII               = umlog10(ArIII_HII_array/ArIV_HII_array)

#Define script name and location
ScriptFolder    = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Total_Grid/'

#Orion nebula has a metallicity of 2/3 solar

Grid_Values = OrderedDict()

Grid_Values['age']      = ['5.0', '5.5', '6.0', '6.5', '7.0', '7.5'] 
Grid_Values['zStars']   = ['-2.4', '-2.1', '-1.7', '-1.31'] 
Grid_Values['zGas']     = ['0.1', '0.31', '0.62'] 
Grid_Values['u']        = ['-3.0', '-2.0', '-1.0'] 

#Dictionary of dictionaries
Grid_frame =  ({k : Series(v) for k, v in Grid_Values.iteritems()})

#Trick to create a frame with different lengths
# Grid_frame = DataFrame({k : Series(v) for k, v in Grid_Values.iteritems()})

#Declare the number of colors
dz.FigConf(n_colors=len(Grid_Values['age'] ))


list_xvalues = array([])
list_yvalues = array([])
list_TSIII = array([])
list_TOIII = array([])

#Fore the each grid point run a cloudy simulation
Model_dict = OrderedDict()
for age in Grid_Values['age']:
    for zStar in Grid_Values['zStars']:
        for zGas in Grid_Values['zGas']:                                        
            for ionization_parameter in Grid_Values['u']: 
                #Format star
                #'S_Ar_Test_' + 'age'+ age + '_zStar' + zStar + '_zGas' + zGas + '_u' + ionization_parameter

                ScriptName = 'S_Ar_Test_' + 'age'+ age + '_zStar' + zStar + '_zGas' + zGas + '_u' + ionization_parameter + '.in'
    
#                 if age not in ['6.5', '7.0', '7.5']:
                if age not in ['9.5']:

#                     #Color by ionization factor
#                     Color = dz.ColorVector[2][Grid_Values['u'].index(ionization_parameter)]
#                     Label = 'u = ' +  ionization_parameter
                    
                    #Color by gas metallicity
#                     Color = dz.ColorVector[2][Grid_Values['zGas'].index(zGas)]
#                     Label = r'$z_{{gas}} = {z}$'.format(z = round(float(zGas) * 0.02 * 2.0/3.0, 3))
                    
                    #Color by age
                    Color = dz.ColorVector[2][Grid_Values['age'].index(age)]
                    Label = r'Stellar population $slog(age) = {age}$'.format(age = age)        
                    
#                     #Color by zStars
#                     Color = dz.ColorVector[2][Grid_Values['zStars'].index(zStar)]
#                     Label = r'Stellar metallicity $log(z_{{Star}}) = {zStar}$'.format(zStar = zStar)                       
                    
                    Line_dict = ct.load_predicted_lines_individual(ScriptName, ScriptFolder)
                    x_values, y_values, TSIII, TOIII = Ar_S_model(Line_dict, diags)
                    dz.data_plot(x_values, y_values, color=Color, label=Label, markerstyle='o')
#                     dz.data_plot(TSIII, y_values, color=Color, label=Label, markerstyle='o')
     
                    list_xvalues    = concatenate([list_xvalues, x_values])
                    list_yvalues    = concatenate([list_yvalues, y_values])
                    list_TSIII      = concatenate([list_TSIII, TSIII])
                    list_TOIII      = concatenate([list_TOIII, TOIII]) 
 
#Purge missing falues or non-physical
Not_infinite = isfinite(list_yvalues)

list_xvalues_clean  = list_xvalues[Not_infinite]
list_yvalues_clean  = list_yvalues[Not_infinite]
list_TSIII_clean    = list_TSIII[Not_infinite]

list_xvalues_Above  = greater(list_xvalues_clean, 0)

list_xvalues_clean_greater  = list_xvalues_clean[list_xvalues_Above]
list_yvalues_clean_greater  = list_yvalues_clean[list_xvalues_Above]
list_TSIII_clean_greater    = list_TSIII_clean[list_xvalues_Above] 


# # #----------------------Plotting temperatures
# #Plot wording
# xtitle  = r'$T[SIII] (K)$'
# ytitle  = r'$log(Ar^{+2}/Ar^{+3})$'
#   
# title   = r'Argon ionic abundance versus $S^{+2}$ temperature in Cloudy models'
# print len(Temps), len(logArII_ArIII)
# dz.data_plot(nominal_values(Temps), nominal_values(logArII_ArIII),  color=dz.ColorVector[1], label='Observations', markerstyle='o', x_error=std_devs(Temps), y_error=std_devs(logArII_ArIII))
# dz.FigWording(xtitle, ytitle, title, axis_Size = 20.0, title_Size = 20.0, legend_size=20.0, legend_loc='upper right')
# 'ArIons_vs_TSIII_Obs'
#   
# #Display figure
# # dz.display_fig()
# dz.savefig(output_address = '/home/vital/Dropbox/Astrophysics/Papers/Elemental_RegressionsSulfur/Cloudy_Models/ArIons_vs_TSIII_Ionization_Obs')
#      
# print 'Data treated'

#----------------------Plotting abundances
#Perform linear regression
zero_vector  = zeros(len(list_xvalues_clean_greater))
m ,n, m_err, n_err, covab = bces(list_xvalues_clean_greater, zero_vector, list_yvalues_clean_greater, zero_vector, zero_vector)
  
x_regresion         = linspace(0, max(list_xvalues_clean_greater), 50)
y_regression        = m[0] * x_regresion + n[0]
    
LinearRegression_Label = r'Linear fitting'.format(n = round(n[0],2) ,nerr = round(n_err[0],2))
dz.data_plot(x_regresion, y_regression, label=LinearRegression_Label, linestyle='--', color=dz.ColorVector[1])
        
logSII_SIII_theo            = m[0] * logArII_ArIII + n[0]
    
dz.data_plot(nominal_values(logArII_ArIII), nominal_values(logSII_SIII_theo),  color=dz.ColorVector[1], label='Observations', markerstyle='o', x_error=std_devs(logArII_ArIII), y_error=std_devs(logSII_SIII_theo))
 
# #Plot fitting formula
formula = r"$log\left(Ar^{{+2}}/Ar^{{+3}}\right) = {m} \cdot log\left(S^{{+2}}/S^{{+3}}\right) + {n}$".format(m='m', n='n')
formula2 = r"$m = {m} \pm {merror}; n = {n} \pm {nerror}$".format(m=round(m[0],3), merror=round(m_err[0],3), n=round(n[0],3), nerror=round(n_err[0],3))
dz.Axis.text(0.50, 0.15, formula, transform=dz.Axis.transAxes, fontsize=20) 
dz.Axis.text(0.50, 0.08, formula2, transform=dz.Axis.transAxes, fontsize=20) 
   
#Plot wording
xtitle  = r'$log(S^{+2}/S^{+3})$'
ytitle  = r'$log(Ar^{+2}/Ar^{+3})$'
title   = 'Argon - Sulfur ionic relation in Cloudy photoionization models'
dz.FigWording(xtitle, ytitle, title, axis_Size = 20.0, title_Size = 20.0, legend_size=20.0, legend_loc='best')
   
#Display figure
dz.display_fig()
# dz.savefig(output_address = '/home/vital/Dropbox/Astrophysics/Papers/Elemental_RegressionsSulfur/Cloudy_Models/ArIons_vs_SIons_Ionization_Obs_allages')
print 'Data treated' 
 
 

 
    
