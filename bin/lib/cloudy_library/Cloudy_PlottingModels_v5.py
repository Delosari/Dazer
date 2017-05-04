from collections                        import OrderedDict

from numpy                              import log10 as nplog10, zeros, min, max, linspace, array, concatenate, isfinite, greater
from uncertainties.unumpy               import uarray, nominal_values, std_devs, log10 as umlog10

from CodeTools.PlottingManager          import myPickle
from ManageFlow                         import DataToTreat
from Math_Libraries.bces_script         import bces
from Plotting_Libraries.dazer_plotter   import Plot_Conf
from cloudy_library.cloudy_methods      import Cloudy_Tools
import pyneb as pn


pn.atomicData.setDataFile('s_iii_coll_HRS12.dat')

def Figure_Legends_Colors(ColorVector):
  
    Model_dict      = OrderedDict()
    Legends_dict    = OrderedDict()
    Colors_dict     = OrderedDict()
  
    Model_dict['_z0.004_age5.0'] = 'log age = 5.0 log z = -2.4, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.004_age5.5'] = 'log age = 5.5 log z = -2.4, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.004_age6.0'] = 'log age = 6.0 log z = -2.4, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.004_age6.5'] = 'log age = 6.5 log z = -2.4, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.004_age7.0'] = 'log age = 7.0 log z = -2.4, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.004_age7.5'] = 'log age = 7.5 log z = -2.4, file="sp-kro_z0001-z05_stellar.mod"'
       
    Model_dict['_z0.008_age5.0'] = 'log age = 5.0 log z = -2.1, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.008_age5.5'] = 'log age = 5.5 log z = -2.1, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.008_age6.0'] = 'log age = 6.0 log z = -2.1, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.008_age6.5'] = 'log age = 6.5 log z = -2.1, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.008_age7.0'] = 'log age = 7.0 log z = -2.1, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.008_age7.5'] = 'log age = 7.5 log z = -2.1, file="sp-kro_z0001-z05_stellar.mod"'
       
    Model_dict['_z0.02_age5.0'] = 'log age = 5.0 log z = -1.7, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.02_age5.5'] = 'log age = 5.5 log z = -1.7, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.02_age6.0'] = 'log age = 6.0 log z = -1.7, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.02_age6.5'] = 'log age = 6.5 log z = -1.7, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.02_age7.0'] = 'log age = 7.0 log z = -1.7, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.02_age7.5'] = 'log age = 7.5 log z = -1.7, file="sp-kro_z0001-z05_stellar.mod"'
      
    Model_dict['_z0.05_age5.0'] = 'log age = 5.0 log z = -1.31, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.05_age5.5'] = 'log age = 5.5 log z = -1.31, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.05_age6.0'] = 'log age = 6.0 log z = -1.31, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.05_age6.5'] = 'log age = 6.5 log z = -1.31, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.05_age7.0'] = 'log age = 7.0 log z = -1.31, file="sp-kro_z0001-z05_stellar.mod"'
    Model_dict['_z0.05_age7.5'] = 'log age = 7.5 log z = -1.31, file="sp-kro_z0001-z05_stellar.mod"'
     
    Legends_dict['_z0.004_age5.0'] = 'log(age) = 5.0'
    Legends_dict['_z0.004_age5.5'] = 'log(age) = 5.5'
    Legends_dict['_z0.004_age6.0'] = 'log(age) = 6.0'
    Legends_dict['_z0.004_age6.5'] = 'log(age) = 6.5'
    Legends_dict['_z0.004_age7.0'] = 'log(age) = 7.0'
    Legends_dict['_z0.004_age7.5'] = 'log(age) = 7.5'
        
    Legends_dict['_z0.008_age5.0'] = 'log(age) = 5.0'
    Legends_dict['_z0.008_age5.5'] = 'log(age) = 5.5'
    Legends_dict['_z0.008_age6.0'] = 'log(age) = 6.0'
    Legends_dict['_z0.008_age6.5'] = 'log(age) = 6.5'
    Legends_dict['_z0.008_age7.0'] = 'log(age) = 7.0'
    Legends_dict['_z0.008_age7.5'] = 'log(age) = 7.5'
        
    Legends_dict['_z0.02_age5.0'] = 'log(age) = 5.0'
    Legends_dict['_z0.02_age5.5'] = 'log(age) = 5.5'
    Legends_dict['_z0.02_age6.0'] = 'log(age) = 6.0'
    Legends_dict['_z0.02_age6.5'] = 'log(age) = 6.5'
    Legends_dict['_z0.02_age7.0'] = 'log(age) = 7.0'
    Legends_dict['_z0.02_age7.5'] = 'log(age) = 7.5'
       
    Legends_dict['_z0.05_age5.0'] = 'log(age) = 5.0'
    Legends_dict['_z0.05_age5.5'] = 'log(age) = 5.5'
    Legends_dict['_z0.05_age6.0'] = 'log(age) = 6.0'
    Legends_dict['_z0.05_age6.5'] = 'log(age) = 6.5'
    Legends_dict['_z0.05_age7.0'] = 'log(age) = 7.0'
    Legends_dict['_z0.05_age7.5'] = 'log(age) = 7.5'
       
    Colors_dict['_z0.004_age5.0'] = ColorVector[2][0]
    Colors_dict['_z0.004_age5.5'] = ColorVector[2][1]
    Colors_dict['_z0.004_age6.0'] = ColorVector[2][2]
    Colors_dict['_z0.004_age6.5'] = ColorVector[2][3]
    Colors_dict['_z0.004_age7.0'] = ColorVector[2][4]
    Colors_dict['_z0.004_age7.5'] = ColorVector[2][5]
        
    Colors_dict['_z0.008_age5.0'] = ColorVector[2][0]
    Colors_dict['_z0.008_age5.5'] = ColorVector[2][1]
    Colors_dict['_z0.008_age6.0'] = ColorVector[2][2]
    Colors_dict['_z0.008_age6.5'] = ColorVector[2][3]
    Colors_dict['_z0.008_age7.0'] = ColorVector[2][4]
    Colors_dict['_z0.008_age7.5'] = ColorVector[2][5]
        
    Colors_dict['_z0.02_age5.0'] = ColorVector[2][0]
    Colors_dict['_z0.02_age5.5'] = ColorVector[2][1]
    Colors_dict['_z0.02_age6.0'] = ColorVector[2][2]
    Colors_dict['_z0.02_age6.5'] = ColorVector[2][3]
    Colors_dict['_z0.02_age7.0'] = ColorVector[2][4]
    Colors_dict['_z0.02_age7.5'] = ColorVector[2][5]
       
    Colors_dict['_z0.05_age5.0'] = ColorVector[2][0]
    Colors_dict['_z0.05_age5.5'] = ColorVector[2][1]
    Colors_dict['_z0.05_age6.0'] = ColorVector[2][2]
    Colors_dict['_z0.05_age6.5'] = ColorVector[2][3]
    Colors_dict['_z0.05_age7.0'] = ColorVector[2][4]
    Colors_dict['_z0.05_age7.5'] = ColorVector[2][5]

#     Legends_dict['_z0.004_age5.0'] = 'z = 0.004'
#     Legends_dict['_z0.004_age5.5'] = 'z = 0.004'
#     Legends_dict['_z0.004_age6.0'] = 'z = 0.004'
#     Legends_dict['_z0.004_age6.5'] = 'z = 0.004'
#     Legends_dict['_z0.004_age7.0'] = 'z = 0.004'
#     Legends_dict['_z0.004_age7.5'] = 'z = 0.004'
#       
#     Legends_dict['_z0.008_age5.0'] = 'z = 0.008'
#     Legends_dict['_z0.008_age5.5'] = 'z = 0.008'
#     Legends_dict['_z0.008_age6.0'] = 'z = 0.008'
#     Legends_dict['_z0.008_age6.5'] = 'z = 0.008'
#     Legends_dict['_z0.008_age7.0'] = 'z = 0.008'
#     Legends_dict['_z0.008_age7.5'] = 'z = 0.008'
#       
#     Legends_dict['_z0.02_age5.0'] = 'z = 0.02'
#     Legends_dict['_z0.02_age5.5'] = 'z = 0.02'
#     Legends_dict['_z0.02_age6.0'] = 'z = 0.02'
#     Legends_dict['_z0.02_age6.5'] = 'z = 0.02'
#     Legends_dict['_z0.02_age7.0'] = 'z = 0.02'
#     Legends_dict['_z0.02_age7.5'] = 'z = 0.02'
#      
#     Legends_dict['_z0.05_age5.0'] = 'z = 0.05'
#     Legends_dict['_z0.05_age5.5'] = 'z = 0.05'
#     Legends_dict['_z0.05_age6.0'] = 'z = 0.05'
#     Legends_dict['_z0.05_age6.5'] = 'z = 0.05'
#     Legends_dict['_z0.05_age7.0'] = 'z = 0.05'
#     Legends_dict['_z0.05_age7.5'] = 'z = 0.05'
#      
#     Colors_dict['_z0.004_age5.0'] = ColorVector[2][0]
#     Colors_dict['_z0.004_age5.5'] = ColorVector[2][0]
#     Colors_dict['_z0.004_age6.0'] = ColorVector[2][0]
#     Colors_dict['_z0.004_age6.5'] = ColorVector[2][0]
#     Colors_dict['_z0.004_age7.0'] = ColorVector[2][0]
#     Colors_dict['_z0.004_age7.5'] = ColorVector[2][0]
#       
#     Colors_dict['_z0.008_age5.0'] = ColorVector[2][1]
#     Colors_dict['_z0.008_age5.5'] = ColorVector[2][1]
#     Colors_dict['_z0.008_age6.0'] = ColorVector[2][1]
#     Colors_dict['_z0.008_age6.5'] = ColorVector[2][1]
#     Colors_dict['_z0.008_age7.0'] = ColorVector[2][1]
#     Colors_dict['_z0.008_age7.5'] = ColorVector[2][1]
#       
#     Colors_dict['_z0.02_age5.0'] = ColorVector[2][2]
#     Colors_dict['_z0.02_age5.5'] = ColorVector[2][2]
#     Colors_dict['_z0.02_age6.0'] = ColorVector[2][2]
#     Colors_dict['_z0.02_age6.5'] = ColorVector[2][2]
#     Colors_dict['_z0.02_age7.0'] = ColorVector[2][2]
#     Colors_dict['_z0.02_age7.5'] = ColorVector[2][2]
#      
#     Colors_dict['_z0.05_age5.0'] = ColorVector[2][3]
#     Colors_dict['_z0.05_age5.5'] = ColorVector[2][3]
#     Colors_dict['_z0.05_age6.0'] = ColorVector[2][3]
#     Colors_dict['_z0.05_age6.5'] = ColorVector[2][3]
#     Colors_dict['_z0.05_age7.0'] = ColorVector[2][3]
#     Colors_dict['_z0.05_age7.5'] = ColorVector[2][3]

    return Model_dict, Legends_dict, Colors_dict

def Ar_S_model_pyneb(Line_dict, diags, Ar3_atom, Ar4_atom, S3_atom, S4_atom):

    TSIII, NSII = diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+', diag_den = '[SII] 6731/6716',
                                            value_tem = Line_dict['6312.06A'] / (Line_dict['9068.62A'] + Line_dict['9532A']),
                                            value_den = Line_dict['S2_6731A']/Line_dict['S2_6716A']) 

    TOIII, NSII_2 = diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+', diag_den = '[SII] 6731/6716',
                                            value_tem = Line_dict['O3_4363A']/(Line_dict['O3_5007A'] + Line_dict['O3_4959A']),
                                            value_den = Line_dict['S2_6731A']/Line_dict['S2_6716A']) 

        
    Ar3         = Ar3_atom.getIonAbundance(int_ratio = Line_dict['Ar3_7751A'] + Line_dict['Ar3_7135A'], tem=TSIII, den=NSII, to_eval = 'L(7751) + L(7136)', Hbeta = Line_dict['H1_4861A'])
    Ar4         = Ar4_atom.getIonAbundance(int_ratio = Line_dict['Ar4_4740A'] + Line_dict['Ar4_4711A'], tem=TOIII, den=NSII, to_eval = 'L(4740) + L(4711)', Hbeta = Line_dict['H1_4861A'])                             
    
    S3          = S3_atom.getIonAbundance(int_ratio = (Line_dict['9068.62A'] + Line_dict['9532A']), tem=TSIII, den=NSII, to_eval = 'L(9069)+L(9531)', Hbeta = Line_dict['H1_4861A'])
    S4          = S3_atom.getIonAbundance(int_ratio = (Line_dict['10.51m']), tem=TOIII, den=NSII, to_eval = 'L(10.51)', Hbeta = Line_dict['H1_4861A'])
    
    x_axis      = nplog10(Ar3) - nplog10(Ar4)
    y_axis      = nplog10(S3) - nplog10(S4)
    
    indexes     = x_axis>0.0

    return x_axis[indexes], y_axis[indexes], TSIII[indexes], TOIII[indexes]

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

pv                              = myPickle()
dz                              = Plot_Conf()
ct                              = Cloudy_Tools()
diags                           = pn.Diagnostics()

#Define data type and location
Catalogue_Dic                   = DataToTreat()
Pattern                         = Catalogue_Dic['Datatype'] + '.fits'

#Define figure format
dz.FigConf(n_colors=6)

#Define script name and location
# ScriptFolder    = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Few_Models/'
ScriptFolder    = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Complete_Model/'
ScriptPrefix    = 'S_Ar_test'

#4 metallicities 0.004, 0.008, 0.02, 0.05
#5 ages 5.0, 5.5, 6.0, 6.5, 7.0, 7.5
Model_dict, Legends_dict, Colors_dict = Figure_Legends_Colors(dz.ColorVector)

list_xvalues = array([])
list_yvalues = array([])
list_TSIII = array([])
list_TOIII = array([])

#----Observations data
FilesList                   = pv.Folder_Explorer(Pattern, Catalogue_Dic['Obj_Folder'], CheckComputer=False)
  
Abundances_Matrix           = import_data_from_objLog_triple(FilesList, pv)
  
Objects                     = Abundances_Matrix[:,0]
ArIII_HII_array             = Abundances_Matrix[:,1]
ArIV_HII_array              = Abundances_Matrix[:,2]
Temps                       = Abundances_Matrix[:,3]    
SIII_HII_array              = Abundances_Matrix[:,4]   
  
print 'longitudes', len(Objects), len(ArIII_HII_array), len(ArIV_HII_array), len(Temps), len(SIII_HII_array)
  
logArII_ArIII               = umlog10(ArIII_HII_array/ArIV_HII_array)
   

for key in Model_dict.keys():
    
    #Declare scripting name
    ScriptName = ScriptPrefix + key + '.in'
    
    #Generate lines dictionary with the output data
    Line_dict = ct.load_predicted_lines(ScriptName, ScriptFolder)

    x_values, y_values, TSIII, TOIII = Ar_S_model(Line_dict, diags)
#     dz.data_plot(x_values, y_values, color=Colors_dict[key], label=Legends_dict[key], markerstyle='o')
    dz.data_plot(TSIII, y_values,  color=Colors_dict[key], label=Legends_dict[key], markerstyle='o')

    list_xvalues    = concatenate([list_xvalues, x_values])
    list_yvalues    = concatenate([list_yvalues, y_values])
    list_TSIII      = concatenate([list_TSIII, TSIII])
    list_TOIII      = concatenate([list_TOIII, TOIII])

Not_infinite = isfinite(list_yvalues)

list_xvalues_clean  = list_xvalues[Not_infinite]
list_yvalues_clean  = list_yvalues[Not_infinite]
list_TSIII_clean    = list_TSIII[Not_infinite]

list_xvalues_Above  = greater(list_xvalues_clean, 0)

list_xvalues_clean_greater  = list_xvalues_clean[list_xvalues_Above]
list_yvalues_clean_greater  = list_yvalues_clean[list_xvalues_Above]
list_TSIII_clean_greater    = list_TSIII_clean[list_xvalues_Above]

# # #----------------------Plotting temperatures
#Plot wording
xtitle  = r'$T[SIII] (K)$'
ytitle  = r'$log(Ar^{+2}/Ar^{+3})$'
 
title   = r'Argon ionic abundance versus $S^{+2}$ temperature in Cloudy models'
print len(Temps), len(logArII_ArIII)
dz.data_plot(nominal_values(Temps), nominal_values(logArII_ArIII),  color=dz.ColorVector[1], label='Observations', markerstyle='o', x_error=std_devs(Temps), y_error=std_devs(logArII_ArIII))
dz.FigWording(xtitle, ytitle, title, axis_Size = 20.0, title_Size = 20.0, legend_size=20.0, legend_loc='upper right')
'ArIons_vs_TSIII_Obs'


dz.Axis.set_xlim(5000,20000)

#Display figure
# dz.display_fig()
dz.savefig(output_address = '/home/vital/Dropbox/Astrophysics/Papers/Elemental_RegressionsSulfur/Cloudy_Models/ArIons_vs_TSIII_Obs')
    
print 'Data treated'

# #----------------------Plotting abundances
# #Perform linear regression
# zero_vector  = zeros(len(list_xvalues_clean_greater))
# m ,n, m_err, n_err, covab = bces(list_xvalues_clean_greater, zero_vector, list_yvalues_clean_greater, zero_vector, zero_vector)
# #     
# x_regresion         = linspace(0, max(list_xvalues_clean_greater), 50)
# y_regression        = m[0] * x_regresion + n[0]
#   
#   
# LinearRegression_Label = r'Linear fitting'.format(n = round(n[0],2) ,nerr = round(n_err[0],2))
# dz.data_plot(x_regresion, y_regression, label=LinearRegression_Label, linestyle='--', color=dz.ColorVector[1])
#      
#    
# logSII_SIII_theo            = m[0] * logArII_ArIII + n[0]
#   
#  
# dz.data_plot(nominal_values(logArII_ArIII), nominal_values(logSII_SIII_theo),  color=dz.ColorVector[1], label='Observations', markerstyle='o', x_error=std_devs(logArII_ArIII), y_error=std_devs(logSII_SIII_theo))
#          
# # #Plot fitting formula
# formula = r"$log\left(Ar^{{+2}}/Ar^{{+3}}\right) = {m} \cdot log\left(S^{{+2}}/S^{{+3}}\right) + {n}$".format(m='m', n='n')
# formula2 = r"$m = {m} \pm {merror}; n = {n} \pm {nerror}$".format(m=round(m[0],3), merror=round(m_err[0],3), n=round(n[0],3), nerror=round(n_err[0],3))
# dz.Axis.text(0.50, 0.15, formula, transform=dz.Axis.transAxes, fontsize=20) 
# dz.Axis.text(0.50, 0.08, formula2, transform=dz.Axis.transAxes, fontsize=20) 
#  
# #Plot wording
# xtitle  = r'$log(S^{+2}/S^{+3})$'
# ytitle  = r'$log(Ar^{+2}/Ar^{+3})$'
# title   = 'Argon - Sulfur ionic relation in Cloudy photoionization models'
# dz.FigWording(xtitle, ytitle, title, axis_Size = 20.0, title_Size = 20.0, legend_size=20.0, legend_loc='best')
# 
# 
# dz.Axis.set_xlim(-0.5, 5)
# dz.Axis.set_ylim(0.5, 7)
# 
# #Display figure
# # dz.display_fig()
# dz.savefig(output_address = '/home/vital/Dropbox/Astrophysics/Papers/Elemental_RegressionsSulfur/Cloudy_Models/ArIons_vs_SIons_Obs')
# print 'Data treated'