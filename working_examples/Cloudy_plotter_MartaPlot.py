import pyneb as pn
from collections                            import OrderedDict
from lmfit.models                           import LinearModel
from numpy                                  import log10 as nplog10, linspace, array, hstack, min as np_min, max as np_max, power, pi, isinf, isnan, vstack, savetxt
from bin.lib.cloudy_library.cloudy_methods  import Cloudy_Tools 
from bin.dazer_methods                      import Dazer

def Ar_S_abundances_model(Line_dict, diags, Ar3_atom, Ar4_atom, S2_atom, S3_atom, S4_atom, threshold = 3):
    
    data_dict = OrderedDict()
    
    #Calculate sulfur properties
    TSIII, NSII = diags.getCrossTemDen(diag_tem = '[SIII] 6312/9200+', diag_den = '[SII] 6731/6716',
                                            value_tem = Line_dict['6312A'] / (Line_dict['9068.62A'] + Line_dict['9532A']),
                                            value_den = Line_dict['6731A']/Line_dict['6716A']) 
    #Calculate Oxygen properties
    TOIII, NSII_2 = diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+', diag_den = '[SII] 6731/6716',
                                            value_tem = Line_dict['4363A']/(Line_dict['5007A'] + Line_dict['4959A']),
                                            value_den = Line_dict['6731A']/Line_dict['6716A']) 
  
    #Determine ionic abundances
    Ar3 = Ar3_atom.getIonAbundance(int_ratio = Line_dict['7751.11A'] + Line_dict['7135A'], tem=TSIII, den=NSII, to_eval = 'L(7751) + L(7136)', Hbeta = Line_dict['4861.36A'])
    Ar4 = Ar4_atom.getIonAbundance(int_ratio = Line_dict['4740.12A'] + Line_dict['4711.26A'], tem=TOIII, den=NSII_2, to_eval = 'L(4740) + L(4711)', Hbeta = Line_dict['4861.36A'])                             
    S2  = S2_atom.getIonAbundance(int_ratio = (Line_dict['6716A'] + Line_dict['6731A']), tem=TSIII, den=NSII, to_eval = 'L(6716)+L(6731)', Hbeta = Line_dict['4861.36A'])
    S3  = S3_atom.getIonAbundance(int_ratio = (Line_dict['9068.62A'] + Line_dict['9532A']), tem=TSIII, den=NSII, to_eval = 'L(9069)+L(9531)', Hbeta = Line_dict['4861.36A'])
    S4  = S4_atom.getIonAbundance(int_ratio = (Line_dict['10.51m']), tem=TOIII, den=NSII_2, wave = 105000., Hbeta = Line_dict['4861.36A'])
     
    #Calculate the logaritmic axis for the plot
    x_axis      = nplog10(Line_dict['6716A'] + Line_dict['6731A'] + Line_dict['9068.62A'] + Line_dict['9532A'])
    y_axis      = 12 + nplog10(S2 + S3)
    
    if (isinf(x_axis) == False) & (isinf(y_axis) == False) & (isnan(x_axis) == False) & (isnan(y_axis) == False):
        
        data_dict['S+']             = S2
        data_dict['S2+']            = S3
        data_dict['S3+']            = S4
        data_dict['Ar2+']           = Ar3
        data_dict['Ar3+']           = Ar4         
        data_dict['Hbeta+']         = Line_dict['4861.36A']
        data_dict['Ar3_7751.11A']   = Line_dict['7751.11A']
        data_dict['Ar3_7135A']      = Line_dict['7135A']
        data_dict['Ar4_4740.12A']   = Line_dict['4740.12A']
        data_dict['Ar4_4711.26A']   = Line_dict['4711.26A']
        data_dict['S2_6716A']       = Line_dict['6716A']
        data_dict['S2_6731A']       = Line_dict['6731A']            
        data_dict['S3_9068.62A']    = Line_dict['9068.62A']        
        data_dict['S3_9532A']       = Line_dict['9532A']
        data_dict['S4_10.51m']      = Line_dict['10.51m']                   
           
        return x_axis, y_axis, data_dict
    else:
        return None, None, None

dz     = Dazer()
ct     = Cloudy_Tools()
diags  = pn.Diagnostics()
 
#Set atomic data and objects
pn.atomicData.setDataFile('s_iii_coll_HRS12.dat')
diags  = pn.Diagnostics()
Ar3 = pn.Atom('Ar', 3)
Ar4 = pn.Atom('Ar', 4)
S2 = pn.Atom('S', 2)
S3 = pn.Atom('S', 3)
S4 = pn.Atom('S', 4)

colors_list = ['#0072B2',
               '#009E73',
               '#D55E00',
               '#CC79A7',
               '#F0E442',
               '#56B4E9',
               '#bcbd22',
               '#7f7f7f',
               '#FFB5B8']

#Declare the number of colors
size_dict = {'axes.labelsize':20, 'legend.fontsize':18, 'font.family':'Times New Roman', 'mathtext.default':'regular', 'xtick.labelsize':18, 'ytick.labelsize':18}
dz.FigConf(plotSize = size_dict)

#Define script name and location
ScriptFolder    = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Total_Q_R_Grid2/'
  
#Load the conditions in the scripts
Grid_Values                 = OrderedDict()
Grid_Values['age']          = ['5.','5.48','5.7','5.85','6.','6.1','6.18','6.24','6.3','6.35','6.4','6.44','6.48','6.51','6.54','6.57','6.6','6.63','6.65','6.68','6.7','6.72']
Grid_Values['clus_mass']    = ['12000.', '20000.', '60000.', '100000.', '200000.'] 
Grid_Values['zGas']         = ['0.0001', '0.0004', '0.004', '0.008']#, '0.02', '0.05'] 
Grid_Values['zStars']       = ['-2.1'] 
  
#Data from popstar
Frame_MetalsEmission = ct.import_popstar_data()
  
#List to store the data
x_linealFitting = array([])
y_linealFitting = array([])
point_data = array([])

#-----------------------------------------ORIGINAL CASE---------------------------------------------
#Fore the each grid point run a cloudy simulation
Model_dict = OrderedDict()
for age in Grid_Values['age']:
    for mass in Grid_Values['clus_mass']:
        for zGas in Grid_Values['zGas']:                                        
            for zStar in Grid_Values['zStars']:
                     
                index                       = (Frame_MetalsEmission["Z"] == float(zGas)) & (Frame_MetalsEmission["M_Msun"] == float(mass)) & (Frame_MetalsEmission["t"] == float(age))
                Model_dict['Name']          = 'S_Ar_Test_' + 'age'+ age + '_zStar' + zStar + '_clusMass' + mass + '_zGas' + zGas
                Model_dict['age']           = age
                Model_dict['zGas']          = zGas
                Model_dict['Metals_frac']   = str(float(zGas) / 0.02)
                Model_dict['zGas']          = zGas
                Model_dict['zStars']        = zStar
      
                Model_dict['Q']             = Frame_MetalsEmission.loc[index, 'Q'].values[0]
                Model_dict['R']             = Frame_MetalsEmission.loc[index, 'logR_cm'].values[0]    
      
                ScriptName                  = Model_dict['Name'] + '.in'
                Line_dict                   = ct.load_predicted_lines_individual(ScriptName, ScriptFolder)
                  
                #Coloring scheme
                parameter_divider = zGas
                color = Grid_Values['zGas'].index(parameter_divider)
                label = r'$Z = {logage}$'.format(logage = Model_dict['zGas'])
                  
                #Calculate the grid point abundances
                #x_values, y_values         = Ar_S_model(Line_dict, threshold = 4, z = float(Model_dict['zGas']))
                x_values, y_values, data_dict = Ar_S_abundances_model(Line_dict, diags, Ar3, Ar4, S2, S3, S4, 3)
                
                if (x_values != None) and (y_values != None):
                    print data_dict.keys()
                    data_list = [float(age), float(mass), float(zGas), float(zStar)] + data_dict.values()
                    print data_list
                    if len(point_data) == 0:
                        point_data      = array(data_list)
                    else:
                        point_data      = vstack((point_data, array(data_list)))
                    
                    x_linealFitting = hstack([x_linealFitting, x_values])
                    y_linealFitting = hstack([y_linealFitting, y_values])
                    dz.data_plot(x_values, y_values,  label=label, markerstyle='o', color=colors_list[color])
  
  
  
#Lineal model
lineal_mod          = LinearModel(prefix='lineal_')
Lineal_parameters   = lineal_mod.guess(y_linealFitting, x=x_linealFitting)
x_lineal            = linspace(np_min(x_linealFitting), np_max(x_linealFitting), 100)
y_lineal            = Lineal_parameters['lineal_slope'].value * x_lineal + Lineal_parameters['lineal_intercept'].value
dz.data_plot(x_lineal, y_lineal, label='Linear fitting', color = 'black', linestyle='-')
formula = r'$12 + log\left(\frac{{S}}{{H}}\right)$ = {m} $\cdot$ $log\left(S_{{23}}\right)$ + {n}'.format(m=round(Lineal_parameters['lineal_slope'].value,3),
                                            n=round(Lineal_parameters['lineal_intercept'].value, 3))
dz.Axis.text(0.35, 0.15, formula, transform=dz.Axis.transAxes, fontsize=20) 
  
#Plot wording
xtitle  =   r'$log\left(S_{23}\right)$'
ytitle  =   r'$12 + log\left(\frac{S}{H}\right)$'
title   =   'Argon - Sulfur ionic abundances\nfor a Z, Mass, log(t) cluster grid'
dz.FigWording(xtitle, ytitle, title, loc='best')

#save the data
savetxt('/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/output_values.txt', point_data, fmt=['%10f'] * 4 + ['%1.4e'] * len(data_dict))

#Display figure
dz.display_fig()

print 'Data treated otro'
