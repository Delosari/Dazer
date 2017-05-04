from numpy                              import linspace, argmin, delete, mean, random
from uncertainties.unumpy               import uarray, nominal_values, std_devs

from Astro_Libraries.Abundances_Class   import Parametrized_Emissivities
from CodeTools.PlottingManager          import myPickle
from Math_Libraries.FittingTools        import Linear_Regressions


def synthetic_abundances():
    
    error_factor = 1.5
    
    N       = 50
    m_true  = 240
    n_true  = 0.24
    
    x_min   = 6.00e-05
    x_max   = 0.0004
    
    sig_x   = 6.5e-6
    sig_y   = 0.015

    x_true  = linspace(x_min, x_max, N)
    y_true  = m_true * x_true + n_true
    
    x_noise = random.normal(0.0, sig_x, N)
    y_noise = random.normal(0.0, sig_y, N)
    
    x_obs   = x_true + x_noise
    y_obs   = y_true + y_noise
    x_error = random.normal(0.0, sig_x, N)
    y_error = random.normal(0.0, sig_y, N)
        
    return x_obs, x_error, y_obs, y_error
    

#Generate dazer object

pv      = myPickle()
lr      = Linear_Regressions()
ch_an   = Parametrized_Emissivities()

#WMAP Primoridal Helium abundance
Y_WMAP  = 0.240

#Load the data from Fabian's thesis
Oxygen_abund, Oxygen_abundError, Nitrogen_abund, Nitrogen_abundError, Ymagnitude, Yerror, y_magnitude, y_error  = ch_an.Abundances_FabianSample()
# Oxygen_abund, Oxygen_abundError, Ymagnitude, Yerror     = synthetic_abundances() 

#Remove Izwicky18 from the sample
Izw18_index = argmin(Ymagnitude)
Oxygen_abund, Ymagnitude, Oxygen_abundError, Yerror, y_magnitude, y_error = delete(Oxygen_abund,  Izw18_index), delete(Ymagnitude,  Izw18_index), delete(Oxygen_abundError,  Izw18_index), delete(Yerror,  Izw18_index), delete(y_magnitude,  Izw18_index), delete(y_error,  Izw18_index)



#Define plot frame and colors
pv.FigFormat_One(ColorConf = 'Day_Elemental')

#WMAP Primoridal Helium abundance
Y_WMAP                          = 0.2470

He_vector   = uarray(y_magnitude,y_error)
O_Vector    = uarray(Oxygen_abund,Oxygen_abundError)

Y_mine      = (4 * He_vector * (1-25*O_Vector)) / (1 + 4*He_vector)

#Plotting the data
pv.DataPloter_One(nominal_values(O_Vector), nominal_values(Y_mine),     'Oxygen abundance',     pv.Color_Vector[2][1],  LineStyle=None,  XError=std_devs(O_Vector), YError=std_devs(Y_mine))
# pv.DataPloter_One(Oxygen_abund, Ymagnitude,     'Oxygen abundance',     pv.Color_Vector[2][1],  LineStyle=None,  XError=Oxygen_abundError, YError=Yerror)
pv.DataPloter_One(0.0,          Y_WMAP,         'True data',            pv.Color_Vector[1],     LineStyle=None)

#Linear regression
x_regresion_range   = linspace(0.0, max(Oxygen_abund)*1.10, 20)

#Load data for regressions
lr.load_obs_data(nominal_values(O_Vector), nominal_values(Y_mine), std_devs(O_Vector), std_devs(Y_mine))


# #Linfit regression
# linfit_dict         = lr.perform_regression(Methodology = 'linfit')
# y_regression_range  = linfit_dict['m'] * x_regresion_range + linfit_dict['n']
# label = linfit_dict['methodology'] + ': ' + str(round(linfit_dict['n'],5)) + r'$\pm$' + str(round(linfit_dict['n_error'],5))
# pv.DataPloter_One(x_regresion_range, y_regression_range, label, pv.Color_Vector[2][2])

#Bces regression    
bces_dict           = lr.perform_regression(Methodology= 'bces')

for i in range(len(bces_dict['m'])):
    y_regression_range  = bces_dict['m'][i] * x_regresion_range + bces_dict['n'][i]
    label = bces_dict['methodology'][i] + ' ' +str(i) + ': ' + str(round(bces_dict['n'][i],5)) + r'$\pm$' + str(round(bces_dict['n_error'][i],5))
    pv.DataPloter_One(x_regresion_range, y_regression_range, label, pv.Color_Vector[2][3+i])
    print 'Methodology: ', bces_dict['methodology'][i] + ' m:' + str(round(bces_dict['m'][i],5)) + r'$\pm$' + str(round(bces_dict['n'][i],5))

# #Scipy regression    
# scipy_dict              = lr.perform_regression(Methodology= 'scipy')
# y_regression_range      = scipy_dict['m'] * x_regresion_range + scipy_dict['n']
# label                   = scipy_dict['methodology'] + ': ' + str(round(scipy_dict['n'],5)) + r'$\pm$' + str(round(scipy_dict['n_error'],5))
# pv.DataPloter_One(x_regresion_range, y_regression_range, label, 'red')
# 
# #kmpfit regression    
# kmpfit_dict             = lr.perform_regression(Methodology= 'kmpfit')
# y_regression_range      = kmpfit_dict['m'] * x_regresion_range + kmpfit_dict['n']
# label                   = kmpfit_dict['methodology'] + ': ' + str(round(kmpfit_dict['n'],5)) + r'$\pm$' + str(round(kmpfit_dict['n_error'],5))
# pv.DataPloter_One(x_regresion_range, y_regression_range, label, 'green')

#kmpfit regression    
# chisqBayes_dict         = lr.perform_regression(Methodology= 'Inference - ChiSq')
# y_regression_range      = chisqBayes_dict['m'] * x_regresion_range + chisqBayes_dict['n']
# label                   = chisqBayes_dict['methodology'] + ': ' + str(round(chisqBayes_dict['n'],5)) + r'$\pm$' + str(round(chisqBayes_dict['n_error'],5))
# pv.DataPloter_One(x_regresion_range, y_regression_range, label, 'black')

# #Kelly bces bayesian regression    
# bcesBayes_dict          = lr.perform_regression(Methodology= 'kelly')
# y_regression_range      = bcesBayes_dict['m'] * x_regresion_range + bcesBayes_dict['n']
# label                   = bcesBayes_dict['methodology'] + ': ' + str(round(bcesBayes_dict['n'],5)) + r'$\pm$' + str(round(bcesBayes_dict['n_error'],5))
# pv.DataPloter_One(x_regresion_range, y_regression_range, label, 'black')

# outliers_dict           = lr.perform_regression(Methodology= 'Inference - Outliers')
# x_values, y_values      = outliers_dict['x_coords_outliers'], outliers_dict['y_coords_outliers']
# label                   = 'Outliers Krough'
# pv.Axis1.scatter(x_values, y_values,  facecolors='#A60628', edgecolors='#7A68A6', label='Outliers', lw=1, s=100, alpha=0.2)


#Labels and wording
pv.Labels_Legends_One(Plot_Title =  'Test Yp linear regressions',  Plot_xlabel = 'Oxygen abundance', Plot_ylabel = 'Helium mass fraction', LegendLocation='best')         

pv.Axis1.set_xlim(0, max(Oxygen_abund)*1.10)
pv.Axis1.set_ylim(0.0, 0.5)

pv.DisplayFigure()