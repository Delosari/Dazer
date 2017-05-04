import numpy as np
from timeit import default_timer as timer

a               = 3.0
Te, ne          = np.empty(1000), np.empty(1000)
Te[:], ne[:]    = np.nan, np.nan

start = timer()
uno = np.isnan(a)
end = timer()
print 'uno', (end - start), uno

start = timer()
cuatro = np.isnan(np.sum(np.isnan(a)))
end = timer()
print 'cuatro', (end - start), cuatro

start = timer()
tres = np.sum(np.isnan(ne))
end = timer()
print 'tres', (end - start), tres


start = timer()
dos = np.isnan(Te)
end = timer()
print 'dos', (end - start), dos





# import numpy as np
# import pandas as pd
# 
# myDict      = {'a':1, 'b':2, 'c':3}
# mySeries    = pd.Series(data = [1,2,3], index = ['a', 'b', 'c']) 
# 
# if myDict.viewkeys() >= {'a', 'b'}:
#     print 'a and b are in dictionary'
# 
# set
# 
# if set(mySeries.index) >= {'a', 'b'}:
#     print 'Bingo'

#'L({})'.format(.x[x.find('_') + 1 :len(x)], coso))
# import pyneb as pn
# import numpy as np
# 
# S2      = pn.Atom('S',2)
# diags   = pn.Diagnostics()
# 
# #abund = S2.getIonAbundance(int_ratio = 0.03, tem=10000, den=100, to_eval = 'L(6731)+L(6716)', Hbeta = 1)   
# 
# abund = S2.getIonAbundance(int_ratio = 0.03, tem=10000, den=100, to_eval = 'L(6731)', Hbeta = 1)   
# print abund
# print S2.getIonAbundance(int_ratio = 0.03, tem=10000, den=100, wave = 6731, Hbeta = 1)

# import numpy as np
# import numexpr as ne
#  
# data = {}
# data['S2_6717A'] = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
# data['S2_6731A'] = np.array([2.0, 2.0, 2.0, 2.0, 2.0])
# 
# coso = data['S2_6717A']
# self.S2_atom.getIonAbundance(int_ratio=self.abunData.R_SII_abund, tem=Te, den=ne, to_eval = 'L(6731)+L(6716)', Hbeta = self.Hbeta_flux)   
# 
# 
# nan_occurences = np.isnan(coso)
# 
# print nan_occurences
# print np.sum(nan_occurences)

# print isinstance(coso, np.nan)
# 
# print np.isnan(coso)
# print np.all(np.isnan(coso))
# 
# expression = 's2_6717a/s2_6731a'
#  
# print ne.evaluate(expression, local_dict=data)

# c = np.array([a,b])
# import numexpr as ne
# In [146]: ne.evaluate('2*a*(b/c)**2', local_dict=var)
# Out[146]: array([ 0.88888889,  0.44444444,  4.        ])
# 
# print np.sum(np.array([a]), axis=0)
# print np.sum(c, axis=0)


# import pandas as pd
# import numpy as np
# from uncertainties.unumpy import uarray, nominal_values, std_devs, log10 as unum_log10, pow as unnumpy_pow
# from uncertainties import ufloat
# 
# data = {'Col1' : [ufloat(3,4),ufloat(6,0.5),np.nan,np.nan], 'Col2' : [10,20,30,30], 'Col3' : [100,50,-30,-50], 'Col4' : ['AAA', 'BBB', 'AAA', 'CCC']}
# df = pd.DataFrame(data=data, index = ['R1','R2','R3','R4'])
# idcs_nan = df['Col1'].isnull()
# original_data = df['Col1']
# 
# print idcs_nan
# 
# print len(idcs_nan)
# 
# df.insert(df.columns.get_loc('Col1') + 1, 'Col1' + '_err', np.full(len(idcs_nan), np.nan))
# 
# print 
# 
# print df
# mag = nominal_values(original_data[~idcs_nan].values)
# errors = std_devs(original_data[~idcs_nan].values)
# 
# df.loc[~idcs_nan, 'Col1']           = mag
# df.loc[~idcs_nan, 'Col1' + '_err']  = errors
# 
# print df

# print df
# idcs_nan = df['Col2'].isnull()
#                              
# print idcs_nan                          
#                              
# df.loc[~idcs_nan, 'Col3'] = uarray(df.loc[~idcs_nan, 'Col1'], df.loc[~idcs_nan, 'Col2'])
# df.loc[idcs_nan, 'Col3'] = df.loc[idcs_nan, 'Col1']
# 
# print df



# columns_toClean = np.array(['Col1','Col2','Col3','Col4', 'Col5'])
# 
# print 'estos', np.intersect1d(df.columns, columns_toClean, assume_unique=True) 
# 
# print df
# 
# np.NaN
# 
# print df.loc['R1']
# 
# df.loc['R1', columns_toClean] = np.NaN
# 
# print df.loc['R1']
# 
# print df


# import pandas as pd
# import numpy as np
# 
# data = {'Col1' : [4,5,6,7], 'Col2' : [10,20,30,40], 'Col3' : [100,50,-30,-50], 'Col4' : ['AAA', 'BBB', 'AAA', 'CCC']}
# df = pd.DataFrame(data=data, index = ['R1','R2','R3','R4'])
# df.loc['R2.5'] = [5.5, 25, 30, 'CCC']
# 
# print df
# target_lines = ['R1', 'R2']
# if df.index.isin(target_lines).sum() == len(target_lines):
#     print 'Bingo'



 
# from numpy import sqrt, pi, exp, linspace, random, empty, array
# from scipy.optimize import curve_fit
# from lmfit import Model
# from timeit import default_timer as timer
# 
# def gaussian(x, amp, cen, wid):
#     return amp * exp(-(x-cen)**2 /wid)
# 
# def gaussian_dicty(x, dicty):
#     return dicty['amp'] * exp(-(x-dicty['cen'])**2 /dicty['wid'])
# 
# init_vals = [1, 0, 1]
# x = linspace(-10,10)
# n_points = len(x) 
# gaussian_true = gaussian(x, 2.33, 0.21, 1.51)
# gmod = Model(gaussian)
# 
# print 'True', 2.33, 0.21, 1.51
# 
# print '--Single fit'
# start = timer()
# best_vals, covar = curve_fit(gaussian, x, gaussian_true + random.normal(0, 0.2, n_points), p0=init_vals)
# end = timer()
# print 'curve fit', best_vals, ' time ', (end - start) 
# start = timer()
# dicty = {'amp':1, 'cen':0, 'wid':1}
# best_vals2, covar = curve_fit(gaussian_dicty, x, gaussian_true + random.normal(0, 0.2, n_points), p0=dicty)
# end = timer()
# print 'curve fit dicty', best_vals2, ' time ', (end - start) 
# 
# start = timer()
# result = gmod.fit(gaussian_true + random.normal(0, 0.2, n_points), x=x, amp=1, cen=0, wid=1)
# end = timer()
# print 'lmfit', array(result.params.valuesdict().values()), (end - start) 


# print '--Single fit'
# start = timer()
# best_vals, covar = curve_fit(gaussian, x, gaussian_true + random.normal(0, 0.2, n_points), p0=init_vals)
# end = timer()
# print 'curve fit', best_vals, ' time ', (end - start) 
# start = timer()
# result = gmod.fit(gaussian_true + random.normal(0, 0.2, n_points), x=x, amp=1, cen=0, wid=1)
# end = timer()
# print 'lmfit', array(result.params.valuesdict().values()), (end - start) 
# 
# 
# print '--Bootstrap'
# results_matrix_curvefit = empty([3,1000])
# results_matrix_lmfit = empty([3,1000])
# start = timer()
# for i in range(1000):    
#     best_vals, covar = curve_fit(gaussian, x, gaussian_true + random.normal(0, 0.2, n_points), p0=init_vals)
#     results_matrix_curvefit[:,i] = best_vals
# end = timer()
# print 'curve fit', results_matrix_curvefit.mean(1), ' time ', (end - start) 
# 
# 
# start = timer()
# for i in range(1000):
#     result = gmod.fit(gaussian_true + random.normal(0, 0.2, n_points), x=x, amp=1, cen=0, wid=1)
#     results_matrix_lmfit[:,i] = array(result.params.valuesdict().values())
# end = timer()
# print 'lmfit', results_matrix_lmfit.mean(1), (end - start), results_matrix_lmfit.std(1)

# from numpy import exp, linspace, random, mean, ones, exp, log, pi, array, empty
# from lmfit import Minimizer, Parameters, minimize, Model
# import matplotlib.pyplot as plt
# from timeit import default_timer as timer
#  
#  
# def gaussian_curve(A, mu, sigma, x):           
#     return A * exp(-(x-mu)*(x-mu)/(2 * sigma * sigma))
#  
# def residual(p, x, y, err):
#     v = p.valuesdict()
#     return (gaussian_curve(v['A'], v['mu'], v['sigma'], x) - y) / err
#  
# def lnprob(p, x, y, err):
#     resid = residual(p, x, y, err)
#     return -0.5 * sum(((resid - y) / err)**2 + log(2 * pi * err**2))
#  
# n_points = 50
# iterations = 500
# wave = linspace(1, 17, n_points)
# random.seed(0)
# gmod = Model(gaussian_curve)
#  
# error_vector = random.normal(0, 0.1, n_points)
# err = mean(error_vector) * ones(n_points)
# flux_true = gaussian_curve(2.33, 8.2, 1.51, wave)
# flux_obs = flux_true + error_vector
#  
# plt.plot(wave, flux_obs, 'o')
#  
# params = Parameters()
# params.add('A', value = 2.0)
# params.add('mu', value = 8.0)
# params.add('sigma', value = 1.0)
#  
#  
# #-- Simple fitting
# fit_simple = minimize(residual, params, args=([wave, flux_obs, err])) #, method='powell'
# start = timer()
# fit_parames_simple = fit_simple.params.valuesdict()
# end = timer()
#  
# plt.plot(wave, gaussian_curve(fit_parames_simple['A'], fit_parames_simple['mu'], fit_parames_simple['sigma'], wave), 'r')
#  
# print '--Fit simple:', array(fit_simple.params.valuesdict().values()), 'time', (end - start)
#  
# #-- Bootstrap
# results_matrix_lmfit = empty([3,iterations])
# start = timer()
# for i in range(iterations):
#     result = minimize(residual, params, args=([wave, flux_obs, err]))    
#     results_matrix_lmfit[:,i] = array(result.params.valuesdict().values())
# end = timer()
# print '--Lmfit', results_matrix_lmfit.mean(1), (end - start)
#  
# #-- Complex fitting
# mini_posterior = Minimizer(lnprob, params, fcn_args = ([wave, flux_obs, err]))
# start = timer()
# fit_mcmc = mini_posterior.emcee(steps=iterations, params=params)
# end = timer()
# fit_parames_mcmc = fit_simple.params.valuesdict()
# output_params = fit_simple.params
# 
# for key in output_params:
#     print key + '_norm', output_params[key].value
#     print key + '_norm_er', output_params[key].stderr    
# 
# plt.plot(wave, gaussian_curve(fit_parames_mcmc['A'], fit_parames_mcmc['mu'], fit_parames_mcmc['sigma'], wave), 'b')
#  
# print '--Fit complex:', array(fit_mcmc.params.valuesdict().values()), 'time', (end - start) 
#  
# plt.show()





# from numpy import exp, linspace, random, sum, log, pi, array
# from lmfit import Parameters, Minimizer
# from timeit import default_timer as timer
#  
# def gaussian_curve(A, mu, sigma, x):           
#     return A * exp(-(x-mu)*(x-mu)/(2 * sigma * sigma))
#  
# def gaussian_residual(params, x, y, err):
#     return (gaussian_curve(params['A'], params['mu'], params['sigma'], x) - y) / err
#  
# def lnprob(p, x, y, err):
#     resid = gaussian_residual(params, x, y, err)
#     return -0.5 * sum((resid)**2 + log(2 * pi * err**2))
#  
# x = linspace(-10,10)
# n_points = len(x) 
#  
# gaussian_true = gaussian_curve(2.33, 0.21, 1.51, x)
# err = random.normal(0, 0.2, n_points)
# y_obs = gaussian_true + err
#  
# params = Parameters()
# params.add('A', value = 1)
# params.add('mu', value = 0)
# params.add('sigma', value = 1)
#          
# mini = Minimizer(lnprob, params, fcn_args = ([x, y_obs, err]))
#  
# out = mini.emcee(steps=1000, params=params)
# 
# print 'Hola'
# print array(out.params.valuesdict().values())


# def gaussian(x, amp, cen, wid):
#     return amp * exp(-(x-cen)**2 /wid)


# print '--Single fit'
# start = timer()
# best_vals, covar = curve_fit(gaussian, x, gaussian_true + random.normal(0, 0.2, n_points), p0=init_vals)
# end = timer()
# print 'curve fit', best_vals, ' time ', (end - start) 
# start = timer()
# result = gmod.fit(gaussian_true + random.normal(0, 0.2, n_points), x=x, amp=1, cen=0, wid=1)
# end = timer()
# print 'lmfit', array(result.params.valuesdict().values()), (end - start) 
#  
# print '--Bootstrap'
# results_matrix_curvefit = empty([3,1000])
# results_matrix_lmfit = empty([3,1000])
# start = timer()
# for i in range(1000):    
#     best_vals, covar = curve_fit(gaussian, x, gaussian_true + random.normal(0, 0.2, n_points), p0=init_vals)
#     results_matrix_curvefit[:,i] = best_vals
# end = timer()
# print 'curve fit', results_matrix_curvefit.mean(1), ' time ', (end - start) 
#  
# start = timer()
# for i in range(1000):
#     result = gmod.fit(gaussian_true + random.normal(0, 0.2, n_points), x=x, amp=1, cen=0, wid=1)
#     results_matrix_lmfit[:,i] = array(result.params.valuesdict().values())
# end = timer()
# print 'lmfit', results_matrix_lmfit.mean(1), (end - start)
# 
# import numpy as np
# import lmfit
# import matplotlib.pyplot as plt
# import corner
# 
# #Calculate a Gaussian
# def gauss(x, a_max, loc, sd):
#     return a_max * np.exp(-((x - loc) / sd)**2)
# 
# #The normalised residual for the data
# def residual(p, just_generative=False):
#     v = p.valuesdict()
#     generative = v['a'] + v['b'] * x
#     M = 0
#     while 'a_max%d' % M in v:
#         generative += gauss(x, v['a_max%d'%M], v['loc%d'%M], v['sd%d'%M])
#         M += 1
# 
#     if just_generative:
#         return generative
#     return (generative - y) / dy
# 
# #Create a Parameter set for the initial guesses
# def initial_peak_params(M):
#     p = lmfit.Parameters()
#     
#     # a and b give a linear background
#     a = np.mean(y)
#     b = 1
#     
#     # a_max, loc and sd are the amplitude, location and SD of each gaussian component
#     a_max = np.max(y)
#     loc = np.mean(x)
#     sd = (np.max(x) - np.min(x)) * 0.5
# 
#     p.add_many(('a', np.mean(y), True, 0, 10), ('b', b, True, 1, 15))
# 
#     for i in range(M):
#         p.add_many(('a_max%d'%i, 0.5 * a_max, True, 10, a_max),
#                    ('loc%d'%i, loc, True, np.min(x), np.max(x)),
#                    ('sd%d'%i, sd, True, 0.1, np.max(x) - np.min(x)))
#     return p
# 
# #This is the log-likelihood probability for the sampling.
# def lnprob(p):
#     resid = residual(p, just_generative=True)
#     return -0.5 * np.sum(((resid - y) / dy)**2 + np.log(2 * np.pi * dy**2))
# 
# 
# x = np.linspace(3, 7, 250)
# np.random.seed(0)
# y = 4 + 10 * x + gauss(x, 200, 5, 0.5) + gauss(x, 60, 5.8, 0.2)
# dy = np.sqrt(y)
# y += dy * np.random.randn(np.size(y))
# 
# plt.errorbar(x, y)
# 
# 
# p1 = initial_peak_params(1)
# mi1 = lmfit.minimize(residual, p1, method='differential_evolution')
# 
# lmfit.printfuncs.report_fit(mi1.params, min_correl=0.5)
# 
# plt.show()


#  
# T = 14000
# t = T/10000
#  
# R_S2 = 1.4528
# R_S2 = 1.35
#  
# a0 = 16.054 - 7.79/t - 11.32*t
# a1 = -22.66 + 11.08/t + 16.02*t
# b0 = -21.61 + 11.89/t + 14.59*t
# b1 = 9.17 - 5.09/t - 6.18 * t
#  
# n = 1000 * (R_S2 * a0 + a1) / (R_S2 * b0 + b1)
#  
# print n

# import pyneb as pn
#  
# diags = pn.Diagnostics()
#  
# lines_dict = {}
#  
# lines_dict['O3_4363A'] = 1.040046e-15
# lines_dict['O3_4959A'] = 1.336816e-14
# lines_dict['O3_5007A'] = 3.974326e-14
# lines_dict['S2_6716A'] = 4.892321e-16
# lines_dict['S2_6731A'] = 2.638651e-16
#  
# print lines_dict['S2_6716A']/lines_dict['S2_6731A']
# 
# S2 = pn.Atom('S', 2)
# 
# S2.getTemDen(int_ratio, tem, den, lev_i1, lev_j1, lev_i2, lev_j2, wave1, wave2, maxError, method, log, start_x, end_x, to_eval, nCut, maxIter)
#  
# tem, den = diags.getCrossTemDen(diag_tem = '[OIII] 4363/5007+', diag_den  = '[SII] 6731/6716',
#                         value_tem = lines_dict['O3_4363A']/(lines_dict['O3_5007A'] + lines_dict['O3_4959A']),
#                         value_den = 0.75)
# 
# print tem
# print den

# df.iloc[2] = series
# 
# print df
# 
# df.iloc[2] = series_inverse
# 
# print df
# 
# df.iloc[2] = series_imcomplete
# 
# print df

# df.iloc[2] = series_imcomplete
# 
# print df

# curve fit [ 2.14292041  0.10294343  1.51744256]  time  0.743386030197
# lmfit [2.1429094932812265, 0.10294902884798313, 1.5174734606419207] 4.23295807838

# 0.109779834747
# 0.0164370536804
# 0.0189270973206
# 
# result for v0.9.1 is
# 
# 20.6462230682
# 1.95246195793
# 1.74573016167
# while result for v0.8.3 is
# 
# 0.400651931763
# 0.0140390396118
# 0.0666930675507

# #!/usr/bin/env python
# 
# """
#   Developed by Hiro Takami (hiro@asiaa.sinica.edu.tw), revised on 3/18/2014.
#   For the public release version I have copied some external files here
#   as classes so that you do not need to download any other files.
# 
#   Run as read_plots.py. The major revisions are summarized below:
# 
#   - Now the figure window is fit in the GUI. 
# 
#   - The GUIs to determine the axes is now organized as the same GUI.
#     When we finish determining the axes, the window is replaced as
#     the one to measure the values in the figure window. 
# 
# """
# 
# ### Orignaly imported ###
# import Tkinter as Tk
# from pylab import *
# from tkFileDialog import askopenfilename,asksaveasfilename
# import glob,pickle,os,tkMessageBox
# from PIL import Image
# ### from tk_window_with_mpl_figure.py ###
# import matplotlib
# matplotlib.use('TkAgg')
# 
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# from matplotlib.figure import Figure
# #import Tkinter as Tk
# #from tkFileDialog import asksaveasfilename
# import sys
# #import mpl_event_handling_sub as event_handling
# 
# ### from mpl_event_handling_sub.py ###
# from pylab import *
# from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
# 
# ##############################################################################################
# ##############################################################################################
# ### Constants ################################################################################
# ##############################################################################################
# ##############################################################################################
# 
# crosshair_is_activated_without_fixed_x_or_y = True  # changed later when activated.
# format_for_values_default                   = "%4.1e"
# 
# # default setting for colors and font
# 
# show_lines_for_recorded_positions=True
# show_dots_for_recorded_positions=True
# 
# color_axes          ='b'
# color_crosshair     ='b'
# color_recorded      ='k'
# linewidth_axes      = 1
# linewidth_crosshair = 1
# linewidth_recorded  = 2
# marker_size_squares = 5
# 
# font=('Hevertica', 12)
# #font=('Times', 12, "bold")
# #font_large=('Times', 20, "bold")
# 
# #color_default         = "systemWindowBody"
# #
# # 3/18/2014
# #
# # It turns out that the above default color does not work with ubuntu,
# # so I have revised the code.
# # 
# root=Tk.Tk()
# color_default=root.cget('bg')
# root.destroy()
# del root
# 
# color_error_entry_box = "#ffCCCC"
# 
# """
#   The classes and modules are described below for the following order:
# 
#   GUIs:
#     Get_values
#     Determine_axis
# 
#   Step 1: determine the lines for the x and y axes
#     Draw_axis
#     Xaxis
#     Yaxis
# 
#   Step 2a: give two ticks to x-axis
#     Xticks
#     Xticks_sub
# 
#   Step 2b: give two ticks to y-axis
#     Yticks
#     Yticks_sub
# 
#   Step 3: give values to ticks
#     input_values_for_ticks
# 
#   Step 4: measure values in the plot
#     Record_positions
#     Show_crosshair
# 
#   Other class and modules:
#     panzoom
#     zoomrect
#     home
#     update_format_for_values
#     draw_again
# 
# """
# 
# ##############################################################################################
# ##############################################################################################
# ### GUIs #####################################################################################
# ##############################################################################################
# ##############################################################################################
# 
# class Get_values:
#   """
#     Widgets and modules to get values. The widgets are described in main_Tk_window.console.
#   """
#   def __init__(self):
# 
#     # Parameters
#     self.constraint_var=Tk.StringVar()
#     self.constraint_var.set("None")
#     self.constraint_last="None"
# 
#     # Menu bar
# 
#     fm=Tk.Frame(main_Tk_window.console, relief=Tk.RAISED, bd=1)
#     #fm.pack(side=Tk.TOP,anchor=Tk.W,fill=Tk.X)
#     fm.pack(fill=Tk.X)
# 
#     mb1=Tk.Menubutton(fm,text='File')
#     mb1.pack(side=Tk.LEFT,anchor=Tk.W)
#     menu1=Tk.Menu(mb1)
#     menu1.add_command(label='Quit',command=main_Tk_window.root.quit)
#     mb1.config(menu=menu1)
# 
#     mb2=Tk.Menubutton(fm,text='Format for values')
#     mb2.pack(side=Tk.LEFT,anchor=Tk.W)
#     menu2=Tk.Menu(mb2)
#     menu2.add_command(label='%.1f',command=self.format_1f)
#     menu2.add_command(label='%.2f',command=self.format_2f)
#     menu2.add_command(label='%.3f',command=self.format_3f)
#     menu2.add_command(label='%.4f',command=self.format_4f)
#     menu2.add_command(label='%.1e',command=self.format_1e)
#     menu2.add_command(label='%.2e',command=self.format_2e)
#     menu2.add_command(label='%.3e',command=self.format_3e)
#     menu2.add_command(label='%.4e',command=self.format_4e)
#     mb2.config(menu=menu2)
# 
#     mb3=Tk.Menubutton(fm,text='Color')
#     mb3.pack(side=Tk.LEFT,anchor=Tk.W)
#     menu3=Tk.Menu(mb3)
#     menu3.add_command(label='blue'   ,command=self.color_b)
#     menu3.add_command(label='green'  ,command=self.color_g)
#     menu3.add_command(label='red'    ,command=self.color_r)
#     menu3.add_command(label='cyan'   ,command=self.color_c)
#     menu3.add_command(label='magenda',command=self.color_m)
#     menu3.add_command(label='yellow' ,command=self.color_y)
#     menu3.add_command(label='black'  ,command=self.color_k)
#     menu3.add_command(label='white'  ,command=self.color_w)
#     mb3.config(menu=menu3)
# 
#     mb4=Tk.Menubutton(fm,text='Line width')
#     mb4.pack(side=Tk.LEFT,anchor=Tk.W)
#     menu4=Tk.Menu(mb4)
#     menu4.add_command(label='1',command=self.linewidth_1)
#     menu4.add_command(label='2',command=self.linewidth_2)
#     menu4.add_command(label='3',command=self.linewidth_3)
#     menu4.add_command(label='4',command=self.linewidth_4)
#     mb4.config(menu=menu4)
# 
#     mb5=Tk.Menubutton(fm,text='Marker')
#     mb5.pack(side=Tk.LEFT,anchor=Tk.W)
#     menu5=Tk.Menu(mb5)
#     menu5.add_command(label='x',command=self.marker_x)
#     menu5.add_command(label='+',command=self.marker_cross)
#     menu5.add_command(label='o',command=self.marker_round)
#     mb5.config(menu=menu5)
# 
#     mb6=Tk.Menubutton(fm,text='Marker size')
#     mb6.pack(side=Tk.LEFT,anchor=Tk.W)
#     menu6=Tk.Menu(mb6)
#     menu6.add_command(label='1',command=self.markersize_1)
#     menu6.add_command(label='3',command=self.markersize_3)
#     menu6.add_command(label='5',command=self.markersize_5)
#     menu6.add_command(label='10',command=self.markersize_10)
#     mb6.config(menu=menu6)
# 
#     # To constrain the movement of the crosshair if you want
#     l1=Tk.Frame(main_Tk_window.console,padx=10,pady=10)
#     l1.pack()
# 
#     Tk.Label(l1,text="Constraint",font=font).pack(side=Tk.LEFT,padx=5)
#     self.optionmenu=Tk.OptionMenu(l1, self.constraint_var, "None","X","Y",command=self.constraint)
#     self.optionmenu.configure(width=10)
#     self.optionmenu.pack(side=Tk.LEFT,padx=5)
#     self.entry=Tk.Entry(l1,width=15,font=font,state=Tk.DISABLED)
#     self.entry.pack(side=Tk.LEFT,padx=5)
#     self.button=Tk.Button(l1,text="Go",command=self.go,state=Tk.DISABLED)
#     self.button.pack(side=Tk.LEFT,padx=5)
# 
#     # table for instructions
#     table=Tk.Frame(main_Tk_window.console,padx=10,pady=10)
#     table.pack()
#     Tk.Label(table,text="Left mouse button"            ,font=font).grid(row=0,column=0,sticky=Tk.W)
#     Tk.Label(table,text="..."                          ,font=font).grid(row=0,column=1)
#     Tk.Label(table,text="record the position"          ,font=font).grid(row=0,column=2,sticky=Tk.W)
#     Tk.Label(table,text="Right mouse button or 'p' key",font=font).grid(row=1,column=0,sticky=Tk.W)
#     Tk.Label(table,text="..."                          ,font=font).grid(row=1,column=1)
#     Tk.Label(table,text="remove the last position"     ,font=font).grid(row=1,column=2,sticky=Tk.W)
# 
#     # buttons to reset or save the record
#     l2=Tk.Frame(main_Tk_window.console,padx=10,pady=10)
#     l2.pack()
#     Tk.Button(l2,text="Save",command=self.save_record).pack(side=Tk.RIGHT,padx=15)
#     Tk.Button(l2,text="Reset the record",command=self.reset_record).pack(side=Tk.RIGHT,padx=15)
# 
#   #-------------------------------
#   def constraint(self,value):
#     if value != self.constraint_last:
#       if value == "None":
#         self.entry.configure(state=Tk.DISABLED)
#         self.button.configure(state=Tk.DISABLED)
# 
#         # added on 10/12/2013
#         crosshair.reset_constraint()
#         global record_positions
#         #record_positions = Record_positions()
#         record_positions.activate()
# 
#       else:
#         self.entry.configure(state=Tk.NORMAL)
#         self.button.configure(state=Tk.NORMAL)
#         #if self.constraint_last == "None":
#           #record_positions.exit()
#           #del record_positions
#           #record_positions.reset_all_positions()
#           #record_positions.deactivate()
# 
#     self.constraint_last=value
# 
#   #-------------------------------
#   def go(self):
#     try:
#       value=float(self.entry.get())
#     except:
#       #print "Input values are wrong"
#       tkMessageBox.showerror("Error","Input the right value in the text box.")
#       return
# 
#     #global record_positions
# 
#     axis=self.constraint_var.get()
#     if axis == "X":
#       crosshair.set_x(value)
#       self.reset_record()
#       #record_positions.reset_all_positions()
#       #record_positions.deactivate()
# 
#     elif axis == "Y":
#       crosshair.set_y(value)
#       self.reset_record()
#       #record_positions.reset_all_positions()
#       #record_positions.deactivate()
# 
#     return
# 
#   #-------------------------------
#   def format_1f(self):
#     update_format_for_values("%.1f")
#   #---
#   def format_2f(self):
#     update_format_for_values("%.2f")
#   #---
#   def format_3f(self):
#     update_format_for_values("%.3f")
#   #---
#   def format_4f(self):
#     update_format_for_values("%.4f")
#   #---
#   def format_1e(self):
#     update_format_for_values("%.1e")
#   #---
#   def format_2e(self):
#     update_format_for_values("%.2e")
#   #---
#   def format_3e(self):
#     update_format_for_values("%.3e")
#   #---
#   def format_4e(self):
#     update_format_for_values("%.4e")
# 
#   #-------------------------------
#   def color_b(self):
#     global record_positions
#     record_positions.color("b")
#   #---
#   def color_g(self):
#     global record_positions
#     record_positions.color("g")
#   #---
#   def color_r(self):
#     global record_positions
#     record_positions.color("r")
#   #---
#   def color_c(self):
#     global record_positions
#     record_positions.color("c")
#   #---
#   def color_m(self):
#     global record_positions
#     record_positions.color("m")
#   #---
#   def color_y(self):
#     global record_positions
#     record_positions.color("y")
#   #---
#   def color_k(self):
#     global record_positions
#     record_positions.color("k")
#   #---
#   def color_w(self):
#     global record_positions
#     record_positions.color("w")
# 
#   #-------------------------------
#   def linewidth_1(self):
#     global record_positions
#     record_positions.linewidth(1)
#   #---
#   def linewidth_2(self):
#     global record_positions
#     record_positions.linewidth(2)
#   #---
#   def linewidth_3(self):
#     global record_positions
#     record_positions.linewidth(3)
#   #---
#   def linewidth_4(self):
#     global record_positions
#     record_positions.linewidth(4)
# 
#   #-------------------------------
#   def marker_round(self):
#     global record_positions
#     record_positions.marker("o")
#   #---
#   def marker_cross(self):
#     global record_positions
#     record_positions.marker("+")
#   #---
#   def marker_x(self):
#     global record_positions
#     record_positions.marker("x")
# 
#   #-------------------------------
#   def markersize_1(self):
#     global record_positions
#     record_positions.markersize(1)
#   #---
#   def markersize_3(self):
#     global record_positions
#     record_positions.markersize(3)
#   #---
#   def markersize_5(self):
#     global record_positions
#     record_positions.markersize(5)
#   #---
#   def markersize_10(self):
#     global record_positions
#     record_positions.markersize(10)
# 
#   #-------------------------------
#   def reset_record(self):
#     global record_positions
#     record_positions.reset_all_positions()
#   #-------------------------------
#   def save_record(self):
#     global record_positions
#     record_positions.save()
# 
# ##############################################################################################
# 
# class Determine_axis:
#   """
#     Widgets and modules to determine the xy axis. The widgets are described in main_Tk_window.console.
#   """
#   def __init__(self):
# 
#     # Parameters for display
#     width_entry           = 10
#     width_message         = 35
#     bg_color_message      = "white"
#     self.fg_color_labels  = "#A0A0A0"
# 
#     # Parameters for analysis
#     self.axis_working_now ='X'  # either 'X' or 'Y'
#     self.status           = 0   # 0 ... setting the axis ; 1 ... setting the ticks and values
# 
#     self.scale_is_linear=Tk.BooleanVar()
#     self.scale_is_linear.set(True)
#     self.entry1_var=Tk.StringVar()
#     self.entry1_var.set("")
#     self.entry2_var=Tk.StringVar()
#     self.entry2_var.set("")
# 
#     # Set the overall window
# 
#     #self.window=Tk.Frame(main_Tk_window.console,relief=Tk.GROOVE,bd=5,padx=5,pady=5)
#     #self.window.place(relx=0.2,rely=0.2,relheight=0.6,relwidth=0.6)
#     self.window=Tk.Frame(main_Tk_window.console,padx=5,pady=5)
#     self.window.pack()
# 
#     l1=Tk.Frame(self.window,pady=0)
#     l1.pack()
#     self.label_message=Tk.Label(l1,font=font,width=width_message,
#                                 text="Set X-axis",
#                                 bg=bg_color_message,
#                                 relief=Tk.RIDGE,pady=3
#                                 )
#     self.label_message.pack(pady=5,padx=10,side=Tk.LEFT)
#     self.button_done=Tk.Button(l1,text='Done',command=self.done)
#     self.button_done.pack(side=Tk.RIGHT)
# 
# 
#     l_table=Tk.Frame(self.window,pady=0)
#     l_table.pack()
# 
#     # Set values
# 
#     self.label_values=Tk.Label(l_table,text='Values',fg=self.fg_color_labels)
#     self.label_values.grid(row=0,column=1)
#     entries=Tk.Frame(l_table,padx=10)
#     entries.grid(row=1,column=1)
#     self.entry1=Tk.Entry(entries,width=width_entry,font=font,state=Tk.DISABLED,
#                          textvariable=self.entry1_var)
#     self.entry1.pack(side=Tk.LEFT)
#     self.entry2=Tk.Entry(entries,width=width_entry,font=font,state=Tk.DISABLED,
#                          textvariable=self.entry2_var)
#     self.entry2.pack(side=Tk.LEFT)
# 
#     # Scale
# 
#     self.label_scale=Tk.Label(l_table,text='Scale',fg=self.fg_color_labels)
#     self.label_scale.grid(row=0,column=2)
# 
#     radiobuttons=Tk.Frame(l_table,padx=10)
#     radiobuttons.grid(row=1,column=2)
#     self.radiobutton1=Tk.Radiobutton(radiobuttons,text = 'Linear', 
#                                      variable = self.scale_is_linear, value = True,
#                                      state=Tk.DISABLED)
#     self.radiobutton1.pack(side=Tk.LEFT)
#     self.radiobutton2=Tk.Radiobutton(radiobuttons,text = 'Log', 
#                                      variable = self.scale_is_linear, value = False,
#                                      state=Tk.DISABLED)
#     self.radiobutton2.pack(side=Tk.LEFT)
# 
#     # Draw x-axis
#     self.draw_xaxis_params=Draw_axis("x")
# 
#   #-------------------------------
#   def done(self):
#     """
#       Executed when the "done" button is pressed.
#     """
# 
#     if self.axis_working_now=='X':
#       # if we finish setting the x axis
#       if self.status == 0:
# 
#         # finish drawing the xaxis if approxriately set
#         #self.draw_xaxis_params.finish()
# 
#         try:
#           self.draw_xaxis_params.axis.get_fitting_params()
#         except:
#           tkMessageBox.showerror("Error","Correctly define the axis line in the figure.")
#           return
# 
#         self.draw_xaxis_params.axis.hide_squares()
#         event_handling.reset()
#         #---
# 
#         self.label_message.configure(text="Set ticks and values for X-axis")
#         xticks.determine()
#         self.status = 1
# 
#         self.label_values.configure(fg="black")
#         self.label_scale.configure(fg="black")
#         self.entry1.configure(state=Tk.NORMAL)
#         self.entry2.configure(state=Tk.NORMAL)
#         self.radiobutton1.configure(state=Tk.NORMAL)
#         self.radiobutton2.configure(state=Tk.NORMAL)
# 
#       # if we finish the ticks and values for the x axis
#       else:
# 
#         if xticks.ticking_completed == False:
#           tkMessageBox.showerror("Error","Set ticks in the axis.")
#           return
# 
#         tmp=self.entry1_var.get()
# 
#         try:
#           self.x1_graph=float(self.entry1_var.get())
#           self.entry1.config(bg=color_default)
#         except:
#           self.entry1.config(bg=color_error_entry_box)
#           tkMessageBox.showerror("Error","Correctly input the red text box.")
#           return
# 
#         try:
#           self.x2_graph=float(self.entry2_var.get())
#           self.entry2.config(bg=color_default)
#         except:
#           self.entry2.config(bg=color_error_entry_box)
#           tkMessageBox.showerror("Error","Correctly input the red text box.")
#           return
# 
#         if self.scale_is_linear.get():
#           self.xscale="linear"
#         else:
#           self.xscale="log"
# 
#         self.label_message.configure(text="Set Y-axis")
#         self.axis_working_now='Y'
#         self.status          = 0
#         self.draw_yaxis_params=Draw_axis("y")
# 
#         self.label_values.configure(fg=self.fg_color_labels)
#         self.label_scale.configure(fg=self.fg_color_labels)
#         self.entry1_var.set("")
#         self.entry2_var.set("")
#         self.entry1.configure(state=Tk.DISABLED)
#         self.entry2.configure(state=Tk.DISABLED)
#         self.radiobutton1.configure(state=Tk.DISABLED)
#         self.radiobutton2.configure(state=Tk.DISABLED)
# 
#     else:
# 
#       # if we finish setting the y axis
#       if self.status == 0:
# 
#         # finish drawing the yaxis if approxriately set
#         #self.draw_yaxis_params.finish()
#         try:
#           self.draw_yaxis_params.axis.get_fitting_params()
#         except:
#           tkMessageBox.showerror("Error","Correctly define the axis line in the figure.")
#           return
# 
#         self.draw_yaxis_params.axis.hide_squares()
#         event_handling.reset()
#         #---
# 
#         self.label_message.configure(text="Set ticks and values for Y-axis")
#         yticks.determine()
#         self.status = 1
# 
#         self.label_values.configure(fg="black")
#         self.label_scale.configure(fg="black")
#         self.entry1.configure(state=Tk.NORMAL)
#         self.entry2.configure(state=Tk.NORMAL)
#         self.radiobutton1.configure(state=Tk.NORMAL)
#         self.radiobutton2.configure(state=Tk.NORMAL)
# 
#       # if we finish the ticks and values for the y axis
#       else:
# 
#         if yticks.ticking_completed == False:
#           tkMessageBox.showerror("Error","Set ticks in the axis.")
#           return
# 
#         try:
#           self.y1_graph=float(self.entry1_var.get())
#           self.entry1.config(bg=color_default)
#         except:
#           self.entry1.config(bg=color_error_entry_box)
#           tkMessageBox.showerror("Error","Correctly input the red text box.")
#           return
# 
#         try:
#           self.y2_graph=float(self.entry2_var.get())
#           self.entry2.config(bg=color_default)
#         except:
#           self.entry2.config(bg=color_error_entry_box)
#           tkMessageBox.showerror("Error","Correctly input the red text box.")
#           return
# 
#         if self.scale_is_linear.get():
#           self.yscale="linear"
#         else:
#           self.yscale="log"
# 
#         event_handling.reset()
#         input_values_for_ticks(self.x1_graph,self.x2_graph,
#                                self.y1_graph,self.y2_graph,
#                                self.xscale,self.yscale)
#         #self.window.place_forget()
#         self.window.pack_forget()
# 
#         # show the crosshair for the next step
#         global crosshair,record_positions
#         #crosshair        = Show_crosshair()
#         #record_positions = Record_positions()
#         crosshair.activate()
#         record_positions.activate()
# 
#         Get_values()
# 
# 
# ##############################################################################################
# ##############################################################################################
# ### Step 1: determine the lines for the x and y axes  ########################################
# ##############################################################################################
# ##############################################################################################
# 
# """
#   The classes below are for determining the axis. These will be used as follows:
# 
#    (main program)
# 
#     xaxis=Xaxis()
#     yaxis=Yaxis()
# 
#    (in the GUI, or commandlines)
# 
#     Draw_axis("x")
#     Draw_axis("y")
# 
# """
# 
# class Draw_axis:
# 
#   """
#   See the above instruction about how the class will be called. Once it is called:
#   1st click ... determine the left or bottom side of the axis
#   2nd click ... determine the right or top side of the axis
#   click edges of the axis ... for minor ajustment
#   press "n" ... completed > now we use a button in the GUI.
# 
#   ### Modules ###
#   
#   We do not have to call any of these modules from outside.
# 
#   self.key_pressed            ... just to be able to press n to notice we have determined the axis.
#   self.determine_1st_position (button_press)
#   move_position               (motion_notify)
#   determine_position          (button_press)
#   pick_edge                   (picker)
#   finish
#   
# 
#   ### Accessible internal parameters (i.e., self.****) ###
# 
#   There are minor so that we do not have to browse them from outside at all.
# 
#   axis                         ... either the parameter xaxis or yaxis is set.
#   key_n_activated (True/False) ... n key is used after we finish determining x or y axis
# 
#   """
#   def __init__(self,axis):  # "axis" is either "x" or "y" 
# 
#  
#     self.key_n_activated = False
# 
#     if   axis=="x": 
#       global xaxis
#       self.axis=xaxis
#     elif axis=="y":
#       global yaxis
#       self.axis=yaxis
#     else:
#       return
#     event_handling.reset()
# 
#     #print 'Click and determine %s-axis. Press "n" if completed.' % (axis)
# 
#     event_handling.activate("button_press_event",self.determine_1st_position)
#     event_handling.activate("key_press_event"   ,self.key_press)
# 
#     return
# 
#   #---
#   def determine_1st_position(self,event):
#     """
#     Draw_axis.determine_1st_position
#     """
#     if not event.inaxes: return
#     elif event.button != 1: return
# 
#     self.axis.update_position_of_edge(0,event.xdata,event.ydata)
#     self.axis.update_position_of_edge(1,event.xdata,event.ydata)
#     self.axis.show_squares()
#     draw_again()
#     event_handling.deactivate("button_press_event")
# 
#     self.position_no=1 # will move and determine the 2nd position
#     event_handling.activate("motion_notify_event",self.move_position)
#     event_handling.activate("button_press_event" ,self.determine_position)
# 
#     # Do not activate the "n" key when the motion notify event is active
#     self.key_n_activated = False
# 
#     return
# 
#   #---
#   def move_position(self,event):
#     """
#     Draw_axis.move_position
#     """
#     if not event.inaxes: return
#     self.axis.update_position_of_edge(self.position_no,event.xdata,event.ydata)
# 
#     return
# 
#   #---
#   def determine_position(self,event):
#     """
#     Draw_axis.determine_position
#     """
#     if not event.inaxes: return
#     elif event.button != 1: return
# 
#     self.axis.update_position_of_edge(self.position_no,event.xdata,event.ydata)
#     event_handling.deactivate("button_press_event")
#     event_handling.deactivate("motion_notify_event")
# 
#     self.key_n_activated = True
# 
#     # still be able to pick two edges of the axis for minor adjustment
#     #if not connect_params.has_key("pick_event"):
#     if event_handling.status.count("pick_event") == 0:
#       self.axis.activate_picker()  # this command is just for displaying small squares at two ends
#       event_handling.activate("pick_event",self.pick_edge)
# 
#     return
# 
#   #---
#   def pick_edge(self,event):
#     """
#     Draw_axis.pick_edge
#     """
#     self.position_no=event.ind
#     event_handling.deactivate("pick_event")
#     self.axis.deactivate_picker()  # this command is just for erasing small squares at two ends
# 
#     event_handling.activate("motion_notify_event",self.move_position)
#     event_handling.activate("button_press_event",self.determine_position)
# 
#     # Do not activate the "n" key when the motion notify event is active
#     self.key_n_activated = False
# 
#   #------------------------------------------------------------------------
#   def key_press(self,event):
#     """
#     Draw_axis.key_press
#     """
# 
#     ### if "n" is pressed ##
#     if self.key_n_activated and event.key=="n":
#       self.finish()
# 
#     return
# 
#   #------------------------------------------------------------------------
#   def finish(self):
# 
#     self.axis.hide_squares()
#     event_handling.reset()
#     self.axis.get_fitting_params()
# 
#     return
# 
# #------------------------------------------------------------------------
# 
# class Xaxis:
#   """
#     Class for the x axis. Called in Read_graph as follows:
#       xaxis=Xaxis()
# 
#     ### Input parameters ###
# 
#     None
# 
#     ### Modules ###
#     self.update_position_of_edge
#        ... as described. Parameters are:
#            i   ... 0 for the left, 1 for the right edge.
#            x,y ... matplotlib coordinate. event.xdata and event.ydata will be
#                    handed.
# 
#     self.get_fitting_params
#        ... get self.a and self.b described below. Will be used when the two edges
#            of the axis line are determined. Then we will be able to
#            describe the axis line in the matplotlib coordinate using these parameters.
# 
#     self.add_parameters_to_determine_x
#        ... get self.d and self.e to convert the matplotlib coordinate to the graph
#            coordinate. Will be used when two ticks are determined along the x axis.
#            The input parameters are:
#              x_mpl,y_mpl     ... matplotlib coodinate (1-D array with two parameters)
#              x_graph,y_graph ... graph coordinate     (1-D array with two parameters)
#              ay              ... inclination for y axis
# 
#     self.get_x
#        ... get the x value in the graph coordinate, with a given (x,y)
#            in the matplotlib coordinate.
# 
#     self.show_squares      ... as described. See each module, and you can immediately
#     self.hide_squares          understand what these are doing.
#     self.activate_picker
#     self.deactivate_picker
#     self.hide
#     self.has
# 
#     ### Output parameters ###
#     self.x_mpl_edge,self.y_mpl_edge ... edges for the line being drawn
#     self.pl_line,self.pl_squares    ... parameters for the plot command
#                                           for the lines and two squares at the end, respectively.
#     self.a,self.b                   ... parameters to decribe the x axis line
#                                           in the matplotlib format
#                                            (y_mpl_edge = a * x_mpl_edge + b)
# 
#   """
# 
#   def __init__(self,params=True): # identical with the class Yaxis
# 
#     # prepare for the plots if the parameters are not defined above.
#     if params == True:
#       self.x_mpl_edge  =[0,0]
#       self.y_mpl_edge  =[0,0]
#       self.pl_line,    =ax.plot([0,0],[0,0],color_axes+'-',linewidth=linewidth_axes,visible=False)
#       self.pl_squares, =ax.plot([0,0],[0,0],color_axes+'s',
#                              markersize=marker_size_squares,visible=False)
# 
#     # otherwise, load the data
#     else:
#       self.a=params['a']
#       self.b=params['b']
#       self.d=params['d']
#       self.e=params['e']
#       self.a_for_yaxis=params['a_for_yaxis']
#       self.x_mpl_edge=params['x_mpl_edge']
#       self.scale     =params['scale']
# 
#   #--- (identical with the class Yaxis) -----------------------------------
#   def update_position_of_edge(self,i,x,y):
#     self.x_mpl_edge[i],self.y_mpl_edge[i]=x,y
# 
#     self.pl_line.set_xdata(self.x_mpl_edge)
#     self.pl_line.set_ydata(self.y_mpl_edge)
#     self.pl_line.set_visible(True)
# 
#     self.pl_squares.set_xdata(self.x_mpl_edge)
#     self.pl_squares.set_ydata(self.y_mpl_edge)
# 
#     draw_again()
# 
#     return
# 
#   #------------------------------------------------------------------------
#   def get_fitting_params(self):
#     """
#       y_mpl_edge = a * x_mpl_edge + b
#     """
#     self.a = (self.y_mpl_edge[1]-self.y_mpl_edge[0])/(self.x_mpl_edge[1]-self.x_mpl_edge[0])
#     self.b =  self.y_mpl_edge[0] - self.a * self.x_mpl_edge[0]
# 
#   #------------------------------------------------------------------------
#   def add_parameters_to_determine_x(self,x_mpl,x_graph,y_mpl,y_graph,ay,
#                                     scale # either "linear" or "log"
#                                     ):
#     if scale != "linear" and scale != "log":
#       print '!!! Error !!! - the scale must be either "linear" or "log".'
# 
#     else: 
#       self.scale       = scale
#       self.a_for_yaxis = ay
# 
#       # x_graph = d * x_mpl + e
# 
#       if self.scale == "linear":
#         self.d = (x_graph[1]-x_graph[0])/(x_mpl[1]-x_mpl[0])
#         self.e =  x_graph[0] - self.d * x_mpl[0]
#       else:
#         self.d = (log10(x_graph[1])-log10(x_graph[0]))/(x_mpl[1]-x_mpl[0])
#         self.e =  log10(x_graph[0]) - self.d * x_mpl[0]
# 
#   #------------------------------------------------------------------------
#   def get_x(self,x_mpl,y_mpl):
#     x=self.d * (x_mpl - self.a_for_yaxis * (y_mpl - self.b)) \
#                     / (1-self.a*self.a_for_yaxis) + self.e
# 
#     if self.scale == "linear":
#       return x
#     else:
#       return 10**x
#   #------------------------------------------------------------------------
#   def params_for_save(self):
#     params={}
#     params['a']=self.a
#     params['b']=self.b
#     params['d']=self.d
#     params['e']=self.e
#     params['a_for_yaxis']=self.a_for_yaxis
#     params['x_mpl_edge']=self.x_mpl_edge
#     params['scale']     =self.scale
# 
#     return params
# 
#   #--- (identical with the class Yaxis) -----------------------------------
#   def show_squares(self):
#     if self.has('pl_squares'):
#       self.pl_squares.set_visible(True)
# 
#   #--- (identical with the class Yaxis) -----------------------------------
#   def hide_squares(self):
#     if self.has('pl_squares'):
#       self.pl_squares.set_visible(False)
# 
#   #--- (identical with the class Yaxis) -----------------------------------
#   def activate_picker(self):
#     if self.has('pl_squares'):
#       self.pl_squares.set_picker(5)
# 
#   #--- (identical with the class Yaxis) -----------------------------------
#   def deactivate_picker(self):
#     if self.has('pl_squares'):
#       self.pl_squares.set_picker(None)
# 
#   #--- (identical with the class Yaxis) -----------------------------------
#   def hide(self):
#     if self.has('pl_line'):
#       self.pl_line.set_visible(False)
#     if self.has('pl_squares'):
#       self.pl_squares.set_visible(False)
# 
#   #--- (identical with the class Yaxis) -----------------------------------
#   def has(self,key): 
#     return dir(self).count(key)
# 
# #------------------------------------------------------------------------
# 
# class Yaxis:
#   """
#     Class for the y axis. Called in Read_graph as follows:
# 
#       yaxis=Yaxis()
# 
#     Almost identical with Xaxis, so see the comments for the other class for details.
#   """
# 
#   def __init__(self,params=True): # identical with the class Xaxis
# 
#     # prepare for the plots if the parameters are not defined above.
#     if params == True:
#       self.x_mpl_edge  =[0,0]
#       self.y_mpl_edge  =[0,0]
#       self.pl_line,    =ax.plot([0,0],[0,0],color_axes+'-',linewidth=linewidth_axes,visible=False)
#       self.pl_squares, =ax.plot([0,0],[0,0],color_axes+'s',
#                              markersize=marker_size_squares,visible=False)
# 
#     # otherwise, load the data
#     else:
#       self.a=params['a']
#       self.b=params['b']
#       self.d=params['d']
#       self.e=params['e']
#       self.a_for_xaxis=params['a_for_xaxis']
#       self.y_mpl_edge=params['y_mpl_edge']
#       self.scale     =params['scale']
# 
#   #--- (identical with the class Xaxis) -----------------------------------
#   def update_position_of_edge(self,i,x,y):
#     self.x_mpl_edge[i],self.y_mpl_edge[i]=x,y
# 
#     self.pl_line.set_xdata(self.x_mpl_edge)
#     self.pl_line.set_ydata(self.y_mpl_edge)
#     self.pl_line.set_visible(True)
# 
#     self.pl_squares.set_xdata(self.x_mpl_edge)
#     self.pl_squares.set_ydata(self.y_mpl_edge)
# 
#     draw_again()
# 
#     return
# 
#   #------------------------------------------------------------------------
#   def get_fitting_params(self):
#     """
#       x_mpl_edge = a * y_mpl_edge + b
#     """
#     self.a = (self.x_mpl_edge[1]- self.x_mpl_edge[0])/(self.y_mpl_edge[1]-self.y_mpl_edge[0])
#     self.b =  self.x_mpl_edge[0] - self.a * self.y_mpl_edge[0]
# 
#   #------------------------------------------------------------------------
#   def add_parameters_to_determine_y(self,x_mpl,x_graph,y_mpl,y_graph,ax,
#                                     scale # either "linear" or "log"
#                                     ):
#     if scale != "linear" and scale != "log":
#       print '!!! Error !!! - the scale must be either "linear" or "log".'
# 
#     else: 
#       self.scale       = scale
#       self.a_for_xaxis=ax
# 
#       # y_graph = d * y_mpl + e
#       if self.scale == "linear":
#         self.d = (y_graph[1]-y_graph[0])/(y_mpl[1]-y_mpl[0])
#         self.e =  y_graph[0] - self.d * y_mpl[0]
#       else:
#         self.d = (log10(y_graph[1])-log10(y_graph[0]))/(y_mpl[1]-y_mpl[0])
#         self.e =  log10(y_graph[0]) - self.d * y_mpl[0]
# 
#   #------------------------------------------------------------------------
#   def get_y(self,x_mpl,y_mpl):
#     y = self.d * (y_mpl - self.a_for_xaxis * (x_mpl - self.b)) \
#                   / (1-self.a*self.a_for_xaxis) + self.e
# 
#     if self.scale == "linear":
#       return y
#     else:
#       return 10**y
# 
#   #------------------------------------------------------------------------
#   def params_for_save(self):
#     params={}
#     params['a']=self.a
#     params['b']=self.b
#     params['d']=self.d
#     params['e']=self.e
#     params['a_for_xaxis']=self.a_for_xaxis
#     params['y_mpl_edge']=self.y_mpl_edge
#     params['scale']     =self.scale
# 
#     return params
# 
#   #--- (identical with the class Xaxis) -----------------------------------
#   def show_squares(self):
#     if self.has('pl_squares'):
#       self.pl_squares.set_visible(True)
# 
#   #--- (identical with the class Xaxis) -----------------------------------
#   def hide_squares(self):
#     if self.has('pl_squares'):
#       self.pl_squares.set_visible(False)
# 
#   #--- (identical with the class Xaxis) -----------------------------------
#   def activate_picker(self):
#     if self.has('pl_squares'):
#       self.pl_squares.set_picker(5)
# 
#   #--- (identical with the class Xaxis) -----------------------------------
#   def deactivate_picker(self):
#     if self.has('pl_squares'):
#       self.pl_squares.set_picker(None)
# 
#   #--- (identical with the class Xaxis) -----------------------------------
#   def hide(self):
#     if self.has('pl_line'):
#       self.pl_line.set_visible(False)
#     if self.has('pl_squares'):
#       self.pl_squares.set_visible(False)
# 
#   #--- (identical with the class Xaxis) -----------------------------------
#   def has(self,key): 
#     return dir(self).count(key)
# 
# ##############################################################################################
# ##############################################################################################
# ### Step 2a: give two ticks to x-axis ########################################################
# ##############################################################################################
# ##############################################################################################
# 
# class Xticks:
#   """
#     To determine the axis. Use as follows:-
# 
#      (main program)
# 
#       xticks=Xticks()
# 
#      (in the GUI or a commandline)
# 
#       xticks.determine()
# 
#     All
# 
#   """
#   def __init__(self):
# 
#     self.ticking_completed=False
# 
#   #---
#   def determine(self):
# 
#     self.tick=[Xticks_sub(xaxis.a,xaxis.b)
#                ,Xticks_sub(xaxis.a,xaxis.b)]
# 
#     self.tick_no=0
#     event_handling.activate("motion_notify_event",self.motion_notify)
#     event_handling.activate("button_press_event" ,self.button_pressed)
#     event_handling.activate("key_press_event"    ,self.key_pressed)
# 
#   #---
#   def motion_notify(self,event):
#     """
#     Xticks.motion_notify
#     """
#     # move the tick
# 
#     if event.xdata > xaxis.x_mpl_edge[0] and event.xdata < xaxis.x_mpl_edge[1]:
#       self.tick[self.tick_no].update(event.xdata)
#     else:
#       self.tick[self.tick_no].hide()
# 
#   #---
#   def button_pressed(self,event):
#     """
#     Xticks.button_pressed
#     """
#     #  read the number at the clicked position
# 
#     if event.xdata < xaxis.x_mpl_edge[0] or event.xdata > xaxis.x_mpl_edge[1]: return
# 
#     self.tick[self.tick_no].update(event.xdata)
#     self.tick[self.tick_no].mpl=event.xdata*1
# 
#     if self.tick_no==0:
#       self.tick_no=1
#       event_handling.activate("motion_notify_event",self.motion_notify)
#       event_handling.activate("button_press_event",self.button_pressed)
#     else:
#       event_handling.reset()
#       del self.tick_no
#       self.ticking_completed=True
# 
#     return
# 
#   #------------------------------------------------------------------------
#   def key_pressed(self,event):
#     """
#     Xticks.key_pressed
#     """
#     # reset the record and start from the beginning, if "r" is pressed.
# 
#     if event.key=="r":
# 
#       self.tick[0].hide()
#       self.tick[1].hide()
#       del self.tick
#       event_handling.reset()
#       self.determine()
#   
# #------------------------------------------------------------------------
# #------------------------------------------------------------------------
# 
# class Xticks_sub:
#   """
#     Subroutine for ticks on the x axis. Called in Read_graph as:
# 
#       self.xtick=[Xtick(xaxis.a,xaxis.b,half_length_of_xtick)
#                  ,Xtick(xaxis.a,xaxis.b,half_length_of_xtick)]
# 
#     a,b are the parameters to describe the line for the x-axis.
#     If we would like to update the position, call, e.g.,
# 
#       self.xtick[0](x)
# 
#     It also has the show/hide modules. (Just change the set_visible parameter.)
# 
#   """
#   def __init__(self,a,b):
#     self.a,self.b = a,b
#     self.pl,=ax.plot([0,0],[0,0],color_axes+'-',linewidth=linewidth_axes)
#     self.hide()
# 
#   #------------------------------------------------------------------------
#   def update(self,x):
#     self.pl.set_xdata([x,x])
#     y0=self.a * x + self.b
#     self.pl.set_ydata([y0-half_length_of_xtick,y0+half_length_of_xtick])
#     self.show()
# 
#   #--- (identical with the class Ytick) -----------------------------------
#   def show(self):
#     self.pl.set_visible(True)
#     draw_again()
# 
#   #--- (identical with the class Ytick) -----------------------------------
#   def hide(self):
#     self.pl.set_visible(False)
#     draw_again()
# 
# 
# ##############################################################################################
# ##############################################################################################
# ### Step 2b: give two ticks to y-axis ########################################################
# ##############################################################################################
# ##############################################################################################
# 
# 
# class Yticks:
#   def __init__(self):
# 
#     self.ticking_completed=False
# 
#   #---
#   def determine(self):
# 
#     self.tick=[Ytick_sub(yaxis.a,yaxis.b)
#                ,Ytick_sub(yaxis.a,yaxis.b)]
#     draw_again()
# 
#     self.tick_no=0
#     event_handling.activate("motion_notify_event",self.motion_notify)
#     event_handling.activate("button_press_event" ,self.button_pressed)
#     event_handling.activate("key_press_event"    ,self.key_pressed)
# 
#   #---
#   def motion_notify(self,event):
#     """
#     Yticks.motion_notify
#     """
#     #  move the tick
# 
#     if event.ydata > yaxis.y_mpl_edge[0] and event.ydata < yaxis.y_mpl_edge[1]:
#       self.tick[self.tick_no].update(event.ydata)
#     else:
#       self.tick[self.tick_no].hide()
# 
#   #---
#   def button_pressed(self,event):
#     """
#     Yticks.button_pressed
#     """
#     # read the number at the clicked position
# 
#     if event.ydata < yaxis.y_mpl_edge[0] or event.ydata > yaxis.y_mpl_edge[1]: return
# 
#     self.tick[self.tick_no].update(event.ydata)
#     self.tick[self.tick_no].mpl=event.ydata*1
# 
#     if self.tick_no==0:
#       self.tick_no=1
#       event_handling.activate("motion_notify_event",self.motion_notify)
#       event_handling.activate("button_press_event",self.button_pressed)
#     else:
#       event_handling.reset()
#       del self.tick_no
#       self.ticking_completed=True
# 
#     return
# 
#  #------------------------------------------------------------------------
#   def key_pressed(self,event):
#     """
#     Yticks.key_pressed
#     """
#     # reset the record and start from the beginning, if "r" is pressed.
# 
#     if event.key=="r":
# 
#       self.tick[0].hide()
#       self.tick[1].hide()
#       del self.tick
#       event_handling.reset()
#       self.determine()
# 
# #------------------------------------------------------------------------
# #------------------------------------------------------------------------
# 
# class Ytick_sub:
#   """
#     Subroutine for ticks on the y axis. See comments for Xtick for details.
#   """
#   def __init__(self,a,b):
#     self.a,self.b = a,b
#     self.pl,=ax.plot([0,0],[0,0],color_axes+'-',linewidth=linewidth_axes)
#     self.hide()
# 
#   #------------------------------------------------------------------------
#   def update(self,y):
#     self.pl.set_ydata([y,y])
#     x0=self.a * y + self.b
#     self.pl.set_xdata([x0-half_length_of_ytick,x0+half_length_of_ytick])
#     self.show()
# 
#   #--- (identical with the class Xtick) -----------------------------------
#   def show(self):
#     self.pl.set_visible(True)
#     draw_again()
# 
#   #--- (identical with the class Xtick) -----------------------------------
#   def hide(self):
#     self.pl.set_visible(False)
#     draw_again()
# 
# ##############################################################################################
# ##############################################################################################
# ### Step 3: give values to ticks, save the conversion factor in the pickle format ############
# ##############################################################################################
# ##############################################################################################
# 
# def input_values_for_ticks(x1_graph,x2_graph,y1_graph,y2_graph,
#                            xscale='linear',yscale='linear' # either 'linear' or 'log'
#                            ):
#   """
#     As described above.
#   """
# 
#   event_handling.reset()
# 
#   x_mpl=[xticks.tick[0].mpl,xticks.tick[1].mpl]
#   y_mpl=[yticks.tick[0].mpl,yticks.tick[1].mpl]
# 
#   xaxis.add_parameters_to_determine_x(x_mpl,[x1_graph,x2_graph],
#                                       y_mpl,[y1_graph,y2_graph],
#                                       yaxis.a,xscale)
# 
#   yaxis.add_parameters_to_determine_y(x_mpl,[x1_graph,x2_graph],
#                                       y_mpl,[y1_graph,y2_graph],
#                                       xaxis.a,yscale)
# 
#   # hide the axes and ticks
#   xaxis.hide()
#   yaxis.hide()
#   xticks.tick[0].hide()
#   xticks.tick[1].hide()
#   yticks.tick[0].hide()
#   yticks.tick[1].hide()
# 
#   # save the fitting
#   params={}
#   params['xaxis']=xaxis.params_for_save()
#   params['yaxis']=yaxis.params_for_save()
# 
#   db=open(filename_for_axes,'wb')
#   pickle.dump(params,db)
#   db.close()
# 
#   return
# 
# ############################################################################################
# ############################################################################################
# ### Step 4: record positions ###############################################################
# ############################################################################################
# ############################################################################################
# 
# class Record_positions:
#   """
#     As described. Use as:
# 
#       a=Record_positions()
#       a.save()
# 
#     Before save:
#       Left click  ... record the new positions 
#       Right click ... remove the last position.
# 
# 
#     ### Input parameters ###
# 
#     None. Just Steps 1-3 should be completed before executing this program.
# 
#     ### Modules ###
# 
#     button_pressed
#     key_pressed
#     remove_last_position
#     update_plot
#     save
# 
#     color
#     linewidth
#     marker
#     markersize
#     
#     ### Accessible internal parameters (self.****) ###
# 
#     positions_graph ... array for positions with graph and matplotlib coordinates.
#     positions_mpl
#  
#     pl_lines
#     pl_dots
# 
#     filename_for_positions
# 
#   """
#   def __init__(self):
# 
#     global ax
# 
#     if not crosshair_is_activated_without_fixed_x_or_y:
#       #print "The recording function is not activated because of your constraint on either x or y."
#       tkMessageBox.showerror("Error","The recording function is not activated because of your constraint on either x or y.")
#       return
# 
#     self.positions_graph =[]
#     self.positions_mpl   =[]
# 
#     # initialize the lines and dots for the recorded points
#     self.pl_lines,=ax.plot([0],[0],color_recorded+"-",
#                                  linewidth=linewidth_recorded,
#                                  visible=show_lines_for_recorded_positions)
#     self.pl_dots,=ax.plot([0],[0],color_recorded+"x",
#                                  markeredgewidth=linewidth_recorded,
#                                  visible=show_dots_for_recorded_positions)
# 
#     self.status=None
#     self.deactivate()
# 
#   #---
#   def activate(self):
# 
#     if self.status != "active":
# 
#       event_handling.activate("button_press_event",self.button_pressed)
#       event_handling.activate("key_press_event",self.key_pressed)
# 
#       self.status="active"
# 
#     return
# 
#   #---
#   def deactivate(self):
# 
#     if self.status != "not active":  
# 
#       event_handling.deactivate("button_press_event")
#       event_handling.deactivate("key_press_event")
# 
#     self.status="not active"
# 
#     return
# 
#  #---
#   def button_pressed(self,event):
#     """
#     Record_positions.button_pressed
#     """
# 
#     if not event.inaxes: return
#     elif event.xdata < xaxis.x_mpl_edge[0] or event.xdata > xaxis.x_mpl_edge[1]: return
#     elif event.ydata < yaxis.y_mpl_edge[0] or event.ydata > yaxis.y_mpl_edge[1]: return
# 
#     # Left click  ... record the new positions 
#     if event.button == 1:
# 
#       x=xaxis.get_x(event.xdata,event.ydata)
#       y=yaxis.get_y(event.xdata,event.ydata)
#       self.positions_graph.append([x,y])
#       self.positions_mpl.append([event.xdata,event.ydata])
#       print format_for_print % (x,y)
#       self.update_plot()
# 
#     # Right click ... remove the last position.
#     if event.button == 3:
#       self.remove_last_position()
# 
#     return
# 
#   #---
#   def key_pressed(self,event):
#     """
#     Record_positions.key_pressed
#     """
#     if   event.key == "r": self.reset_all_positions()
#     elif event.key == "p": self.remove_last_position()
#     elif event.key == "w": self.save()
# 
#     return
# 
#   #---
#   def reset_all_positions(self):
#     self.positions_graph=[]
#     self.positions_mpl=[]
#     self.update_plot()
# 
#     return
# 
#   #---
#   def remove_last_position(self):
#     if self.positions_graph != []:
#       self.positions_graph.pop()
#       self.positions_mpl.pop()
#       self.update_plot()
# 
#     return
# 
#   #---
#   def update_plot(self):
#     # if there is no data point
#     if self.positions_mpl == []:
#       self.pl_lines.set_visible(False)
#       self.pl_dots.set_visible(False)
#       #draw_again()
#       #main_Tk_window.canvas_fig.draw()
# 
#     # if there is only data point
#     elif len(self.positions_mpl) == 1:
#       self.pl_lines.set_visible(False)
#       if show_dots_for_recorded_positions:
#         self.pl_dots.set_xdata([self.positions_mpl[0][0],self.positions_mpl[0][0]])
#         self.pl_dots.set_ydata([self.positions_mpl[0][1],self.positions_mpl[0][1]])
#         self.pl_dots.set_visible(True)
#         #draw_again()
#         #main_Tk_window.canvas_fig.draw()
# 
#     # if there are two more more data points, and if you would like to display them
#     elif show_lines_for_recorded_positions or show_dots_for_recorded_positions: 
#       positions=array(self.positions_mpl)
# 
#       if show_lines_for_recorded_positions:
#         self.pl_lines.set_xdata(positions[:,0])
#         self.pl_lines.set_ydata(positions[:,1])
#         self.pl_lines.set_visible(True)
#       else:
#         self.pl_lines.set_visible(False)
# 
#       if show_dots_for_recorded_positions:
#         self.pl_dots.set_xdata(positions[:,0])
#         self.pl_dots.set_ydata(positions[:,1])
#         self.pl_dots.set_visible(True)
#       else:
#         self.pl_dots.set_visible(False)
# 
#       #draw_again()
#     main_Tk_window.canvas_fig.draw()
# 
#     return
# 
#   #---
#   def save(self):
#     try:
#       filename_tmp=filename_for_graph.split(".")
#       filename_for_positions_init = filename_for_graph.replace("."+filename_tmp[-1],".txt")
#       self.filename_for_positions=asksaveasfilename(
#         initialfile=filename_for_positions_init,filetypes=[('txt files','*.txt')])
# 
#       savetxt(self.filename_for_positions,self.positions_graph,
#               fmt=format_for_save
#               )
#     except:
#       tkMessageBox.showinfo("Info","The values have not been saved to a text file.")
# 
#     return
# 
#   #---------------------------------------------------------------------------
#   #  The modules below are to change the color, linewidth etc. for the plots.
#   #---------------------------------------------------------------------------
# 
#   def color(self,color):
#     self.pl_lines.set_color(color)
#     marker=self.pl_dots.get_marker()
#     if marker == "x" or marker == "+":
#       self.pl_dots.set_markeredgecolor(color)
#     else:
#       self.pl_dots.set_markerfacecolor(color)
#       self.pl_dots.set_markeredgecolor(color)
#     #draw_again()
#     main_Tk_window.canvas_fig.draw()
# 
#     return
# 
#   #---
#   def linewidth(self,linewidth):
#     self.pl_lines.set_linewidth(linewidth)
#     self.pl_dots.set_markeredgewidth(linewidth)
#     #draw_again()
#     main_Tk_window.canvas_fig.draw()
# 
#     return
# 
#   #---
#   def marker(self,marker):
#     self.pl_dots.set_marker(marker)
#     #draw_again()
#     main_Tk_window.canvas_fig.draw()
# 
#     return
# 
#   #---
#   def markersize(self,markersize):
#     self.pl_dots.set_markersize(markersize)
#     #draw_again()
#     main_Tk_window.canvas_fig.draw()
# 
#     return
# 
#   #---
#   def alpha(self,alpha):
#     self.pl_lines.set_alpha(alpha)
#     self.pl_dots.set_alpha(alpha)
#     main_Tk_window.canvas_fig.draw()
#     #draw_again()
# 
#     return
# 
# ############################################################################################
# ############################################################################################
# ### other class and modules ################################################################
# ############################################################################################
# ############################################################################################
# 
# class Show_crosshair:
# 
#   """
#     As described. Used as:
# 
#       crosshair_params=Show_crosshair()
# 
#     or just
# 
#       Show_crosshair()
# 
#     ### Input parameters ###
# 
#     x,y (optional) ... If we set x or y, the crosshair is move along the specified x or y axis.
#                        Otherwise, we can move the crosshair in any positon within two axes.
# 
#     ### Internal parameters ###
# 
#     self.crosshair_x_graph ... constraint of the position (if one of them is defined) 
#     self.crosshair_y_graph
# 
#   """
# 
#   def __init__(self,x=None,y=None):
# 
#     global ax
# 
#     if x != None and y!= None: return
# 
#     self.crosshair_x_graph=x;self.crosshair_y_graph=y
# 
#     # initialize the crosshair
#     self.pl_crosshair_h,=ax.plot([0,0],[0,0],color_crosshair+"-",
#                       linewidth=linewidth_crosshair)
#     self.pl_crosshair_v,=ax.plot([0,0],[0,0],color_crosshair+"-",
#                       linewidth=linewidth_crosshair)
# 
#     if x == y == None:
#       crosshair_is_activated_without_fixed_x_or_y=True
#     else:
#       crosshair_is_activated_without_fixed_x_or_y=False
# 
#     self.status=None
#     self.deactivate()
# 
#   #---
#   def activate(self):
# 
#     if self.status != "active":
# 
#       self.pl_crosshair_h.set_visible(True)
#       self.pl_crosshair_v.set_visible(True)
#       #range_for_plot=[ax.get_xlim()[0],ax.get_xlim()[1],ax.get_ylim()[0],ax.get_ylim()[1]]
#       #print range_for_plot
#       #draw_again()
#       main_Tk_window.canvas_fig.draw()
#       event_handling.activate("motion_notify_event",self.motion_notify)
# 
#       self.status="active"
# 
#     return
# 
#   #---
#   def deactivate(self):
# 
#     #print self.status,self.status,self.status
#     if self.status != "not active":
# 
#       self.pl_crosshair_h.set_visible(False)
#       self.pl_crosshair_v.set_visible(False)
#       main_Tk_window.canvas_fig.draw()
#       event_handling.deactivate("motion_notify_event")
# 
#       self.status="not active"
# 
#     return
# 
#   #---
#   def motion_notify(self,event):
#     """
#     Show_crosshair.motion_notify
#     """
#     if not event.inaxes: return
# 
#     if event.xdata < xaxis.x_mpl_edge[0] or event.xdata > xaxis.x_mpl_edge[1]: return
#     if event.ydata < yaxis.y_mpl_edge[0] or event.ydata > yaxis.y_mpl_edge[1]: return
# 
#     # if the x axis of the crosshair is determined
#     if self.crosshair_x_graph != None:
#       x_mpl_on_axis=(self.crosshair_x_graph-xaxis.e)/xaxis.d
#       y_mpl_on_axis=xaxis.a * x_mpl_on_axis + xaxis.b
#       x_mpl        =x_mpl_on_axis+(event.ydata-y_mpl_on_axis)*yaxis.a
#       self.motion_notify_sub(x_mpl,event.ydata)
# 
#     # if the y axis of the crosshair is determined
#     elif self.crosshair_y_graph != None:
#       y_mpl_on_axis=(self.crosshair_y_graph-yaxis.e)/yaxis.d
#       x_mpl_on_axis=yaxis.a * y_mpl_on_axis + yaxis.b
#       y_mpl        =y_mpl_on_axis+(event.xdata-x_mpl_on_axis)*xaxis.a
#       self.motion_notify_sub(event.xdata,y_mpl)
# 
#     # if no constraint for moving the crosshair
#     else:
#       self.motion_notify_sub(event.xdata,event.ydata)
# 
#     return
# 
#   #---
#   def motion_notify_sub(self,x_mpl,y_mpl):
# 
#     global ax
# 
#     # show crosshair
#     self.pl_crosshair_h.set_xdata(xaxis.x_mpl_edge)
#     self.pl_crosshair_h.set_ydata(xaxis.a * (array(xaxis.x_mpl_edge) - x_mpl)
#                                   + y_mpl)
#     self.pl_crosshair_v.set_ydata(yaxis.y_mpl_edge)
#     self.pl_crosshair_v.set_xdata(yaxis.a * (array(yaxis.y_mpl_edge) - y_mpl)
#                                   + x_mpl)
# 
#     # show the coordinate
#     ax.set_title(format_for_title % (
#            xaxis.get_x(x_mpl,y_mpl),yaxis.get_y(x_mpl,y_mpl)))
# 
#     #draw_again()
#     main_Tk_window.canvas_fig.draw()
# 
#     return
# 
#   #---
#   def set_x(self,x):
#     if xaxis.scale=="linear":
#       self.crosshair_x_graph=x
#     else:
#       self.crosshair_x_graph=log10(x)
#     self.crosshair_y_graph=None
#     crosshair_is_activated_without_fixed_x_or_y=False
# 
#     return
# 
#   #---
#   def set_y(self,y):
#     self.crosshair_x_graph=None
#     if yaxis.scale=="linear":
#       self.crosshair_y_graph=y
#     else:
#       self.crosshair_y_graph=log10(y)
#     crosshair_is_activated_without_fixed_x_or_y=False
# 
#     return
# 
#   #---
#   def reset_constraint(self):
#     self.crosshair_x_graph=None;self.crosshair_y_graph=None
#     crosshair_is_activated_without_fixed_x_or_y=True
# 
#     return
# 
# #------------------------------------------------------------------------
# def panzoom(event):
#   """
#     Action for the pan/zoom button
#   """
#   event_handling.panzoom_default_event(event)
#   if event_handling.status == "PAN":
#     crosshair.deactivate()
#     record_positions.deactivate()
#   else:
#     crosshair.activate()
#     record_positions.activate()
# 
#   return
# 
# #------------------------------------------------------------------------
# def zoomrect(event):
#   """
#     Action for the zoom rect button
#   """
# 
#   event_handling.zoomrect_default_event(event)
#   if event_handling.status == "ZOOM":
#     crosshair.deactivate()
#     record_positions.deactivate()
#   else:
#     crosshair.activate()
#     record_positions.activate()
# 
#   return
# 
# #------------------------------------------------------------------------
# def home(event):
#   """
#     Action for the home button
#   """
# 
#   event_handling.home_default_event(event)
# 
#   return
# 
# #------------------------------------------------------------------------
# 
# def update_format_for_values(txt):
#   """
#     As described. Input should be like "%.2e" or "%.4f".
#   """
# 
#   global format_for_title,format_for_print,format_for_save
# 
#   format_for_title = "x="+txt+", y="+txt
#   format_for_print = "["+txt+","+txt+"] is recorded."
#   format_for_save  = txt + "   " + txt
# 
#   return
# 
# #------------------------------------------------------------------------
# 
# def draw_again(): 
# 
#   #global ax,main_Tk_window
# 
#   ax.set_xlim(range_for_plot[0],range_for_plot[1])
#   ax.set_ylim(range_for_plot[2],range_for_plot[3])
#   main_Tk_window.canvas_fig.draw()
# 
# ##############################################################################################
# ##############################################################################################
# ### Copy of mpl_event_handling_sub ###########################################################
# ###     (copied on 2/6/2014)       ###########################################################
# ##############################################################################################
# ##############################################################################################
# 
# #---
# #from pylab import *
# #from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
# 
# class Mpl_event_handling_sub:
# 
#   def __init__(self):
#     self._connect_params={}
# 
#     self.status=None
# 
#     self.debugging=False
# 
#     self._panzoom_default =NavigationToolbar2TkAgg.pan
#     self._zoomrect_default=NavigationToolbar2TkAgg.zoom
#     self._home_default    =NavigationToolbar2TkAgg.home
# 
#   ### Modules ###
# 
#   def activate(self,func,module):
#     #global fig,status
#     if self.status != "PAN" and self.status != "ZOOM":
#       self.deactivate(func)
#       self._connect_params[func]=self.fig.canvas.mpl_connect(func,module)
#       self.status=self._connect_params.keys()
#     return
# 
#   #------------------------------------------------------------------------
#   def deactivate(self,func):
#     #global fig,status
#     if self.status != "PAN" and self.status != "ZOOM" and self._connect_params.has_key(func):
#       self.fig.canvas.mpl_disconnect(self._connect_params[func])
#       del self._connect_params[func]
#       self.status=self._connect_params.keys()
#     return
# 
#   #------------------------------------------------------------------------
#   def reset(self):
#     #global status
#     if self.status != "PAN" and self.status != "ZOOM" and self.status != None:
#       for func in self.status:
#         self.deactivate(func)
#       self.status=None
#     return
# 
#   #------------------------------------------------------------------------
#   def panzoom_default(self):
#     """
#       Equivalent to the pan/zoom button for the figure window.
#     """
#     #global fig,status
#     if self.status != "PAN":
#       self.reset()
#     self.fig.canvas.toolbar.pan()
#     self.status=fig.canvas.toolbar._active
#     return
# 
#   #------------------------------------------------------------------------
#   def panzoom_default_event(self,event):
#     """
#       Same as the above but used with the toolbar
#     """
#     #global fig,status
#     if self.debugging: print "Pan!"
#     if self.status != "PAN":
#       self.reset()
#     self._panzoom_default(event)
#     self.status=self.fig.canvas.toolbar._active
#     if self.debugging: print self.status
#     return
# 
#   #------------------------------------------------------------------------
#   def zoomrect_default(self):
#     """
#       Equivalent to the zoom rect button for the figure window.
#     """
#     #global fig,status
#     if self.status != "ZOOM":
#       self.reset()
#     self.fig.canvas.toolbar.zoom()
#     self.status=self.fig.canvas.toolbar._active
#     return
# 
#   #------------------------------------------------------------------------
#   def zoomrect_default_event(self,event):
#     """
#       Same as the above but used with the toolbar
#     """
#     #global fig,status
#     if self.debugging: print "Zoom!"
#     if self.status != "ZOOM":
#       self.reset()
#     self._zoomrect_default(event)
#     self.status=self.fig.canvas.toolbar._active
#     if self.debugging: print self.status
#     return
# 
#   #------------------------------------------------------------------------
#   def home_default(self):
#     """
#       Equivalent to the home button for the figure window.
#     """
#     #global fig,status
#     self.fig.canvas.toolbar.home()
#     self.status=self.fig.canvas.toolbar._active
#     return
# 
#   #------------------------------------------------------------------------
#   def home_default_event(self,event):
#     """
#       Same as the above but used with the toolbar
#     """
#     #global fig,status
#     if self.debugging: print "Home!"
#     self._home_default(event)
#     self.status=self.fig.canvas.toolbar._active
#     if self.debugging: print self.status
#     return
# 
# #import mpl_event_handling_sub as event_handling
# event_handling=Mpl_event_handling_sub()
# 
# 
# ##############################################################################################
# ##############################################################################################
# ### copy of tk_window_with_mpl_figure ########################################################
# ###       (copied on 2/6/2014)        ########################################################
# ##############################################################################################
# ##############################################################################################
# 
# 
# #---
# #import matplotlib
# #matplotlib.use('TkAgg')
# 
# #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# #from matplotlib.figure import Figure
# #import Tkinter as Tk
# #from tkFileDialog import asksaveasfilename
# #from sys import exit
# #import mpl_event_handling_sub as event_handling
# 
# class Tk_window_with_mpl_figure:
# 
#   ###########################################################################
#   ### Modules ###############################################################
#   ###########################################################################
# 
#   def open_window(self,figsize=None,
#                   dpi=None,
#                   console_width=None,
#                   console_height=None,
#                   console_side="bottom",
#                   show_toolbar=True
#                   ):
#     """
#         (the following parameters are not used anymore)
#         reset_button_function     ... function for indivdual buttons for the toolbar I made.
#         panzoom_button_function       The default is the same as the default toolbar for the mpl
#         zoomrect_button_function      figure window, but without interfering event handing
#         save_button_function          we separately sef 
# 
#     """
# 
#     #global fig,root,canvas_fig,console,toolbar       
# 
#     if console_side == "right":
#       side_f=Tk.LEFT
#       side_c=Tk.RIGHT
#     else:
#       side_f=Tk.TOP
#       side_c=Tk.BOTTOM
# 
#     self.fig = Figure(figsize=figsize,dpi=dpi)
#     event_handling.fig = self.fig
# 
#     console_outer=Tk.Frame(self.root,padx=2,pady=2)
#     console_outer.pack(side=side_c, fill=Tk.BOTH, expand=1)
#     self.console=Tk.Frame(console_outer,relief=Tk.RIDGE,bd=3,width=None,height=None)
#     self.console.pack(fill=Tk.BOTH, expand=1)
# 
#     self.canvas_fig = FigureCanvasTkAgg(self.fig, master=self.root)
#     self.canvas_fig.show()
# 
#     self.toolbar = NavigationToolbar2TkAgg( self.canvas_fig, self.root )
#     self.toolbar.update()
#     self.canvas_fig._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH)
#     self.canvas_fig.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
# 
#     if show_toolbar==False: self.fig.canvas.toolbar.pack_forget()
# 
#   #-------------------------------
#   def title(self,text):
#     self.root.wm_title(text)
#     return
# 
#   #-------------------------------
#   def panzoom_func(self,func=event_handling.panzoom_default_event):
#     NavigationToolbar2TkAgg.pan  = func
#     return
#   #-------------------------------
#   def zoomrect_func(self,func=event_handling.zoomrect_default_event):
#     NavigationToolbar2TkAgg.zoom = func
#     return
#   #-------------------------------
#   def home_func(self,func=event_handling.home_default_event):
#     NavigationToolbar2TkAgg.home = func
#     return
# 
#   #-------------------------------
#   def exit_event(self,event):
#     sys.exit()
#     return
# 
#   #-------------------------------
#   def quit_event(self,event):
#     sys.exit()
#     return
# 
#   #-------------------------------
#   def quit(self):
#     sys.exit(self)
#     return
# 
#   ###########################################################################
#   ### Main ##################################################################
#   ###########################################################################
# 
#   def __init__(self):
# 
#     self.root = Tk.Tk()
#     self.root.bind('<Control-q>',self.exit_event)
# 
#     ### Set the default function for the buttons in the tool bar. ###
# 
#     self.panzoom_func()
#     self.zoomrect_func()
#     self.home_func()
# 
#     return
# 
# #import tk_window_with_mpl_figure as main_Tk_window
# main_Tk_window=Tk_window_with_mpl_figure()
# 
# ##############################################################################################
# ##############################################################################################
# ### Begin ####################################################################################
# ##############################################################################################
# ##############################################################################################
# 
# if __name__ == '__main__':
# 
#   """
#     ### Parameters ###
# 
#     im                ... parameter for the Image command         
#     ax
#     range_for_plot
#     half_length_of_xtick
#     half_length_of_ytick
#     
#     xaxis,yaxis      ... for the classes Xaxis  and Yaxis
#     xticks,yticks    ... for the classes Xticks and Yticks
# 
#   """
# 
#   ### Set the GUI window ###
# 
#   main_Tk_window.title('read_plots.py')
# 
#   # define the function corresponding to the button
#   main_Tk_window.panzoom_func(panzoom)
#   main_Tk_window.zoomrect_func(zoomrect)
#   main_Tk_window.home_func(home)
# 
#   main_Tk_window.open_window(console_width=708,console_height=196)
# 
#   ### read the filename, and initialize parameters ###
# 
#   #filename_for_graph="Cotera01_f7.png"
#   filename_for_graph=askopenfilename(initialdir=os.getcwd(),
#                                     filetypes=[('png/jpg/gif/tiff','*.png *.jpg *.gif *.tiff')])
# 
#   update_format_for_values(format_for_values_default)
# 
#   # define the header of the file name
#   filename_tmp      = filename_for_graph.split(".")
#   filename_header   = filename_for_graph.replace("."+filename_tmp[-1],"")
#   filename_for_axes = filename_header+".axes"
#   del filename_tmp,filename_header
# 
#   # specify the file name if necessary
#   if filename_for_graph == "":
#     filename_for_graph=askopenfilename(filetypes=[('png','*.png'),('jpg','*.jpg'),('gif','*.gif')])
# 
#   ### show the original image specified by the file name
#   im=Image.open(filename_for_graph)
#   ax=main_Tk_window.fig.add_axes([0,0,1,0.93])
#   ax.imshow(im,origin="lower",interpolation="nearest")
# 
#   for tick in ax.xaxis.get_major_ticks():
#     tick.tick1On=False
#     tick.tick2On=False
# 
#   for tick in ax.yaxis.get_major_ticks():
#     tick.tick1On=False
#     tick.tick2On=False
# 
#   ax.set_xticklabels("")
#   ax.set_yticklabels("")
# 
#   # keep the present xy range
# 
#   range_for_plot=[ax.get_xlim()[0],ax.get_xlim()[1],
#                   ax.get_ylim()[0],ax.get_ylim()[1]]
# 
#   # define crosshair and record_positions
# 
#   crosshair = Show_crosshair()
#   record_positions = Record_positions()
# 
#   # read the parameters for axes if we have saved before
#   if glob.glob(filename_for_axes) != []:
# 
#     db=open(filename_for_axes,'rb')
#     params=pickle.load(db)
#     db.close()
#     xaxis = Xaxis(params['xaxis'])
#     yaxis = Yaxis(params['yaxis'])
#     del db,params
# 
#     crosshair.activate()
#     record_positions.activate()
#     
#     # Organize the widgets to measure values
#     Get_values()
# 
#   # otherwise, prepare for defining the axes
#   else:
# 
#     # set the half length of the tick
#     half_length_of_xtick=(range_for_plot[3]-range_for_plot[2])*0.04
#     half_length_of_ytick=(range_for_plot[1]-range_for_plot[0])*0.02
# 
#     # set the Axis vector
# 
#     xaxis=Xaxis()
#     yaxis=Yaxis()
#     xticks=Xticks()
#     yticks=Yticks()
# 
#     # organize the widgets for determining the axes
# 
#     Determine_axis()
# 
#   draw_again()
# 
#   Tk.mainloop()
# 
# 
# # # import pandas as pd
# # # import numpy as np
# # # 
# # # data = {'Col1' : [1,2,3,4], 'Col2' : [10,20,30,40], 'Col3' : [100,50,-30,-50], 'Col4' : ['AAA', 'BBB', 'AAA', 'CCC']}
# # # 
# # # df = pd.DataFrame(data=data, index = ['R1','R2','R3','R4'])
# # # 
# # # print df
# # # 
# # # 
# # # for row_index, row in df.iterrows():
# # #     print row
# # #     row.Col1 =5
# # #     print row
# # #     
# # #     
# # # print df
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # import pickle
# # from dazer_methods  import Dazer
# # from numpy import array
# #  
# # dz = Dazer()
# #   
# # # wavelength, Flux_array, Header_0 = dz.get_spectra_data('/home/vital/Astrodata/WHT_2016_04/Night1/objects/MAR2071_Blue_cr_f_t_w_e_fglobal.fits', ext=0)
# #   
# # dz.FigConf()
# # 
# # 
# # pickle_file = '/home/vital/Dropbox/Astrophysics/Data/WHT_observations/objects/3/A3_3_RedShiftCorrection.dz_pickle'
# # 
# # 
# # with open(pickle_file, 'rb') as fid:
# #     dz.Fig = pickle.load(fid)
# # 
# # dz.display_fig()
# # 
# # # 
# # # wavelength, Flux_array = array([0,1,2,3,4,5]), array([0,1,2,3,4,5])
# # # 
# # # 
# # # dz.data_plot(wavelength, Flux_array, 'spectrum1')
# # # dz.data_plot(wavelength, Flux_array*1.10, 'spectrum2')
# # # dz.data_plot(wavelength, Flux_array*1.20, 'spectrum3')
# # # dz.data_plot(wavelength, Flux_array*1.30, 'spectrum4')
# # # dz.data_plot(wavelength, Flux_array*1.40, 'spectrum5')
# # # dz.data_plot(wavelength, Flux_array*1.50, 'spectrum6')
# # # dz.data_plot(wavelength, Flux_array*1.60, 'spectrum7')
# # # dz.data_plot(wavelength, Flux_array*1.70, 'spectrum8')
# # # dz.data_plot(wavelength, Flux_array*1.80, 'spectrum9')
# # # dz.data_plot(wavelength, Flux_array*1.60, 'spectrum7', color='#bcbd22')
# # # dz.data_plot(wavelength, Flux_array*1.70, 'spectrum8', color='#7f7f7f')
# # # dz.data_plot(wavelength, Flux_array*1.80, 'spectrum9', color='#FFB5B8')
# # # # dz.data_plot(wavelength, Flux_array*1.90, 'spectrum10')
# # # 
# # # dz.FigWording(r'Wavelength $(\AA)$', 'Flux ' + r'$(erg\,cm^{-2} s^{-1} \AA^{-1})$', 'Plotting example')
# # #  
# # # dz.display_fig()
# # # 
# # # 
# # # 
# # # # from matplotlib import pyplot as plt, style
# # # # style.use('dark_background')
# # # # 
# # # # Fig = plt.figure(figsize = (12, 8))  
# # # # Axis = Fig.add_subplot(111)
# # # # 
# # # # 
# # # # Axis.plot([1,2,3,4], [1,2,3,4], label = 'coso')
# # # # 
# # # # Axis.legend()
# # # # 
# # # # plt.show()