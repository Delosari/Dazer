import numpy as np
from dazer_methods import Dazer
from lib.Astro_Libraries.Nebular_Continuum import NebularContinuumCalculator
from timeit import default_timer as timer
import timeit

dz      = Dazer()
neb     = NebularContinuumCalculator()
for_dic = {'figure.figsize':(14,20)}
dz.FigConf(for_dic)

Te          = 10000.0
wave        = np.arange(2000, 7000)
HeI,HeII    = 0.01, 0.001
Flux_recomb = 1.5e-13

wave_ryd = 3

neb.PropertiesfromUser(Te, Flux_recomb, HeI, HeII, wave)
 
neb.load_neb_constants()
wave_ryd            = (neb.const['h'] * neb.const['c_Angs']) / (neb.const['Ryd2erg'] * wave)
Gamma_FF_HI_fast    = neb.free_free_gCont(wave, Te)
Gamma_2q_HI_fast    = neb.two_photon_gCont(wave, Te)
Gamma_FB_HI_fast    = neb.free_bound_gCont(wave_ryd, Te, neb.HI_fb_dict)
Gamma_FB_HeI_fast   = neb.free_bound_gCont(wave_ryd, Te, neb.HeI_fb_dict)
Gamma_FB_HeII_fast  = neb.free_bound_gCont(wave_ryd, Te, neb.HeII_fb_dict)
Gamma_FB_HI         = neb.FreeBoundContinuum_EP("HeII")

Gamma_Total, Gamma_lambda, Gamma_FB_HI, Gamma_FB_HeI, Gamma_FB_HeII, Gamma_2q, Gamma_FF = neb.Calculate_Nebular_gamma()
Flux_Total = neb.Zanstra_Calibration('Halpha', Flux_recomb, Gamma_Total)
dz.data_plot(wave, Flux_Total,    label = 'Total old')

neb_gCont = neb.calculate_neb_gCont(wave, Te, HeI, HeII)
neb_fCont = neb.gCont_calibration(wave, Te, Flux_recomb, neb_gCont)
dz.data_plot(wave, neb_fCont,  label = 'Total new')




# dz.data_plot(wave, Gamma_FF_HI_fast,    label = 'H Free-Free')
# dz.data_plot(wave, 4 * Gamma_FF_HI_fast,label = 'He Free-Free')
# dz.data_plot(wave, Gamma_2q_HI_fast,    label = 'H Two-photon')
# dz.data_plot(wave, Gamma_FB_HI_fast,    label = 'H Free-Bound')
# dz.data_plot(wave, Gamma_FB_HeI_fast,   label = r'$He$ Free-Bound')
# dz.data_plot(wave, Gamma_FB_HeII_fast,  label = r'$He^{+}$ Free-Bound')
# dz.data_plot(wave, Gamma_FB_HI,    label = 'H Free-Bound EP')

dz.Axis.set_yscale('log')    
dz.FigWording(r'Wavelength $(\AA)$', r'$\gamma$ $\left(erg\,cm^{3}\,s^{-1}\,Hz^{-1}\right)$', title = 'Hydrogen and helium nebular components')
# dz.savefig('/home/vital/Dropbox/Astrophysics/Seminars/Cloudy School 2017/nebular_components')
dz.display_fig()
# start = timer()   
# Gamma_FF_HI = neb.FreeFreeContinuum("HI")
# end = timer()
# print 'FF', ' time ', (end - start)
# 
# start = timer()   
# Gamma_2q    = neb.TwoPhotonContinuum()
# end = timer()
# print '2q', ' time ', (end - start) 
# 
# start = timer()   
# Gamma_FB_HI = neb.FreeBoundContinuum_EP("HI")
# end = timer()
# print 'FB_HI', ' time ', (end - start)
# 
# start = timer()   
# Gamma_FF_HI_fast = neb.free_free_gCont(wave, Te)
# end = timer()
# print 'FF', ' time ', (end - start) 
# 
# start = timer()   
# Gamma_2q_HI_fast = neb.two_photon_gCont(wave, Te)
# end = timer()
# print '2q', ' time ', (end - start) 
# 
# start = timer()   
# Gamma_FB_HI_fast = neb.free_bound_gCont(wave, Te, 'HI')
# end = timer()
# print 'FB_HI', ' time ', (end - start)

# start = timer()   
# Gamma_FB_HI_fast = neb.free_bound_gCont(wave, Te, 'HeI')
# end = timer()
# print 'FB_HeI', ' time ', (end - start)
# 
# start = timer()   
# Gamma_FB_HI_fast = neb.free_bound_gCont(wave, Te, 'HeII')
# end = timer()
# print 'FB_HeII', ' time ', (end - start)

# import numpy
# from scipy import interpolate
# x = numpy.array([0.0, 0.60, 1.0])
# y = numpy.array([0.0, 0.25, 0.80, 1.0])
# z = numpy.array([ 
#    [ 1.4 ,  6.5 ,  1.5 ,  1.8 ],
#    [ 8.9 ,  7.3 ,  1.1 ,  1.09],
#    [ 4.5 ,  9.2 ,  1.8 ,  1.2 ]])
# # you have to set kx and ky small for this small example dataset
# # 3 is more usual and is the default
# # s=0 will ensure this interpolates.  s>0 will smooth the data
# # you can also specify a bounding box outside the data limits
# # if you want to extrapolate
# sp = interpolate.RectBivariateSpline(x, y, z, kx=2, ky=2, s=0)
# 
# print sp([0.60], [0.25])  # array([[ 7.3]])
# print sp([0.25], [0.60])  # array([[ 2.66427408]])
# 
# print sp([0.60], y)  # array([[ 7.3]])
# 
# print sp(x, 0.8)  # array([[ 7.3]])



# 
# x_array = np.arange(0, 15)
# y_delta = np.zeros(len(x_array))
# y_delta[3], y_delta[7], y_delta[13] = 1, 2, 3
# 
# step_function = foo2(y_delta)
# 
# x_high_resolution = np.linspace(0, 15, 30)
# print 'bicho', x_high_resolution
# 
# 
# # dz.data_plot(x_array, y_delta, label='delta array', markerstyle='o')
# # dz.Axis.step(x_array, step_function, label='Step function', color='red')
# 
# delta_interpolated  = np.interp(x_high_resolution, x_array, y_delta)
# step_interpolated   = np.interp(x_high_resolution, x_array, step_function)
# 
# dz.data_plot(x_high_resolution, delta_interpolated, label='delta interpolation', markerstyle='o')
# dz.Axis.step(x_high_resolution, step_interpolated, label='step interpolation', linestyle='--', color='red')
# 
# dz.FigWording(r'X', r'Y', title = 'Increasing resolution in delta and step functions')
#  
# dz.display_fig()
