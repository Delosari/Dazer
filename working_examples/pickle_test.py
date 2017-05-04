import numpy as np
import lmfit
import matplotlib.pyplot as plt

def gauss_curve(p, x):
    v = p.valuesdict()
    return v['a1'] * np.exp(-x / v['t1']) + v['a2'] * np.exp(-(x - 0.1) / v['t2'])

def residual(p, x, y):
    v = p.valuesdict()
    return v['a1'] * np.exp(-x / v['t1']) + v['a2'] * np.exp(-(x - 0.1) / v['t2']) - y

def lnprob(p, x, y):
    resid = residual(p, x, y)
    s = p['f']
    resid *= 1 / s
    resid *= resid
    resid += np.log(2 * np.pi * s**2)
    return -0.5 * np.sum(resid)


wave = np.linspace(1, 10, 250)
np.random.seed(0)
flux = 3.0 * np.exp(-wave / 2) - 5.0 * np.exp(-(wave - 0.1) / 10.) + 0.1 * np.random.randn(len(wave))
plt.plot(wave, flux)

p = lmfit.Parameters()
p.add_many(('a1', 4.), ('a2', 4.), ('t1', 3.), ('t2', 3., True))

mi = lmfit.minimize(residual, p, args=([wave, flux]), method='powell')

lmfit.printfuncs.report_fit(mi.params, min_correl=0.5)

plt.plot(wave, flux)
plt.plot(wave, gauss_curve(mi.params, wave), 'r')

mi.params.add('f', value=1, min=0.001, max=2)
 
mini = lmfit.Minimizer(lnprob, mi.params,  fcn_args=([wave, flux]))
 
res = mini.emcee(burn=300, steps=600, thin=10, params=mi.params)
 
lmfit.report_fit(res.params)
 
plt.plot(wave, gauss_curve(res.params, wave), 'b')

plt.show()












# import pyfits 
# import pysao 
# ds9 = pysao.ds9() 
# f = pyfits.open('/home/vital/Astrodata/WHT_2011_11/Night1/objects/9_Red_cr_f_t_w.fits') 
# ds9.view(f[0]) 
#
#
# from dazer_methods import Dazer
# import pickle
# 
# #Generate dazer object
# dz = Dazer()
#   
# # #Choose plots configuration
# # dz.FigConf()
# #   
# # #Import catalogue
# # Catalogue = dz.import_catalogue()
# #   
# # #Perform operations
# #   
# # x = [1,2,3,4,5,6]
# # y = [1,2,3,4,5,6]
# #   
# # #Plot the data
# # dz.data_plot(x, y, markerstyle = 'o')
# # 
# # pickle.dump(dz.Fig, file('/home/vital/Desktop/sinus.pickle','wb'))
# 
# # Load figure from disk and display
# dz.Fig = pickle.load(open('/home/vital/Desktop/sinus.pickle','rb'))
# dz.display_fig()
