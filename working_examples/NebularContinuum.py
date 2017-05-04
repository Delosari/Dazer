#!/usr/bin/env python

import pymc
import pyneb as pn
import numpy as np
from dazer_methods import Dazer
from scipy.interpolate import interp1d
from libraries.Astro_Libraries.Nebular_Continuum import Nebular_Bayesian
from lmfit.models import LinearModel

H1 = pn.RecAtom('H', 1)

print H1.getEmissivity(tem=10000, den=100, wave=6563) / H1.getEmissivity(tem=10000, den=100, wave=4861)


#Declare code classes
dz = Dazer()
nb = Nebular_Bayesian()
 
#Define plot frame and colors
dz.FigConf(n_colors=5)
 
nb.ObjectData_dict = {}
nb.ObjectData_dict['nHeII_HII']             = 0.075
nb.ObjectData_dict['TOIII']                 = 10000
nb.ObjectData_dict['TOIII_error']           = 500.0
nb.ObjectData_dict['Flux_Hbeta_Normalize']  = 2.86  * 1e-14  / 1e-14
nb.ObjectData_dict['Error_Hbeta_Normalize'] = 0.05   * 1e-14  / 1e-14
 
wave = np.linspace(3500, 7000, 7000 - 3500 - 1)
nebular_flux_n = nb.Calculate_Nebular_gamma(nb.ObjectData_dict['TOIII'], nb.ObjectData_dict['Flux_Hbeta_Normalize'], nb.ObjectData_dict['nHeII_HII'] , 0.00, wave)
dz.data_plot(wave, nebular_flux_n, 'Nebular flux normalized', dz.ColorVector[2][1])
 
 
#---------------Bayesian mod-----------------
    
#Declare the MCMC dictionary
MCMC_dict = nb.model_parameters(wave, nebular_flux_n, 0.05)
           
#Run MCMC with MAP
MAP_Model = pymc.MAP(MCMC_dict)
MAP_Model.fit(method = 'fmin_powell') 
MAP_Model.revert_to_max()
      
version = '_real_2000_100'
Folder_database     = '/home/vital/Desktop/db_Testing/'
Db_name             = 'db_Continuum' + version
Db_global_file_name = 'db_Continuum_Global' + version
      
M = pymc.MCMC(MAP_Model.variables, db = 'pickle', dbname = Folder_database + Db_name)
M.sample(iter=2000, burn=100) 
M.write_csv(Folder_database + Db_global_file_name, variables=['He_abud', 'Te', 'Flux_Recomb'])
M.db.close() 
  
dz.FigWording(r'Wavelength $(\AA)$', 'Flux' + r'$(erg\,cm^{-2} s^{-1} \AA^{-1})$', 'Fitting nebular continuum')   
  
               
#-----------------------------------------------------------------------------------------------------
dz.display_fig()

print "All plots generated"






# import pymc as mc
# import numpy as np
# import pylab as pl
# 
# def GaussFunc(x, amplitude, centroid, sigma):
#     return amplitude * np.exp(-0.5 * ((x - centroid) / sigma)**2)
# 
# 
# 
# wavelength = np.arange(5000, 5050, 0.02)
# 
# # Profile 1
# centroid_one = 5025.0
# sigma_one = 2.2
# height_one = 0.8
# profile1 = GaussFunc(wavelength, height_one, centroid_one, sigma_one, )
# 
# # Profile 2
# centroid_two = 5027.0
# sigma_two = 1.2
# height_two = 0.5
# profile2 = GaussFunc(wavelength, height_two, centroid_two, sigma_two, )
# 
# 
# est_centroid_one = mc.Uniform("est_centroid_one", 5000, 5050 )
# est_centroid_two = mc.Uniform("est_centroid_two", 5000, 5050 )
# 
# est_sigma_one = mc.Uniform( "est_sigma_one", 0, 5 )
# est_sigma_two = mc.Uniform( "est_sigma_two", 0, 5 )
# 
# est_height_one = mc.Uniform( "est_height_one", 0, 5 ) 
# est_height_two = mc.Uniform( "est_height_two", 0, 5 ) 
# 
# #std deviation of the noise, converted to precision by tau = 1/sigma**2
# precision= 1./mc.Uniform("std", 0, 1)**2
# 
# #Set up the model's relationships.
# 
# @mc.deterministic( trace = False) 
# def est_profile_1(x = wavelength, centroid = est_centroid_one, sigma = est_sigma_one, height= est_height_one):
#     return GaussFunc( x, height, centroid, sigma )
# 
# 
# @mc.deterministic( trace = False) 
# def est_profile_2(x = wavelength, centroid = est_centroid_two, sigma = est_sigma_two, height= est_height_two):
#     return GaussFunc( x, height, centroid, sigma )
# 
# 
# @mc.deterministic( trace = False )
# def mean( profile_1 = est_profile_1, profile_2 = est_profile_2 ):
#     return profile_1 + profile_2
# 
# # Measured values
# noise = np.random.normal(0.0, 0.02, len(wavelength))
# combined = profile1 + profile2 + noise
# 
# 
# observations = mc.Normal("obs", mean, precision, value = combined, observed = True)
# 
# 
# model = mc.Model([est_centroid_one, 
#               est_centroid_two, 
#                 est_height_one,
#                 est_height_two,
#                 est_sigma_one,
#                 est_sigma_two,
#                 precision])
# 
# #always a good idea to MAP it prior to MCMC, so as to start with good initial values
# map_ = mc.MAP( model )
# map_.fit()
# 
# mcmc = mc.MCMC( model )
# mcmc.sample( 50000,40000 ) 
# 
# # Measured values
# noise = np.random.normal(0.0, 0.02, len(wavelength))
# combined = profile1 + profile2 + noise
# 
# # If you want to plot what this looks like
# pl.plot(wavelength, combined, label="Measured")
# pl.plot(wavelength, profile1, color='red', linestyle='dashed', label="1")
# pl.plot(wavelength, profile2, color='green', linestyle='dashed', label="2")
# pl.title("Feature One and Two")
# pl.legend()
# 
# pl.show()
# 
# import pymc
# import numpy as np
# from dazer_methods import Dazer
# from libraries.Astro_Libraries.Nebular_Continuum import NebularContinuumCalculator
# from libraries.Astro_Libraries.Nebular_Continuum import Nebular_Bayesian
# 
# def Bayesian_ContinuumFit():
# 
#     y_plus      = pymc.Uniform('He_abud', 0.050, 0.15)
#     Te          = pymc.Uniform('Te', 10000, 14000)
#     Flux_Hbeta  = pymc.Uniform('Flux_Hbeta', 0.5, 1.5)
#      
#     #Calculate nebular continuum
#     @pymc.deterministic
#     def Calculate_Continuum(y_plus=y_plus, Te=Te, Flux_Hbeta=Flux_Hbeta):
#          
#         nebCalc.PropertiesfromUser(Te, Flux_Hbeta * 1e-14, y_plus, 0.0, Wave, Calibration = 'Zanstra')
#          
#         Gamma_Total, Gamma_lambda, Gamma_FB_HI, Gamma_FB_HeI, Gamma_FB_HeII, Gamma_2q, Gamma_FF = nebCalc.Calculate_Nebular_gamma()
#          
#         nebCalc.Zanstra_Calibration('Hbeta', Flux_Hbeta * 1e-14, Gamma_Total)
#          
#         return nebCalc.Zanstra_Calibration('Hbeta', Flux_Hbeta, Gamma_Total)
#      
#     #Likelihood
#     @pymc.stochastic(observed=True)
#     def Likelihood_model(value=NebularInt_Hbeta, nebInt_cont = Calculate_Continuum, sigmaLines = 0.05):
#         chi_F           = np.sum(np.square(nebInt_cont - value) / np.square(sigmaLines))
#         return - chi_F / 2
# 
#     return locals()
# 
# #Declare objects
# 
# dz      = Dazer()
# nb      = Nebular_Bayesian()
# nebCalc = NebularContinuumCalculator()
#  
# #Define operation
# Catalogue_Dic       = dz.import_catalogue()
# pattern             = Catalogue_Dic['Datatype'] + '.fits'
# Lineslog_extension  = '_' + Catalogue_Dic['Datatype'] + '_LinesLog_v3.txt' 
#  
# #Locate files on hard drive
# FilesList           = dz.Folder_Explorer(pattern, Catalogue_Dic['Obj_Folder'], CheckComputer=False)
# 
# #Define plot frame and colors
# dz.FigConf(n_colors=5)
#   
# #Extract the data from fits files
# Wave = np.linspace(3000, 9000, 9000 - 3000)
#  
# Flux_Hbeta  = 1e-14
# Te          = 12000
# nHeII_HII   = 0.1
# nHeIII_HII  = 0.0
# 
# #-- Calculate nebular continuum
# nebCalc.PropertiesfromUser(Te, Flux_Hbeta, nHeII_HII, nHeIII_HII, Wave, Calibration = 'Zanstra')
# nb.calculate_wavelength_intervals(Wave)
#  
# #-- Calculate continuous emissino coefficients:
# Gamma_Total, Gamma_lambda, Gamma_FB_HI, Gamma_FB_HeI, Gamma_FB_HeII, Gamma_2q, Gamma_FF = nebCalc.Calculate_Nebular_gamma()
#  
# #-- Caculate nebular flux with different calibration methods
# NebularInt_Hbeta = nebCalc.Zanstra_Calibration('Hbeta', Flux_Hbeta, Gamma_Total)
# NebularFlux_lambda = nb.Calculate_Nebular_gamma(Te, Flux_Hbeta, nHeII_HII, nHeIII_HII)
# 
# noise = np.random.normal(0,  0.05, size = len(Wave))
# 
# Flux_Observed = NebularInt_Hbeta + NebularInt_Hbeta * noise
# 
# #Plotting the data
# 
# dz.data_plot(Wave, NebularInt_Hbeta, 'Nebular continuum - Hbeta calibration', dz.ColorVector[2][1])
# # dz.data_plot(Wave, Flux_Observed, 'Noisy continuum', dz.ColorVector[2][2])
# dz.data_plot(Wave, NebularFlux_lambda, 'Nebular continuum NEW', dz.ColorVector[2][3])
# 
# #Format the graphs
# PlotTitle = r'Nebular continuum calculation'
# dz.FigWording(r'Wavelength $(\AA)$', 'Flux' + r'$(erg\,cm^{-2} s^{-1} \AA^{-1})$', PlotTitle)       
# 
# 
# # #---------------Bayesian mod-----------------
# #  
# # #Declare the MCMC dictionary
# # MCMC_dict = Bayesian_ContinuumFit()
# #        
# # #Run MCMC with MAP
# # MAP_Model = pymc.MAP(MCMC_dict)
# # MAP_Model.fit(method = 'fmin_powell') 
# # MAP_Model.revert_to_max()
# # 
# # version = '_v5000_1000'
# # Folder_database     = '/home/vital/Desktop/db_Testing/'
# # Db_name             = 'db_Continuum' + version
# # Db_global_file_name = 'db_Continuum_Global' + version
# # 
# # M = pymc.MCMC(MAP_Model.variables, db = 'pickle', dbname = Folder_database + Db_name)
# # M.sample(iter=5000, burn=1000) 
# # M.write_csv(Folder_database + Db_global_file_name, variables=['He_abud', 'Te', 'Flux_Hbeta'])
# # M.db.close() 
# 
# dz.display_fig()
#  
# print 'All data treated', dz.display_errors()