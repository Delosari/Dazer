'''
Created on Oct 19, 2016
 
@author: vital
'''

import numpy            as np
from sklearn.mixture    import GMM, GaussianMixture
import matplotlib.pyplot as plt
from scipy.stats        import norm

# import pandas           as pd
# from astropy.io.fits    import getheader, getdata
# from lmfit.models       import GaussianModel, LinearModel
# from scipy.integrate    import simps
# from scipy.signal       import argrelextrema
# from collections        import OrderedDict

# 
# def get_object_files(reference_object):
#     
#     object_files = {}
# 
#     object_files['log'] = '/home/vital/Dropbox/Astrophysics/Data/WHT_Catalogue_SulfurRegression/Objects/{code_name}/{code_name}_WHT_LinesLog_v3.txt'.format(code_name = reference_object) 
#     object_files['fits'] = '/home/vital/Dropbox/Astrophysics/Data/WHT_Catalogue_SulfurRegression/Objects/{code_name}/obj{code_name}_WHT.fits'.format(code_name = reference_object) 
# 
#     return object_files
# 
# def load_spectrum(file_address, ext=1):
#     
#     header    = getheader(file_address, ext=1)
#     data_page = getdata(file_address, ext=ext)
# 
#     #In the case a fits file I made
#     if ('WHTJOINW' in header) or ("STALIGHT" in header) or ("NEBUSPEC" in header):
#         wavelength = data_page['Wave']
#         flux_array = data_page['Int']
# 
#     elif "COEFF0" in header:
#         dw              = 10.0**header['COEFF1']                # dw = 0.862936 INDEF (Wavelength interval per pixel)
#         Wmin            = 10.0**header['COEFF0']
#         pixels          = header['NAXIS1']                      # nw = 3801 number of output pixels
#         Wmax            = Wmin + dw * pixels
#         wavelength      = np.linspace(Wmin,Wmax,pixels,endpoint=False)
#         flux_array      = data_page['Int']
#         
#     elif "LTV1" in header:
#         StartingPix     = -1 * header['LTV1']                   # LTV1 = -261. 
#         Wmin_CCD        = header['CRVAL1']
#         dw              = header['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
#         pixels          = header['NAXIS1']                      # nw = 3801 number of output pixels
#         Wmin            = Wmin_CCD + dw * StartingPix
#         Wmax            = Wmin + dw * pixels
#         wavelength      = np.linspace(Wmin,Wmax,pixels,endpoint=False)
#         flux_array      = data_page['Int']
#     
#     else:
#         Wmin            = header['CRVAL1']
#         dw              = header['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
#         pixels          = header['NAXIS1']                      # nw = 3801 number of output pixels
#         Wmax            = Wmin + dw * pixels
#         wavelength      = np.linspace(Wmin,Wmax,pixels,endpoint=False)
#         flux_array      = data_page['Int']
# 
#     return wavelength, flux_array, header
# 
# def load_log(line_reference, lines_log_address):
#         
#     lines_frame         = pd.read_csv(lines_log_address, skiprows = [0], delim_whitespace = True, header = 0, index_col = 0)
#       
#     line_row           = lines_frame.loc[(line_reference), :]
#       
#     return lines_frame, line_row
# 
# def get_line_region(wave, flux, zone_ids):
#     
#     blue_wave, blue_flux                    = wave[zone_ids[0]:zone_ids[1]], flux[zone_ids[0]:zone_ids[1]]
#     line_wave, line_flux                    = wave[zone_ids[2]:zone_ids[3]], flux[zone_ids[2]:zone_ids[3]]
#     red_wave, red_flux                      = wave[zone_ids[4]:zone_ids[5]], flux[zone_ids[4]:zone_ids[5]]
# 
#     Continuum_wave, Continuum_flux          = np.hstack([blue_wave, red_wave]), np.hstack([blue_flux, red_flux])
#     
#     lineal_mod                              = LinearModel(prefix='lineal_') #This should go in the init
#     Lineal_parameters                       = lineal_mod.guess(Continuum_flux, x=Continuum_wave)
#     line_zerolev                            = Lineal_parameters['lineal_slope'].value * line_wave + Lineal_parameters['lineal_intercept'].value
#     err_continuum                           = np.std(Lineal_parameters['lineal_slope'].value * Continuum_wave + Lineal_parameters['lineal_intercept'].value - Continuum_flux)
#     
#     return line_wave, line_flux, Continuum_wave, Continuum_flux, line_zerolev, err_continuum
# 
# def normalize_line_region(line_wave, line_flux, Continuum_wave, Continuum_flux, line_zerolev, err_continuum):
#     
#     idx_maxFlux         = np.argmax(line_flux)
#     maxWave, maxFlux    = line_wave[idx_maxFlux], line_flux[idx_maxFlux]
#     
#     #Move the line to zero wavelength with respect to peak 
#     line_wave_n, Continuum_wave_n = line_wave - maxWave, Continuum_wave - maxWave
# 
#     #Normalize peak intensity with respect curve area (without continuum)     
#     area_line = simps(line_flux, line_wave_n) - simps(line_zerolev, line_wave_n)
#     line_flux_n, Continuum_flux_n, line_zerolev_n, err_continuum_n = line_flux/area_line, Continuum_flux/area_line, line_zerolev/area_line, err_continuum/area_line
# 
#     return line_wave_n, line_flux_n, Continuum_wave_n, Continuum_flux_n, line_zerolev_n, err_continuum_n

# current_line            = 'H1_6563A'
# 
# obj_files               = get_object_files('11')
# lines_df, l_series      = load_log(current_line, obj_files['log'])
# wave, flux, header      = load_spectrum(obj_files['fits'])
#  
# wave_line               = l_series['TheoWavelength']
# 
# #We should self this
# regions_indeces = wave.searchsorted((l_series.Wave1, l_series.Wave2, l_series.Wave3, l_series.Wave4, l_series.Wave5, l_series.Wave6)) 
# 
# #Slice the spectrum region of interest
# line_wave, line_flux, Continuum_wave, Continuum_flux, line_zerolev, err_continuum = get_line_region(wave, flux, regions_indeces)
# 
# #Normalize the data
# line_wave_n, line_flux_n, Continuum_wave_n, Continuum_flux_n, line_zerolev_n, err_continuum_n = normalize_line_region(line_wave, line_flux, Continuum_wave, Continuum_flux, line_zerolev, err_continuum)
# 
# #Geometry of the line
# max_idces   = argrelextrema(line_flux, np.greater)[0]
# peak_fluxes = line_flux_n[max_idces]
# peak_waves  = line_wave_n[max_idces]
# 
# for j in range(len(line_wave)):
#     print '[{wave_v}, {flux_v}],'.format(wave_v = line_wave[j], flux_v = line_flux[j])




#Raw data
data = np.array([[6535.62597656, 7.24362260936e-17],
        [6536.45898438, 6.28683338273e-17],
        [6537.29248047, 5.84596729207e-17],
        [6538.12548828, 8.13193914837e-17],
        [6538.95849609, 6.70583742068e-17],
        [6539.79199219, 7.8511483881e-17],
        [6540.625, 9.22121293063e-17],
        [6541.45800781, 7.81353615478e-17],
        [6542.29150391, 8.58095991639e-17],
        [6543.12451172, 9.30569784967e-17],
        [6543.95800781, 9.92541957936e-17],
        [6544.79101562, 1.1682282379e-16],
        [6545.62402344, 1.21238102142e-16],
        [6546.45751953, 1.51062780724e-16],
        [6547.29052734, 1.92193416858e-16],
        [6548.12402344, 2.12669644265e-16],
        [6548.95703125, 1.89356624109e-16],
        [6549.79003906, 1.62571112976e-16],
        [6550.62353516, 1.73262984876e-16],
        [6551.45654297, 1.79300635724e-16],
        [6552.29003906, 1.93990357551e-16],
        [6553.12304688, 2.15530881856e-16],
        [6553.95605469, 2.13273711105e-16],
        [6554.78955078, 3.03175829363e-16],
        [6555.62255859, 3.17610250166e-16],
        [6556.45556641, 3.75917668914e-16],
        [6557.2890625, 4.64631505826e-16],
        [6558.12207031, 6.9828152092e-16],
        [6558.95556641, 1.19680535606e-15],
        [6559.78857422, 2.18677945421e-15],
        [6560.62158203, 4.07692754678e-15],
        [6561.45507812, 5.89089137849e-15],
        [6562.28808594, 7.48005986578e-15],
        [6563.12158203, 7.49293900174e-15],
        [6563.95458984, 4.59418727426e-15],
        [6564.78759766, 2.25848015792e-15],
        [6565.62109375, 1.04438093017e-15],
        [6566.45410156, 6.61019482779e-16],
        [6567.28759766, 4.45881319808e-16],
        [6568.12060547, 4.1486649376e-16],
        [6568.95361328, 3.69435405178e-16],
        [6569.78710938, 2.63747028003e-16],
        [6570.62011719, 2.58619514057e-16],
        [6571.453125, 2.28424298265e-16],
        [6572.28662109, 1.85772271843e-16],
        [6573.11962891, 1.90082094593e-16],
        [6573.953125, 1.80158097764e-16],
        [6574.78613281, 1.61992695352e-16],
        [6575.61914062, 1.44038495311e-16],
        [6576.45263672, 1.6536593789e-16],
        [6577.28564453, 1.48634721076e-16],
        [6578.11914062, 1.28145245545e-16],
        [6578.95214844, 1.30889102898e-16],
        [6579.78515625, 1.42521644591e-16],
        [6580.61865234, 1.6919170778e-16],
        [6581.45166016, 2.35394744146e-16],
        [6582.28515625, 2.75400454352e-16],
        [6583.11816406, 3.42150435774e-16],
        [6583.95117188, 3.06301301529e-16],
        [6584.78466797, 2.01059337187e-16],
        [6585.61767578, 1.36484708427e-16],
        [6586.45068359, 1.26422274651e-16],
        [6587.28417969, 9.79250952203e-17],
        [6588.1171875, 8.77299287344e-17],
        [6588.95068359, 6.6478752208e-17],
        [6589.78369141, 4.95864370066e-17]])



#Get the data
obs_wave, obs_flux = data[:,0], data[:,1]

#Center the x data in zero and normalized the y data to the area of the curve
n_wave = obs_wave - obs_wave[np.argmax(obs_flux)]
n_flux = obs_flux / sum(obs_flux) 

#Generate a distribution of points matcthing the curve
line_distribution   = np.random.choice(a = n_wave, size = 100000, p = n_flux)
number_points       = len(line_distribution)

#Run the fit
gmm = GaussianMixture(n_components = 4)
gmm.fit(np.reshape(line_distribution, (number_points, 1)))
gauss_mixt = np.array([p * norm.pdf(n_wave, mu, sd) for mu, sd, p in zip(gmm.means_.flatten(), np.sqrt(gmm.covariances_.flatten()), gmm.weights_)])
gauss_mixt_t = np.sum(gauss_mixt, axis = 0)  

#Plot the data
fig, axis = plt.subplots(1, 1, figsize=(10, 12))
axis.plot(n_wave, n_flux, label = 'Normalized observed flux')
axis.plot(n_wave, gauss_mixt_t, label = '4 components fit')
  
for i in range(len(gauss_mixt)):
    axis.plot(n_wave, gauss_mixt[i], label = 'Gaussian '+str(i), linestyle = '--')

axis.set_xlabel('normalized wavelength', fontsize = 15)
axis.set_ylabel('normalized flux', fontsize = 15)
axis.set_title('Sklearn GM fit', fontsize = 15)

axis.legend(fontsize = 15)
plt.show()


# gmm         = GMM(n_components = 3, min_covar = 1e-12)
#  
# gmm.fit(np.reshape(line_distribution,(number_points, 1)))
#  
# gauss_mixt      = np.array([p * norm.pdf(line_wave_n, mu, sd) for mu, sd, p in zip(gmm.means_.flatten(), np.sqrt(gmm.covars_.flatten()), gmm.weights_)])
# gauss_mixt_t    = np.sum(gauss_mixt, axis = 0)  
#  
# # gmm.means_ = np.array([[-15], [0], [0], [20]])
# # gmm.covars_ = np.array([[1], [1], [4], [1]]) ** 2
#  
# #Figure conf
# fig, axis = plt.subplots(1, 1, figsize=(10, 12))
# axis.plot(line_wave_n, Gaussian_Mixture)
# axis.plot(Continuum_wave_n, Continuum_flux_n)
# axis.plot(line_wave_n, line_zerolev_n)
#  
# # axis.hist(line_distribution, bins = len(line_wave_n[0:-1]), color = 'grey', edgecolor = 'k', label = "Empirical")
#   
# for i in range(len(gauss_mixt)):
#     axis.plot(line_wave_n, gauss_mixt[i], label='Gaussian '+str(i))
#  
# axis.legend()
#  
# plt.show()




