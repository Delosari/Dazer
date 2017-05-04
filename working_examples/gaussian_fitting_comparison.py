import  numpy as np
from    matplotlib import pyplot as plt, rc
from    lmfit.models import GaussianModel, ExponentialModel, LinearModel
from    lmfit import Parameters, minimize, report_fit
from    collections import OrderedDict
from    scipy.signal import argrelextrema
from    pyneb import Atom
from    scipy.integrate    import simps
from    scipy.integrate import trapz
from    uncertainties import ufloat, unumpy
from    itertools import izip

def Emission_Threshold(LineLoc, TotalWavelen, TotalInten, BoxSize = 70):
    
    #Use this method to determine the box and location of the emission lines
    Bot                 = LineLoc - BoxSize
    Top                 = LineLoc + BoxSize
    
    indmin, indmax      = np.searchsorted(TotalWavelen, (Bot, Top))
    if indmax > (len(TotalWavelen)-1):
        indmax = len(TotalWavelen)-1
    
    PartialWavelength   = TotalWavelen[indmin:indmax]
    PartialIntensity    = TotalInten[indmin:indmax]
    
    Bot                 = LineLoc - 2
    Top                 = LineLoc + 2
    
    indmin, indmax      = np.searchsorted(PartialWavelength, (Bot, Top))
    
    LineHeight          = np.max(PartialIntensity[indmin:indmax])
    LineExpLoc          = np.median(PartialWavelength[np.where(PartialIntensity == LineHeight)])
          
    return PartialWavelength, PartialIntensity, LineHeight, LineExpLoc

def ufloatDict_nominal(ufloat_dict):
    'This gives us a dictionary of nominal values from a dictionary of uncertainties'
    return OrderedDict(izip(ufloat_dict.keys(), map(lambda x: x.nominal_value, ufloat_dict.values())))
 
def ufloatDict_stdev(ufloat_dict):
    'This gives us a dictionary of nominal values from a dictionary of uncertainties'
    return OrderedDict(izip(ufloat_dict.keys(), map(lambda x: x.std_dev, ufloat_dict.values())))
 
def Observational_data_blended():
     
    #Since we know the emision line is Halpha we may declare the their locations
    mu_Emissions = np.array([6548.0, 6563.0, 6584.0])
     
    #Wavelengths marking the blue continuum, line region and red continuum
    Wavelength_regions = np.array([6490.55, 6524.84, 6539.56, 6592.46, 6594.82, 6615.0])
 
    #Halpha line
    Spectrum = np.array([[ 6478.64 , 8.44735e-17 ],
                        [ 6480.38 , 8.72528e-17 ],
                        [ 6482.11 , 8.64812e-17 ],
                        [ 6483.85 , 8.56387e-17 ],
                        [ 6485.58 , 8.30909e-17 ],
                        [ 6487.32 , 8.13269e-17 ],
                        [ 6489.05 , 7.95901e-17 ],
                        [ 6490.79 , 8.54554e-17 ],
                        [ 6492.52 , 8.04148e-17 ],
                        [ 6494.26 , 7.956e-17 ],
                        [ 6495.99 , 8.30647e-17 ],
                        [ 6497.73 , 7.7904e-17 ],
                        [ 6499.46 , 8.87215e-17 ],
                        [ 6501.2 , 8.64422e-17 ],
                        [ 6502.93 , 8.16061e-17 ],
                        [ 6504.67 , 8.65604e-17 ],
                        [ 6506.4 , 8.31527e-17 ],
                        [ 6508.14 , 8.30319e-17 ],
                        [ 6509.87 , 8.63514e-17 ],
                        [ 6511.61 , 8.71934e-17 ],
                        [ 6513.34 , 9.19315e-17 ],
                        [ 6515.08 , 8.56164e-17 ],
                        [ 6516.81 , 8.59196e-17 ],
                        [ 6518.55 , 9.06681e-17 ],
                        [ 6520.28 , 8.68873e-17 ],
                        [ 6522.02 , 8.20436e-17 ],
                        [ 6523.75 , 8.25962e-17 ],
                        [ 6525.49 , 9.21053e-17 ],
                        [ 6527.22 , 9.26558e-17 ],
                        [ 6528.96 , 9.78202e-17 ],
                        [ 6530.69 , 8.98786e-17 ],
                        [ 6532.43 , 9.69294e-17 ],
                        [ 6534.16 , 9.45276e-17 ],
                        [ 6535.9 , 9.68097e-17 ],
                        [ 6537.63 , 9.49528e-17 ],
                        [ 6539.37 , 9.91305e-17 ],
                        [ 6541.1 , 1.08216e-16 ],
                        [ 6542.84 , 1.29092e-16 ],
                        [ 6544.57 , 2.15399e-16 ],
                        [ 6546.31 , 3.5338e-16 ],
                        [ 6548.04 , 3.83513e-16 ],
                        [ 6549.78 , 2.93025e-16 ],
                        [ 6551.52 , 2.42892e-16 ],
                        [ 6553.25 , 2.6806e-16 ],
                        [ 6554.99 , 4.25984e-16 ],
                        [ 6556.72 , 8.69919e-16 ],
                        [ 6558.46 , 3.12376e-15 ],
                        [ 6560.19 , 1.05061e-14 ],
                        [ 6561.93 , 1.72276e-14 ],
                        [ 6563.66 , 1.45295e-14 ],
                        [ 6565.4 , 7.38004e-15 ],
                        [ 6567.13 , 2.50114e-15 ],
                        [ 6568.87 , 8.43756e-16 ],
                        [ 6570.6 , 3.74166e-16 ],
                        [ 6572.34 , 2.25341e-16 ],
                        [ 6574.07 , 1.77545e-16 ],
                        [ 6575.81 , 1.53118e-16 ],
                        [ 6577.54 , 1.72322e-16 ],
                        [ 6579.28 , 3.17761e-16 ],
                        [ 6581.01 , 7.27118e-16 ],
                        [ 6582.75 , 9.71497e-16 ],
                        [ 6584.48 , 7.37653e-16 ],
                        [ 6586.22 , 4.04425e-16 ],
                        [ 6587.95 , 1.87279e-16 ],
                        [ 6589.69 , 1.24133e-16 ],
                        [ 6591.42 , 1.11212e-16 ],
                        [ 6593.16 , 9.87069e-17 ],
                        [ 6594.89 , 9.51901e-17 ],
                        [ 6596.63 , 8.85843e-17 ],
                        [ 6598.36 , 9.26498e-17 ],
                        [ 6600.1 , 9.4736e-17 ],
                        [ 6601.83 , 8.44449e-17 ],
                        [ 6603.57 , 9.11648e-17 ],
                        [ 6605.3 , 8.95719e-17 ],
                        [ 6607.04 , 8.81881e-17 ],
                        [ 6608.77 , 8.77689e-17 ],
                        [ 6610.51 , 8.31259e-17 ],
                        [ 6612.24 , 8.37382e-17 ],
                        [ 6613.98 , 7.99306e-17 ],
                        [ 6615.71 , 9.14558e-17 ],
                        [ 6617.45 , 8.54267e-17 ]])
 
    #The method returns the wavelength, flux, and spectrum regions        
    return Spectrum[:,0], Spectrum[:,1], Wavelength_regions, mu_Emissions
 
def Observational_data_single():
     
    #Since we know the emision line is Halpha we may declare the their locations
    mu_Emissions = np.array([5875.624])
     
    #Wavelengths marking the blue continuum, line region and red continuum
    Wavelength_regions = np.array([5815.48, 5835.232, 5866.90, 5883.494, 5896.930, 5927.779])
 
    #HeI line
    Spectrum = np.array([[5807.15966797, 9.68819344411e-17],
                        [5808.89501953, 9.63200074899e-17],
                        [5810.62988281, 1.03990076116e-16],
                        [5812.36523438, 1.03326908875e-16],
                        [5814.10009766, 9.72317656657e-17],
                        [5815.83544922, 9.83536144822e-17],
                        [5817.5703125, 9.42696186573e-17],
                        [5819.30566406, 1.0271498051e-16],
                        [5821.04052734, 1.03455492447e-16],
                        [5822.77587891, 9.65017953187e-17],
                        [5824.51074219, 1.03579238667e-16],
                        [5826.24609375, 9.95591938291e-17],
                        [5827.98095703, 1.02687498262e-16],
                        [5829.71630859, 9.84125957686e-17],
                        [5831.45117188, 9.43101041852e-17],
                        [5833.18603516, 1.02751813209e-16],
                        [5834.92138672, 1.00710510128e-16],
                        [5836.65625, 1.00407828198e-16],
                        [5838.39160156, 9.27246636487e-17],
                        [5840.12646484, 1.02179457164e-16],
                        [5841.86181641, 9.57215853301e-17],
                        [5843.59667969, 1.07000140043e-16],
                        [5845.33203125, 9.61419717523e-17],
                        [5847.06689453, 1.01047516744e-16],
                        [5848.80224609, 9.74931017997e-17],
                        [5850.53710938, 9.64077415744e-17],
                        [5852.27246094, 9.61284721647e-17],
                        [5854.00732422, 9.65234277461e-17],
                        [5855.74267578, 9.80389615947e-17],
                        [5857.47753906, 1.04452781098e-16],
                        [5859.21289062, 1.0083067631e-16],
                        [5860.94775391, 9.77636560345e-17],
                        [5862.68310547, 9.93312360872e-17],
                        [5864.41796875, 9.9555904959e-17],
                        [5866.15283203, 1.03696764488e-16],
                        [5867.88818359, 1.07675536322e-16],
                        [5869.62304688, 1.29282149049e-16],
                        [5871.35839844, 2.46428910001e-16],
                        [5873.09326172, 5.60902676573e-16],
                        [5874.82861328, 8.12969305522e-16],
                        [5876.56347656, 6.70266594404e-16],
                        [5878.29882812, 3.76000598734e-16],
                        [5880.03369141, 1.91817916565e-16],
                        [5881.76904297, 1.32638292109e-16],
                        [5883.50390625, 1.11387810414e-16],
                        [5885.23925781, 1.04052716849e-16],
                        [5886.97412109, 8.94234255287e-17],
                        [5888.70947266, 9.85098788261e-17],
                        [5890.44433594, 9.48601726751e-17],
                        [5892.1796875, 8.5499181441e-17],
                        [5893.91455078, 9.35529956967e-17],
                        [5895.64941406, 9.48250869822e-17],
                        [5897.38476562, 8.97595784947e-17],
                        [5899.11962891, 9.23991449166e-17],
                        [5900.85498047, 9.7875100424e-17],
                        [5902.58984375, 9.59906572572e-17],
                        [5904.32519531, 9.40512429756e-17],
                        [5906.06005859, 9.35506994433e-17],
                        [5907.79541016, 9.68653246544e-17],
                        [5909.53027344, 9.80831528917e-17],
                        [5911.265625, 9.80379755954e-17],
                        [5913.00048828, 9.09356903594e-17],
                        [5914.73583984, 9.27302421548e-17],
                        [5916.47070312, 9.34611587963e-17],
                        [5918.20605469, 8.8935117652e-17],
                        [5919.94091797, 9.18060299476e-17],
                        [5921.67626953, 9.70743498864e-17],
                        [5923.41113281, 9.41662012284e-17],
                        [5925.14599609, 9.15681394209e-17],
                        [5926.88134766, 9.50579614857e-17],
                        [5928.61621094, 9.89626311713e-17],
                        [5930.3515625, 9.61828080048e-17],
                        [5932.08642578, 9.07040334658e-17],
                        [5933.82177734, 9.42495677992e-17],
                        [5935.55664062, 9.7377594299e-17],
                        [5937.29199219, 9.01313002271e-17],
                        [5939.02685547, 9.48850741202e-17],
                        [5940.76220703, 9.72331222419e-17],
                        [5942.49707031, 9.14327928204e-17],
                        [5944.23242188, 9.07852890717e-17]])
 
    #The method returns the wavelength, flux, and spectrum regions        
    return Spectrum[:,0], Spectrum[:,1], Wavelength_regions, mu_Emissions
 
def Gaussian_Curve(A, mu, sigma, x, zerolev):
     
    return A * np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + zerolev
 
def CompositeModel_zerolev(params, x, zerolev, Ncomps):
 
    y_model = 0.0
    for i in range(Ncomps):
        index   = str(i)
        A       = params['A' + index] 
        mu      = params['mu' + index] 
        sigma   = params['sigma' + index] 
        y_model = y_model + A * np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))
         
    return y_model + zerolev
 
def CompResid_zerolev(params, x, y, zerolev, Ncomps, err):
 
    return (CompositeModel_zerolev(params.valuesdict(), x, zerolev, Ncomps) - y) / err
 
def Load_lmfit_parameters_blended(N_comps, Initial_guesses_dic, wide_component = False, mu_precission = 1):
 
    N2 = Atom('N', 2)
    N2_6548A = N2.getEmissivity(tem=10000, den=100, wave=6548)
    N2_6584A = N2.getEmissivity(tem=10000, den=100, wave=6584)
     
    print 'Radio', N2_6584A / N2_6548A
     
    params = Parameters()
 
    for i in range(N_comps):
        index = str(i)
        params.add('A'      + index, value = Initial_guesses_dic['A'][i])
        params.add('mu'     + index, value = Initial_guesses_dic['mu'][i], min = Initial_guesses_dic['mu'][i] - mu_precission, max = Initial_guesses_dic['mu'][i] + mu_precission) #One angstrom tolerance for min and max value of mu
        params.add('sigma'  + index, value = Initial_guesses_dic['sigma'][i], min = Initial_guesses_dic['min_sigma'][i], max = 5.0)
        params.add('fwhm'   + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
        params.add('FWHM'   + index, expr = '({fwhm}/{mu}) * 2.99792458e5'.format(fwhm = 'fwhm' + index, mu = 'mu' + index))
        params.add('Flux'   + index, expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + index, fwhm = 'fwhm' + index))
 
    #Set all sigmas to the same value
    for i in range(1, N_comps):
        index = str(i)
        params['sigma'      + index].expr='sigma0'
# 
    #Set flux in N2_6548A as 1/3 N2_6584A
    params.add('Flux0', expr = 'Flux2 / {ratio}'.format(ratio = N2_6584A / N2_6548A))
 
    if wide_component:
        index = str(N_comps)
        params.add('A'      + index, value = Initial_guesses_dic['A'][N_comps])
        params.add('mu'     + index, expr = 'mu1')
        params.add('sigma'  + index, value = Initial_guesses_dic['sigma'][N_comps], min = Initial_guesses_dic['min_sigma'][N_comps], max = 8.0)
        params.add('fwhm'   + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
        params.add('FWHM'   + index, expr = '({fwhm}/{mu}) * 2.99792458e5'.format(fwhm = 'fwhm' + index, mu = 'mu' + index))
        params.add('Flux'   + index, expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + index, fwhm = 'fwhm' + index))
 
#     if wide_component:
#         index = str(N_comps)
#         print 'Y la europea', 'A'      + index
#         params.add('A'      + index, value = Initial_guesses_dic['A'][N_comps], min=0.04748900 * 0.95, max = 0.04748900 * 1.05)
#         params.add('mu'     + index, value = Initial_guesses_dic['mu'][N_comps], min = Initial_guesses_dic['mu'][N_comps] - 2, max = Initial_guesses_dic['mu'][N_comps] + 2)
#         params.add('sigma'  + index, value = Initial_guesses_dic['sigma'][N_comps], min = Initial_guesses_dic['min_sigma'][N_comps], max = 8.0)
#         params.add('fwhm'   + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
#         params.add('FWHM'   + index, expr = '({fwhm}/{mu}) * 2.99792458e5'.format(fwhm = 'fwhm' + index, mu = 'mu' + index))
#         params.add('Flux'   + index, expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + index, fwhm = 'fwhm' + index))
  
#     if wide_component:
#         index = str(N_comps + 1)
#         params.add('sigm_tol', value = 2.5, min = 2, max = 4)
#         params.add('A'      + index, value = Initial_guesses_dic['A'][N_comps])
#         params.add('mu'     + index, expr = 'mu2')
#         params.add('sigma'  + index, expr = 'sigm_tol*sigma2')
#         params.add('sigma'  + index, value = Initial_guesses_dic['sigma'][N_comps], min = Initial_guesses_dic['min_sigma'][N_comps], max = 8.0)
#         params.add('fwhm'   + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
#         params.add('FWHM'   + index, expr = '({fwhm}/{mu}) * 2.99792458e5'.format(fwhm = 'fwhm' + index, mu = 'mu' + index))
#         params.add('Flux'   + index, expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + index, fwhm = 'fwhm' + index))
 
    return params
 
def Load_lmfit_parameters_single(N_comps, Initial_guesses_dic, wide_component = False, mu_precission = 1):
     
    params = Parameters()
 
    for i in range(N_comps):
        index = str(i)
        params.add('A'      + index, value = Initial_guesses_dic['A'][i])
        params.add('mu'     + index, value = Initial_guesses_dic['mu'][i]) #One angstrom tolerance for min and max value of mu
        params.add('sigma'  + index, value = Initial_guesses_dic['sigma'][i], min = Initial_guesses_dic['min_sigma'][i], max = 5.0)
        params.add('fwhm'   + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
        params.add('FWHM'   + index, expr = '({fwhm}/{mu}) * 2.99792458e5'.format(fwhm = 'fwhm' + index, mu = 'mu' + index))
        params.add('Flux'   + index, expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + index, fwhm = 'fwhm' + index))
 
    if wide_component:
        index = str(N_comps + 1)
        params.add('A'      + index, value = Initial_guesses_dic['A'][N_comps], min=0.04748900 * 0.95, max = 0.04748900 * 1.05)
        params.add('mu'     + index, value = Initial_guesses_dic['mu'][N_comps], min = Initial_guesses_dic['mu'][N_comps] - 2, max = Initial_guesses_dic['mu'][N_comps] + 2)
        params.add('sigma'  + index, value = Initial_guesses_dic['sigma'][N_comps], min = Initial_guesses_dic['min_sigma'][N_comps], max = 8.0)
        params.add('fwhm'   + index, expr = '2.354820045 * {sigma}'.format(sigma = 'sigma'  + index))
        params.add('FWHM'   + index, expr = '({fwhm}/{mu}) * 2.99792458e5'.format(fwhm = 'fwhm' + index, mu = 'mu' + index))
        params.add('Flux'   + index, expr = '({A}*{fwhm})/(2.35*0.3989)'.format(A = 'A'  + index, fwhm = 'fwhm' + index))
 
    return params
 
def rescale_parameters(fit_output, x_scale, y_scale, N_comps, wide_component = False):
     
#     scale_params = OrderedDict()
    scale_params = Parameters()
    scale_params = OrderedDict()
    if wide_component:
        N_comps = N_comps + 1
 
    print 'antes'
    print fit_output['A0'].value, fit_output['A0'].stderr
    print fit_output['mu0'].value, fit_output['mu0'].stderr
    print fit_output['sigma0'].value, fit_output['sigma0'].stderr
    print fit_output['Flux0'].value, fit_output['Flux0'].stderr
 
    for i in range(N_comps):
        index = str(i)
        scale_params.add('A'    + index, value = fit_output['A'+ index] * y_scale)
        scale_params.add('mu'   + index, value = fit_output['mu'+ index] + x_scale)
        scale_params.add('sigma'+ index, value = fit_output['sigma'+ index])
        scale_params.add('fwhm' + index, value = fit_output['fwhm'+ index])
        scale_params.add('FWHM' + index, value = scale_params['fwhm'+ index] / scale_params['mu'+ index] * 2.99792458e5)
        scale_params.add('Flux' + index, value = (scale_params['A'+index] * scale_params['fwhm'+index]) / (2.35 * 0.3989))
 
    print 'despues'
    print fit_output['A0'].value, fit_output['A0'].stderr
    print fit_output['mu0'].value, fit_output['mu0'].stderr
    print fit_output['sigma0'].value, fit_output['sigma0'].stderr
    print fit_output['Flux0'].value, fit_output['Flux0'].stderr
     
    return scale_params
 
def rescale_parameters2(fit_output, x_scale, y_scale, N_comps, wide_component = False):
     
    scale_params = OrderedDict()
 
    for i in range(N_comps + 1):
        index = str(i)
        scale_params['A'+ index]     = ufloat(fit_output['A'+ index].value, fit_output['A'+ index].stderr) * y_scale
        scale_params['mu'+ index]    = ufloat(fit_output['mu'+ index].value, fit_output['mu'+ index].stderr) + x_scale
        scale_params['sigma'+ index] = ufloat(fit_output['sigma'+ index].value, fit_output['sigma'+ index].stderr)
        scale_params['fwhm'+ index]  = ufloat(fit_output['fwhm'+ index].value, fit_output['fwhm'+ index].stderr)
        scale_params['FWHM'+ index]  = ufloat(fit_output['FWHM'+ index].value, fit_output['FWHM'+ index].stderr)
        scale_params['Flux'+ index]  = ufloat(fit_output['Flux'+ index].value, fit_output['Flux'+ index].stderr) * y_scale
        
    return scale_params
 
# # #Import data
# Wavelength, Flux, Spectrum_regions, mu_Emissions = Observational_data_single() 
#  
# #Normalize the data
# Max_Flux        = np.max(Flux)
# Max_wavelength  = Wavelength[np.argmax(Flux)] 
# WavelengthN, FluxN, Spectrum_regionsN, mu_EmissionsN = Wavelength - Max_wavelength, Flux/Max_Flux, Spectrum_regions - Max_wavelength, mu_Emissions 
# 
# line_region             = (Wavelength >= Spectrum_regions[2]) & (Wavelength <= Spectrum_regions[3])
# blue_region, red_region = (Wavelength >= Spectrum_regions[0]) & (Wavelength <= Spectrum_regions[1]), (Wavelength >= Spectrum_regions[4]) & (Wavelength <= Spectrum_regions[5])
# line_wave, line_flux    = Wavelength[line_region], Flux[line_region]
# blue_wave, blue_flux    = Wavelength[blue_region], Flux[blue_region]
# red_wave, red_flux      = Wavelength[red_region], Flux[red_region]
# 
# line_regionN              = (WavelengthN >= Spectrum_regionsN[2]) & (WavelengthN <= Spectrum_regionsN[3])
# blue_regionN, red_regionN = (WavelengthN >= Spectrum_regionsN[0]) & (WavelengthN <= Spectrum_regionsN[1]), (WavelengthN >= Spectrum_regionsN[4]) & (WavelengthN <= Spectrum_regionsN[5])
# line_waveN, line_fluxN    = WavelengthN[line_region], FluxN[line_region]
# blue_waveN, blue_fluxN    = WavelengthN[blue_region], FluxN[blue_region]
# red_waveN, red_fluxN      = WavelengthN[red_region], FluxN[red_region]
# 
# # #Declaring a linear continuum uppon which the line is located
# # lineal_mod                      = LinearModel(prefix='lineal_')
# # Continuum_wave, Continuum_flux  = np.hstack([blue_wave, red_wave]), np.hstack([blue_flux, red_flux])
# # Lineal_parameters               = lineal_mod.guess(Continuum_flux, x=Continuum_wave)
# # lineal_zerolev                  = Lineal_parameters['lineal_slope'].value * line_wave + Lineal_parameters['lineal_intercept'].value
# # err_continuum                   = np.std(Lineal_parameters['lineal_slope'].value * Continuum_wave + Lineal_parameters['lineal_intercept'].value - Continuum_flux)
# 
# #Declaring a linear continuum uppon which the line is located
# lineal_modN                      = LinearModel(prefix='lineal_')
# Continuum_waveN, Continuum_fluxN = np.hstack([blue_waveN, red_waveN]), np.hstack([blue_fluxN, red_fluxN])
# Lineal_parametersN               = lineal_modN.guess(Continuum_fluxN, x=Continuum_waveN)
# lineal_zerolevN                  = Lineal_parametersN['lineal_slope'].value * line_waveN + Lineal_parametersN['lineal_intercept'].value
# err_continuumN                   = np.std(Lineal_parametersN['lineal_slope'].value * Continuum_waveN + Lineal_parametersN['lineal_intercept'].value - Continuum_fluxN)
# 
# print 'Trapz', (trapz(line_fluxN, line_waveN) - trapz(lineal_zerolevN, line_waveN)) * Max_Flux
# print 'simps', (simps(line_fluxN, line_waveN) - simps(lineal_zerolevN, line_waveN)) * Max_Flux
# 
# print 'Trapz', (trapz(line_fluxN) - trapz(lineal_zerolevN)) * Max_Flux
# print 'simps', (simps(line_fluxN) - simps(lineal_zerolevN)) * Max_Flux

plt.style.use('ggplot')
font = {'size'   : 16}
  
rc('font', **font)
  
#Import data
Wavelength, Flux, Spectrum_regions, mu_Emissions = Observational_data_blended()
  
#Normalize the data
Max_Flux        = np.max(Flux)
Max_wavelength  = Wavelength[np.argmax(Flux)] 
Wavelength, Flux, Spectrum_regions, mu_Emissions = Wavelength - Max_wavelength, Flux/Max_Flux, Spectrum_regions - Max_wavelength, mu_Emissions 
  
#Variables to improve the initial guesses
Halpha_Peak     = np.max(Flux)
mu_precission   = 3 #Maximun tolerance on the lines center (Angtroms)
  
#We stablish the continuum blue and red regions (for the linear component calculation) and we stablish the line region (for the gaussian components calculation)
line_region             = (Wavelength >= Spectrum_regions[2]) & (Wavelength <= Spectrum_regions[3])
blue_region, red_region = (Wavelength >= Spectrum_regions[0]) & (Wavelength <= Spectrum_regions[1]), (Wavelength >= Spectrum_regions[4]) & (Wavelength <= Spectrum_regions[5])
line_wave, line_flux    = Wavelength[line_region], Flux[line_region]
blue_wave, blue_flux    = Wavelength[blue_region], Flux[blue_region]
red_wave, red_flux      = Wavelength[red_region], Flux[red_region]

# for local maxima
max_index = argrelextrema(line_flux, np.greater)[0]
peak_fluxes = line_flux[max_index]
peak_waves  = line_wave[max_index]
  
# for local minima
min_index = argrelextrema(line_flux, np.less)[0]
Minima_fluxes = line_flux[min_index]
Minima_waves  = line_wave[min_index]
  
# for local minima
min_index = argrelextrema(line_flux, np.less)[0]
Minima_fluxes = line_flux[min_index]
Minima_waves  = line_wave[min_index]
 
#Lmfit parameters
Ncomps = 3
Initial_guesses_dic                     = OrderedDict()
Initial_guesses_dic['A']                = np.array([peak_fluxes[0], peak_fluxes[1], peak_fluxes[2], peak_fluxes[1]/20])
Initial_guesses_dic['mu']               = np.array([peak_waves[0], peak_waves[1], peak_waves[2], peak_waves[1]])
Initial_guesses_dic['sigma']            = np.array([1.0, 1.0, 1.0, 5.0])
Initial_guesses_dic['min_sigma']        = np.zeros(Ncomps + 1)
params = Load_lmfit_parameters_blended(Ncomps, Initial_guesses_dic, wide_component = True, mu_precission = mu_precission)
    
#Declaring a linear continuum uppon which the line is located
lineal_mod                      = LinearModel(prefix='lineal_')
Continuum_wave, Continuum_flux  = np.hstack([blue_wave, red_wave]), np.hstack([blue_flux, red_flux])
Lineal_parameters               = lineal_mod.guess(Continuum_flux, x=Continuum_wave)
lineal_zerolev                  = Lineal_parameters['lineal_slope'].value * line_wave + Lineal_parameters['lineal_intercept'].value
err_continuum                   = np.std(Lineal_parameters['lineal_slope'].value * Continuum_wave + Lineal_parameters['lineal_intercept'].value - Continuum_flux)
  
Lineal_parameters_scale         = lineal_mod.guess(Continuum_flux * Max_Flux, x=Continuum_wave + Max_wavelength)
lineal_zerolev_scale            = Lineal_parameters['lineal_slope'].value * line_wave + Lineal_parameters['lineal_intercept'].value


#Make the fitting
out = minimize(CompResid_zerolev, params, args=(line_wave, line_flux, lineal_zerolev, Ncomps, err_continuum))
  
#Rescale the parameters
scale_params        = rescale_parameters2(out.params, Max_wavelength, Max_Flux, Ncomps, wide_component = False)
scale_params_values = ufloatDict_nominal(scale_params)
  
#Make the plots
x_resample      = np.linspace(line_wave[0], line_wave[-1], 100) + Max_wavelength
lineal_resample = Lineal_parameters_scale['lineal_slope'].value * x_resample + Lineal_parameters_scale['lineal_intercept'].value 
Gauss_mix_curve = CompositeModel_zerolev(scale_params_values, x_resample, lineal_resample, Ncomps)

plt.step(Wavelength + Max_wavelength, Flux * Max_Flux, color= 'black', label = "Complete spectrum",  where='mid')
plt.plot(x_resample, Gauss_mix_curve, 'r-', label = 'Fitted line')
plt.plot(line_wave[max_index] + Max_wavelength, line_flux[max_index] * Max_Flux, 'o', color = 'purple', label = 'Maxima')
plt.plot(line_wave[min_index] + Max_wavelength, line_flux[min_index] * Max_Flux, 'o', color = 'orange', label = 'Minima')
plt.plot(line_wave + Max_wavelength, lineal_zerolev * Max_Flux, 'g:', label = 'Fitted continuum')

SubWave, SubInt, EM_Height, EM_ExpLoc = Emission_Threshold(6563.0, Wavelength + Max_wavelength, Flux * Max_Flux)
        
OverallMedian = np.median(Flux * Max_Flux)

plt.xlim([SubWave[0],SubWave[-1]])
plt.ylim(OverallMedian/10, EM_Height*2)

Gaussian_line = None
for i in range(Ncomps + 1):       
    index   = str(i)
    A       = scale_params_values['A' + index]
    mu      = scale_params_values['mu' + index] 
    sigma   = scale_params_values['sigma' + index] 
    Gaussian_line = Gaussian_Curve(A, mu, sigma, x_resample, lineal_resample)
    plt.plot(x_resample, Gaussian_line, '--', color = 'blue', label = 'Gaussian components' if i == 0 else "")
          
plt.legend(loc = 'best',numpoints=1)
plt.title(r'$H\alpha$ gaussian mixture fitting')
plt.xlabel(r'Wavelength $(\AA)$')
plt.ylabel(r'Flux $(erg/s/cm^2)$')
   
print '--Area simps integration'
line_waveS, line_fluxS = line_wave + Max_wavelength, line_flux * Max_Flux
y_continuum     = Lineal_parameters_scale['lineal_slope'].value * line_wave + Lineal_parameters_scale['lineal_intercept'].value
y_continuumS    = Lineal_parameters_scale['lineal_slope'].value * line_waveS + Lineal_parameters_scale['lineal_intercept'].value
 
print 'Points',             simps(line_fluxS, line_waveS) - simps(y_continuumS, line_waveS)
print 'Gaussian points',    simps(Gauss_mix_curve, x_resample) - simps(lineal_resample, x_resample)
print 'Fit flux',           (scale_params['Flux0'] + scale_params['Flux1'] + scale_params['Flux2']).nominal_value

plt.show()


#---Single line
# plt.style.use('ggplot')
# font = {'size'   : 16}
#   
# rc('font', **font)
#   
# #Import data
# Wavelength, Flux, Spectrum_regions, mu_Emissions = Observational_data_single()
#   
# #Normalize the data
# Max_Flux        = np.max(Flux)
# Max_wavelength  = Wavelength[np.argmax(Flux)] 
# Wavelength, Flux, Spectrum_regions, mu_Emissions = Wavelength - Max_wavelength, Flux/Max_Flux, Spectrum_regions - Max_wavelength, mu_Emissions 
#   
# #Variables to improve the initial guesses
# Halpha_Peak     = np.max(Flux)
# mu_precission   = 2 #Maximun tolerance on the lines center (Angtroms)
#   
# #We stablish the continuum blue and red regions (for the linear component calculation) and we stablish the line region (for the gaussian components calculation)
# line_region             = (Wavelength >= Spectrum_regions[2]) & (Wavelength <= Spectrum_regions[3])
# blue_region, red_region = (Wavelength >= Spectrum_regions[0]) & (Wavelength <= Spectrum_regions[1]), (Wavelength >= Spectrum_regions[4]) & (Wavelength <= Spectrum_regions[5])
# line_wave, line_flux    = Wavelength[line_region], Flux[line_region]
# blue_wave, blue_flux    = Wavelength[blue_region], Flux[blue_region]
# red_wave, red_flux      = Wavelength[red_region], Flux[red_region]
#   
# # for local maxima
# max_index = argrelextrema(line_flux, np.greater)[0]
# peak_fluxes = line_flux[max_index]
# peak_waves  = line_wave[max_index]
#   
# # for local minima
# min_index = argrelextrema(line_flux, np.less)[0]
# Minima_fluxes = line_flux[min_index]
# Minima_waves  = line_wave[min_index]
#   
# # for local minima
# min_index = argrelextrema(line_flux, np.less)[0]
# Minima_fluxes = line_flux[min_index]
# Minima_waves  = line_wave[min_index]
#  
# # #Lmfit parameters
# # Ncomps = 3
# # Initial_guesses_dic                     = OrderedDict()
# # Initial_guesses_dic['A']                = np.array([peak_fluxes[0], peak_fluxes[1], peak_fluxes[2], peak_fluxes[1]/20])
# # Initial_guesses_dic['mu']               = np.array([peak_waves[0], peak_waves[1], peak_waves[2], peak_waves[1]])
# # Initial_guesses_dic['sigma']            = np.array([1.0, 1.0, 1.0, 5.0])
# # Initial_guesses_dic['min_sigma']        = np.zeros(Ncomps + 1)
# # params = Load_lmfit_parameters(Ncomps, Initial_guesses_dic, wide_component = True, mu_precission = mu_precission)
#   
# #Lmfit parameters
# Ncomps = 1
# Initial_guesses_dic                     = OrderedDict()
# Initial_guesses_dic['A']                = np.array([peak_fluxes[0]])
# Initial_guesses_dic['mu']               = np.array([peak_waves[0]])
# Initial_guesses_dic['sigma']            = np.array([1.0])
# Initial_guesses_dic['min_sigma']        = np.zeros(Ncomps)
#   
# params = Load_lmfit_parameters_single(Ncomps, Initial_guesses_dic, wide_component = False, mu_precission = mu_precission)
# 
# print 'Estos parametros', params.values()
# 
# #Declaring a linear continuum uppon which the line is located
# lineal_mod                      = LinearModel(prefix='lineal_')
# Continuum_wave, Continuum_flux  = np.hstack([blue_wave, red_wave]), np.hstack([blue_flux, red_flux])
# Lineal_parameters               = lineal_mod.guess(Continuum_flux, x=Continuum_wave)
# lineal_zerolev                  = Lineal_parameters['lineal_slope'].value * line_wave + Lineal_parameters['lineal_intercept'].value
# err_continuum                   = np.std(Lineal_parameters['lineal_slope'].value * Continuum_wave + Lineal_parameters['lineal_intercept'].value - Continuum_flux)
#   
# Lineal_parameters_scale         = lineal_mod.guess(Continuum_flux * Max_Flux, x=Continuum_wave + Max_wavelength)
# lineal_zerolev_scale            = Lineal_parameters['lineal_slope'].value * line_wave + Lineal_parameters['lineal_intercept'].value
#   
# #Make the fitting
# out = minimize(CompResid_zerolev, params, args=(line_wave, line_flux, lineal_zerolev, Ncomps, err_continuum))
# # report_fit(out.params)
#   
# #Rescale the parameters
# scale_params        = rescale_parameters2(out.params, Max_wavelength, Max_Flux, Ncomps, wide_component = False)
# scale_params_values = ufloatDict_nominal(scale_params)
#  
# #Make the plots
# x_resample      = np.linspace(line_wave[0], line_wave[-1], 100) + Max_wavelength
# plt.plot(Wavelength + Max_wavelength, Flux * Max_Flux,  'o', color= 'black', label = "Complete spectrum")
#  
# lineal_resample = Lineal_parameters_scale['lineal_slope'].value * x_resample + Lineal_parameters_scale['lineal_intercept'].value 
# plt.plot(x_resample, CompositeModel_zerolev(scale_params_values, x_resample, lineal_resample, Ncomps), 'r-', label = 'Fitted line')
# plt.plot(line_wave[max_index] + Max_wavelength, line_flux[max_index] * Max_Flux, 'o', color = 'purple', label = 'Maxima')
# plt.plot(line_wave[min_index] + Max_wavelength, line_flux[min_index] * Max_Flux, 'o', color = 'orange', label = 'Minima')
# plt.plot(line_wave + Max_wavelength, lineal_zerolev * Max_Flux, 'g:', label = 'Fitted continuum')
#   
# Gaussian_line = None
# for i in range(Ncomps):    
#     index   = str(i)
#     A       = scale_params_values['A' + index]
#     mu      = scale_params_values['mu' + index] 
#     sigma   = scale_params_values['sigma' + index] 
#     Gaussian_line = Gaussian_Curve(A, mu, sigma, x_resample, lineal_resample)
#     plt.plot(x_resample, Gaussian_line, '--', color = 'blue', label = 'Gaussian components' if i == 0 else "")
#           
# plt.legend(loc = 'best',numpoints=1)
# plt.title(r'$H\alpha$ gaussian mixture fitting')
# plt.xlabel(r'Wavelength $(\AA)$')
# plt.ylabel(r'Flux $(erg/s/cm^2)$')
#    
# print '--Area simps integration'
# line_waveS, line_fluxS = line_wave + Max_wavelength, line_flux * Max_Flux
# y_continuum     = Lineal_parameters_scale['lineal_slope'].value * line_wave + Lineal_parameters_scale['lineal_intercept'].value
# y_continuumS    = Lineal_parameters_scale['lineal_slope'].value * line_waveS + Lineal_parameters_scale['lineal_intercept'].value
#  
# print 'Points',             simps(line_fluxS, line_waveS) - simps(y_continuumS, line_waveS)
# print 'Gaussian points',    simps(Gaussian_line, x_resample) - simps(lineal_resample, x_resample)
# print 'Fit flux',           repr(scale_params['Flux0'])
# 
# # print '--Area trapz integration'
# # print 'Points',     trapz(line_fluxS, line_waveS) - trapz(y_continuumS, line_waveS)
# # print 'Gaussian',   trapz(Gaussian_line, x_resample) - trapz(lineal_resample, x_resample)
# #  
# # print '--Integration scaled and rescale'
# # print 'Trapz',      (trapz(line_flux, line_wave) - trapz(y_continuum, line_wave)) * Max_Flux
# # print 'simpz',      (simps(line_flux, line_wave) - simps(y_continuum, line_wave)) * Max_Flux
# #  
# # print '-- Lmfit output values'
# # print 'A0'     ,    repr(scale_params['A0']) 
# # print 'mu0'    ,    repr(scale_params['mu0'])    
# # print 'sigma0' ,    repr(scale_params['sigma0']) 
# # print 'Flux0' ,     repr(scale_params['Flux0'])
# #   
# # #Display graphs
# plt.show()
