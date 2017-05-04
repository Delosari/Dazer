from scipy.integrate    import simps
from dazer_methods import Dazer
from numpy import median, random, log as np_ln

#Declare object
dz = Dazer()

#Define operation
Catalogue_Dic   = dz.import_catalogue()
Pattern         = '51959-092_WHT.fits'

#Locate files on hard drive
FilesList = dz.Folder_Explorer(Pattern, Catalogue_Dic['Obj_Folder'], CheckComputer=False)

#Generate plot frame and colors
dz.FigConf(n_colors=10)

for i in range(len(FilesList)):

    #Load the file
    CodeName, FileName, FileFolder  = dz.Analyze_Address(FilesList[i])
    Wave, Flux, ExtraData           = dz.File_to_data(FileFolder, FileName)

    #Get the emission line location
    Line_label, line_wavelength     = 'N2_6548A', 6548.050
    Line_selections                 = dz.RangesR([], line_wavelength, FileFolder + CodeName + '_WHT_LinesLog_v3.txt')
    print 'selections', Line_selections
      
    Line_selections                 = [6485.44, 6515.86, 6535.03, 6593.01, 6601.26, 6616.07]
    subWave, subFlux, line_heith, line_exploc = dz.Emission_Threshold(5016.0, Wave, Flux)
    
    #Fit the line
    Fitting_single = dz.command_LineMesuring(subWave, subFlux, Line_selections, line_wavelength, Line_label, Measuring_Method = 'lmfit', store_data = False)
    
    #Plotting data
    dz.data_plot(Wave, Flux, label = CodeName + ' spectrum', color=dz.ColorVector[2][0], linestyle='step')
    dz.data_plot(Wave[dz.ind3:dz.ind4], Fitting_single['ContinuumFlux'][dz.ind3:dz.ind4], label = CodeName + ' spectrum', color=dz.ColorVector[2][0], linestyle='step')
    dz.data_plot(Fitting_single['x_resample'], Fitting_single['y_resample'][0], label = 'Simple fitting', color=dz.ColorVector[2][1])
    dz.data_plot(Fitting_single['x_resample'], Fitting_single['zerolev_resample'], label = 'simple zero level', linestyle= ':', color=dz.ColorVector[2][1])
    dz.data_plot(subWave[dz.ind3:dz.ind4][dz.Fitting_dict['wide mask']], subFlux[dz.ind3:dz.ind4][dz.Fitting_dict['wide mask']], label = CodeName + ' original', color='red', markerstyle='o')
    dz.data_plot(dz.Fitting_dict['Maxima_Waves'], dz.Fitting_dict['Maxima_peaks'], label = 'Maxima peaks = ' + str(len(dz.Fitting_dict['Maxima_peaks'])), color='orange', markerstyle='o')
#     dz.data_plot(dz.Fitting_dict['Maxima_Waves_pu'], dz.Fitting_dict['Maxima_peaks_pu'], label = 'Maxima peaks = ' + str(len(dz.Fitting_dict['Maxima_Waves_pu'])), color='black', markerstyle='o')
    #dz.data_plot(subWave[dz.ind3:dz.ind4], dz.Fitting_dict['emis_limpio'] * dz.Fitting_dict['y_scaler'], label = CodeName + ' original', color='green', markerstyle='o')

#     #Plot the masks regions
#     sigmas_limits = [1.5,4,3]
#     for k in [0,1,2]:
#         index_str = str(k)
#         limit = dz.Fitting_dict['sigma'+index_str].nominal_value * sigmas_limits[k]
#         dz.Axis.axvspan(dz.Fitting_dict['mu'+index_str].nominal_value - limit, dz.Fitting_dict['mu'+index_str].nominal_value + limit, facecolor = 'yellow', alpha=0.5)
# 
#     for z in [1,2,3,4]:
#         dz.data_plot(Fitting_single['x_resample'], Fitting_single['y_resample'][z], label = 'Simple fitting', color='red', linestyle='--')
    
    suma = 0  
    for i in range(dz.Fitting_dict['blended number']):
        index = str(i)
        print 'FluxG' + index, dz.Fitting_dict['FluxG' + index]
        suma += dz.Fitting_dict['FluxG' + index]
     
    print 'Total gaussian', suma
    print 'integrated',  dz.Fitting_dict['FluxI']
#     print 'Special relation', dz.Fitting_dict['FluxG2'] / dz.Fitting_dict['FluxG0']
#     print 'For N2_6548A: ', dz.Fitting_dict['A0'].nominal_value * dz.Fitting_dict['sigma0'].nominal_value * dz.sqrt2pi ,dz.Fitting_dict['A0'].nominal_value, dz.Fitting_dict['mu0'].nominal_value, dz.Fitting_dict['sigma0'].nominal_value
#     print 'For N2_6584A: ', dz.Fitting_dict['A2'].nominal_value * dz.Fitting_dict['sigma2'].nominal_value * dz.sqrt2pi ,dz.Fitting_dict['A2'].nominal_value, dz.Fitting_dict['mu2'].nominal_value, dz.Fitting_dict['sigma2'].nominal_value
#     print 'Obs, ratio', dz.Fitting_dict['A2'].nominal_value * dz.Fitting_dict['sigma2'].nominal_value / (dz.Fitting_dict['A0'].nominal_value * dz.Fitting_dict['sigma0'].nominal_value)
#     

    #Repeating the integration
    FluxI   = simps(subFlux[dz.ind3:dz.ind4], subWave[dz.ind3:dz.ind4]) - simps(Fitting_single['ContinuumFlux'][dz.ind3:dz.ind4], subWave[dz.ind3:dz.ind4])
    FluxG   = simps(Fitting_single['y_resample'][0], Fitting_single['x_resample']) - simps(Fitting_single['zerolev_resample'], Fitting_single['x_resample'])
       
    #Plot format
    dz.Axis.set_xlim(subWave[0], subWave[-1])
    dz.Axis.set_ylim(median(subFlux)/10, line_heith * 2)
    dz.FigWording(r'Wavelength $(\AA)$', 'Flux ' + r'$(erg\,cm^{-2} s^{-1} \AA^{-1})$', 'Gaussian fitting')
    dz.display_fig()
    dz.save_manager('/home/vital/Desktop/Testing_Halpha/' + CodeName + 'O3_5007')

print 'data treated'