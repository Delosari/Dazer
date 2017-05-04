from DZ_DataExplorer                import Plots_Manager


# We declare the folder and log file to drop the lines data
dz = Plots_Manager()
  
# Forcing the remake of new log file
dz.RemakeFiles = True
    
#Identify spectrum (It assumes it is on the same folder of the code)
FileFolder, FileName = '', 'objSHOC579_WHT.fits'
  
#Load fits file data
Wave, Int, Header = dz.File2Data(FileFolder, FileName)

#Six wavelengths which mark the blue continuum, the line, and the red continuum (example for OIII 5007Angstroms)
Wavelength_points = [4967.758991, 4984.263585, 4999.928963, 5012.657082, 5030.420501, 5046.505486]

#Using a direct gausshermite
results_dict = dz.Command_LineMesuring(Wave, Int, Wavelength_points, Measuring_Method = 'kmpfit_GaussHermite')

#Using a MC gauss-hermite
results_dict_MCMC = dz.Command_LineMesuring(Wave, Int, Wavelength_points, Measuring_Method = 'kmpfit_GaussHermite_MCMC')

#Dictionary with the results from the analysis
print '\nResults example1'
for key in results_dict:
    print key, results_dict[key]

print '\nResults example2'
#Dictionary with the results from the analysis
for key in results_dict_MCMC:
    print key, results_dict_MCMC[key]