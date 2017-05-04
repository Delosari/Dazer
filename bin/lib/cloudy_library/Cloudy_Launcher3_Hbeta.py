from collections                        import OrderedDict

from numpy                              import log10 as nplog10, zeros, min, max, linspace, array, concatenate, isfinite, greater
from pandas                             import Series

from Math_Libraries.bces_script         import bces
from Plotting_Libraries.dazer_plotter   import Plot_Conf
from cloudy_library.cloudy_methods      import Cloudy_Tools


dz              = Plot_Conf()
ct              = Cloudy_Tools()

#Define script name and location
ScriptFolder    = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Ionization_Models_Hbeta_trans/'

#Orion nebula has a metallicity of 2/3 solar

Grid_Values = OrderedDict()
Grid_Values['age']      = ['5.0', '5.5', '6.0', '6.5', '7.0', '7.5'] 
Grid_Values['zStars']   = ['-2.4', '-2.1', '-1.7', '-1.31'] 
Grid_Values['zGas']     = ['0.1', '0.31', '0.62'] 
Grid_Values['u']        = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5'] 

#Dictionary of dictionaries
Grid_frame =  ({k : Series(v) for k, v in Grid_Values.iteritems()})

#Trick to create a frame with different lengths
# Grid_frame = DataFrame({k : Series(v) for k, v in Grid_Values.iteritems()})

#Generate the scripts with the lines we want to print the flux
ct.lines_to_extract(ScriptFolder)

#Fore the each grid point run a cloudy simulation
Model_dict = OrderedDict()
for age in Grid_Values['age']:
    for zStar in Grid_Values['zStars']:
        for zGas in Grid_Values['zGas']:                                        
            for ionization_parameter in Grid_Values['u']: 
                
                Model_dict['Name']      = 'S_Ar_Test_' + 'age'+ age + '_zStar' + zStar + '_zGas' + zGas + '_u' + ionization_parameter
                Model_dict['age']       = age
                Model_dict['zStars']    = zStar
                Model_dict['zGas']      = zGas
                Model_dict['u']         = ionization_parameter
    
                ScriptName              = Model_dict['Name'] + '.in'
    
                #Generate the script
                ct.popstar_S_Ar_Grid_Hbeta(ScriptFolder, Model_dict)
                  
                #Run the cloudy script
                ct.run_script(ScriptName, ScriptFolder)
  
print 'Data treated'