from numpy import log10 as nplog10
from collections                        import OrderedDict
from pandas                             import Series
from cloudy_library.cloudy_methods      import Cloudy_Tools

ct              = Cloudy_Tools()

#Define script name and location
ScriptFolder    = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Total_Grid_SAL2/'

#Orion nebula has a metallicity of 2/3 solar
Grid_Values         = OrderedDict()

Solar_abundances    = OrderedDict()
Solar_abundances['Li']  = -8.69
Solar_abundances['Be']  = -10.58
Solar_abundances['B']   = -9.21
Solar_abundances['C']   = -3.61
Solar_abundances['N']   = -4.07
Solar_abundances['O']   = -3.31
Solar_abundances['Ne']  = -4.00
Solar_abundances['F']   = -7.52
Solar_abundances['Ne']  = -4.0
Solar_abundances['Na']  = -5.67
Solar_abundances['Mg']  = -4.46
Solar_abundances['Al']  = -5.53
Solar_abundances['Si']  = -4.46
Solar_abundances['P']   = -6.50
Solar_abundances['S']   = -4.74
Solar_abundances['Cl']  = -6.72
Solar_abundances['Ar']  = -5.60
Solar_abundances['K']   = -6.88
Solar_abundances['Ca']  = -5.64
Solar_abundances['Sc']  = -8.83
Solar_abundances['Ti']  = -6.98
Solar_abundances['V']   = -8.00
Solar_abundances['Cr']  = -6.33
Solar_abundances['Mn']  = -6.54
Solar_abundances['Fe']  = -4.55
Solar_abundances['Ni']  = -5.75
Solar_abundances['Cu']  = -7.79
Solar_abundances['Zn']  = -7.40

Grid_Values['age']      = ['5.0', '5.25', '5.5', '5.75', '6.0', '6.25', '6.5', '6.75'] 
Grid_Values['zStars']   = ['-2.1'] 
Grid_Values['zGas']     = ['0.05', '0.2', '0.4'] 
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

                Model_dict['abundances_command'] = 'abundances H=0.0 He=-1.0'
                for key in Solar_abundances:
                    Model_dict['abundances_command'] = '{previous_line} {new_element}={abundance}'.format(previous_line = Model_dict['abundances_command'], new_element=key, abundance=round(Solar_abundances[key] + nplog10(float(zGas)), 3))
                
                ScriptName = Model_dict['Name'] + '.in'
     
                #Generate the script
                ct.popstar_S_Ar_Grid(ScriptFolder, Model_dict)
#                 ct.popstar_S_Ar_Grid_absAbundances(ScriptFolder, Model_dict)
                   
                #Run the cloudy script
                ct.run_script(ScriptName, ScriptFolder)
   
print 'Data treated'