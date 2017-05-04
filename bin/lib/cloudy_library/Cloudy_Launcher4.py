import pandas as pd
from numpy                              import log10 as nplog10, pi, power
from collections                        import OrderedDict
from cloudy_library.cloudy_methods      import Cloudy_Tools

def import_popstar_data():
    
    FilesFolder                 = '/home/vital/Dropbox/Astrophysics/Lore/PopStar/' 
    TableAddressn100            = 'mnras0403-2012-SD2_clusterTable2.txt'
    
    Frame_MetalsEmission        = pd.read_csv(FilesFolder + TableAddressn100, delim_whitespace = True)
    
    nH          = 10.0 #cm^-3
    c           = 29979245800.0 #cm/s
    pc_to_cm    = 3.0856776e18 #cm/pc
    
    Frame_MetalsEmission['logR_cm'] = nplog10(Frame_MetalsEmission['logR'] * pc_to_cm)
    Frame_MetalsEmission['Q']       = nplog10(power(10, Frame_MetalsEmission['logU']) * 4 * pi * c * nH * power(Frame_MetalsEmission['logR'] * pc_to_cm, 2))

    return Frame_MetalsEmission

ct              = Cloudy_Tools()

#Define script name and location
ScriptFolder = '/home/vital/Dropbox/Astrophysics/Tools/Cloudy/S_Ar_test/Total_Q_R_Grid2/'

Grid_Values                 = OrderedDict()

Grid_Values['age']          = ['5.','5.48','5.7','5.85','6.','6.1','6.18','6.24','6.3','6.35','6.4','6.44','6.48','6.51','6.54','6.57','6.6','6.63','6.65','6.68','6.7','6.72']
Grid_Values['clus_mass']    = ['12000.', '20000.', '60000.', '100000.', '200000.'] 
Grid_Values['zGas']         = ['0.0001', '0.0004', '0.004', '0.008', '0.02', '0.05'] 
Grid_Values['zStars']       = ['-2.1'] 

#Data from popstar
Frame_MetalsEmission = import_popstar_data()              
                         
#Generate the scripts with the lines we want to print the flux
ct.lines_to_extract(ScriptFolder)
 
#Fore the each grid point run a cloudy simulation
Model_dict = OrderedDict()
for age in Grid_Values['age']:
    for mass in Grid_Values['clus_mass']:
        for zGas in Grid_Values['zGas']:                                        
            for zStar in Grid_Values['zStars']:
                
                Model_dict['Name']          = 'S_Ar_Test_' + 'age'+ age + '_zStar' + zStar + '_clusMass' + mass + '_zGas' + zGas
                Model_dict['age']           = age
                Model_dict['zGas']          = zGas
                Model_dict['Metals_frac']   = str(float(zGas) / 0.02)
                Model_dict['zGas']          = zGas
                Model_dict['zStars']        = zStar
                
                index = (Frame_MetalsEmission["Z"] == float(zGas)) & (Frame_MetalsEmission["M_Msun"] == float(mass)) & (Frame_MetalsEmission["t"] == float(age))
                
                print '--Going for conditions: ', age, zGas, zStar, mass
                
                Model_dict['Q']             = Frame_MetalsEmission.loc[index, 'Q'].values[0]
                Model_dict['R']             = Frame_MetalsEmission.loc[index, 'logR_cm'].values[0]
    
                ScriptName = Model_dict['Name'] + '.in'
                print ScriptName
                 
                #Generate the script
                ct.popstar_S_Ar_QR_grid(ScriptFolder, Model_dict)
                     
                #Run the cloudy script
                ct.run_script(ScriptName, ScriptFolder)
                
print 'Data treated'