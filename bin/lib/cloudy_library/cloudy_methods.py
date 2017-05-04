from collections    import OrderedDict
from os             import environ, chdir, system
from os.path        import isfile
from subprocess     import PIPE
from subprocess     import Popen
from subprocess     import STDOUT

from numpy          import loadtxt, pi, log10 as nplog10, power
from pandas         import Series, DataFrame, read_csv


class cloudy_scripts():
    
    def __init__(self):
        
        self.default_script_folder = None

    def select_script(self, simulation_type, ScriptName, ScriptFolder,key = None, data_dict = None):
        
        #Quick start basic example
        if simulation_type == 'pn':
            self.pn_example(ScriptFolder, ScriptName)
        
        #Popstar Sulfur argon model
        if simulation_type == 'Argon_Sulfur_LinesModel':
            self.popstar_S_Ar_relation(ScriptFolder, ScriptName)
        
        #Popstar Sulfur argon model
        if simulation_type == 'Argon_Sulfur_LinesModel_Ages':
            self.popstar_S_Ar_ages_relation(ScriptFolder, ScriptName, key, data_dict)        

        #Popstar Sulfur argon model
        if simulation_type == 'Argon_Sulfur_LinesModel_Ionization':
            self.popstar_S_Ar_ages_IonizationFactor(ScriptFolder, ScriptName, key, data_dict)        
 
        #Popstar Sulfur argon model
        if simulation_type == 'Argon_Sulfur_LinesModel_Hbeta':
            self.popstar_S_Ar_ages_Hbeta(ScriptFolder, ScriptName, key, data_dict)   
              
        return

    def pn_example(self, ScriptFolder, ScriptName):
        
        script_Lines = [
                        'blackbody, T=1e5 K',
                        'luminosity total 38',
                        'radius 18',
                        'hden 5',
                        'sphere',
                        'abundances planetary nebula',
                        'save overview ".ovr"',
                        'save continuum ".con" units microns'
                        ]
        self.save_script(ScriptFolder + ScriptName, script_Lines)

        return

    def absolute_abundances(self, chemical_model = None):
        
        abundances_dict         = OrderedDict()

        if chemical_model == None:
            
            abundances_dict['Li']   = -8.69
            abundances_dict['Be']   = -10.58
            abundances_dict['B']    = -9.21
            abundances_dict['C']    = -3.61
            abundances_dict['N']    = -4.07
            abundances_dict['O']    = -3.31
            abundances_dict['Ne']   = -4.00
            abundances_dict['F']    = -7.52
            abundances_dict['Ne']   = -4.0
            abundances_dict['Na']   = -5.67
            abundances_dict['Mg']   = -4.46
            abundances_dict['Al']   = -5.53
            abundances_dict['Si']   = -4.46
            abundances_dict['P']    = -6.50
            abundances_dict['S']    = -4.74
            abundances_dict['Cl']   = -6.72
            abundances_dict['Ar']   = -5.60
            abundances_dict['K']    = -6.88
            abundances_dict['Ca']   = -5.64
            abundances_dict['Sc']   = -8.83
            abundances_dict['Ti']   = -6.98
            abundances_dict['V']    = -8.00
            abundances_dict['Cr']   = -6.33
            abundances_dict['Mn']   = -6.54
            abundances_dict['Fe']   = -4.55
            abundances_dict['Ni']   = -5.75
            abundances_dict['Cu']   = -7.79
            abundances_dict['Zn']   = -7.40
            
        return
    
    def S_Ar_star(self, ScriptFolder, ScriptName, data_dict):
    
        script_Lines = [
                                
                        'hden = 2',
                        'constant density',                        
                        'blackbody, T={temperature} K'.format(temperature=data_dict['temp']),
                        'ionization parameter ' + data_dict['u'],
                        'radius 18',
                        'sphere',
                        'CMB',
                        'abundances OLD SOLAR 84',
                        'metals ' + data_dict['zGas'],
                        'save overview ".ovr"',
                        'save continuum ".con"' + ' units angstroms',
                        'save line list ".linesPrediction"' + ' "S_Ar_test.lineslist" last no hash',
                        ]
        

        self.save_script(ScriptFolder + ScriptName, script_Lines) 
        
        return   
               
    def lines_to_extract(self, ScriptFolder):

        #Argon and Sulfur lines to store
        emission_lines_list = [
                        'H  1 4861.36A',
                        'H  1 6563A',
                        'TOTL 3727A',                     
                        'O  3  5007A',
                        'O  3  4959A',
                        'TOTL  4363A',
                        'S  3 9068.62A',
                        'H  1 4.051m',
                        'S  3 18.67m',
                        'S  3 9532A',
                        'S  3 6312A',
                        'S II 6731A',
                        'S II 6716A',                        
                        'S  4 10.51m',
                        'AR 3 7135A',
                        'AR 3 7751.11A',
                        'AR 4 4711.26A',
                        'AR 4 4740.12A',                        
                             ]         
         
        self.save_script(ScriptFolder + "S_Ar_test.lineslist", emission_lines_list)

        return

    def popstar_S_Ar_Grid(self, ScriptFolder, data_dict):
                
        #Cloudy script
        script_Lines = [
                        'hden = 2',
                        'constant density',  
                        'ionization parameter ' + data_dict['u'],
                        'radius = 19.5',
                        'table starburst log age = {age} log z = {zStar}, file="sp-sal2.mod"'.format(age = data_dict['age'], zStar = data_dict['zStars']),
                        'sphere static',
                        'diffuse field OTS',
#                         'set save resolving power 10000',
                        'CMB',
                        'abundances OLD SOLAR 84',
                        'metals ' + data_dict['zGas'],
#                         'save transmitted continuum file = "{File_name}" last units angstroms'.format(File_name = '.transContinuum'),
#                         'save incident continuum file = "{File_name}" last units angstroms'.format(File_name = '.inciContinuum'),                        
                        'save overview ".ovr"',
                        'save continuum ".con"' + ' units angstroms',
                        'save line list ".linesPrediction"' + ' "S_Ar_test.lineslist" last no hash',
                        ]
        
        self.save_script(ScriptFolder + data_dict['Name'] + '.in', script_Lines)
           
        return  

    def popstar_S_Ar_QR_grid(self, ScriptFolder, data_dict):
                
        #Cloudy script
        script_Lines = [
                        'hden = 2',
                        'constant density',
                        'Q(H)= {Q_H}'.format(Q_H = round(data_dict['Q'], 4)),
                        'radius= {radious_pc}'.format(radious_pc = round(data_dict['R'], 4)),                                                                 
                        'table starburst log age = {age} log z = {zStar}, file="sp-sal2.mod"'.format(age = data_dict['age'], zStar = data_dict['zStars']),
                        'sphere static',
                        'CMB',
                        'abundances OLD SOLAR 84',
                        'iterations=2',
                        'metals ' + data_dict['Metals_frac'],             
                        'save overview ".ovr"',
                        'save continuum ".con"' + ' units angstroms',
                        'save line list ".linesPrediction"' + ' "S_Ar_test.lineslist" last no hash',
                        ]
        
        self.save_script(ScriptFolder + data_dict['Name'] + '.in', script_Lines)
           
        return

    def popstar_S_Ar_Grid_absAbundances(self, ScriptFolder, data_dict):
                
        #Cloudy script
        script_Lines = [
                        'hden = 2',
                        'constant density',  
                        'ionization parameter ' + data_dict['u'],
                        'radius = 19.5',
                        'table starburst log age = {age} log z = {zStar}, file="sp-sal2.mod"'.format(age = data_dict['age'], zStar = data_dict['zStars']),
                        'sphere static',
                        'diffuse field OTS',
                        'stop temp = 4000',
                        'CMB',
                        data_dict['abundances_command'],
                        'save transmitted continuum file = "{File_name}" last units angstroms'.format(File_name = '.transContinuum'),
                        'save incident continuum file = "{File_name}" last units angstroms'.format(File_name = '.inciContinuum'),                        
                        'save overview ".ovr"',
                        'save continuum ".con"' + ' units angstroms',
                        'save line list ".linesPrediction"' + ' "S_Ar_test.lineslist" last no hash',
                        'save grid ".grd"'                          
                             ]
        
        self.save_script(ScriptFolder + data_dict['Name'] + '.in', script_Lines)
           
        return

    def popstar_S_Ar_Grid_Hbeta(self, ScriptFolder, data_dict):
                
        #Cloudy script
        script_Lines = [
                        'hden = 2',
                        'ionization parameter ' + data_dict['u'],
                        'radius = 18',
                        'table starburst log age = {age} log z = {zStar}, file="sp-kro_z0001-z05_stellar.mod"'.format(age = data_dict['age'], zStar = data_dict['zStars']),
                        'sphere static',
                        'CMB',
                        'abundances hii region',
                        'metals ' + data_dict['zGas'],
                        'normalize to "Inci" 4860A scale factor = 4861',
                        'save overview ".ovr"',
                        'save continuum ".con"' + ' units microns',
                        'save line list ".linesPrediction"' + ' "S_Ar_test.lineslist" last no hash',
                        'save grid ".grd"'                          
                             ]
        
        self.save_script(ScriptFolder + data_dict['Name'] + '.in', script_Lines)
        
        return     

    def save_script(self, scriptAddress, lines_list):
        with open(scriptAddress, 'w') as f:
            for line in lines_list:
                f.write(line + '\n')
                
        return

class cloudy_plots():
    
    def __init__(self):
        
        self.plotting_vector = None
        
    def Incident_Transmitted_Radiation(self, OuputFile, OutputFolder, plotvect):
        
        nu, incident, trans, total = loadtxt(OutputFolder + OuputFile, skiprows = 0, usecols = (0, 1, 2, 4), unpack = True)

        return

class Cloudy_Tools(cloudy_scripts):
    
    def __init__(self):
        
        #Load scripts class
        cloudy_scripts.__init__(self)
        
        #Adding the cloudy path to the environment
        self.my_env = environ
        self.my_env["PATH"] = "/usr/sbin:/sbin:/home/vital/.my_bin:" + self.my_env["PATH"] #This variable should be adapted to the computer
        
    def run_script(self, ScriptName, ScriptFolder):
        
        #Move to the script folder
        chdir(ScriptFolder)
         
        #Preparing the command
        Command = 'cloudy ' + ScriptName

        print "\n--Launching command:"
        print "  ", Command, '@', ScriptFolder, '\n'
        
        #Run the command
        p = Popen(Command, shell=True, stdout=PIPE, stderr=STDOUT, env=self.my_env)

        #Terminal output in terminal
        if len(p.stdout.readlines()) > 0:
            print '-- Code output wording\n'
            for line in p.stdout.readlines():
                print line
                
        return
        
    def open_output(self, ScriptName, ScriptFolder):
        
        OutputName = ScriptName.replace('.in', '.out')
        
        system('gedit '+ ScriptFolder + OutputName)
        
        proc = Popen(['gedit', ScriptFolder + OutputName])

        proc.wait()
                
        return

    def import_popstar_data(self):
         
        FilesFolder                 = '/home/vital/Dropbox/Astrophysics/Lore/PopStar/' 
        TableAddressn100            = 'mnras0403-2012-SD2_clusterTable2.txt'
         
        Frame_MetalsEmission        = read_csv(FilesFolder + TableAddressn100, delim_whitespace = True)
         
        nH          = 10.0 #cm^-3
        c           = 29979245800.0 #cm/s
        pc_to_cm    = 3.0856776e18 #cm/pc
         
        Frame_MetalsEmission['logR_cm'] = nplog10(Frame_MetalsEmission['logR'] * pc_to_cm)
        Frame_MetalsEmission['Q']       = nplog10(power(10, Frame_MetalsEmission['logU']) * 4 * pi * c * nH * power(Frame_MetalsEmission['logR'] * pc_to_cm, 2))
     
        return Frame_MetalsEmission

    def load_predicted_lines(self, ScriptName, ScriptFolder):
        
#         Lines_file  = ScriptName.replace('.in', '.lineslist')
        Lines_file  = 'S_Ar_test.lineslist'
        Lines_list  = loadtxt(ScriptFolder + Lines_file, dtype=str, usecols=[2])
        
        
        Fluxes_file = ScriptName.replace('.in', '.linesPrediction')

        Fluxes_list = loadtxt(ScriptFolder + Fluxes_file, usecols=range(2, len(Lines_list)+2), unpack=False, comments='#')
        
        Line_dict   = OrderedDict()
        
        for i in range(len(Lines_list)):
            Line            = Lines_list[i]
            Line_dict[Line] = Fluxes_list[:,i]
        
        return Line_dict

    def load_predicted_lines_individual(self, ScriptName, ScriptFolder):
        
#         Lines_file  = ScriptName.replace('.in', '.lineslist')
        Lines_file  = 'S_Ar_test.lineslist'
        Lines_list  = loadtxt(ScriptFolder + Lines_file, dtype=str, usecols=[2])
        
        Fluxes_file = ScriptName.replace('.in', '.linesPrediction')

        Fluxes_list = loadtxt(ScriptFolder + Fluxes_file, usecols=range(2, len(Lines_list)+2), unpack=False, comments='#')
        
        Line_dict   = OrderedDict()
        
        for i in range(len(Lines_list)):
            Line            = Lines_list[i]
            Line_dict[Line] = Fluxes_list[i]
        
        return Line_dict

    
# def Launch_Command(ScriptName, ScriptFolder, Variable_Saving = False):
#     
#     #Load the path into the enviroment
#     my_env = environ
#     my_env["PATH"] = "/usr/sbin:/sbin:/home/vital/.my_bin:" + my_env["PATH"]
# 
#     Command = 'cloudy ' + ScriptFolder + ScriptName
#     
#     if Variable_Saving == True:
#         OutPutName = CloudySaveManager(Seed = ScriptName)
#         Command = 'cloudy ' + ScriptName + ' ' + OutPutName
# 
#     print "-- Launching command:"
#     print "  ", Command, '\n'
#     
#     p = Popen(Command, shell=True, stdout=PIPE, stderr=STDOUT, env=self.my_env)
#          
#     for line in p.stdout.readlines():
#         print line
#         
#     retval = p.wait()
#     
#     return
# 
# def CloudySaveManager(Seed, OuputFolder = '/home/vital/Workspace/Astro/Cloudy/Output_Files/'):
#     
#     i = 0
#     NewName = None
#     FileExists = True
#     
#     while FileExists:
#         NewName = Seed + '_' + str(i)
#         print 'Proposed Name:', OuputFolder + NewName
#         if isfile(OuputFolder + NewName):
#             i = i + 1
#         else:
#             FileExists = False
#     
#     return NewName