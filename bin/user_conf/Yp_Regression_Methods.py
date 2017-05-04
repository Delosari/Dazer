from collections                        import OrderedDict 

from numpy                              import array, intersect1d, random, zeros
from pymc                               import database
from uncertainties                      import ufloat
from uncertainties.unumpy               import uarray, nominal_values, std_devs

from Astro_Libraries.Abundances_Class   import Chemical_Analysis


def import_plots_wording(pv):

    Titles_dict             = OrderedDict()
    Colors_dict             = OrderedDict()

    Titles_dict['O_Regression']                     = ('Oxygen abundance',                                  'Primordial Helium regression: Oxygen metallicity tracer',    r'$\frac{O}{H}$', r'$Y_{\frac{O}{H}}$', "Oxigen",   1e-6)
    Titles_dict['N_Regression']                     = ('Nitrogen abundance',                                'Primordial Helium regression: Nitrogen metallicity tracer',  r'$\frac{N}{H}$', r'$Y_{\frac{N}{H}}$', "Nitrogen", 1e-5)
    Titles_dict['S_Regression']                     = ('Sulfur abundance',                                  'Primordial Helium regression: Sulfur metallicity tracer',    r'$\frac{S}{H}$', r'$Y_{\frac{S}{H}}$', "Sulfur",   1e-5)
    Titles_dict['S_ArCorr_Regression']              = ('Sulfur abundance with argon correction',            'Primordial Helium regression: Sulfur metallicity tracer',    r'$\frac{S}{H}$', r'$Y_{\frac{S}{H}}$', "Sulfur",   1e-5)
    Titles_dict['O_Regression_Inference']           = ('Oxygen abundance inference',                        'Primordial Helium regression: Oxygen metallicity tracer',    r'$\frac{O}{H}$', r'$Y_{\frac{O}{H}}$', "Oxigen",   1e-6)
    Titles_dict['N_Regression_Inference']           = ('Nitrogen abundance inference',                      'Primordial Helium regression: Nitrogen metallicity tracer',  r'$\frac{N}{H}$', r'$Y_{\frac{N}{H}}$', "Nitrogen", 1e-5)
    Titles_dict['S_Regression_Inference']           = ('Sulfur abundance inference',                        'Primordial Helium regression: Sulfur metallicity tracer',    r'$\frac{S}{H}$', r'$Y_{\frac{S}{H}}$', "Sulfur",   1e-5)
    Titles_dict['S_ArCorr_Regression_Inference']    = ('Sulfur abundance with argon correction inference',  'Primordial Helium regression: Sulfur metallicity tracer',    r'$\frac{S}{H}$', r'$Y_{\frac{S}{H}}$', "Sulfur",   1e-5)

    Colors_dict['O_Regression']                     = 1
    Colors_dict['N_Regression']                     = 2
    Colors_dict['S_Regression']                     = 3
    Colors_dict['S_ArCorr_Regression']              = 3
    Colors_dict['O_Regression_Inference']           = 1
    Colors_dict['N_Regression_Inference']           = 2
    Colors_dict['S_Regression_Inference']           = 3
    Colors_dict['S_ArCorr_Regression_Inference']    = 3
             
    return Titles_dict, Colors_dict
    
def import_data_from_objLog(FilesList, Objects_Include, pv):
    
    
    List_Abundances     = ['OI_HI', 'NI_HI', 'SI_HI', 'SI_HI_ArCorr', 'Y_Mass_O', 'Y_Mass_S', 'Y_Inference_O', 'Y_Inference_S']
    #List_Abundances    = ['OI_HI', 'NI_HI', 'SI_HI', 'SI_HI_ArCorr', 'Y_Mass_O', 'Y_Mass_S', 'Y_inf_O', 'Y_inf_S']

    #Dictionary of dictionaries to store object abundances
    Abund_dict = OrderedDict()
    for abund in List_Abundances:
        Abund_dict[abund] = OrderedDict()

    #Loop through files
    for i in range(len(FilesList)):
        #Analyze file address
        CodeName, FileName, FileFolder  = pv.Analyze_Address(FilesList[i])  
      
        if CodeName in Objects_Include:
            #Loop through abundances in the log
            for abund in List_Abundances:
                Abund_Mag = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter = abund, Assumption = 'float')
                #If the abundance was measure store it 
                if Abund_Mag != None:
                    Abund_dict[abund][CodeName] = Abund_Mag
        
    #Dictionary to store objects with abundances pairs for regressions. 
    #As an initial value for the keys we define the abundances we want to use for the regression
    Abundances_Pairs_dict = OrderedDict()
    Abundances_Pairs_dict['O_Regression']                   = ('OI_HI','Y_Mass_O')      
    Abundances_Pairs_dict['N_Regression']                   = ('NI_HI','Y_Mass_O')      
    Abundances_Pairs_dict['S_Regression']                   = ('SI_HI','Y_Mass_S')      
    Abundances_Pairs_dict['S_ArCorr_Regression']            = ('SI_HI_ArCorr','Y_Mass_S')
    Abundances_Pairs_dict['O_Regression_Inference']         = ('OI_HI','Y_Inference_O')      
    Abundances_Pairs_dict['N_Regression_Inference']         = ('NI_HI','Y_Inference_O')      
    Abundances_Pairs_dict['S_Regression_Inference']         = ('SI_HI','Y_Inference_S')      
    Abundances_Pairs_dict['S_ArCorr_Regression_Inference']  = ('SI_HI_ArCorr','Y_Inference_S') 
        
    #Loop through the regression lists and get objects with both abundances observed
    for j in range(len(Abundances_Pairs_dict)):
        
        #Get the elements keys for the regression
        Vector, Elem_X, Elem_Y = Abundances_Pairs_dict.keys()[j], Abundances_Pairs_dict.values()[j][0], Abundances_Pairs_dict.values()[j][1]
        
        #Determine objects with both abundances observed
        Obj_vector  = intersect1d(Abund_dict[Elem_X].keys(), Abund_dict[Elem_Y].keys(), assume_unique = True)
        X_vector    = zeros(len(Obj_vector))
        Y_vector    = zeros(len(Obj_vector))
        X_vector_E  = zeros(len(Obj_vector))
        Y_vector_E  = zeros(len(Obj_vector))
                        
        #Generate abundances vectors
        for z in range(len(Obj_vector)):  
            X_vector[z] = nominal_values(Abund_dict[Elem_X][Obj_vector[z]])
            X_vector_E[z] = std_devs(Abund_dict[Elem_X][Obj_vector[z]])            
            Y_vector[z] = nominal_values(Abund_dict[Elem_Y][Obj_vector[z]])
            Y_vector_E[z] = std_devs(Abund_dict[Elem_Y][Obj_vector[z]])
    
        Abundances_Pairs_dict[Vector] = (list(Obj_vector), uarray(X_vector, X_vector_E), uarray(Y_vector, Y_vector_E))
        
    return Abundances_Pairs_dict

def import_data_from_objLog_pyneb(FilesList, Objects_Include, pv):
    
    List_Abundances     = ['OI_HI_pn', 'NI_HI_pn', 'SI_HI_pn', 'SI_HI_ArCorr_pn', 'Y_Mass_O_pn', 'Y_Mass_S_pn', 'Y_Inference_O_pn', 'Y_Inference_S_pn']
    #List_Abundances    = ['OI_HI', 'NI_HI', 'SI_HI', 'SI_HI_ArCorr', 'Y_Mass_O', 'Y_Mass_S', 'Y_inf_O', 'Y_inf_S']

    #Dictionary of dictionaries to store object abundances
    Abund_dict = OrderedDict()
    for abund in List_Abundances:
        Abund_dict[abund] = OrderedDict()

    #Loop through files
    for i in range(len(FilesList)):
        #Analyze file address
        CodeName, FileName, FileFolder  = pv.Analyze_Address(FilesList[i])  
      
        if CodeName in Objects_Include:
            #Loop through abundances in the log
            for abund in List_Abundances:
                Abund_Mag = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter = abund, Assumption = 'float')
                #If the abundance was measure store it 
                if Abund_Mag != None:
                    Abund_dict[abund][CodeName] = Abund_Mag
        
    #Dictionary to store objects with abundances pairs for regressions. 
    #As an initial value for the keys we define the abundances we want to use for the regression
    Abundances_Pairs_dict = OrderedDict()
    Abundances_Pairs_dict['O_Regression']                   = ('OI_HI_pn','Y_Mass_O_pn')      
    Abundances_Pairs_dict['N_Regression']                   = ('NI_HI_pn','Y_Mass_O_pn')      
    Abundances_Pairs_dict['S_Regression']                   = ('SI_HI_pn','Y_Mass_S_pn')      
    Abundances_Pairs_dict['S_ArCorr_Regression']            = ('SI_HI_ArCorr_pn','Y_Mass_S_pn')
    Abundances_Pairs_dict['O_Regression_Inference']         = ('OI_HI_pn','Y_Inference_O_pn')      
    Abundances_Pairs_dict['N_Regression_Inference']         = ('NI_HI_pn','Y_Inference_O_pn')      
    Abundances_Pairs_dict['S_Regression_Inference']         = ('SI_HI_pn','Y_Inference_S_pn')      
    Abundances_Pairs_dict['S_ArCorr_Regression_Inference']  = ('SI_HI_ArCorr_pn','Y_Inference_S_pn') 
        
    #Loop through the regression lists and get objects with both abundances observed
    for j in range(len(Abundances_Pairs_dict)):
        
        #Get the elements keys for the regression
        Vector, Elem_X, Elem_Y = Abundances_Pairs_dict.keys()[j], Abundances_Pairs_dict.values()[j][0], Abundances_Pairs_dict.values()[j][1]
        
        #Determine objects with both abundances observed
        Obj_vector  = intersect1d(Abund_dict[Elem_X].keys(), Abund_dict[Elem_Y].keys(), assume_unique = True)
        X_vector    = zeros(len(Obj_vector))
        Y_vector    = zeros(len(Obj_vector))
        X_vector_E  = zeros(len(Obj_vector))
        Y_vector_E  = zeros(len(Obj_vector))
                        
        #Generate abundances vectors
        for z in range(len(Obj_vector)):  
            X_vector[z] = nominal_values(Abund_dict[Elem_X][Obj_vector[z]])
            X_vector_E[z] = std_devs(Abund_dict[Elem_X][Obj_vector[z]])            
            Y_vector[z] = nominal_values(Abund_dict[Elem_Y][Obj_vector[z]])
            Y_vector_E[z] = std_devs(Abund_dict[Elem_Y][Obj_vector[z]])
    
        Abundances_Pairs_dict[Vector] = (list(Obj_vector), uarray(X_vector, X_vector_E), uarray(Y_vector, Y_vector_E))
        
    return Abundances_Pairs_dict


def import_data_from_objLog_Triple(FilesList, Objects_Include, pv):
    
    List_Abundances     = ['OI_HI', 'NI_HI', 'SI_HI', 'SI_HI_ArCorr', 'Y_Mass_O', 'Y_Mass_S', 'Y_Inference_O', 'Y_Inference_S']
    #List_Abundances    = ['OI_HI', 'NI_HI', 'SI_HI', 'SI_HI_ArCorr', 'Y_Mass_O', 'Y_Mass_S', 'Y_inf_O', 'Y_inf_S']

    #Dictionary of dictionaries to store object abundances
    Abund_dict = OrderedDict()
    for abund in List_Abundances:
        Abund_dict[abund] = OrderedDict()

    #Loop through files
    for i in range(len(FilesList)):
        #Analyze file address
        CodeName, FileName, FileFolder  = pv.Analyze_Address(FilesList[i])  
      
        if CodeName in Objects_Include:
            #Loop through abundances in the log
            for abund in List_Abundances:
                Abund_Mag = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter = abund, Assumption = 'float')
                #If the abundance was measure store it 
                if Abund_Mag != None:
                    Abund_dict[abund][CodeName] = Abund_Mag
        
    #Dictionary to store objects with abundances pairs for regressions. 
    #As an initial value for the keys we define the abundances we want to use for the regression
    Abundances_Pairs_dict = OrderedDict()
    Abundances_Pairs_dict['O_Regression_Inference']         = ('OI_HI','Y_Inference_O')      
    Abundances_Pairs_dict['N_Regression_Inference']         = ('NI_HI','Y_Inference_O')      
    Abundances_Pairs_dict['S_ArCorr_Regression_Inference']  = ('SI_HI_ArCorr','Y_Inference_S') 
        
    #Loop through the regression lists and get objects with both abundances observed
    for j in range(len(Abundances_Pairs_dict)):
        
        #Get the elements keys for the regression
        Vector, Elem_X, Elem_Y = Abundances_Pairs_dict.keys()[j], Abundances_Pairs_dict.values()[j][0], Abundances_Pairs_dict.values()[j][1]
        
        #Determine objects with both abundances observed
        Obj_vector  = intersect1d(Abund_dict[Elem_X].keys(), Abund_dict[Elem_Y].keys(), assume_unique = True)
        X_vector    = zeros(len(Obj_vector))
        Y_vector    = zeros(len(Obj_vector))
        X_vector_E  = zeros(len(Obj_vector))
        Y_vector_E  = zeros(len(Obj_vector))
                        
        #Generate abundances vectors
        for z in range(len(Obj_vector)):  
            X_vector[z] = nominal_values(Abund_dict[Elem_X][Obj_vector[z]])
            X_vector_E[z] = std_devs(Abund_dict[Elem_X][Obj_vector[z]])            
            Y_vector[z] = nominal_values(Abund_dict[Elem_Y][Obj_vector[z]])
            Y_vector_E[z] = std_devs(Abund_dict[Elem_Y][Obj_vector[z]])
    
        Abundances_Pairs_dict[Vector] = (list(Obj_vector), uarray(X_vector, X_vector_E), uarray(Y_vector, Y_vector_E))
        
    return Abundances_Pairs_dict

def Galaxy_sample():
    
#     Objects_Include = ['70', 'SDSS2', '06', '08']
#     Objects_Include  = ['03', '04v2', '06', '08', '09', '10', '14', '24', '27', '71', 'SDSS1', '70', 'SDSS3', '51959-092', '51991-224', '52235-602', '52319-521', 'J2225', 'SHOC575v1', 'SHOC579'] 
#     Objects_Include  = ['03', '04v2', '06', '08', '09', '10', '14', '24', '27', '71', 'SDSS1', '70', 'SDSS3', '51959-092', '51991-224', '52235-602', '52319-521', 'J2225', 'SHOC575v2', 'SHOC579'] 
    Objects_Include  = ['03', 
#                         '04v1', #Repeated
                        '04v2', #Repeated
                        '06',
                        '08',
                        '09',
#                         '10',
#                         '11',
                        '14',
                        '24',
                        '27',
#                         '70', #RepeatedSDSS2 #It is poorly loaded but why it is not loaded
                        '71',
                        'SDSS1', 
                        'SDSS2', 
                        'SDSS3', 
#                         '51959-092',#Repeated11
                        '51991-224',
                        '52235-602',
                        '52319-521',
                        '52703-612v2', #Repeated24 
                        'J2225',
#                         'SHOC036',#Repeated06 No TSIII measurement
#                         'SHOC575v1', #Repeated Poorly starlight fitting
#                         'SHOC575v2',
#                         'SHOC579',
                        'SHOC588', #RepeatedSDSS3 Helium abundance well fitted but sulfur worse SDSS3
#                         'SHOC593', R_S3 = 6
                        ]
        
        
    return Objects_Include

def WMAP_Coordinates():
        
    WMAP_coordinates = array([ufloat(0.0, 0.0), ufloat(0.24709, 0.00025)])
    
    return WMAP_coordinates

def Get_Traces(pv, Regressiontype, Abundances_dict, catalogue_folder, database_extension):
    
    Y_dist_dict         = OrderedDict()
    Metal_Dist_dict     = OrderedDict()
    
    Candidate_Objects   = Abundances_dict[Regressiontype][0]
    Element_abundance   = Abundances_dict[Regressiontype][1]
    Y_abundance         = Abundances_dict[Regressiontype][2]
    
    ch_an                   = Chemical_Analysis()
    
    for i in range(len(Candidate_Objects)):
        
        CodeName        = Candidate_Objects[i]
        Metallic_elem   = Element_abundance[i]
        Element         = Regressiontype[0:Regressiontype.find('_')]

        FileFolder      = catalogue_folder + CodeName + '/'
        db_Address      = FileFolder + 'J0_' +  CodeName + database_extension

        #SI_HI           = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter='SI_HI_ArCorr', Assumption='float')
        HeIII_HII       = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter='HeIII_HII', Assumption='float')

        #Get y_plus trace
        pymc_database   = database.pickle.load(db_Address)
        HeII_HII_trace  = pymc_database.trace('He_abud')[:]
        pymc_database.close()
        
        #Generate the 5000 array from the distribution
        HeII_HII_dist   = random.choice(HeII_HII_trace, size=len(HeII_HII_trace)*2) 
        Elemental_dist  = random.normal(Metallic_elem.nominal_value, Metallic_elem.std_dev, size = len(HeII_HII_trace)*2)
                
        #Set to zero the HeIII ionic abundance if not observed
        if HeIII_HII != None:
            HeIII_HII_dist = random.normal(HeIII_HII.nominal_value, HeIII_HII.std_dev, size = len(HeII_HII_trace)*2)
        else:
            HeIII_HII_dist = zeros(len(HeII_HII_trace)*2)
            
        #Calculate the HeI/HI distribution
        HeI_dist = HeII_HII_dist + HeIII_HII_dist
        
        #Calculate the Y distribution for the Hydrogen and helium abundance
        if Regressiontype in ['S_Regression', 'S_ArCorr_Regression', 'S_Regression_Inference', 'S_ArCorr_Regression_Inference']:
            OI_SI_dist                  = random.normal(ch_an.OI_SI.nominal_value, ch_an.OI_SI.std_dev, size = len(HeII_HII_trace)*2)
            Y_dist_dict[CodeName]       = (4 * HeI_dist * (1 - 20 * OI_SI_dist * Elemental_dist)) / (1 + 4 * HeI_dist) #WARNING WE ARE NOT TAKING INTO CONSIDERATION THE INCREASE IN ERROR DUE TO THE UNCERTAINTY IN OI_SI
            Metal_Dist_dict[CodeName]   = Elemental_dist
            
        elif Regressiontype in ['O_Regression', 'O_Regression_Inference']:
            Y_dist_dict[CodeName]       = (4 * HeI_dist * (1 - 20 * Elemental_dist)) / (1 + 4 * HeI_dist) 
            Metal_Dist_dict[CodeName]   = Elemental_dist
            
        elif Regressiontype in ['N_Regression', 'N_Regression_Inference']:
            OI_HI                       = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter='OI_HI', Assumption='float')
            OI_HI_dist                  = random.normal(OI_HI.nominal_value, OI_HI.std_dev, size = len(HeII_HII_trace)*2)
            Y_dist_dict[CodeName]       = (4 * HeI_dist * (1 - 20 * OI_HI_dist)) / (1 + 4 * HeI_dist) #WARNING WE ARE NOT TAKING INTO CONSIDERATION THE INCREASE IN ERROR DUE TO THE UNCERTAINTY IN OI_SI
            Metal_Dist_dict[CodeName]   = Elemental_dist            
            
    return Metal_Dist_dict, Y_dist_dict

def Get_Traces_pn(pv, Regressiontype, Abundances_dict, catalogue_folder, database_extension):
    
    Y_dist_dict         = OrderedDict()
    Metal_Dist_dict     = OrderedDict()
    
    Candidate_Objects   = Abundances_dict[Regressiontype][0]
    Element_abundance   = Abundances_dict[Regressiontype][1]
    Y_abundance         = Abundances_dict[Regressiontype][2]
    
    ch_an                   = Chemical_Analysis()
    
    for i in range(len(Candidate_Objects)):
        
        CodeName        = Candidate_Objects[i]
        Metallic_elem   = Element_abundance[i]
        Element         = Regressiontype[0:Regressiontype.find('_')]

        FileFolder      = catalogue_folder + CodeName + '/'
        db_Address      = FileFolder + 'J0_' +  CodeName + database_extension

        #SI_HI           = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter='SI_HI_ArCorr', Assumption='float')
        HeIII_HII       = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter='HeIII_HII', Assumption='float')

        #Get y_plus trace
        pymc_database   = database.pickle.load(db_Address)
        HeII_HII_trace  = pymc_database.trace('He_abud')[:]
        pymc_database.close()
        
        #Generate the 5000 array from the distribution
        HeII_HII_dist   = random.choice(HeII_HII_trace, size=len(HeII_HII_trace)*2) 
        Elemental_dist  = random.normal(Metallic_elem.nominal_value, Metallic_elem.std_dev, size = len(HeII_HII_trace)*2)
                
        #Set to zero the HeIII ionic abundance if not observed
        if HeIII_HII != None:
            HeIII_HII_dist = random.normal(HeIII_HII.nominal_value, HeIII_HII.std_dev, size = len(HeII_HII_trace)*2)
        else:
            HeIII_HII_dist = zeros(len(HeII_HII_trace)*2)
            
        #Calculate the HeI/HI distribution
        HeI_dist = HeII_HII_dist + HeIII_HII_dist
        
        #Calculate the Y distribution for the Hydrogen and helium abundance
        if Regressiontype in ['S_Regression', 'S_ArCorr_Regression', 'S_Regression_Inference', 'S_ArCorr_Regression_Inference']:
            OI_SI_dist                  = random.normal(ch_an.OI_SI.nominal_value, ch_an.OI_SI.std_dev, size = len(HeII_HII_trace)*2)
            Y_dist_dict[CodeName]       = (4 * HeI_dist * (1 - 20 * OI_SI_dist * Elemental_dist)) / (1 + 4 * HeI_dist) #WARNING WE ARE NOT TAKING INTO CONSIDERATION THE INCREASE IN ERROR DUE TO THE UNCERTAINTY IN OI_SI
            Metal_Dist_dict[CodeName]   = Elemental_dist
            
        elif Regressiontype in ['O_Regression', 'O_Regression_Inference']:
            Y_dist_dict[CodeName]       = (4 * HeI_dist * (1 - 20 * Elemental_dist)) / (1 + 4 * HeI_dist) 
            Metal_Dist_dict[CodeName]   = Elemental_dist
            
        elif Regressiontype in ['N_Regression', 'N_Regression_Inference']:
            OI_HI                       = pv.GetParameter_ObjLog(CodeName, FileFolder, Parameter='OI_HI_pn', Assumption='float')
            OI_HI_dist                  = random.normal(OI_HI.nominal_value, OI_HI.std_dev, size = len(HeII_HII_trace)*2)
            Y_dist_dict[CodeName]       = (4 * HeI_dist * (1 - 20 * OI_HI_dist)) / (1 + 4 * HeI_dist) #WARNING WE ARE NOT TAKING INTO CONSIDERATION THE INCREASE IN ERROR DUE TO THE UNCERTAINTY IN OI_SI
            Metal_Dist_dict[CodeName]   = Elemental_dist            
            
    return Metal_Dist_dict, Y_dist_dict
