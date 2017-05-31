#!/usr/bin/python

import  pymc
from    CodeTools.PlottingManager                               import myPickle
from    ManageFlow                                              import DataToTreat
from    Astro_Libraries.Abundances_InferenceModel_Helium_v12    import HeAbundance_InferenceStructure

#Import classes
pv = myPickle()
bm = HeAbundance_InferenceStructure()

#Define data type and location
Catalogue_Dic           = DataToTreat()
Pattern                 = Catalogue_Dic['Datatype'] + '_EmissionRedden2nd.fits'                            #First batch process for untreated spectra
AbundancesFileExtension = '_'+ Catalogue_Dic['Datatype'] + '_EmissionFlux_LinesLog_v3.txt'  
database_extension      = '_extandar_30000_5000_10_Revision3'
globalfile_extension    = '_global_30000_5000_10_Revision3'

#Locate files on hard drive
FilesList               = pv.Folder_Explorer(Pattern, Catalogue_Dic['Obj_Folder'], CheckComputer=False)

#Define plot frame and colors
pv.FigFormat_One(ColorConf = 'Night1')

objects_que_petaron = ['09','10'] 
# objects_que_petaron = ['SHOC579']

for i in range(len(FilesList)):
    
    try:
                
        #Analyze file address
        CodeName, FileName, FileFolder  = pv.Analyze_Address(FilesList[i])
        
        if CodeName in objects_que_petaron:
        
            print 'Going for', CodeName, 'number', i
            #Make sure the object logs ready
            pv.SetLogFile(CodeName + pv.ObjectLog_extension, FileFolder)
            
            #Import the lines log data
            bm.Import_Object_Data(FileFolder, CodeName, AbundancesFileExtension)
        
            #Check the observational data for Hydrogen and Helium lines which can be used for the analysis:
            bm.Check_EmissionLinesObserved(bm.ObjectData_dict)
            
            #Load the data into vectors for the MCMC
            bm.Preload_Data(Hbeta_Normalized = False)
            
            #Declare the MCMC dictionary
            MCMC_dict = bm.Bayesian_HeliumAbundance_Analysis()
         
            #Run MCMC with MAP
            MAP_Model = pymc.MAP(MCMC_dict)
            MAP_Model.fit(method = 'fmin_powell') 
            
            M = pymc.MCMC(MCMC_dict, db = 'pickle', dbname =  FileFolder + pv.ScriptCode + '_' + CodeName + database_extension)
            M.sample(iter=30000, burn=5000, thin=10)
            M.write_csv(FileFolder + pv.ScriptCode + '_' + CodeName + globalfile_extension, variables=['ChiSq', 'He_abud', 'T_e', 'n_e', 'Tau', 'c_Hbeta', 'Xi'])
            M.db.close() 
                    
    except:
        pv.log_error(CodeName) 
    
    print i+1, '/' ,len(FilesList) 
    
print '\nAll data treated\n'
pv.display_errors()
