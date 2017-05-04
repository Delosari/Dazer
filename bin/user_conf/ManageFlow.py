def DataToTreat(Catalogue = 'WHT_observations'):
    
    Catalogue_Dictionary = {}
    
    if Catalogue == 'WHT_observations':
        Catalogue_Dictionary['Folder']      = '/home/vital/Dropbox/Astrophysics/Data/WHT_observations/'
        Catalogue_Dictionary['Datatype']    = 'WHT'
        Catalogue_Dictionary['Obj_Folder']  = '/home/vital/Dropbox/Astrophysics/Data/WHT_observations/' + 'objects/'
        Catalogue_Dictionary['Data_Folder'] = '/home/vital/Dropbox/Astrophysics/Data/WHT_observations/' + 'data/'
        Catalogue_Dictionary['dataframe']   = '/home/vital/Dropbox/Astrophysics/Data/WHT_observations/catalogue_df'
    
    if Catalogue == 'WHT_HII_Galaxies':
        Catalogue_Dictionary['Folder']      = '/home/vital/Dropbox/Astrophysics/Data/WHT_Catalogue_SulfurRegression/'
        Catalogue_Dictionary['Datatype']    = 'WHT'
        Catalogue_Dictionary['Obj_Folder']  = '/home/vital/Dropbox/Astrophysics/Data/WHT_Catalogue_SulfurRegression/' + 'Objects/SHOC579/'
        Catalogue_Dictionary['Data_Folder'] = '/home/vital/Dropbox/Astrophysics/Data/WHT_Catalogue_SulfurRegression/' + 'Data/'
        
    if Catalogue == 'WHT_CandiatesObjects':
        Catalogue_Dictionary['Folder']      = "Dropbox/Astrophysics/Data/WHT_CandiatesObjects/"
        Catalogue_Dictionary['Datatype']    = "dr10"
        Catalogue_Dictionary['Obj_Folder']  = "Dropbox/Astrophysics/Data/WHT_CandiatesObjects/"

    if Catalogue == 'WHT_CandiatesObjectsFabian':
        Catalogue_Dictionary['Folder']      = '/home/vital/Dropbox/Astrophysics/Data/Fabian_Catalogue/'
        Catalogue_Dictionary['Datatype']    = "dr10"
        Catalogue_Dictionary['Obj_Folder']  = '/home/vital/Dropbox/Astrophysics/Data/Fabian_Catalogue/'

    if Catalogue == 'Marta_Catalogue':
        Catalogue_Dictionary['Folder']      = "/home/vital/Dropbox/Astrophysics/Data/WHT_MartaCandidates_2016/"
        Catalogue_Dictionary['Datatype']    = "dr10" 
        Catalogue_Dictionary['Obj_Folder']  = "/home/vital/Dropbox/Astrophysics/Data/WHT_MartaCandidates_2016/Objects/"

    if Catalogue == 'SDSS_Catalogue':
        Catalogue_Dictionary['Folder']      = "Dropbox/Astrophysics/Data/Fabian_Catalogue/"
        Catalogue_Dictionary['Datatype']    = "dr10" 
        Catalogue_Dictionary['Obj_Folder']  = "Dropbox/Astrophysics/Data/Fabian_Catalogue/"

    if Catalogue == 'Testing_Pypeline':
        Catalogue_Dictionary['Folder']      = "Dropbox/Astrophysics/Data/ToCompare/"
        Catalogue_Dictionary['Datatype']    = "dr10" 
        Catalogue_Dictionary['Obj_Folder']  = "Dropbox/Astrophysics/Data/ToCompare/"

    return Catalogue_Dictionary








