from lib.Astro_Libraries.Abundances_Class import Chemical_Analysis
from lib.Astro_Libraries.Abundances_Class import Parametrized_Emissivities


#Declare code classes
Article_data                = Parametrized_Emissivities()
Abund_Calculator            = Chemical_Analysis()

#Declare the parameters to calculate
EmLine_List                 = ['R_SII', 'R_SIIprime', 'R_SIII', 'R_NII', 'R_OII', 'R_OIII']
Den_List                    = ['nSII']
Temp_List                   = ['TOIII', 'TOII', 'TSII', 'TSIII','TNII', 'TOII_approx_TOIII', 'TSIII_approx_TOIII', 'TOIII_approx_TSIII', 'TNII_approx_TOIII']
IonicAbund_List             = ['SII_HII', 'SIII_HII', 'SIV_HII', 'OII_HII', 'OII_HII_3279A', 'OII_HII_7319A', 'NII_HII', 'HeII_HII', 'HeIII_HII', 'ArIII_HII', 'ArIV_HII']
Element_List                = ['SI_HI', 'SI_HI_ArCorr', 'OI_HI', 'NI_OI', 'NI_HI', 'HeI_HI']
Abund_Calculator.Properties_dict  = dict.fromkeys(Element_List + IonicAbund_List + Den_List + Temp_List + EmLine_List)

#Data dictionaries
para_dict                   = {}
data_dict                   = Article_data.Haegele2008_Fluxes()

keys, values = data_dict.keys(), data_dict.values()

for i in range(len(keys)):
    print keys[i], values[i]


#Sulfur Density
para_dict['TSIII']              = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'S3',     Parameter='TSIII')
data_dict['Temp']               = para_dict['TSIII']
para_dict['nSII']               = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'S3',     Parameter='nSII')
data_dict['Den']                = para_dict['nSII'] 

#Oxygen
para_dict['TOIII']              = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'O3',     Parameter='TOIII')
data_dict['Temp']               = para_dict['TOIII']
para_dict['TOII_approx_TOIII']  = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'O2',     Parameter='TOII_approx_TOIII')
para_dict['TOII']               = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Haegele2008', Ion = 'O2',     Parameter='TOII')
para_dict['OIII_HII']           = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion ='O3',      Parameter='OIII_HII')
data_dict['Temp']               = para_dict['TOII_approx_TOIII']
para_dict['OII_HII_3279A']      = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion ='O2',      Parameter='OII_HII_3279A')
para_dict['OII_HII_7319A']      = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Fabian2006',  Ion ='O2',      Parameter='OII_HII_7319A')

print '\nOxygen\n'

print 'Using density:',         para_dict['nSII'].nominal_value, para_dict['nSII'].std_dev
print 'TOIII',                  para_dict['TOIII'].nominal_value, para_dict['TOIII'].std_dev
print 'TOII_approx_TOIII',      para_dict['TOII_approx_TOIII'].nominal_value, para_dict['TOII_approx_TOIII'].std_dev
print 'TOII',                   para_dict['TOII'].nominal_value, para_dict['TOII'].std_dev
print 'OIII_HII',               para_dict['OIII_HII'].nominal_value, para_dict['OIII_HII'].std_dev
print 'OII_HII_3279A',          para_dict['OII_HII_3279A'].nominal_value, para_dict['OII_HII_3279A'].std_dev


#Sulfur
para_dict['TSIII']              = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'S3',     Parameter='TSIII')
para_dict['TSII']               = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Haegele2008', Ion = 'S2',     Parameter='TSII')
data_dict['Temp']               = para_dict['TSIII']
para_dict['nSII']               = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'S3',     Parameter='nSII')
data_dict['Den']                = para_dict['nSII'] 
para_dict['SIII_HII']           = Abund_Calculator.empiric_formulae(data_dict, methodology ='Haegele2008',  Ion='S3',       Parameter='SIII_HII')
data_dict['Temp']               = para_dict['TSII']
para_dict['SII_HII']            = Abund_Calculator.empiric_formulae(data_dict, methodology ='Haegele2008',  Ion='S2',       Parameter='SII_HII')

print '\nSulfur\n'
 
print 'T_SIII',                 para_dict['TSIII'].nominal_value, para_dict['TSIII'].std_dev
print 'T_SII',                  para_dict['TSII'].nominal_value, para_dict['TSII'].std_dev
print 'n_SII',                  para_dict['nSII'].nominal_value, para_dict['nSII'].std_dev
print 'SIII_HII',               para_dict['SIII_HII'].nominal_value, para_dict['SIII_HII'].std_dev
print 'SII_HII',                para_dict['SII_HII'].nominal_value, para_dict['SII_HII'].std_dev

#Argon
data_dict['Temp']               = para_dict['TSIII']
para_dict['ArIII_HII']          = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Haegele2008', Ion = 'Ar3',    Parameter='ArIII_HII')
data_dict['Temp']               = para_dict['TOIII']
para_dict['ArIV_HII']           = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Haegele2008', Ion = 'Ar4',    Parameter='ArIV_HII')

print '\nArgon\n'

print 'ArIII_HII',              para_dict['ArIII_HII'].nominal_value, para_dict['ArIII_HII'].std_dev
print 'ArIV_HII',               para_dict['ArIV_HII'].nominal_value, para_dict['ArIV_HII'].std_dev

print '\nSulfur IV\n'

data_dict['ArIII_HII']          = para_dict['ArIII_HII']
data_dict['ArIV_HII']           = para_dict['ArIV_HII']
data_dict['SIII_HII']           = para_dict['SIII_HII']
para_dict['SIV_HII']            = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Angeles2015', Ion = 'S4',    Parameter='SIV_HII')
print 'SIV_HII',                para_dict['SIV_HII'].nominal_value, para_dict['SIV_HII'].std_dev

#Nitrogen
print '\nNitrogen\n'

para_dict['TNII']               = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'N2',     Parameter='TNII')
data_dict['Temp']               = para_dict['TOIII']
para_dict['TNII_approx_TOIII']  = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'N2',     Parameter='TNII_approx_TOIII')
data_dict['Temp']               = para_dict['TNII_approx_TOIII']
para_dict['NII_HII']            = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Epm2014',     Ion = 'N2',     Parameter='NII_HII')

print 'T_NII',                  para_dict['TNII'].nominal_value, para_dict['TNII'].std_dev
print 'TNII_approx_TOIII',      para_dict['TNII_approx_TOIII'].nominal_value, para_dict['TNII_approx_TOIII'].std_dev
print 'NII_HII',                para_dict['NII_HII'].nominal_value, para_dict['NII_HII'].std_dev
print 'NI_OI',                  para_dict['NI_OI'].nominal_value, para_dict['NI_OI'].std_dev
print 'NI_HI',                  para_dict['NI_HI'].nominal_value, para_dict['NI_HI'].std_dev

print '\nHelium\n'

#Helium
data_dict['Temp']               = para_dict['TOIII']
data_dict['Den']                = para_dict['nSII'] 
para_dict['HeII_HII']           = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Fabian2006',  Ion = 'HeI',    Parameter='HeII_HII')
para_dict['HeIII_HII']          = Abund_Calculator.empiric_formulae(data_dict, methodology = 'Fabian2006',  Ion = 'HeII',   Parameter='HeIII_HII')

# print '\nHelium\n'
# 
print 'HeII_HII',               para_dict['HeII_HII'].nominal_value, para_dict['HeII_HII'].std_dev
print 'HeIII_HII',              para_dict['HeIII_HII'].nominal_value, para_dict['HeIII_HII'].std_dev

