import  pyneb as pn
import  matplotlib.pyplot as plt
from    numpy import linspace, array
from    bin.dazer_methods import Dazer

rc = pn.RedCorr()
red_curve, Rv = 'G03', 2.5
rc.R_V = Rv
rc.law = 'G03 LMC'

dz = Dazer()

f, ax = plt.subplots(figsize = (10,10))
rc.plot(laws='G03 LMC', ax=ax)

wave    = linspace(3000, 7000, 1000)
Xx      = dz.reddening_Xx(wave, red_curve, Rv)
print 'Mine', Xx
print 'pn', rc.X(wave)
Xx_beta = dz.reddening_Xx(array([4861.331]), red_curve, Rv) 
print 'bicho', Xx_beta
ax.plot(wave,Xx,label='mine')
ax.plot(wave,Xx/Xx_beta,label='mine by Hbeta')
ax.plot(wave,Xx/2,label='mine by Hbeta')

plt.show()

from dazer_methods  import Dazer
from numpy          import nan as np_nan, array, power
from uncertainties.unumpy import pow as unum_pow
import pyneb as pn

#Generate dazer object
dz = Dazer()
dz.load_elements()

#Load catalogue dataframe
catalogue_dict  = dz.import_catalogue()
catalogue_df    = dz.load_excel_DF('/home/vital/Dropbox/Astrophysics/Data/WHT_observations/WHT_Galaxies_properties.xlsx')

#Declare data for the analisis
AbundancesFileExtension = '_' + catalogue_dict['Datatype'] + '_linesLog_emission_2nd.txt'
 
#Reddening properties
R_v         = 3.4
red_curve   = 'G03'
cHbeta_type = 'cHbeta_emis'

wavelengths_test = array([4363.0, 5007.0, 9069.0])

#Loop through objects:
for objName in catalogue_df.index:
        
#         try:
        if objName == '8':
        
            print '\n-Treating: {} {}/{}'.format(objName, catalogue_df.index.get_loc(objName), len(catalogue_df))
                        
            ouput_folder = '{}{}/'.format(catalogue_dict['Obj_Folder'], objName) 
            lineslog_address = '{objfolder}{codeName}{lineslog_extension}'.format(objfolder = ouput_folder, codeName=objName, lineslog_extension=AbundancesFileExtension)
                                       
            #Load lines frame
            lineslog_frame = dz.load_lineslog_frame(lineslog_address)
           
            #Perform the reddening correction
            cHbeta = catalogue_df.loc[objName, cHbeta_type]
            dz.deredden_lines(lineslog_frame, reddening_curve=red_curve, cHbeta=cHbeta, R_v=R_v)
            
            rc          = pn.RedCorr()
            rc.R_V      = 3.4
            rc.law      = 'G03 LMC'
            rc.cHbeta   = cHbeta.nominal_value
            corr_red    = rc.getCorr(wavelengths_test)
            print 'Estas seguro', rc.E_BV
            
            Xx_test     = dz.reddening_Xx(wavelengths_test, 'G03', R_v)
            Ebv_test    = dz.Ebv_from_cHbeta(cHbeta.nominal_value, 'G03', R_v)
            print 'Estas cosas'
            print 'pyneb',              rc.X(wavelengths_test)
            print 'mia',                Xx_test
            print 'E(B-V) mio',         Ebv_test
            print 'E(B-V) pyneb',       rc.EbvFromCHbeta(cHbeta.nominal_value)
            print 'Correction mia',     power(10, 0.4 *  Xx_test * Ebv_test)
            print 'Correction mia2',    10 ** (0.4 *  Xx_test * Ebv_test)
            print 'Correction mia3',    unum_pow(10,  0.4 * Xx_test * Ebv_test)
            print 'Correction pn',      corr_red
            
            print lineslog_frame.loc['O3_4363A', 'line_Flux'], '->', lineslog_frame.loc['O3_4363A', 'line_Int'], 'change', round((lineslog_frame.loc['O3_4363A', 'line_Int']/lineslog_frame.loc['O3_4363A', 'line_Flux']).nominal_value - 1, 3) * 100
            print lineslog_frame.loc['O3_5007A', 'line_Flux'], '->', lineslog_frame.loc['O3_5007A', 'line_Int'], 'change', round((lineslog_frame.loc['O3_5007A', 'line_Int']/lineslog_frame.loc['O3_5007A', 'line_Flux']).nominal_value - 1, 3) * 100
            print lineslog_frame.loc['S3_9069A', 'line_Flux'], '->', lineslog_frame.loc['S3_9069A', 'line_Int'], 'change', round((lineslog_frame.loc['S3_9069A', 'line_Int']/lineslog_frame.loc['S3_9069A', 'line_Flux']).nominal_value - 1, 3) * 100
