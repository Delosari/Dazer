

from string import ascii_uppercase

new_row = list(ascii_uppercase)

for i in range(3):
    for letter in ascii_uppercase:
        new_row.append(ascii_uppercase[i] + letter)
        
print new_row

# import pyneb as pn
# from dazer_methods import Dazer
# from numpy import array
# from uncertainties.unumpy           import uarray, exp as unum_exp, log10 as unum_log10, pow as unum_pow
# from uncertainties.umath            import pow as umath_pow, log10 as umath_log10, exp as umath_exp
# 
# Rv = 3.4
# 
# rc_cardelli = pn.RedCorr(R_V = Rv, law='CCM89')
# rc_cardelli2 = pn.RedCorr(law='CCM89')
# 
# # rc_cardelli2.E_BV = 1.2
# # rc_cardelli2.R_V = 3.2
# # 
# # print rc_cardelli.cHbeta
# # print rc_cardelli2.cHbeta
# # print rc_gordon.cHbeta
# # 
# # print rc_cardelli.getCorr(4861.)
# # print rc_cardelli2.getCorr(4861.)
# # print rc_gordon.getCorr(4861.)
# # 
# # print rc_cardelli.getCorr(3726.)
# # print rc_cardelli2.getCorr(3726.)
# # print rc_gordon.getCorr(3726.)
# # 
# # print rc_cardelli.getCorr(6563.)
# # print rc_cardelli2.getCorr(6563.)
# # print rc_gordon.getCorr(6563.)
# # print '\n\n'
# 
# 
# #Generate dazer object
# dz = Dazer()
# dz.FigConf()
# 
# #Load catalogue dataframe
# catalogue_dict = dz.import_catalogue()
# catalogue_df = dz.load_dataframe(catalogue_dict['dataframe'])
# 
# #Declare data for the analisis
# AbundancesFileExtension = '_' + catalogue_dict['Datatype'] + '_linesLog_reduc.txt'
# cHbeta_type = 'cHbeta_reduc'
# 
# #Locate the objects
# objName = '8'
# ouput_folder = '{}{}/'.format(catalogue_dict['Obj_Folder'], objName) 
# lineslog_address = '{objfolder}{codeName}{lineslog_extension}'.format(objfolder = ouput_folder, codeName=objName, lineslog_extension=AbundancesFileExtension)
# 
# #Load lines frame
# lineslog_frame = dz.load_lineslog_frame(lineslog_address)
# 
# #Perform the reddening correction
# cHbeta      = catalogue_df.loc[objName, cHbeta_type]
# rc_gordon   = pn.RedCorr(cHbeta = cHbeta.nominal_value, R_V = Rv, law='G03 LMC')
# 
# dz.deredden_lines(lineslog_frame, reddening_curve = 'G03', cHbeta = cHbeta, R_v = Rv)
# print 'cHbeta', cHbeta
# 
# for label in ['O3_5007A', 'H1_4861A']:
#     print label
#     print lineslog_frame.loc[label, 'line_Flux']
#     print lineslog_frame.loc[label, 'line_Int'], '\n'
# 
# print 'my Xx', lineslog_frame.loc['O3_5007A', 'line_Xx'], 'pn Xx', rc_gordon.X(5007.)
# print 'mine', lineslog_frame.loc['O3_5007A', 'line_Int']/lineslog_frame.loc['O3_5007A', 'line_Flux'], 'pn', rc_gordon.getCorr(5007.), '\n'
# 
# print 'E(B-V) pn', rc_gordon.E_BV         
# print 'E(B-V) mine', dz.Ebv_from_cHbeta(cHbeta, reddening_curve= 'G03', R_v=3.4)
# 
# lineslog_frame['line_f'] = dz.flambda_from_Xx(lineslog_frame['line_Xx'], reddening_curve = 'G03', R_v=3.4)
# 
# print 'f_5007', lineslog_frame.loc['O3_5007A', 'line_f']
# 
# print 'Correction E(B-V)', lineslog_frame.loc['O3_5007A', 'line_Int'], lineslog_frame.loc['O3_5007A', 'line_Int']/lineslog_frame.loc['O3_5007A', 'line_Flux']
# 
# print 'Factoruco', umath_pow(10, cHbeta * (1 - lineslog_frame.loc['O3_5007A', 'line_f']))
# 
# print 'Correction f(l)', lineslog_frame.loc['O3_5007A', 'line_Flux'] * umath_pow(10, (1 - lineslog_frame.loc['O3_5007A', 'line_f']) * cHbeta)
# 
# 
# 
#         
#         