from os                             import chdir
from shutil                         import copy
from subprocess                     import Popen, PIPE, STDOUT

from matplotlib                     import pyplot as plt
from matplotlib.cm                  import hot as cm_hot
from numpy                          import loadtxt, absolute, log10, zeros, floor, arange
from scipy.interpolate              import interp1d
from uncertainties                  import ufloat
from uncertainties.unumpy           import exp as unum_exp
from uncertainties.unumpy           import log10 as unum_log10

# from CodeTools.vitools              import FamilyOfItemsInArray
# from CodeTools.vitools              import File2Lines, LineFinder
# from CodeTools.vitools              import GetParameterFromDigitalLines
# from CodeTools.vitools              import Log_2_Parameter, Parameter_2_Log

# class AbundanceCalculator():

#     def __init__(self):
#         
#         # # Emision lines used to calculate abundances
#         # Hydrogen_Lines  =    ["H1_4861A"]
#         # Helium_Lines    =    ["He1_4472A","He2_4686A","He1_5876A","He1_6678A"]
#         # Oxygen_Lines    =    ["O3_4363A","O3_4959A","O3_5007A","O2_7319A"]
#         # Nitrogen_Lines  =    ["N2_6548A","N2_6584A"]
#         # Sulphur_Lines   =    ["S3_6312A","S2_6716A","S2_6731A","S3_9069A","S3_9531A"]
#         
#         #Files (THESE ATRIBUTES MIGHT BECOME CONFLICTED IF WE NEST THIS FUNCTION INSIDE ANOTHER ONE
#         self.FileFolder     = None
#         self.FileName       = None
#         self.CodeName       = None
#         
#         self.Flux_Header         = 'FluxGauss'
#         self.FluxError_Header    = 'ErrorEL_MCMC'        
#                      
#     def DefineObject(self, FileFolder, FileName, CodeName):
#         
#         self.FileFolder = FileFolder
#         self.FileName   = FileName
#         self.CodeName   = CodeName
#         
#         #Elemental Abundances
#         self.OI_HI              = None
#         self.NI_HI              = None    
#         self.SI_HI              = None    
#         self.HeI_HI             = None    
#         self.Y_Mass             = None
#         self.SI_HIAr            = None 
#        
#         #Ionic Abundances
#         self.OII_HII            = None
#         self.OIII_HII           = None
#         
#         self.SII_HII            = None
#         self.SIII_HII           = None
#         self.SIV_HII            = None
#         
#         self.HeII_HII           = None
#         self.HeIII_HII          = None
#         
#         self.NII_HII            = None
#              
#         #Line Parameters
#         self.R_OIII             = None           #(O3_4959A + O3_5007A) / O3_4363A
# 
#         self.R_SII              = None           #S2_6716A / S2_6731A
#         self.R_SIII             = None           #(S3_9069A + S3_9531A) / S3_6312A
#         self.R_SII_Prime        = None
# 
#         #Physical Parameters
#         self.T_OIII             = None
#         self.T_OII              = None
#         self.n_SII              = None
#         
#         #Checks
#         self.Total_MissingLines = []
#         self.MissingLines       = []
#         self.DisplayWarnings    = True
#         self.Methodology        = None
#         
#         self.ArIII_HII          = None
#         self.ArIV_HII           = None
#         
#     def Den_elec(self, Ion, Methodology = 'Fabian2006', Store = True, Temp0 = None):
# 
#         n_e                     = None
#         self.Methodology        = Methodology
# 
#         ReasonableParameters    = False
#         self.Methodology        = Methodology
# 
#         #Ion for which we calculate the parameter
#         if Ion == 'SII':
#             
#             #Scheme used to calculate the parameter
#             if Methodology == 'Fabian2006':
# 
#                 emLine                      = ["O3_4363A", "O3_4959A", "O3_5007A", "S2_6716A", "S2_6731A"]
#                 Flux, CalculationPosible    = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
# 
#                 if (self.T_OIII != None) and (self.T_OIII > 0):
#                     ReasonableParameters = True
#                 else:
#                     self.Errors_DeviantValues([self.T_OIII], ['T_OIII'])
#                 
#                 #Calculate the parameter
#                 if CalculationPosible and ReasonableParameters:
#                 
#                     self.R_OIII   = (Flux[emLine.index('O3_4959A')] + Flux[emLine.index('O3_5007A')]) / Flux[emLine.index('O3_4363A')]
#                                                 
#                     self.R_SII    = Flux[emLine.index('S2_6716A')] / Flux[emLine.index('S2_6731A')]
#                     t = self.T_OIII / 10000
#                     
#                     a_0     =  2.21 - 1.3  / t - 1.25 * t + 0.23 * t*t
#                     a_1     = -3.35 + 1.94 / t + 1.93 * t - 0.36 * t*t
#                     b_0     = -4.33 + 2.33 / t + 2.72 * t - 0.57 * t*t
#                     b_1     =  1.84 - 1.0  / t - 1.14 * t + 0.24 * t*t
#                                         
#                     n_SII   = (10.0**3) * (self.R_SII*a_0 + a_1) / (self.R_SII*b_0 + b_1)
#                     
#                     n_e     =  n_SII
#                  
#                 self.Errors_Density(n_e)
# 
#             if Methodology == 'Epm2014':
# 
#                 emLine                      = ["O3_4363A", "O3_4959A", "O3_5007A", "S2_6716A", "S2_6731A"]
#                 Flux, CalculationPosible    = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
# 
#                 if (self.T_OIII != None) and (self.T_OIII > 0):
#                     ReasonableParameters = True
#                 else:
#                     self.Errors_DeviantValues([self.T_OIII], ['T_OIII'])
#                 
#                 #Calculate the parameter
#                 if CalculationPosible and ReasonableParameters:
#                                                                 
#                     self.R_SII    = Flux[emLine.index('S2_6716A')] / Flux[emLine.index('S2_6731A')]
#                     
#                     t = self.T_OIII / 10000.0
#                     
#                     a_0     =   16.054 - 7.79/t - 11.32 * t
#                     a_1     =   -22.66 + 11.08/t + 16.02 * t
#                     b_0     =   -21.61 + 11.89/t + 14.59 * t
#                     b_1     =   9.17 - 5.09/t - 6.18 * t
#                                         
#                     n_SII   = (10.0**3) * (self.R_SII*a_0 + a_1) / (self.R_SII*b_0 + b_1)
#                     
#                     n_e     =  n_SII
#                  
#                 self.Errors_Density(n_e)
#    
#         return n_e
#      
#     def Tem_elec(self, Ion, Methodology, Store = True, den = None):
#        
#         T_e                     =  None
#         self.Methodology        = Methodology
#         
#         ReasonableParameters    = False
#                 
#         if Ion == 'OIII':
#             if Methodology == 'Fabian2006':
#                 
#                 #Declare the emission lines we are going to need
#                 emLine                      = ["O3_4363A", "O3_4959A", "O3_5007A"]
#                 Flux, CalculationPosible    = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
#                 
#                 #Calculate the parameter
#                 if CalculationPosible:
#                     self.R_OIII      = (Flux[emLine.index('O3_4959A')] + Flux[emLine.index('O3_5007A')]) / Flux[emLine.index('O3_4363A')]
#                     T_OIII          = 0.8254 - 0.0002415 * self.R_OIII + (47.77 / self.R_OIII)
#                                         
#                     T_e = T_OIII * 10000
#                     
#         elif Ion == 'OII':
#             if Methodology == 'Fabian2006':
#                 
#                 emLine                      = ["O3_4363A", "O3_4959A", "O3_5007A"]
#                 Flux, CalculationPosible    = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)   
#                      
#                 #Security check against extreme conditions
#                 if (self.T_OIII != None) and (self.n_SII > 0):  
#                     ReasonableParameters = True
#                 else:
#                     self.Errors_DeviantValues([self.T_OIII, self.n_SII], ['T_OIII', 'n_SII'])
#                 
#                 #Calculate the parameter
#                 if CalculationPosible and ReasonableParameters:
#                     T_OII   = (1.2 + 0.002 * self.n_SII + 4.2/self.n_SII) / (1/(self.T_OIII/10000) + 0.08 + 0.003 * self.n_SII + 2.5/self.n_SII)
#                     T_e     = T_OII * 10000
#     
#         elif Ion == 'SII':
#             if Methodology == 'Haegele2008':
#         
#                 emLine                      =  ["S2_6716A","S2_6731A"] + ["S2_4069A","S2_4076A"]
#                 Flux, CalculationPosible    = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)    
#                 Flux, CalculationPosible    = self.SulfurII_MissingLines(emLine, Flux)                                 #Special case when we have one of the 4000 SII lines
#                                                       
#                 #Security check against extreme conditions
#                 if (self.n_SII > 0):  
#                     ReasonableParameters = True
#                 else:
#                     self.Errors_DeviantValues([self.n_SII], ['n_SII'])
#                         
#                 #Calculate the parameter
#                 if ReasonableParameters:
#                     self.R_SII       = Flux[emLine.index('S2_6716A')] / Flux[emLine.index('S2_6731A')]
#                     self.R_SII_Prime = (Flux[emLine.index('S2_6716A')] + Flux[emLine.index('S2_6731A')]) / (Flux[emLine.index('S2_4069A')] + (Flux[emLine.index('S2_4076A')])) 
#                                         
#                     T_SII   = 1.92 - 0.0375 * self.R_SII_Prime - (14.5 / self.R_SII_Prime) + (105.64 / (self.R_SII_Prime**2))    # Need to add a missing function here         1    
#                     
#                     T_e     = T_SII * 10000
#                     
#         elif Ion == 'SIII':
#             if Methodology == 'Haegele2008':
#                 emLine                      = ["S3_9069A", "S3_9531A", "S3_6312A"]
#                 Flux, CalculationPosible    = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
#                 Flux, CalculationPosible    = self.SulfurIII_MissingLines(emLine, Flux)                                      #Special case when we have one of the 9000 SIII lines
#             print 'Possible?', CalculationPosible
#             if CalculationPosible:
#                 self.R_SIII = (Flux[emLine.index('S3_9069A')] + Flux[emLine.index('S3_9531A')]) / Flux[emLine.index('S3_6312A')]
#                 T_SIII      = (self.R_SIII + 36.4) / (1.8 * self.R_SIII - 3.01)
#                 
#                 T_e = T_SIII * 10000
#     
#     
#         return T_e
#      
#     def Oxygen_Abundance(self, CodeName, FileName, FileFolder, Methodology):
#                 
#         if Methodology == 'Fabian2006':    
#             
#             #Calculate physical properties
#             self.T_OIII     = self.Tem_elec(Ion = 'OIII', Methodology = "Fabian2006")
#             self.n_SII      = self.Den_elec(Ion = 'SII', Methodology = "Fabian2006")
#             self.T_OII      = self.Tem_elec(Ion = 'OII', Methodology = "Fabian2006")
# 
#             #Get line intensities corrected from intrinsic reddening
#             emLine                          = ["O3_4959A","O3_5007A","H1_4861A"]
#             Flux, CalculationPosible        = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
# 
#             #Check for reasonable parameters
#             ReasonableParameters = False
#             if (self.T_OIII != None) and (self.n_SII > 0) and (self.T_OII != None): 
#                 ReasonableParameters = True
#             else:
#                 self.Errors_DeviantValues([self.T_OIII, self.n_SII, self.T_OII], ['T_OIII', 'n_SII', 'T_OII'])
#             
#             #Calculate OIII_HII
#             if (ReasonableParameters and CalculationPosible):
#                 t_3 = self.T_OIII / 10000
#                 logOIII_logHII  = -12 + unum_log10((Flux[emLine.index('O3_4959A')] + Flux[emLine.index('O3_5007A')]) / Flux[emLine.index('H1_4861A')] ) + 6.144 + 1.251 / t_3 - 0.55 * unum_log10(t_3)
#                 self.OIII_HII   = 10**(logOIII_logHII)
#                         
#             #Calculating OII_HII
#             emLine                          = ["O2_7319A","H1_4861A"]
#             Flux, CalculationPosible        = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)           
# 
#             ReasonableParameters = False
#             if (self.n_SII > 0) and (self.T_OII != None): 
#                 ReasonableParameters = True
#                 
#             if (ReasonableParameters and CalculationPosible):
#                 logOII_logHII = -12 + unum_log10(Flux[emLine.index('O2_7319A')] / Flux[emLine.index('H1_4861A')]) + 6.895 + 2.44/(self.T_OII/10000) - 0.58*unum_log10(self.T_OII/10000) - unum_log10(1.0 + 0.0047 * self.n_SII)
#                 self.OII_HII = 10**(logOII_logHII)
#                 
#             #Calculating OI_HI
#             if (self.OIII_HII != None) and (self.OII_HII != None):
#                 self.OI_HI = self.OII_HII + self.OIII_HII
# 
#         
#         return self.OI_HI
#         
#     def Nitrogen_Abundance(self, CodeName, FileName, FileFolder, Methodology):
#         
#         self.NI_HI      = None
#         self.NII_HII    = None
#     
#         if Methodology == 'Fabian2006':
#             
#             self.T_OIII  = self.Tem_elec(Ion = 'OIII', Methodology = Methodology)
#             self.n_SII  = self.Den_elec(Ion = 'SII', Methodology = Methodology)
#             self.T_OII  = self.Tem_elec(Ion = 'OII', Methodology = Methodology)
#             
#             #Calculating NII_HII
#             emLine                          = ["N2_6548A","N2_6584A", "H1_4861A"]
#             Flux, CalculationPosible        = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
#             
#             #Check reasonable conditions
#             ReasonableParameters = False
#             if (self.T_OII != None) and (self.OII_HII != None) and (self.OI_HI != None): 
#                 ReasonableParameters = True   
#             else:
#                 self.Errors_DeviantValues([self.T_OII, self.OII_HII, self.OI_HI], ['T_OII', 'OII_HII', 'OI_HI'])   
#           
#             if (ReasonableParameters and CalculationPosible):
#                 t_2 = self.T_OII / 10000
#                 logNII_logHII = -12 + unum_log10((Flux[emLine.index('N2_6548A')] + Flux[emLine.index('N2_6584A')])/ Flux[emLine.index('H1_4861A')]) + 6.273 + 0.894/t_2 - 0.592*unum_log10(t_2)
#                 self.NII_HII  = 10**(logNII_logHII)
#         
#                 #Calculating NI_HI:
#                 if (self.OI_HI != None) and (self.OII_HII != None):
#                     
#                     NI_OI = (self.NII_HII) / (self.OII_HII) 
#                     
#                     self.NI_HI = NI_OI * self.OI_HI 
#             
#         return self.NI_HI
#         
#     def Argon_Abundance(self, CodeName, FileName, FileFolder, Methodology):
#         
#         self.ArIII_HII  = None
#         self.ArIV_HII   = None
#         
#         if Methodology == "EMontero":
#             self.T_OIII     = self.Tem_elec('OIII', 'Fabian2006')
#             self.n_SII      = self.Den_elec('SII', 'Fabian2006')
#             self.T_SIII     = self.Tem_elec('SIII', 'Haegele2008')
# 
#             emLine                           = ["Ar4_4740A","Ar3_7136A"] + ["S3_9069A","S3_9531A","H1_4861A", "S3_6312A"] + ["S2_6716A","S2_6731A"] + ["S2_4069A","S2_4076A"] 
#             Flux, Missing_Lines_Check        = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
#             Flux, Missing_Lines_Check        = self.ArgonCheck_MissingLines(emLine, Flux)  #Special case when we have one of the 9000 SIII line            
#         
#             ReasonableParameters = False
#             if (self.T_SIII != None) and (self.n_SII != None) and (self.T_OIII != None):
#                 ReasonableParameters = True
#             else:
#                 self.Errors_DeviantValues([self.T_SIII, self.n_SII], ['T_SIII', 'n_SII'])           
#             
#             if Missing_Lines_Check and ReasonableParameters:
#                 t_o3    = self.T_OIII / 10000
#                 t_s3    = self.T_SIII / 10000
#                 logArIII_HII    = -12 + unum_log10(Flux[emLine.index('Ar3_7136A')] / Flux[emLine.index('H1_4861A')]) + 6.157 + 0.808/t_s3 - 0.508 * unum_log10(t_s3) 
#                 self.ArIII_HII  = 10**logArIII_HII
#                 logArIV_HII     = -12 + unum_log10(Flux[emLine.index('Ar4_4740A')] / Flux[emLine.index('H1_4861A')]) + 4.705 + 1.246/t_o3 - 0.156 * unum_log10(t_o3)
#                 self.ArIV_HII   = 10**logArIV_HII
#                 
#         return 
#         
#     def Sulfur_Abundance(self, CodeName, FileName, FileFolder, Methodology):
#         self.SI_HI      = None
#         self.SI_HIAr    = None 
#         self.SII_HII    = None
#         self.SIII_HII   = None
#         self.SIV_HII    = None
#         
#         if Methodology == 'Haegele2008':
#             
#             #Calculating SIII_HII
#             self.T_OIII     = self.Tem_elec('OIII', 'Fabian2006')
#             self.n_SII      = self.Den_elec('SII', 'Fabian2006')
# 
#             self.T_SIII     = self.Tem_elec('SIII', 'Haegele2008')
#             self.T_SII      = self.Tem_elec('SII', 'Haegele2008')
# 
#             emLine                           = ["S3_9069A","S3_9531A","H1_4861A", "S3_6312A"] + ["S2_6716A","S2_6731A"] + ["S2_4069A","S2_4076A"] 
#             Flux, Missing_Lines_Check        = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
#             Flux, Missing_Lines_Check        = self.SulfurII_and_SulfurIII_MissingLines(emLine, Flux)  #Special case when we have one of the 9000 SIII lines
#                         
#             ReasonableParameters = False
#             if (self.T_SIII != None) and (self.T_SII != None) and (self.n_SII != None):
#                 ReasonableParameters = True
#             else:
#                 self.Errors_DeviantValues([self.T_SIII, self.T_SII, self.n_SII], ['T_SIII', 'T_SII', 'n_SII'])   
#                         
#             if Missing_Lines_Check and ReasonableParameters:
#                
#                 t_3 = self.T_SIII / 10000
#                 t_2 = self.T_SII / 10000
#                                 
#                 logSIII_logHII  = -12 + unum_log10((Flux[emLine.index('S3_9069A')] + Flux[emLine.index('S3_9531A')]) / Flux[emLine.index('H1_4861A')]) + 5.8 + 0.77/t_3 - 0.22 * unum_log10(t_3) 
#                 self.SIII_HII   = 10**logSIII_logHII
# 
#                 logSII_logHII   = -12 + unum_log10((Flux[emLine.index('S2_6716A')] + Flux[emLine.index('S2_6731A')]) / Flux[emLine.index('H1_4861A')]) + 5.423 + 0.929/t_2 - 0.280 * unum_log10(t_2) + unum_log10(1 + 0.0001 * self.n_SII)
#                 self.SII_HII    = 10**logSII_logHII
#                 
#             #Calculating SII_HII
#             if (self.SII_HII != None) and (self.SIII_HII != None):
#                 self.SI_HI = self.SII_HII + self.SIII_HII
# 
#         if Methodology == 'Angeles2014':
#             
#             #Calculating SIII_HII
#             self.T_OIII     = self.Tem_elec('OIII', 'Fabian2006')
#             self.n_SII      = self.Den_elec('SII', 'Fabian2006')
# 
#             self.T_SIII     = self.Tem_elec('SIII', 'Haegele2008')
#             self.T_SII      = self.T_SIII 
# 
#             emLine                          = ["S3_9069A","S3_9531A","H1_4861A", "S3_6312A"] + ["S2_6716A","S2_6731A"]
#             Flux, Missing_Lines_Check        = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
#             Flux, Missing_Lines_Check        = self.SulfurIII_MissingLines(emLine, Flux)                                      #Special case when we have one of the 9000 SIII lines
#             
#             ReasonableParameters = False
#             if (self.T_SIII != None) and (self.T_SII != None) and (self.n_SII != None):
#                 ReasonableParameters = True
#             else:
#                 self.Errors_DeviantValues([self.T_SIII, self.T_SII, self.n_SII], ['T_SIII', 'T_SII', 'n_SII'])   
#                 
#             if Missing_Lines_Check and ReasonableParameters:
#                 t_3 = self.T_SIII / 10000
#                 t_2 = self.T_SII / 10000
#                                 
#                 logSIII_logHII  = -12 + unum_log10((Flux[emLine.index('S3_9069A')] + Flux[emLine.index('S3_9531A')]) / Flux[emLine.index('H1_4861A')]) + 5.8 + 0.77/t_3 - 0.22 * unum_log10(t_3) 
#                 self.SIII_HII   = 10**logSIII_logHII
# 
#                 logSII_logHII   = -12 + unum_log10((Flux[emLine.index('S2_6716A')] + Flux[emLine.index('S2_6731A')]) / Flux[emLine.index('H1_4861A')]) + 5.423 + 0.929/t_2 - 0.280 * unum_log10(t_2) + unum_log10(1 + 0.0001 * self.n_SII)
#                 self.SII_HII    = 10**logSII_logHII
#                 
#             #Calculating SII_HII
#             if (self.SII_HII != None) and (self.SIII_HII != None):
#                 self.SI_HI = self.SII_HII + self.SIII_HII        
#                 
#                 if self.ArIII_HII != None and self.ArIV_HII != None:
#                     
#                     y_Ar = unum_log10(self.ArIII_HII/self.ArIV_HII)
#                     
#                     x_S = (y_Ar - 0.126) / 1.192
#                     
#                     self.SIV_HII = self.SIII_HII / 10**(x_S)
#                     
#                     self.SI_HIAr = self.SII_HII + self.SIII_HII + self.SIV_HII
#                     
#                 else:
#                                         
#                     self.SI_HIAr = self.SII_HII + self.SIII_HII
#                     
#         return self.SI_HI
#         
#     def Helium_Abundance(self,  CodeName, FileName, FileFolder, Methodology):
#     
#         if Methodology == 'Fabian2006':
# 
#             self.T_OIII  = self.Tem_elec(Ion = 'OIII', Methodology = Methodology)
#             self.n_SII  = self.Den_elec(Ion = 'SII', Methodology = Methodology)
# 
#             emLine                          = ["He1_4472A","He1_5876A","He1_6678A", "He2_4686A", "H1_4861A"]
#             Flux, CalculationPosible        = self.GetFluxes(emLine, self.FileFolder, self.FileName, self.CodeName)
#             CalculationPosible              = self.HeII_MissingCheck(Flux, emLine)
# 
#             ReasonableParameters = False
#             if (self.T_OIII != None) and (self.n_SII > 0): 
#                 ReasonableParameters = True
# 
#             if (ReasonableParameters and CalculationPosible):
#                 
#                 t = self.T_OIII / 10000
#                                                
#                 y0_I4471 = 2.04   * (t**0.13) * Flux[emLine.index('He1_4472A')]  / Flux[emLine.index('H1_4861A')]
#                 y0_I5876 = 0.738  * (t**0.23) * Flux[emLine.index('He1_5876A')]  / Flux[emLine.index('H1_4861A')]
#                 y0_I6678 = 2.58   * (t**0.25) * Flux[emLine.index('He1_6678A')]  / Flux[emLine.index('H1_4861A')]
#                 
#                 D = 1.0 + 3110.0  * (t**-0.51) * (1.0/self.n_SII)
#                 g_I4471 = 6.11    * (t**0.02)  * (unum_exp(-4.544)) / D
#                 g_I5876 = (7.12   * (t**0.14)  * (unum_exp(-3.776/t)) + 1.47 * (t**-0.28) * (unum_exp(-4.544/t)))/ D
#                 g_I6678 = (3.27   * (t**-0.41) * (unum_exp(-3.777/t)) + 0.49 * (t**-0.52) * (unum_exp(-4.544/t)))/ D
#                 
#                 y_I4471 = y0_I4471 / (1.0 + g_I4471)
#                 y_I5876 = y0_I5876 / (1.0 + g_I5876)
#                 y_I6678 = y0_I6678 / (1.0 + g_I6678)
#                 
#                 #Can we calculate with only one ionization?
#                 self.HeII_HII   = (3.0/5.0)*(y_I5876 + (1.0/3.0)*y_I4471 + (1.0/3.0)*y_I6678)
#                 
#                 if Flux[emLine.index('He2_4686A')] != None:
#                     y_PPI4686 = 0.084 * (t**0.14) * Flux[emLine.index('He2_4686A')]  / Flux[emLine.index('H1_4861A')]
#                     self.HeIII_HII  = y_PPI4686
#                 else:
#                     self.HeIII_HII = 0.0
#                  
#                 self.HeI_HI     = self.HeII_HII  + self.HeIII_HII
#                         
#         return self.HeI_HI        
# 
#     def Y_Calculation(self):
#         
#         if (self.HeI_HI != None) and (self.OI_HI != None): 
#             self.Y_Mass = (4.0 * self.HeI_HI*(1.0 - 20.0 * self.OI_HI )) / (1.0 + 4.0 * self.HeI_HI)
# 
#         else:
#             self.Y_Mass = None
#             
#         return self.Y_Mass
# 
#     def GetFluxes(self, emLine, FileFolder, FileName, CodeName):
# 
#         self.MissingLines = []
#         
#         if isinstance(emLine, list):
#             Flux    = [0] * len(emLine)
#         
#         else:
#             emLine = [emLine]
#             Flux    = [0] * len(emLine)
#         
#         #Check if we have all the lines to calculate this parameter
#         CalculationPosible = True
#         for i in range(len(emLine)):
#             Flux[i] = self.GetLineFlux(emLine[i], FileFolder, FileName, CodeName)
#             
#             if Flux[i] == None:
#                 self.MissingLines.append(emLine[i])
#                 CalculationPosible = False
#                                 
#         self.Errors_MissingLines(emLine, Flux)
#        
#         return Flux, CalculationPosible
#           
#     def GetLineFlux(self, Line_Label, FolderName, FileName, CodeName):
#         # 'NOT GETTING THE ERROR FROM CHBETA!!!'
#         # NEGATIVE values error??
#         R_V = 3.2
#                 
#         Valid_Flag = True
#     
#         LogFileName = FileName
#         LogLines    = File2Lines(FolderName, LogFileName)
#         
#         Flux_Brute  = GetParameterFromDigitalLines(LogLines,Line_Label,self.Flux_Header,2)
#         Flux_Error  = GetParameterFromDigitalLines(LogLines,Line_Label,self.FluxError_Header,2)
#         
#         Wavelength  = float(GetParameterFromDigitalLines(LogLines,Line_Label,"TheoWavelength",2))
#         Wave_Error  = 1.0
# 
# 
#         if (Flux_Brute != None) and (Flux_Error != None):
#             
#             Flux_V      = ufloat(float(Flux_Brute), absolute(float((Flux_Error))))
#             Wave_V      = ufloat(float(Wavelength), float(Wave_Error))
#              
#             f_HBeta     = Reddening_f_Cal(4861.333, R_V)
#             f_lambda    = Reddening_f_Cal(Wavelength, R_V)
# 
#             cHBeta      = Log_2_Parameter(FolderName, CodeName, "cHBeta")
#             
#             if cHBeta != '-':
#                 cHBeta = float(cHBeta)
#             
#                 if cHBeta < 0:
#                     Flux_Deredd = Flux_V
#                 else:
#                     Flux_Deredd = Flux_V * 10**(cHBeta*(f_lambda-f_HBeta))
#                     
#             else:
#                 Flux_Deredd = Flux_V
#             
#             LineData = Flux_Deredd
#                         
#         else:
#             Valid_Flag = False
#             LineData = None
#             
#         return LineData        
#         
#     def Errors_MissingLines(self, emLine, Flux, DisplayErrorMessages = True):
#         
#         if self.DisplayWarnings == True:
#             if len(emLine) != len(Flux):
#                 print '--Methodology ' + "'" + self.Methodology + "'" +' not possible:'
#                 print '--- Missing lines:'
#                 for line in self.MissingLines:
#                         print '   ', line
# 
#     def Errors_Density(self, n_e):
#         if n_e < 100:
#             print '--Warning: Density very low'
#             if n_e < 0:
#                 print '--ZERO value'   
#      
#     def Errors_DeviantValues(self, mylist, list_Names = None):
#         print '--Deviant values encountered'
# 
#         if list_Names != None:
#             for i in range(len(mylist)):
#                 print '   ', list_Names[i], mylist[i]
#         else:
#             for i in mylist:
#                 print '   ', i
# 
#     def SulfurII_and_SulfurIII_MissingLines(self, emLine, Flux):
# 
#         Check1 = False
#         if (('S3_9069A' in self.MissingLines) and ('S3_9531A' not in self.MissingLines)):
#             Flux[emLine.index('S3_9069A')] = Flux[emLine.index('S3_9531A')] / 2.4688
#             self.MissingLines.remove('S3_9069A') 
#             Check1 = True
#             
#         elif (('S3_9069A' not in self.MissingLines) and ('S3_9531A' in self.MissingLines)):
#             Flux[emLine.index('S3_9531A')] = Flux[emLine.index('S3_9069A')] * 2.4688
#             self.MissingLines.remove('S3_9531A') 
#             Check1 = True
#             
#         elif (('S3_9531A' not in self.MissingLines) and ('S3_9531A' not in self.MissingLines)):
#             Check1 = True
# 
#         Check2 = False
#         if (('S2_4069A' in self.MissingLines) and ('S2_4076A' not in self.MissingLines)):
#             Flux[emLine.index('S2_4069A')] = Flux[emLine.index('S2_4076A')] * 3.11
#             self.MissingLines.remove('S2_4069A') 
#             Check2 = True
#             
#         elif (('S2_4069A' not in self.MissingLines) and ('S2_4076A' in self.MissingLines)):
#             Flux[emLine.index('S2_4076A')] = Flux[emLine.index('S2_4069A')] / 3.11
#             self.MissingLines.remove('S2_4076A') 
#             Check2 = True      
#             
#         elif ('S2_4069A' not in self.MissingLines) and ('S2_4076A' not in self.MissingLines):
#             Check2 = True
# 
#         Check = True
#         if Check1 and Check2:
#             for i in range(len(emLine)):
#                 if Flux[i] == 0:
#                     Check = False 
#         else:
#             Check = False
#             
#         return Flux, Check
#     
#     def ArgonCheck_MissingLines(self, emLine, Flux):
# 
#         Check1 = False
#         if (('S3_9069A' in self.MissingLines) and ('S3_9531A' not in self.MissingLines)):
#             Flux[emLine.index('S3_9069A')] = Flux[emLine.index('S3_9531A')] / 2.4688
#             self.MissingLines.remove('S3_9069A') 
#             Check1 = True
#             
#         elif (('S3_9069A' not in self.MissingLines) and ('S3_9531A' in self.MissingLines)):
#             Flux[emLine.index('S3_9531A')] = Flux[emLine.index('S3_9069A')] * 2.4688
#             self.MissingLines.remove('S3_9531A') 
#             Check1 = True
#             
#         elif (('S3_9531A' not in self.MissingLines) and ('S3_9531A' not in self.MissingLines)):
#             Check1 = True
# 
#         Check2 = False
#         if (('S2_4069A' in self.MissingLines) and ('S2_4076A' not in self.MissingLines)):
#             Flux[emLine.index('S2_4069A')] = Flux[emLine.index('S2_4076A')] * 3.11
#             self.MissingLines.remove('S2_4069A') 
#             Check2 = True
#             
#         elif (('S2_4069A' not in self.MissingLines) and ('S2_4076A' in self.MissingLines)):
#             Flux[emLine.index('S2_4076A')] = Flux[emLine.index('S2_4069A')] / 3.11
#             self.MissingLines.remove('S2_4076A') 
#             Check2 = True      
#             
#         elif ('S2_4069A' not in self.MissingLines) and ('S2_4076A' not in self.MissingLines):
#             Check2 = True
# 
#         Check = True
#         if Check1 and Check2:
#             for i in range(len(emLine)):
#                 if Flux[i] == 0 or Flux[i] == 'None' or Flux[i] == None:
#                     Check = False 
#         else:
#             Check = False
#             
#         return Flux, Check
#     
#     def SulfurIII_MissingLines(self, emLine, Flux):
#             
#         Check1 = False
#         if (('S3_9069A' in self.MissingLines) and ('S3_9531A' not in self.MissingLines)):
#             Flux[emLine.index('S3_9069A')] = Flux[emLine.index('S3_9531A')] / 2.4688
#             self.MissingLines.remove('S3_9069A') 
#             Check1 = True
#             
#         elif (('S3_9069A' not in self.MissingLines) and ('S3_9531A' in self.MissingLines)):
#             Flux[emLine.index('S3_9531A')] = Flux[emLine.index('S3_9069A')] * 2.4688
#             self.MissingLines.remove('S3_9531A') 
#             Check1 = True
#             
#         elif (('S3_9531A' not in self.MissingLines) and ('S3_9531A' not in self.MissingLines)):
#             Check1 = True
#                         
#         return Flux, Check1
# 
#     def SulfurII_MissingLines(self, emLine, Flux):
#             
#         Check2 = False
#         if (('S2_4069A' in self.MissingLines) and ('S2_4076A' not in self.MissingLines)):
#             Flux[emLine.index('S2_4069A')] = Flux[emLine.index('S2_4076A')] * 3.11
#             self.MissingLines.remove('S2_4069A') 
#             Check2 = True
#             
#         elif (('S2_4069A' not in self.MissingLines) and ('S2_4076A' in self.MissingLines)):
#             Flux[emLine.index('S2_4076A')] = Flux[emLine.index('S2_4069A')] / 3.11
#             self.MissingLines.remove('S2_4076A') 
#             Check2 = True      
#             
#         elif ('S2_4069A' not in self.MissingLines) and ('S2_4076A' not in self.MissingLines):
#             Check2 = True
#             
#         return Flux, Check2
#    
#     def HeII_MissingCheck(self, Flux, emLine):
#         
#         Check = True
#         
#         for i in range(len(emLine)):
#             if Flux[i] == None and emLine != "He2_4686A":
#                 Check = False
#         
#         return Check
#         
#     def StoreValues(self, FileFolder, CodeName,  DisplayOnScreen = True, StoreData = False):
#         
#         if DisplayOnScreen:
#             print "\n", "\n", "\n", "\n", "\n", "\n"
#             print '-Physical properties'
#             print 'T_OII', self.T_OII
#             print 'T_OIII', self.T_OIII
#             print 'T_SII', self.T_SII
#             print 'T_SIII', self.T_SIII
#             print 'n_SII', self.n_SII,'\n'
#             
#             print '-Element abundances'
#             print 'OI_HI', self.OI_HI.nominal_value
#             print 'NI_HI', self.NI_HI
#             print 'SI_HI', self.SI_HI
#             print 'SI_HI', self.SI_HI
#             print 'SI_HI_Ar', self.SI_HIAr
#             print 'Y Mass', self.Y_Mass, '\n'
#             
#             print '-Ionic abundances'
#             print 'HeII_HII', self.HeII_HII
#             print 'HeIII_HII', self.HeIII_HII            
#             print 'OII_HII', self.OII_HII
#             print 'OIII_HII', self.OIII_HII
#             print 'NII_HII', self.NII_HII
#             print 'SII_HII', self.SII_HII
#             print 'SIII_HII', self.SIII_HII
#             print 'SIV_HII', self.SIV_HII
#             print 'ArIII_HII', self.ArIII_HII
#             print 'ArIV_HII', self.ArIV_HII, '\n'
#             
#         if StoreData:
#             Parameter_2_Log(FileFolder, CodeName, "te_OII",     PrepareToStore(self.T_OII))
#             Parameter_2_Log(FileFolder, CodeName, "te_OIII",    PrepareToStore(self.T_OIII))
#             Parameter_2_Log(FileFolder, CodeName, "te_SII",     PrepareToStore(self.T_SII))
#             Parameter_2_Log(FileFolder, CodeName, "te_SII",     PrepareToStore(self.T_SIII))
#             Parameter_2_Log(FileFolder, CodeName, "ne_SII",     PrepareToStore(self.n_SII))
#             
#             Parameter_2_Log(FileFolder, CodeName, "OI_HI",      PrepareToStore(self.OI_HI))
#             Parameter_2_Log(FileFolder, CodeName, "NI_HI",      PrepareToStore(self.NI_HI))
#             Parameter_2_Log(FileFolder, CodeName, "SI_HI",      PrepareToStore(self.SI_HI))
#             Parameter_2_Log(FileFolder, CodeName, "HeI_HI",     PrepareToStore(self.HeI_HI))
#             Parameter_2_Log(FileFolder, CodeName, "Y_Mass",     PrepareToStore(self.Y_Mass))   
#             Parameter_2_Log(FileFolder, CodeName, "SI_HI_Ar",   PrepareToStore(self.SI_HIAr))
# 
#             Parameter_2_Log(FileFolder, CodeName, "HeII_HII",   PrepareToStore(self.HeII_HII))
#             Parameter_2_Log(FileFolder, CodeName, "HeIII_HII",  PrepareToStore(self.HeIII_HII))            
#             Parameter_2_Log(FileFolder, CodeName, "OII_HII",    PrepareToStore(self.OII_HII))
#             Parameter_2_Log(FileFolder, CodeName, "OIII_HII",   PrepareToStore(self.OIII_HII))
#             Parameter_2_Log(FileFolder, CodeName, "NII_HII",    PrepareToStore(self.NII_HII))
#             Parameter_2_Log(FileFolder, CodeName, "SII_HII",    PrepareToStore(self.SII_HII))
#             Parameter_2_Log(FileFolder, CodeName, "SIII_HII",   PrepareToStore(self.SIII_HII))
#             Parameter_2_Log(FileFolder, CodeName, "SIV_HII",    PrepareToStore(self.SIV_HII))


def PrepareToStore(Parameter):
        
    if Parameter != None:
        Data = str(Parameter.nominal_value) + '+/-' + str(Parameter.std_dev)
    else:
        Data = str(Parameter)

    return Data
                
def GenerateStarlightFiles(FileFolder, FileName, CodeName, X, Y, ExtraData = None, Velocity_Vector = ['FIT',  '0.0',  "10.0"], ComputerRoot = '/home/vital/', EmLinesFileExtension = "LickIndexes.txt", TwoArmsMode = False):
    #Velocity_Vector =                                                                               [Kinematics calculation type,    v0_start,    vd_start]
    
    Starlight_Folder    = ComputerRoot + 'Starlight/'
    Default_InputFolder = ComputerRoot + 'Starlight/Obs/'
    Default_MaskFolder  = ComputerRoot + 'Starlight/Masks/'
    Default_BasesFolder = ComputerRoot + 'Starlight/Bases/'
    Default_OutputFoler = ComputerRoot + 'Starlight/Output/'

    if TwoArmsMode == True:
        Color = None
        if 'Blue' in FileName:
            Color = 'Blue'
        elif 'Red' in FileName:
            Color = 'Red'
        
    #-----------------------     Generating Base File    ----------------------------------------------
    BaseDataFile         = 'Dani_Bases_296.txt'

    #-----------------------     Generating Configuration File    -------------------------------------
    ConfigurationFile   = 'Sl_Config_v1.txt'
    
    #-----------------------Generating input spectra Textfile---------------------------------
    Interpolation       = interp1d(X, Y, kind = 'slinear')
    Wmin                = int(round(X[0],0))
    Wmax                = int(round(X[-1],0))
    
#   Interpolate the new spectra to one angstrom per pixel resolution
    X_1Angs = range(Wmin+1,Wmax-1,1)
    Y_1Angs = Interpolation(X_1Angs)
    
    Sl_Input_Filename = FileName.replace(".fits", ".slInput")
    SaveSpectra_2_Starlight([X_1Angs, Y_1Angs], Default_InputFolder + Sl_Input_Filename)
    print '-- Starlight File:', Default_InputFolder + Sl_Input_Filename
    #-----------------------     Generating Mask File    ----------------------------------------------
    
    if 'WHT' in FileName and TwoArmsMode == False:
        SpectraMeet             = Log_2_Parameter(FileFolder, CodeName, "Spectra_Meet")
        WHT_Wmatch_Wavelength   = float(Log_2_Parameter(FileFolder, CodeName, "WHT_Wmatch"))
        JoiningRegion_Begining  = WHT_Wmatch_Wavelength - 100
        JoiningRegion_End       = WHT_Wmatch_Wavelength + 100
        Not_MiddleRegionEncountered = False                #Boolean to check when we pass the match region... 
        
    
    else:
        Not_MiddleRegionEncountered = False                #Boolean to check when we pass the match region... 
        JoiningRegion_Begining = 0.0
        JoiningRegion_End = 0.0


    #Block Initial region
    Masks = []
    EdgeBoundary = 100
    Masks.append([Wmin, Wmin + EdgeBoundary, 'Lower_Edge'])
    
    #Import emision line location from lick indexes file
    LogName = CodeName + '_LickIndexes.txt'                             #Should change this to the Leak indexes plot  
    Labels_List = loadtxt(FileFolder + LogName, dtype=str, skiprows = 1, usecols = [0])
    IniWave_List, FinWave_List = loadtxt(FileFolder + LogName, dtype = float, skiprows = 1, usecols = (6,7), unpack = True)
    
    #In this first import we only load the masks within the spectra wavelength range
    for i in range(Labels_List.size):
        if (Wmin < IniWave_List[i])  and  (FinWave_List[i] < Wmax):
            Masks.append([IniWave_List[i],FinWave_List[i],Labels_List[i]])
    
    #Block final region
    Masks.append([Wmax - EdgeBoundary, Wmax , 'Upper_Edge'])
    MaskVector = [[],[],[]]
    
    #Generate the masks      
    for j in range(len(Masks)):
        Label               = Masks[j][2]
        Begining_point      = Masks[j][0]
        Ending_point        = Masks[j][1]
        

        if (Label == 'H1_6563A') or (Label == 'H1_4340A') or (Label == 'N2_6584A') or (Label == 'O3_5007A') or (Label == 'S3_9531A'):
            Increment = round((Masks[j][1] - Masks[j][0]) / 0.7, 0)
        elif (Label == 'Lower_Edge') or (Label == 'Upper_Edge'):  
            Increment = 0
        else:
            Increment = round((Masks[j][1] - Masks[j][0]) / 4, 0)
        
        IniWave = Masks[j][0] - Increment
        FinWave = Masks[j][1] + Increment
        
        if j > 0:
            PrevWave = MaskVector[1][-1]
            if IniWave <= PrevWave:
                if FinWave <= PrevWave: 
                    MaskVector[2][-1] = MaskVector[2][-1] + ' ' +Label
                else:    
                    MaskVector[1][-1] = FinWave
                    MaskVector[2][-1] = MaskVector[2][-1] + ' ' +Label
    
            else:
                MaskVector[0].append(IniWave)
                MaskVector[1].append(FinWave)
                MaskVector[2].append(Label)
        else:
            MaskVector[0].append(IniWave)
            MaskVector[1].append(FinWave)
            MaskVector[2].append(Label)
            
        Case_Inside     = ((MaskVector[0][-1] >= JoiningRegion_Begining) and (MaskVector[1][-1] <= JoiningRegion_End))
        Case_WeMissedIt = ((MaskVector[0][-1] >= JoiningRegion_Begining) and (MaskVector[1][-1] >= JoiningRegion_End) and (Not_MiddleRegionEncountered == True))
        
        if Case_Inside:
            if Not_MiddleRegionEncountered == True:
                MaskVector[0][-1] = JoiningRegion_Begining
                MaskVector[1][-1] = JoiningRegion_End
                MaskVector[2][-1] = 'Joining region'
                Not_MiddleRegionEncountered = False
            else:
                del MaskVector[0][-1]
                del MaskVector[1][-1]
                del MaskVector[2][-1]
        
        if Case_WeMissedIt:
            Ini_0   = MaskVector[0][-1]
            Fin_1   = MaskVector[1][-1]
            Lab     = MaskVector[2][-1]
            MaskVector[0][-1] = JoiningRegion_Begining
            MaskVector[1][-1] = JoiningRegion_End
            MaskVector[2][-1] = 'Joining region'
            MaskVector[0].append(Ini_0)
            MaskVector[1].append(Fin_1)
            MaskVector[2].append(Lab)            
            Not_MiddleRegionEncountered = False
    
    
    
    if TwoArmsMode == True:
        MaskFileName = CodeName + "_" + Color + '_Mask.lineslog'
    else:
        MaskFileName = CodeName + '_Mask.lineslog'
    
    #Esto como que jode el invento de antes....
    if SpectraMeet == 'True':
        MaskVector[0].append(JoiningRegion_Begining)
        MaskVector[1].append(JoiningRegion_End)
        MaskVector[2].append('Spectra meeting region')        
        
    
    File = open(FileFolder + MaskFileName, "w")
    File.write(str(len(MaskVector[0])) + '\n')
    for k in range(len(MaskVector[0])):
        Line = str(MaskVector[0][k]) + '  ' + str(MaskVector[1][k]) + '  0.0  ' + str(MaskVector[2][k]) + '\n'
        File.write(Line)
    
    File.close()

    copy(FileFolder + MaskFileName, Default_MaskFolder + MaskFileName)
    print '-- Mask File:', Default_MaskFolder + MaskFileName
    
    #-----------------------     Generating output files    -------------------------------------------

    if ".fits" in FileName:
        Sl_Output_Filename = FileName.replace(".fits", ".slOutput")
    else:
        Sl_Output_Filename = FileName.replace(FileName[FileName.rfind("."):len(FileName)],'.slOutput')

    Sl_OutputFolder = Default_OutputFoler
    print '-- Output address:', Sl_OutputFolder + Sl_Output_Filename

    #-----------------------Generating Grid file---------------------------------
     
    GridLines = []
    GridLines.append("1")                               #"[Number of fits to run]"])
    GridLines.append(Default_BasesFolder)               #"[base_dir]"])
    GridLines.append(Default_InputFolder)               #"[obs_dir]"])
    GridLines.append(Default_MaskFolder)                #"[mask_dir]"])
    GridLines.append(Default_OutputFoler)               #"[out_dir]"])
    GridLines.append("-652338184")                      #"[your phone number]"])
    GridLines.append("4500.0 ")                         #"[llow_SN]   lower-lambda of S/N window"])
    GridLines.append("4550.0")                          #"[lupp_SN]   upper-lambda of S/N window"])
    GridLines.append("3400.0")                          #"[Olsyn_fin] upper-lambda for fit"])
    GridLines.append("12000.0")                         #"[Olsyn_fin] upper-lambda for fit"])
    GridLines.append("1.0")                             #"[Odlsyn]    delta-lambda for fit"])
    GridLines.append("1.0")                             #"[fscale_chi2] fudge-factor for chi2"])
    GridLines.append(Velocity_Vector[0])                #"[FIT/FXK] Fit or Fix kinematics"])
    GridLines.append("0")                               #"[IsErrSpecAvailable]  1/0 = Yes/No"])
    GridLines.append("0")                               #"[IsFlagSpecAvailable] 1/0 = Yes/No"])

    Redlaw = 'CCM'
    v0_start = Velocity_Vector[1]
    vd_start = Velocity_Vector[2]

    GridLines.append([Sl_Input_Filename, ConfigurationFile, BaseDataFile, MaskFileName, Redlaw, v0_start, vd_start, Sl_Output_Filename])     
    Grid_FileName = FileName.replace(FileName[FileName.rfind("."):len(FileName)],'.slGrid')
            
    File = open(Starlight_Folder + Grid_FileName,"w")        
    
    print '-- Grid File:', Starlight_Folder + Grid_FileName
    
    for i in range(len(GridLines) - 1):
        Parameter = GridLines[i]
        Element = str(Parameter) + "\n"
        File.write(Element)
    
    Element = "  ".join(GridLines[-1])+'\n'
    File.write(Element)
    File.close()
        
    return Grid_FileName, Sl_Output_Filename, Sl_OutputFolder, X_1Angs, Y_1Angs

def Starlight_Launcher(Grid_FileName, ComputerRoot = '/home/vital/'):
    
    chdir(ComputerRoot + 'Starlight/')
    Command = './StarlightChains_v04.exe < ' + Grid_FileName
    print "Launch command:", Command
#     Command = './Starlight_v04_Mac.exe < grid_example1.in'  
    
    p = Popen(Command, shell=True, stdout=PIPE, stderr=STDOUT)
         
    for line in p.stdout.readlines():
        print line,
        
    retval = p.wait()
    
    return
    
def Sl_CumulativeLine_Plotter(Axis, index, x_j, Mini_j, Mcor_j, age_j, Z_j, LbyM, Mstars):

    Age = []
    PerMass = []
    PerMassPrev = []
    zValues = FamilyOfItemsInArray(Z_j)   
    legend = map(str, zValues)
    width = 0.07
    
    zPast = 0.0   
       
    for z in range(len(zValues)):              
        zTarget = zValues[z]
        
        # Get index of stars groups with same metallicity
        z_indez = []    
        for i in range(len(Z_j)):
            zFound = float(Z_j[i])
            if zFound == zTarget:
                z_indez.append(i)

        #Get mass fraction of stars meeting the metallicity criteria         
        PerMass = []
        for i in range(len(z_indez)):
            PerMass.append(x_j[z_indez[i]])
        
        #Get mass fraction of stars meeting the metallicity criteria
        Age = []
        for i in range(len(z_indez)):
            Age.append(log10(age_j[z_indez[i]]))
        
        #Get the mass which should go below each bar
        PerMassPrev = [] 
        for i in range(len(PerMass)):
            Bottom = 0.0
            for j in range(len(x_j)):
                AgePast = log10(round(age_j[j], 2-int(floor(log10(age_j[j])))-1))
                if (Z_j[j]==zPast) and (Age[i]==AgePast):
                    Bottom = x_j[j]
            PerMassPrev.append(Bottom)
   
        myColor = cm_hot(z/float(len(zValues)),1)
        plt.bar(Age,PerMass,width,bottom=PerMassPrev,color=myColor,align='center',label="z = "+legend[z])
        plt.legend(loc = 2)   
        zPast = zTarget

    return      
        
def SaveSpectra_2_Starlight(TableOfValues, FileAddress):

    File = open(FileAddress,"w")
    
    for i in range(len(TableOfValues[0])):
        Sentence = ''
        for j in range(len(TableOfValues)):
            Sentence = Sentence + ' ' + str(TableOfValues[j][i])
            
        File.write(Sentence + '\n')
    
    File.close()

    return       
        
def SlOutput_2_StarData(FileFolder, FileName):
    
    Sl_Data = File2Lines(FileFolder, FileName)
    
    BasesLine = LineFinder(Sl_Data, "[N_base]")                                                          #Location of my normalization flux in starlight output
    Bases = int(Sl_Data[BasesLine].split()[0])
      
    Sl_DataHeader = LineFinder(Sl_Data, "# j     x_j(%)      Mini_j(%)     Mcor_j(%)     age_j(yr)")        #Location of my normalization flux in starlight output
    Ind_i = Sl_DataHeader + 1
    Ind_f = Sl_DataHeader + Bases                                    
    
    index = []
    x_j = []
    Mini_j = []
    Mcor_j = []
    age_j = []
    Z_j = []
    LbyM =[]
    Mstars = []

    for j in range(Ind_i, Ind_f + 1): 
        myDataLine = Sl_Data[j].split()
        index.append(float(myDataLine[0]))
        x_j.append(float(myDataLine[1]))
        Mini_j.append(float(myDataLine[2]))
        Mcor_j.append(float(myDataLine[3]))
        age_j.append(float(myDataLine[4]))
        Z_j.append(float(myDataLine[5]))
        LbyM.append(float(myDataLine[6]))
        Mstars.append(float(myDataLine[7]))

    return index, x_j, Mini_j, Mcor_j, age_j, Z_j, LbyM, Mstars     

def Sl_CumulativeHistogram(Folder, File, Type, Color_Vector):
        
    OutPutData = []
    
    index, x_j, Mini_j, Mcor_j, age_j, Z_j, LbyM, Mstars = SlOutput_2_StarData(Folder, File)

    if Type == 'x_j':
        ParameterVector = x_j
        
    if Type == 'Mcor_j':
        ParameterVector = Mcor_j
        
    zValues = FamilyOfItemsInArray(Z_j)   
    
    AgeBins_HW = 0.20
    
    AgeBins = arange(4.80, 10, AgeBins_HW)
    DataVector = []
                    
    for i in range(len(zValues)):
        z               = zValues[i]
        MassSumVector   = zeros(len(AgeBins))
        LightMassVector = zeros(len(AgeBins))
        for j in range(len(Z_j)):
            if Z_j[j] == z:
                for k in range(len(AgeBins)):
                    Bin = AgeBins[k]
                    if (Bin - AgeBins_HW/2) <= log10(age_j[j]) < (Bin + AgeBins_HW/2):
                        MassSumVector[k] = MassSumVector[k] + ParameterVector[j]

        DataVector.append([z, Color_Vector[2][i], MassSumVector, LightMassVector])
        
    for i in range(len(AgeBins)):
        Total = 0  
        ListMasses = []
        for j in range(len(zValues)):
            ListMasses.append(DataVector[j][2][i])
            Total = Total +DataVector[j][2][i] 
                      
        Indexes_Small2Big = [n[0] for n in sorted(enumerate(ListMasses), key=lambda x:x[1])]
        Indexes_Small2Big.reverse()
 
        for Index in Indexes_Small2Big:
            if DataVector[Index][2][i] > 0:
#                         self.Histogram_One([AgeBins[i]- AgeBins_HW/4], [DataVector[Index][2][i]], str(DataVector[Index][0]), DataVector[Index][1], 1.0, Width = AgeBins_HW/2, Edgecolor = DataVector[Index][1], Linewidth = 1)

                X_Coord         = [AgeBins[i]- AgeBins_HW/4]
                Y_Coord         = [DataVector[Index][2][i]]
                Label           = str(DataVector[Index][0])
                Color           = DataVector[Index][1]
                Transparency    = 1.0
                Width           = AgeBins_HW/2
                Fill            = True
                Edgecolor       = DataVector[Index][1]
                Linestyle       = 'solid'
                Linewidth       = 1
                Log             = False
                
                OutPutData.append([X_Coord, Y_Coord, Label, Color, Transparency, Width, Fill, Edgecolor, Linestyle, Linewidth, Log])
                
                
        if Total > 0:

            X_Coord         = [AgeBins[i]- AgeBins_HW/2]
            Y_Coord         = [Total]
            Label           =  'Bin_total'
            Color           = Color_Vector[0]
            Transparency    = 1
            Width           = AgeBins_HW
            Fill            = False
            Edgecolor       = Color_Vector[1]
            Linestyle       = 'dotted'
            Linewidth       = 1
            Log             = False
            
            OutPutData.append([X_Coord, Y_Coord, Label, Color, Transparency, Width, Fill, Edgecolor, Linestyle, Linewidth, Log])
    
    return OutPutData