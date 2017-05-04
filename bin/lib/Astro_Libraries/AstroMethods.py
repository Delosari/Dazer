'''
Created on Jan 29, 2014

@author: vital
'''
import pyfits

import CodeTools.vitools as vit
import numpy as np
import scipy as sp


def Fits2Data(FolderName,FileName):
    
    #Get the fits file main data and header
    Data, Header_0 = pyfits.getdata(FolderName + FileName, 1, header=True)

    if ('WHTJOINW' in Header_0) or ("STALIGHT" in Header_0) or ("NEBUSPEC" in Header_0):
        x = Data['Wave']
        y = Data['Int']

        return y, x, [Header_0]
    
    #Procedure for dr10 spectra 
    #It returns the spectra redshift corrected and in absolute units
    elif ("COEFF0" in Header_0 and "dr10" in FileName):
        FitsFile = pyfits.open(FolderName + FileName) 
        Spectra = FitsFile[1].data
        Header_2  = FitsFile[2].data
        Header_3 = FitsFile[3].data
        
        Int = Spectra['flux']
        Int_abs = Int / 1e17
        
        WavelenRange = 10.0**Spectra['loglam']
        SDSS_z = float(Header_2["z"][0] + 1)
        Wavelength_z = WavelenRange / SDSS_z

        Headers = (Header_0,Header_2,Header_3)

        return Int_abs, Wavelength_z, Headers
        
    else:
        if Data.ndim == 1: 
            Int = Data
        else:
            Int = Data[0]
            
        if "COEFF0" in Header_0:
            dw = 10.0**Header_0['COEFF1']               # dw = 0.862936 INDEF (Wavelength interval per pixel)
            Wmin = 10.0**Header_0['COEFF0']
            pixels = Header_0['NAXIS1']                 # nw = 3801 number of output pixels
            Wmax     = Wmin + dw * pixels
            WavelenRange = sp.linspace(Wmin,Wmax,pixels,endpoint=False)
            return Int, WavelenRange, [Header_0]
        
        elif "LTV1" in Header_0:
            StartingPix     = -1 * Header_0['LTV1']                   # LTV1 = -261. 
            Wmin_CCD        = Header_0['CRVAL1']
            dw              = Header_0['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
            pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
            Wmin            = Wmin_CCD + dw * StartingPix
            Wmax            = Wmin + dw * pixels
            WavelenRange    = np.linspace(Wmin,Wmax,pixels,endpoint=False)
            
            return Int, WavelenRange, [Header_0]
        
        else:
            Wmin            = Header_0['CRVAL1']
            dw              = Header_0['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
            pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
            Wmax            = Wmin + dw * pixels
            WavelenRange    = np.linspace(Wmin,Wmax,pixels,endpoint=False)
                            
            return WavelenRange, Int, [Header_0]


def Fits2Data_Old(FolderName,FileName):
    FitsFile = pyfits.open(FolderName + FileName) 
    TotalInt = FitsFile[0].data
    if TotalInt.ndim == 1: 
        Int = TotalInt
        MaxInten = max(Int)
    else:
        Int = TotalInt[0,0,:]
        MaxInten = max(Int)
    N_Inten = Int / MaxInten

    if "COEFF0" in FitsFile[0].header:
        dw = 10.0**FitsFile[0].header['COEFF1']               # dw = 0.862936 INDEF (Wavelength interval per pixel)
        Wmin = 10.0**FitsFile[0].header['COEFF0'] + dw
        pixels = FitsFile[0].header['NAXIS1']                 # nw = 3801 number of output pixels
        Wmax = Wmin + dw * pixels
        WavelenRange = sp.linspace(Wmin,Wmax,pixels,endpoint=False)
        FitsFile.close()
        return Int, N_Inten, WavelenRange, Wmin, Wmax 
    
    else:
        if "LTV1" in FitsFile[0].header:
            StartingPix = -1 * FitsFile[0].header['LTV1']     # LTV1 = -261. 
            Wmin_CCD = FitsFile[0].header['CRVAL1']
            dw = FitsFile[0].header['CD1_1']                  # dw = 0.862936 INDEF (Wavelength interval per pixel)
            pixels = FitsFile[0].header['NAXIS1']             # nw = 3801 number of output pixels
            Wmin = Wmin_CCD + dw * StartingPix
            Wmax = Wmin + dw * pixels
            WavelenRange = sp.linspace(Wmin,Wmax,pixels,endpoint=False)
            FitsFile.close()
            return Int, N_Inten, WavelenRange, Wmin, Wmax
        else:
            Wmin = FitsFile[0].header['CRVAL1']
            dw = FitsFile[0].header['CD1_1']                  # dw = 0.862936 INDEF (Wavelength interval per pixel)
            pixels = FitsFile[0].header['NAXIS1']             # nw = 3801 number of output pixels
            Wmax = Wmin + dw * pixels
            WavelenRange = sp.linspace(Wmin,Wmax,pixels,endpoint=False)
            FitsFile.close()
            return Int, N_Inten, WavelenRange, Wmin, Wmax 
       
def StarlightFileManager(FileFolder, FileName):
    DataFile = open(FileFolder + FileName,"r")
    
    StarlightOutput = DataFile.readlines()
    
    DataFile.close()
        
    ## Synthesis Results - Best model ##
    Chi2Line = vit.LineFinder(StarlightOutput, "[chi2/Nl_eff]")
    AdevLine = vit.LineFinder(StarlightOutput, "[adev (%)]")
    SumXdevLine = vit.LineFinder(StarlightOutput, "[sum-of-x (%)]")
    v0_min_Line = vit.LineFinder(StarlightOutput, "[v0_min  (km/s)]")
    vd_min_Line = vit.LineFinder(StarlightOutput, "[vd_min  (km/s)]")
    Av_min_Line = vit.LineFinder(StarlightOutput, "[AV_min  (mag)]")

    Nl_eff_line = vit.LineFinder(StarlightOutput, "[Nl_eff]")
    
    
    SignalToNoise_Line = vit.LineFinder(StarlightOutput, "## S/N")
        
    l_norm_Line = vit.LineFinder(StarlightOutput, "## Normalization info") + 1 
    llow_norm_Line = vit.LineFinder(StarlightOutput, "## Normalization info") + 2 
    lupp_norm_Line = vit.LineFinder(StarlightOutput, "## Normalization info") + 3     
    NormFlux_Line = vit.LineFinder(StarlightOutput, "## Normalization info") + 4 
    
    SpecLine = vit.LineFinder(StarlightOutput, "## Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei")     #Location of my Spectrum in starlight output
   
    #Quality of fit                                         
    Chi2 =  float(StarlightOutput[Chi2Line].split()[0])                                        
    Adev = float(StarlightOutput[AdevLine].split()[0])
    SumXdev = float(StarlightOutput[SumXdevLine].split()[0])
    Nl_eff = float(StarlightOutput[Nl_eff_line].split()[0])
    v0_min = float(StarlightOutput[v0_min_Line].split()[0])
    vd_min = float(StarlightOutput[vd_min_Line].split()[0])
    Av_min = float(StarlightOutput[Av_min_Line].split()[0])
 
    #Signal to noise configuration                                           
    SignalToNoise_lowWave = float(StarlightOutput[SignalToNoise_Line + 1].split()[0]) 
    SignalToNoise_upWave = float(StarlightOutput[SignalToNoise_Line + 2].split()[0]) 
    SignalToNoise_magnitudeWave = float(StarlightOutput[SignalToNoise_Line + 3].split()[0]) 
    
    #Flux normailzation parameters                                          
    l_norm = float(StarlightOutput[l_norm_Line].split()[0])
    llow_norm = float(StarlightOutput[llow_norm_Line].split()[0])
    lupp_norm = float(StarlightOutput[lupp_norm_Line].split()[0])
    FluxNorm = float(StarlightOutput[NormFlux_Line].split()[0])

    Parameters = (Chi2, Adev, SumXdev, Nl_eff, v0_min, vd_min, Av_min, SignalToNoise_lowWave, SignalToNoise_upWave, SignalToNoise_magnitudeWave, l_norm, llow_norm, lupp_norm)
    
    #Spectra pixels location
    Pixels_Number = int(StarlightOutput[SpecLine+1].split()[0])                     #Number of pixels in the spectra
    Ind_i = SpecLine+2                                                              #First pixel location      
    Ind_f = Ind_i + Pixels_Number                                                   #Final pixel location
    
    Input_Wavelength = np.zeros(Pixels_Number)
    Input_Flux = np.zeros(Pixels_Number)
    Output_Flux = np.zeros(Pixels_Number)
    Output_Mask = np.zeros(Pixels_Number)
    
    for i in range(Ind_i, Ind_f): 
        Index = i - Ind_i
        Line = StarlightOutput[i].split()
        Input_Wavelength[Index] = float(Line[0])
        Input_Flux[Index] = float(Line[1])*FluxNorm
        Output_Flux[Index] = float(Line[2])*FluxNorm
        Output_Mask[Index] = float(Line[3])
    
    MaskPixels = [[],[]]                     #The 0 tag
    ClippedPixels = [[],[]]                  #The -1 tag
    FlagPixels = [[],[]]                     #The -2 tag
    
    for j in range(len(Output_Mask)):
        PixelTag = Output_Mask[j]
        Wave = Input_Wavelength[j]
        if PixelTag == 0:
            MaskPixels[0].append(Wave)
            MaskPixels[1].append(Input_Flux[j])
        if PixelTag == -1:
            ClippedPixels[0].append(Wave)
            ClippedPixels[1].append(Input_Flux[j])            
        if PixelTag == -2:
            FlagPixels[0].append(Wave)
            FlagPixels[1].append(Input_Flux[j])            
  
    return Input_Wavelength, Input_Flux, Output_Flux, MaskPixels, ClippedPixels, FlagPixels, Parameters
    
def Fnu_to_FmAB(Flux_nu):
    
    Flux_mAB = -2.5 * np.log10(Flux_nu) - 48.6
    
    return Flux_mAB

def Flam_to_Fnu(Flux_lam, Wavelength):
    
    c_AperS = 2.99792458e18
    
    Flux_nu = (np.power(Wavelength,2)/c_AperS) * Flux_lam
    
    return Flux_nu

def Fnu_to_Flam(Flux_nu, Wavelength):

    c_AperS = 2.99792458e18
    
    Flux_lam = (c_AperS / np.power(Wavelength, 2)) * Flux_nu

    return Flux_lam

def mAb_to_Fnu(Flux_mAB):

    fnu_zero = 3.68e-20

    Flux_nu = fnu_zero * np.power(10, -0.4 * (Flux_mAB))
    
    return Flux_nu

def Gauss_Function(x,*Funct_Coeff):
    A, Mu, Sigma = Funct_Coeff
    return A*np.exp(-(x-Mu)**2/(2.*Sigma**2))

# def TheoHBetaCoefficients():
#     
#     #Theoretical coefficients 
#     TextLines = vit.File2Lines("/home/vital/Dropbox/Astrophysics/Data/","cHBeta_TheoCoefficients.csv")
#     
#     TheoCHBetaTable = []
#     
#     for i in range(len(TextLines)):
#         row = TextLines[i].split()
#         TheoCHBetaTable.append(row)
#         
#     return TheoCHBetaTable

def Reddening_CCM98(Flux,Wave,EBV,R_V,):
    
    x = 1 / (Wave / 10000)
    y = x - 1.82

    y_coeffs = np.array([np.ones(len(y)), y, np.power(y, 2), np.power(y, 3), np.power(y, 4), np.power(y, 5), np.power(y, 6), np.power(y, 7)])
    a_coeffs = np.array([1, 0.17699,    -0.50447,   -0.02427,   0.72085,    0.01979,    -0.77530,   0.32999])
    b_coeffs = np.array([0, 1.41338,    2.28305,   1.07233,   -5.38434,    -0.62251,    5.30260,   -2.09002])
    
    a_x =  np.dot(a_coeffs,y_coeffs)
    b_x =  np.dot(b_coeffs,y_coeffs)
    
    A_V = EBV * R_V
    
    A_lambda = (a_x + b_x / R_V) * A_V
    
    Flux_New = Flux * np.power(10, 0.4 * A_lambda)
    
#     Fig2=plt.figure(figsize=(16,10))
#     Axes = Fig2.add_subplot(111)
#     Axes.plot(Wave,Flux_New)
#     Axes.plot(Wave,Flux)
#     plt.show()
    
    return Flux_New

def Reddening_CCM98_err(Flux,Wave,EBV,R_V,FolderName,FileLog):
    
    x = 1 / (Wave / 10000)
    y = x - 1.82

    y_coeffs = np.array([np.ones(len(y)), y, np.power(y, 2), np.power(y, 3), np.power(y, 4), np.power(y, 5), np.power(y, 6), np.power(y, 7)])
    a_coeffs = np.array([1, 0.17699,    -0.50447,   -0.02427,   0.72085,    0.01979,    -0.77530,   0.32999])
    b_coeffs = np.array([0, 1.41338,    2.28305,   1.07233,   -5.38434,    -0.62251,    5.30260,   -2.09002])
    
    a_x =  np.dot(a_coeffs,y_coeffs)
    b_x =  np.dot(b_coeffs,y_coeffs)
    
    A_V = EBV * R_V
    
    A_lambda = (a_x + b_x / R_V) * A_V
    
    Flux_New = Flux * np.power(10, 0.4 * A_lambda)
    
#     Fig2=plt.figure(figsize=(16,10))
#     Axes = Fig2.add_subplot(111)
#     Axes.plot(Wave,Flux_New)
#     Axes.plot(Wave,Flux)
#     plt.show()
    
    return Flux_New
    
def Reddening_SingleLine_CCM98(Flux,Wave,EBV,R_V,):
    
    x = 1 / (Wave / 10000)
    y = x - 1.82

    y_coeffs = np.array([np.ones(len(y)), y, np.power(y, 2), np.power(y, 3), np.power(y, 4), np.power(y, 5), np.power(y, 6), np.power(y, 7)])
    a_coeffs = np.array([1, 0.17699,    -0.50447,   -0.02427,   0.72085,    0.01979,    -0.77530,   0.32999])
    b_coeffs = np.array([0, 1.41338,    2.28305,   1.07233,   -5.38434,    -0.62251,    5.30260,   -2.09002])
    
    a_x =  np.dot(a_coeffs,y_coeffs)
    b_x =  np.dot(b_coeffs,y_coeffs)
    
    A_V = EBV * R_V
    
    A_lambda = (a_x + b_x / R_V) * A_V
    
    Flux_New = Flux * np.power(10, 0.4 * A_lambda)
        
    return Flux_New
    
def Reddening_f_Cal(wave,R_v):
    
    y=(1/(wave/10000))
    y_beta = (1/(4861.333/10000))
    dm_beta = 0.74*y_beta-0.34+0.341*R_v-1.014
    
    if y <= 2.29:
        dm_lam = 0.74*y-0.34+0.341*R_v-1.014
         
    else:
        dm_lam = 0.43*y+0.37+0.341*R_v-1.014
    
    f = dm_lam/dm_beta
    
    return f

def Log_2_Parameter(FileFolder, ObjectCode, ParameterToFind):
    
    LogExtension = "_log.txt"
     
    Log_File_Address = FileFolder + ObjectCode + LogExtension
    
    Obj_Log = vit.File2Table(Log_File_Address, "")
    
    Index = "None"
    
    for i in range(len(Obj_Log)):
        Item = Obj_Log[i][0]
        if ParameterToFind == Item:
            Index = i
    
    if Index == "None":
        print "WARNING: The parameter cannot be found in log"
        print ParameterToFind
        print Obj_Log
        print FileFolder + ObjectCode + LogExtension
        
    Magnitude = Obj_Log[Index][1]
                
    return Magnitude

def PyNeb_ValidDiagnostics(Diags,Obs,MisDiagnosticos):
    
    ErrorCount = 0
    Data = []
    
    print "--------------------------"
    
    for i in range(len(MisDiagnosticos.keys())):
        for j in range(len(MisDiagnosticos.keys())):
                if i != j:
                    try:
                        Te, Ne = Diags.getCrossTemDen(MisDiagnosticos.keys()[i] ,MisDiagnosticos.keys()[j], obs=Obs)
                        if str(Te) != "nan" and str(Ne) != "nan":
                            print "Combining: "
                            print str(MisDiagnosticos.keys()[i])
                            print str(MisDiagnosticos.keys()[j])
                            print Te, Ne
                            #print type(Te)
                            Data.append([MisDiagnosticos.keys()[i],MisDiagnosticos.keys()[j],Te,Ne])
                    except:
                        print ErrorCount + 1
    
    print "--------------------------"
                    
    return Data

def PyNeb_EmLineNormalizer(EmList,ValuesList):
    
    HBeta_Value = "None"
    HBetaLabel = "H1_4861A"
    
    for i in range(len(EmList)):
        if EmList[i] == HBetaLabel:
            HBeta_Value = float(ValuesList[i])
    
    if HBeta_Value == "None":
        print "WARNING: HBeta value not found"
    
    print "The HBeta value is: "
    print HBeta_Value
    
    n_ValuesList = []
    n_ValuesList.append(ValuesList[0])
    
    for i in range(1,len(ValuesList)):
        Flux = float(ValuesList[i])
        n_Flux = Flux / HBeta_Value          
        n_ValuesList.append(str(n_Flux))
        
    return n_ValuesList, str(HBeta_Value)

def EP_Constant(a_b_List,Sigma_x_y_List):
    Sigma_Y = np.sqrt(np.sum(np.power(a_b_List*Sigma_x_y_List,2)))
    return Sigma_Y

def EP_Sum(Sigma_x_y_List):
    Sigma_Y = np.sqrt(np.sum(np.power(Sigma_x_y_List,2)))
    return Sigma_Y

def EP_MulDiv(z,x_y_List,Sigma_x_y_List):
    Sigma_Z = z * np.sqrt(np.power(Sigma_x_y_List/x_y_List,2))
    return Sigma_Z

def EP_Powers(z,x_y_List,Sigma_x_y_List,a_b_List):
    Sigma_Z = z * np.sqrt(np.power(a_b_List*Sigma_x_y_List/x_y_List,2))
    return Sigma_Z

def EP_log10(x,Sigma_X,m):
    Sigma_Z = m * Sigma_X / (2.303 * x)
    return Sigma_Z

def HMS2deg(ra='', dec=''):
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    
    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)
    
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC 

def Te_ne_Diagnostics():
    Temperature_Diags = ["[ArIII] 5192/7136","[ArIII] 5192/7300+","[NII] 5755/6584","[NII] 5755/6584+" ,"[OIII] 4363/5007","[OIII] 4363/5007+","[SIII] 6312/9069","[SIII] 6312/9200+"]
    Density_Diags = ["[ArIV] 4740/4711","[ClIII] 5538/5518","[FeIII] 4987/4659","[FeIII] 4987/4703","[FeIII] 4987/4882","[OII] 3726/3729","[SII] 6731/6716"]
        
    return Temperature_Diags, Density_Diags

def Bilinear_Interpolation(X,Y,Q_11,Q_12,Q_21,Q_22):
    
    if len(Q_11) != 3 or len(Q_12) != 3 or len(Q_21) != 3 or len(Q_22) != 3:
        print "WARNING: Elements do not have the right size"

    x_1 = Q_11[1]
    x_2 = Q_22[1]

    y_1 = Q_11[2]
    y_2 = Q_22[2]
        
    Q = 1 / ((x_2 - x_1) * (y_2 - y_1)) * (
                                          Q_11[0] * (x_2 - X) * (y_2 - Y)
                                        + Q_21[0] * (X - x_1) * (y_2 - Y)
                                        + Q_12[0] * (x_2 - X) * (Y - y_1)
                                        + Q_22[0] * (X - x_1) * (Y - y_1)  
                                          )
             
    return Q




# def AB_to_Erg(Flujo):
#     
#     a = np.array([1,2,3])
#     
#     if type(Flujo) == "list" or type(Flujo) == type(a):
#         
#         Erg_list = []
#         
#         for i in range(len(Flujo)):
#             AB = Flujo[i]
#             
#             if type(AB) != "float":
#                 AB = float(AB)
#                         
#             Erg = math.pow(10,-0.4*(AB+48.6))
#             
#             Erg_list.append(Erg)
#         
#         return Erg_list
#         
#     if type(Flujo) == float:
#         
#         AB = Flujo
#         
#         Erg = math.pow(10,-0.4*(AB+48.6))
#         
#         return Erg
#       
# def Erg_to_AB(Flujo):
#     
#     a = np.array([1,2,3])
#     
#     if type(Flujo) == "list" or type(Flujo) == type(a):
#         print "Echo"
#         
#         AB_list = []
#         
#         for i in range(len(Flujo)):
#             Erg = Flujo[i]
#             
#             if type(Erg) != "float":
#                 Erg = float(Erg)
#             
#             AB = -2.5 * math.log10(Erg) - 48.6
#             
#             AB_list.append(Erg)
#         
#         return AB_list
#         
#     if type(Flujo) == float:
#         
#         Erg = Flujo
#         
#         AB = - 2.5 * math.log10(Erg) - 48.6
#         
#         return AB


# def AB_to_Erg(Flujo):
#     
#     a = np.array([1,2,3])
#     
#     if type(Flujo) == "list" or type(Flujo) == type(a):
#         
#         Erg_list = []
#         
#         for i in range(len(Flujo)):
#             AB = Flujo[i]
#             
#             if type(AB) != "float":
#                 AB = float(AB)
#                         
#             Erg = math.pow(10,-0.4*(AB+48.6))
#             
#             Erg_list.append(Erg)
#         
#         return Erg_list
#         
#     if type(Flujo) == float:
#         
#         AB = Flujo
#         
#         Erg = math.pow(10,-0.4*(AB+48.6))
#         
#         return Erg
#       
# def Erg_to_AB(Flujo):
#     
#     a = np.array([1,2,3])
#     
#     if type(Flujo) == "list" or type(Flujo) == type(a):
#         print "Echo"
#         
#         AB_list = []
#         
#         for i in range(len(Flujo)):
#             Erg = Flujo[i]
#             
#             if type(Erg) != "float":
#                 Erg = float(Erg)
#             
#             AB = -2.5 * math.log10(Erg) - 48.6
#             
#             AB_list.append(Erg)
#         
#         return AB_list
#         
#     if type(Flujo) == float:
#         
#         Erg = Flujo
#         
#         AB = - 2.5 * math.log10(Erg) - 48.6
#         
#         return AB
