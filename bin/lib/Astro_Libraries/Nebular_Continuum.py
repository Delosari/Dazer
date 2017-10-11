'''
Created on Jul 24, 2015

@author: vital
'''
import pymc
from bisect         import bisect_left, bisect_right
from string         import split
from lmfit.models   import LinearModel
from numpy          import sum as np_sum,square,  mean, array, power, log10, linspace, exp, zeros, ceil, interp, isnan, searchsorted, median, argmin, concatenate, invert, argwhere, diff, r_, std, empty,\
    loadtxt, nonzero, flatnonzero, copy, where, arange, invert, logical_not, digitize
from lib.Math_Libraries.fitting_methods import Python_linfit
from scipy import interpolate

class NebularContinuumCalculator():

    def __init__(self):
        
        self.Calibration        = None
        
        self.Te                 = None
        self.Flux_Recomb        = None

        self.HII_HI             = 1
        self.HeII_HII           = None
        self.HeIII_HII          = None
        
        self.Wavelength_Range   = None
        
        self.DataRoot           = None
        
        self.Data_grid          = None

    def load_neb_constants(self, data_folder):
        
        #Dictionary with the calculation constants
        self.const = {}
         
        #Browns and Seaton FF methodology
        self.const['h']             = 6.626068e-27                            #erg s,  cm2 g / s  ==  erg s
        self.const['c_CGS']         = 2.99792458e10                           #cm / s
        self.const['c_Angs']        = 2.99792458e18                           #cm / s
        self.const['eV2erg']        = 1.602177e-12
        self.const['pi']            = 3.141592
        self.const['masseCGS']      = 9.1096e-28
        self.const['e_proton']      = 4.80320425e-10                          #statCoulomb = 1 erg^1/2 cm^1/2     # Eperez definition    electronCGS = 1.60217646e-19 * 3.e9 # Coulomb  # 1eV = 1.60217646e-19 Jul
        self.const['k']             = 1.3806503e-16                           #erg / K
        self.const['H0_ion_Energy'] = 13.6057                                 #eV
        self.const['nu_0']          = self.const['H0_ion_Energy'] * self.const['eV2erg'] / self.const['h'] #Hz
        self.const['H_Ryd_Energy']  = self.const['h'] * self.const['nu_0']

        #Coefficients for calculating A_2q The total radiative probability 2s -> 1s (s^-1)
        self.const['alpha_A']       = 0.88                                                                          
        self.const['beta_A']        = 1.53
        self.const['gamma_A']       = 0.8
        self.const['lambda_2q']     = 1215.7   # Angstroms
        self.const['C_A']           = 202.0    # (s^-1)
        self.const['A2q']           = 8.2249   # (s^-1) Transition probability at lambda = 1215.7                                                              

        #Free Bound constants
        self.const['Ryd2erg']       = 2.1798723e-11 # Rydberg to erg   # (s^-1) Transition probability at lambda = 1215.7                                                              
                
        #Load files
        self.HI_fb_dict             = self.read_Ercolano_FB_files(data_folder + 't3_elec.ascii')
        self.HeI_fb_dict            = self.read_Ercolano_FB_files(data_folder + 't5_elec.ascii')
        self.HeII_fb_dict           = self.read_Ercolano_FB_files(data_folder + 't4_elec.ascii')

        return

    def read_Ercolano_FB_files(self, file_address):
        
        dict_ion = {}
        
        #Reading the text files
        with open(file_address, 'r') as f:
            
            a = f.readlines()
            
            dict_ion['nTe']     = int(split(a[0])[0])               # number of Te columns
            dict_ion['nEner']   = int(split(a[0])[1])               # number of energy points rows
            dict_ion['skip']    = int(1+ceil(dict_ion['nTe']/8.))   # 8 es el numero de valores de Te por fila.
            dict_ion['temps']   = zeros(dict_ion['nTe'])
            dict_ion['lines']   = a
            
            #Storing temperature range            
            for i in range(1, dict_ion['skip']) :
                tt = split(a[i])
                for j in range(0,len(tt)) :
                    dict_ion['temps'][8*(i-1)+j] = tt[j]
        
            #Storing gamma_cross grids
            dict_ion['matrix'] = loadtxt(file_address, skiprows=dict_ion['skip'])
        
        #Get wavelengths corresponding to table threshold and remove zero entries
        wave_thres  = dict_ion['matrix'][:,0] * dict_ion['matrix'][:,1]        
        idx_zero = (wave_thres == 0)
        dict_ion['wave_thres'] = wave_thres[~idx_zero]
                
        return dict_ion

    def interp_Ercolano_grid(self,  wave_ryd, Ion_data_dict):

        #Bin the data into thresholds interval
        thres_idxbin        = digitize(wave_ryd, Ion_data_dict['wave_thres'])        
        interpol_grid       = interpolate.interp2d(Ion_data_dict['temps'], Ion_data_dict['matrix'][:,1], Ion_data_dict['matrix'][:,2:], kind='linear')
        lower_ener_thres    = Ion_data_dict['wave_thres'][thres_idxbin-1] #WARNING: This one could be an issue for the 
        #WARNING: This - 1 could be an issue for edge cases.. could it be done automatically-        
        
        return interpol_grid, lower_ener_thres
  
    def PropertiesfromUser(self, Te, Flux_Recomb, HeII_HII, HeIII_HII, Wavelength_Range = None, Calibration = 'Zanstra', Data_Root = '/home/vital/Dropbox/Astrophysics/Lore/NebularContinuum/'):

        self.Te             = Te
        self.Flux_Recomb    = Flux_Recomb

        self.HeII_HII       = HeII_HII
        self.HeIII_HII      = HeIII_HII

        self.Calibration    = Calibration
        
        if Wavelength_Range is None:
            self.Wavelength_Range = linspace(912, 30000, 20000-912, endpoint=True)
        else:
            self.Wavelength_Range = Wavelength_Range
        
        if Data_Root == None:
            self.DataRoot = '/home/vital/Workspace/Dazer/Astro_Libraries/'

        else:
            self.DataRoot = Data_Root
        
    def Calculate_Nebular_gamma(self):
        
        Frac            = 1 + self.HeII_HII*4 + self.HeIII_HII*4
        
        Gamma_2q        = self.TwoPhotonContinuum()
             
        Gamma_FF_HI     = self.FreeFreeContinuum("HI")
        Gamma_FF        = Frac * Gamma_FF_HI
            
        Gamma_FB_HI     = self.FreeBoundContinuum_EP("HI")
        Gamma_FB_HeI    = self.FreeBoundContinuum_EP("HeI")
        Gamma_FB_HeII   = self.FreeBoundContinuum_EP("HeII")
        Gamma_FB        = self.HII_HI * Gamma_FB_HI + self.HeII_HII * Gamma_FB_HeI + self.HeIII_HII * Gamma_FB_HeII
    
        Gamma_Total     = Gamma_FB + Gamma_FF + Gamma_2q
        
        c_AperS         = 2.99792458e18 # A / s
        Gamma_lambda    = Gamma_Total * (c_AperS / power(self.Wavelength_Range,2))
        
        return Gamma_Total, Gamma_lambda, Gamma_FB_HI, Gamma_FB_HeI, Gamma_FB_HeII, Gamma_2q, Gamma_FF

    def calculate_neb_gCont(self, wave, Te, HeII_HII, HeIII_HII):
        
        H_He_frac = 1 + HeII_HII * 4 + HeIII_HII * 4
        
        gamma_2q        = self.two_photon_gCont(wave, Te)
        
        gamma_ff        = H_He_frac * self.free_free_gCont(wave, Te, Z_ion = 1.0)
        
        #Get the wavelength range in ryddbergs for the Ercolano grids
        wave_ryd        = (self.const['h'] * self.const['c_Angs']) / (self.const['Ryd2erg'] * wave)
        gamma_fb_HI     = self.free_bound_gCont(wave_ryd, Te, self.HI_fb_dict)
        gamma_fb_HeI    = self.free_bound_gCont(wave_ryd, Te, self.HeI_fb_dict)
        gamma_fb_HeII   = self.free_bound_gCont(wave_ryd, Te, self.HeII_fb_dict)
        gamma_fb        = gamma_fb_HI + HeII_HII * gamma_fb_HeI + HeIII_HII * gamma_fb_HeII
        
        gamma_neb_nu    = gamma_2q + gamma_ff + gamma_fb
        #gamma_neb_lam   = gamma_neb_nu * self.const['nu2lamb_fact']
               
        return gamma_neb_nu
    
    def TwoPhotonContinuum(self):
        
        Gamma_2q = []
        
        #Coefficients for calculating A_2q The total radiative probability 2s -> 1s (s^-1)
        alpha_A = 0.88                                                                          
        beta_A = 1.53
        gamma_A = 0.8
        C_A = 202.0                                                     # (s^-1)
        c_CGS = 2.99792458e10                                           #cm / s  
        h = 6.626068e-27                                                #erg s
        H0_Ionization_Energy = 13.6                                     #eV
        eV2erg = 1.602177e-12
        nu_0 = H0_Ionization_Energy * eV2erg / h                        #Hz
        
        A2q = 8.2249                                                    #(s^-1) Transition probability at lambda = 1215.7                                                              
        alpha_eff_2q = 6.5346e-11 * self.Te**-0.72315                   #(cm^3 s^-1) Effective Recombination coefficient (Este es mi ajuste. El de Molla es: 0.647e-10 * Te**-0.722)
        q2 = 5.92e-4 - 6.1e-9 * self.Te                                 #(cm^3 s^-1) Collisional transition rate coefficient for protons and electrons
        Lambda_0 =  1215.7
        
        for i in range(len(self.Wavelength_Range)):
            Lambda = float(self.Wavelength_Range[i])
            if Lambda > 1215.7:
                nu = c_CGS / (Lambda * 1e-8)
                nu_0 = 3.e18/Lambda_0 
                y = nu / nu_0
                
                A_y = C_A * (
                             y*(1-y) * (1-(4*y*(1-y))**gamma_A)
                             + alpha_A*(y*(1-y))**beta_A * (4*y*(1-y))**gamma_A
                             )
                
                g_nu = h * nu / nu_0 / A2q * A_y                                                    # (erg Hz^-1)
                Gamma_2q.append(alpha_eff_2q * g_nu / (1+q2/A2q))                                # (erg cm^3 s^-1) We take the density of protons and electrons = n_e
    #             Gamma_2q.append(alpha_eff_2q * g_nu / (1+n_e*q2/A2q))                          # (erg cm^3 s^-1) We take the density of protons and electrons = n_e
            else:
                Gamma_2q.append(0.0)
    
        Np_Gamma_2q = array(Gamma_2q)    
        Np_Gamma_2q[isnan(Np_Gamma_2q)] = 0.0    
        
        return Np_Gamma_2q   
    
    def two_photon_gCont(self, wave, Te):
        
        #Prepare arrays
        idx_limit       = (wave > self.const['lambda_2q'])
        gamma_array     = zeros(wave.shape)

        #Params   
        alpha_eff_2q    = 6.5346e-11 * Te**-0.72315 #(cm^3 s^-1) Effective Recombination coefficient (Este es mi ajuste. El de Molla es: 0.647e-10 * Te**-0.722)
        q2              = 5.92e-4 - 6.1e-9 * Te     #(cm^3 s^-1) Collisional transition rate coefficient for protons and electrons
        
        nu_array        = self.const['c_Angs'] / wave[idx_limit]
        nu_limit        = self.const['c_Angs'] / self.const['lambda_2q']
        y               = nu_array / nu_limit
        
        A_y             = self.const['C_A'] * (y*(1-y) * (1-(4*y*(1-y))**self.const['gamma_A']) + self.const['alpha_A']*(y*(1-y))**self.const['beta_A'] * (4*y*(1-y))**self.const['gamma_A'])
        
        g_nu1           = self.const['h'] * nu_array / nu_limit / self.const['A2q'] * A_y 
        g_nu2           = alpha_eff_2q * g_nu1 / (1+q2/self.const['A2q'])
        
        gamma_array[idx_limit] = g_nu2[:]    
        #Np_Gamma_2q = array(Gamma_2q)    
        #Np_Gamma_2q[isnan(Np_Gamma_2q)] = 0.0    
                
        return gamma_array
    
    def free_free_gCont(self, wave, Te, Z_ion = 1.0):
                
        cte_A       = (32 * (Z_ion**2) * (self.const['e_proton']**4)*self.const['h']) / (3*(self.const['masseCGS']**2)*(self.const['c_CGS']**3))
        cte_B       = ((self.const['pi'] * self.const['H_Ryd_Energy'] / (3 * self.const['k'] * Te))**0.5)
        cte_Total   = cte_A * cte_B 
        
        nu_array    = self.const['c_Angs'] / wave
        
        gamma_Comp1 = exp(((-1 * self.const['h'] * nu_array) / (self.const['k'] * Te)))
        gamma_Comp2 = (self.const['h'] * nu_array) / ((Z_ion**2) * self.const['e_proton'] * 13.6057)
        gamma_Comp3 = self.const['k'] * Te / (self.const['h'] * nu_array)
               
        gff         = 1 + 0.1728*power(gamma_Comp2, 0.33333) * (1+2*gamma_Comp3) - 0.0496*power(gamma_Comp2, 0.66667) * (1 + 0.66667*gamma_Comp3 + 1.33333*power(gamma_Comp3, 2))
        
        gamma_array = cte_Total * gamma_Comp1 * gff
             
        return gamma_array

    def free_bound_gCont(self, wave_ryd, Te, data_dict):
        
        #Temperature entry
        t4      = Te / 10000.0
        logTe   = log10(Te)
                
        #Interpolating the grid for the right wavelength        
        thres_idxbin        = digitize(wave_ryd, data_dict['wave_thres'])        
        interpol_grid       = interpolate.interp2d(data_dict['temps'], data_dict['matrix'][:,1], data_dict['matrix'][:,2:], kind='linear')
        ener_low            = data_dict['wave_thres'][thres_idxbin-1] #WARNING: This one could be an issue for wavelength range table limits 
        
        #Interpolate table for the right temperature
        #gamma_cross_tableWave   = data_dict['interpol_grid'](logTe, data_dict['matrix'][:,1])[:, 0]
        gamma_inter_Te      = interpol_grid(logTe, data_dict['matrix'][:,1])[:, 0] 
        gamma_inter_Te_Ryd  = interp(wave_ryd, data_dict['matrix'][:,1], gamma_inter_Te)

        Gamma_fb_f          = gamma_inter_Te_Ryd * 1e-40 * t4**-1.5 * exp(-15.7887 * (wave_ryd - ener_low) / t4)
                    
        return Gamma_fb_f

    def gCont_calibration(self, wave, Te, flux_Emline, gNeb_cont_nu, lambda_EmLine = 6562.819):

        t4 = Te / 10000.0
        
        alfa_eff_alpha      = 2.708e-13 * t4**-0.648 / (1 + 1.315*t4**0.523)
        #alfa_eff_beta      = 0.668e-13 * t4**-0.507 / (1 + 1.221*t4**0.653)
 
        fNeb_cont_lambda    = gNeb_cont_nu * lambda_EmLine * flux_Emline / (alfa_eff_alpha * self.const['h'] * wave * wave)
        
        return fNeb_cont_lambda

    def FreeFreeContinuum(self, Ion):
        
        #Browns and Seaton methodology
        h = 6.626068e-27                                                #erg s
        c_CGS = 2.99792458e10                                           #cm / s   
        eV2erg = 1.602177e-12
        pi = 3.141592
        masseCGS = 9.1096e-28
        e_proton = 4.80320425e-10                                       #statCoulomb = 1 erg^1/2 cm^1/2     # Eperez definition    electronCGS = 1.60217646e-19 * 3.e9 # Coulomb  # 1eV = 1.60217646e-19 Jul
        k = 1.3806503e-16                                               #erg / K
        H0_Ionization_Energy = 13.6057                                  #eV
        nu_0 = H0_Ionization_Energy * eV2erg / h                        #Hz
        H_RydbergEnergy = h * nu_0
    
        Flux_List = []
        
        if Ion == "HI":
            Z = 1
         
        CteTotal = (32 * (Z**2) * (e_proton**4)*h) / (3*(masseCGS**2)*(c_CGS**3)) * ((pi * H_RydbergEnergy / (3 * k * self.Te))**0.5) 
                    
        for i in range(len(self.Wavelength_Range)):
            Wave = float(self.Wavelength_Range[i])
            nu = 2.99792458e18 / Wave  
            
            Flux_Comp1 = exp(((-1*h*nu)/(k*self.Te)))                                                                    #No unit    
            Flux_Comp2 = h*nu/(Z**2*e_proton*13.6057)
            Flux_Comp3 = k*self.Te/(h*nu)
            
            gff_mio = 1 + 0.1728*(Flux_Comp2)**0.33333*(1+2*Flux_Comp3) - 0.0496*(Flux_Comp2)**0.66667*(1+0.66667*Flux_Comp3+1.33333*(Flux_Comp3)**2)
            
            FluxValue = CteTotal * Flux_Comp1 * gff_mio                 # This is gamma_nu
            Flux_List.append(FluxValue)
         
        return array(Flux_List)

    def FreeBoundContinuum_EP(self, Ion):
        
        planckCGS = 6.626068e-27 # cm2 g / s  ==  erg s
        cspeedA = 3.e18 # A / s
        rydbergerg = 2.1798741e-11 # Rydberg to erg
        t4 = self.Te/1.e4
    
        PixelNumber = len(self.Wavelength_Range)
        
        nu          = zeros(PixelNumber)
        eryd        = zeros(PixelNumber)
        eryd_low    = zeros(PixelNumber)
        
        for i in range(len(self.Wavelength_Range)): 
            wave = self.Wavelength_Range[i]                             # wavelength en Angstroms
            nu[i] = cspeedA/wave                                        # frecuencia en Hz
            eryd[i] = planckCGS*cspeedA/(rydbergerg*wave)               # energia en Rydbergs
        
        if Ion == "HI":
            gammafile = self.DataRoot + "t3_elec.ascii"
        if Ion == "HeI":
            gammafile = self.DataRoot + "t5_elec.ascii"
        if Ion == "HeII":
            gammafile = self.DataRoot + "t4_elec.ascii"
    
        f = open(gammafile,'r')
        a = f.readlines()
        f.close()
        nTe = int(split(a[0])[0])                                       # number of Te columns
        nener = int(split(a[0])[1])                                     # number of energy points rows 
        skip = int(1+ceil(nTe/8.))                                      # 8 es el numero de valores de Te por fila.
        
        temp = zeros(nTe)
        for i in range(1,skip) :
            tt = split(a[i])
            for j in range(0,len(tt)) :
                temp[8*(i-1)+j] = tt[j]
    #    tindex = 2 + np.nonzero(temp==np.log10(Te))[0]
        
        rte2 = 0
        if bisect_right(temp,log10(self.Te)) - bisect_left(temp,log10(self.Te)) == 1 : 
            rte1 = bisect_left(temp,log10(self.Te))  # Te existe
        elif bisect_right(temp,log10(self.Te)) - bisect_left(temp,log10(self.Te)) == 0 :
            rte1 = bisect_left(temp,log10(self.Te))-1 # interpola Te
            rte2 = bisect_right(temp,log10(self.Te))  # 
        else :
            print 'ERROR IN Te COLUMN DETECTION FOR Gamma_fb'
    
        Gamma_fb_f = zeros(nener)
        ener = zeros(nener)
        gamma = zeros(nener)
        thresh = zeros(nener) 
        ener_low = zeros(nener) 
           
        for i in range(skip,skip+nener) :
            thresh[i-skip] = int(split(a[i])[0])         # indicador de threshold (1) o de punto nodal (0)
            ener[i-skip] = float(split(a[i])[1])         # energia en Ryd
            gamma[i-skip] = float(split(a[i])[2+rte1])   # [22] es para Te=1e4
            if rte2 > 0 : 
                gamma[i-skip] = float(split(a[i])[2+rte1]) + (float(split(a[i])[2+rte2])-float(split(a[i])[2+rte1]))/(10**temp[rte2]-10**temp[rte1]) * (self.Te-10**temp[rte1])
            
            ener_low = thresh*ener
            atmp=thresh.nonzero()[0]
            for i in range(1,len(atmp),2) : ener_low[atmp[i]] = ener_low[atmp[i-1]]
            for i in range(0,len(ener_low)) : 
                if ener_low[i] == 0 : ener_low[i] = ener_low[i-1]
            
            eryd_low = interp(eryd, ener, ener_low)
            Gamma_fb_f = interp(eryd, ener, gamma)  # linear interpolation at grid of interest
            Gamma_fb_f = Gamma_fb_f * 1e-40*t4**-1.5*exp(-15.7887*(eryd-eryd_low)/t4)
    
        return Gamma_fb_f      

    def Zanstra_Calibration(self, EmLine, Flux_EmLine, Gamma_nu):
        
        # Pequignot et al. 1991
        t4 = self.Te / 10000
        
        if EmLine == 'Hbeta':
            alfa_eff = 0.668e-13 * t4**-0.507 / (1 + 1.221*t4**0.653)
            lambda_EmLine = 4862.683
        elif EmLine == 'Halpha':
            alfa_eff= 2.708e-13 * t4**-0.648 / (1 + 1.315*t4**0.523)
            lambda_EmLine = 6562.819
    
        planckCGS = 6.626068e-27  #erg s-1 cm-2
                     
        NebularFlux = Gamma_nu * lambda_EmLine * Flux_EmLine / (alfa_eff * planckCGS * self.Wavelength_Range**2)                        # erg s-1 [cm-2] A-1 
        
        return NebularFlux 

    def Estimate_BalmerJumpFlux(self, Int):
        
        #This method is use for plotthing a trendline from the green region to the location of the balmer jump
        
        RegionsList = [[4500,4600],
                       [4600,4675],
                       [4750,4800],
                       [4800,4850],
                       [4875,4940],
                       [5040,5100],
                       [5200,5270],
                       [5300,5370],
                       [5600,5670]]
        
        Wave_Regr   = zeros(len(RegionsList))
        Int_Regre   = zeros(len(RegionsList))
            
        for i in range(len(RegionsList)):
            
            Region          = RegionsList[i]
            
            indmin, indmax  = searchsorted(self.Wavelength_Range, (Region[0], Region[1]))
            indmax          = min(len(self.Wavelength_Range)-1, indmax)
            
            Median          = median(Int[indmin:indmax])
            Loc             = (Region[0] + Region[-1]) / 2
            
            Wave_Regr[i]    = Loc
            Int_Regre[i]    = Median
        
        
        
        slope, intercept    = Python_linfit(Wave_Regr, Int_Regre, y_err=None, errors_output = False)
        
        x_trendline         = array(range(3500,5000,1))
        y_trendline                   = slope * x_trendline + intercept
        BalmerJump_Int      = slope * 3646 + intercept
        
        if self.Wavelength_Range[0] < 3646:
            
            Index_blue, Index_red = searchsorted(self.Wavelength_Range, 3646 - 25),  searchsorted(self.Wavelength_Range, 3646 + 25)
            BalmerJump_Int = median(Int[Index_blue:Index_red])
    
        return x_trendline, y_trendline, BalmerJump_Int
    
    def Extend_WavelengthRange(self, Wave):
    
        if Wave[0] >= 3500:
        
            FirstPoint = int(round(Wave[0],0)-1)
            Interval = array(range(3500,FirstPoint,1))
            
            Wave_Neb = concatenate((Interval, Wave))
                        
        else:
            Wave_Neb    = Wave
            
        return Wave_Neb
