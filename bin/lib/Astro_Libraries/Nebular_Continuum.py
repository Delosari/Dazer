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
from bin.lib.Math_Libraries.fitting_methods import Python_linfit
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

    def load_neb_constants(self, wave, data_folder = '/home/vital/Dropbox/Astrophysics/Lore/NebularContinuum/'):
        
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
        self.const['wave_Ryd']      = (self.const['h'] * self.const['c_Angs']) / (self.const['Ryd2erg'] * wave)   # Wavelength range in Rydbergs
                
        #Load files
        self.HI_fb_dict             = self.read_Ercolano_FB_files(data_folder + 't3_elec.ascii', wave)
        self.HeI_fb_dict            = self.read_Ercolano_FB_files(data_folder + 't5_elec.ascii', wave)
        self.HeII_fb_dict           = self.read_Ercolano_FB_files(data_folder + 't4_elec.ascii', wave)

        return

    def read_Ercolano_FB_files(self, file_address, wave):
        
        dict_ion = {}
        
        #Reading the text files
        with open(file_address, 'r') as f:
            
            a = f.readlines()
            
            dict_ion['nTe']     = int(split(a[0])[0])              # number of Te columns
            dict_ion['nEner']   = int(split(a[0])[1])             # number of energy points rows
            dict_ion['skip']    = int(1+ceil(dict_ion['nTe']/8.)) # 8 es el numero de valores de Te por fila.
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
        wave_thres = dict_ion['matrix'][:,0] * dict_ion['matrix'][:,1]        
        idx_zero = (wave_thres == 0)
        dict_ion['wave_thres'] = wave_thres[~idx_zero]
        
        #Bin the data into thresholds interval
        dict_ion['thres_idxbin']        = digitize(self.const['wave_Ryd'], dict_ion['wave_thres'])        
        dict_ion['interpol_grid']       = interpolate.interp2d(dict_ion['temps'], dict_ion['matrix'][:,1], dict_ion['matrix'][:,2:], kind='linear')
        dict_ion['lower_ener_thres']    = dict_ion['wave_thres'][dict_ion['thres_idxbin']-1] #WARNING: This one could be an issue for the 
        #WARNING: This - 1 could be an issue for edge cases.. could it be done automatically-
        
        return dict_ion
     
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
        
        c_AperS         = 2.99792458e18                                             # A / s
        Gamma_lambda    = Gamma_Total * (c_AperS / power(self.Wavelength_Range,2))
        
        return Gamma_Total, Gamma_lambda, Gamma_FB_HI, Gamma_FB_HeI, Gamma_FB_HeII, Gamma_2q, Gamma_FF
    
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
            
#         Np_Gamma_2q = array(Gamma_2q)    
#         Np_Gamma_2q[isnan(Np_Gamma_2q)] = 0.0    
                
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

    def free_bound_gCont(self, wave, Te, data_dict):
               
        t4                  = Te / 10000.0
        logTe               = log10(Te)
        eryd                = self.const['wave_Ryd']   # energy in Rydbergs
        ener_low            = data_dict['lower_ener_thres']
                                  
        #Interpolate table for the right temperature
        #gamma_cross_tableWave   = data_dict['interpol_grid'](logTe, data_dict['matrix'][:,1])[:, 0]
        gamma_inter_Te      = data_dict['interpol_grid'](logTe, data_dict['matrix'][:,1])[:, 0] 
        gamma_inter_Te_Ryd  = interp(self.const['wave_Ryd'], data_dict['matrix'][:,1], gamma_inter_Te)

        Gamma_fb_f          = gamma_inter_Te_Ryd * 1e-40 * t4**-1.5 * exp(-15.7887 * (eryd - ener_low) / t4)
                            
        return Gamma_fb_f

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
    
class Nebular_Bayesian():

    def __init__(self):
        
        self.Calibration        = None
        
        self.c_AperS            = 2.99792458e18 # A / s
        self.planckCGS          = 6.626068e-27  #erg s-1 cm-2
        
        self.DataRoot   = '/home/vital/git/Dazer/Dazer/dazer/lib/Astro_Libraries/'
        self.lineal_mod = LinearModel(prefix='lineal_')

        gammafile = self.DataRoot + "HI_t3_elec.ascii"
        f = open(gammafile,'r')
        self.a_HI = f.readlines()
        f.close()
        
        gammafile = self.DataRoot + "HeI_t5_elec.ascii"
        f = open(gammafile,'r')
        self.a_HeI = f.readlines()
        f.close()
        
        gammafile = self.DataRoot + "HeII_t4_elec.ascii"
        f = open(gammafile,'r')
        self.a_HeII = f.readlines()
        f.close()
    
    def calculate_wavelength_intervals(self, total_wavelength, total_flux, lineslog_frame):

        #Another value
        empty_array = zeros(len(total_wavelength), dtype = bool)
        for j in range(len(lineslog_frame.index)):
    
            blue_limit  = lineslog_frame.iloc[j]['Wave3']  
            red_limit   = lineslog_frame.iloc[j]['Wave4']
            
            if lineslog_frame.index.values[j] in ['O2_3726A', 'O3_4959A', 'O3_5007A', 'H1_6563A']:
                tolerance = 15
            else:
                tolerance = 4
                
            indeces_line    = (total_wavelength >= blue_limit - tolerance) & (total_wavelength <= red_limit + tolerance)   
            empty_array     = empty_array + indeces_line
                            
        #Extra for weird lines
        for interval in [[3630,3640],[3655, 3677], [3816, 3824], [3537, 3542], [6990,7010], [5238, 5330]]:
    
            blue_limit  = interval[0] 
            red_limit   = interval[1]
            tolerance   = 0

            indeces_line = (total_wavelength >= blue_limit - tolerance) & (total_wavelength <= red_limit + tolerance)   
            empty_array = empty_array + indeces_line
            
        self.clean_indeces = invert(empty_array)
        
        #Calculating the error on the continuum
        self.linear_indeces = argwhere(diff(r_[False, self.clean_indeces, False])).reshape(-1, 2)
        self.linear_indeces[:, 1] -= 1

        std_arrays = []
        for interval in self.linear_indeces:
            #dz.data_plot(wave_obs[interval[0]:interval[1]], flux_obs[interval[0]:interval[1]], 'initial', 'blue', markerstyle='o')
            #dz.data_plot(wave_obs[interval[0]], flux_obs[interval[0]], 'initial', 'blue', markerstyle='o')
            #dz.data_plot(wave_obs[interval[1]], flux_obs[interval[1]], 'final', 'red', markerstyle='o')
            try:
                if (interval[1] - interval[0]) > 1:
                    x_region    = total_wavelength[interval[0]:interval[1]]
                    y_region    = total_flux[interval[0]:interval[1]]
                    lin_param   = self.lineal_mod.guess(y_region, x=x_region)
                    y_lineal    = lin_param['lineal_slope'].value * x_region + lin_param['lineal_intercept'].value
                    std_arrays.append(std(y_lineal - y_region))
                    #dz.data_plot(x_region, y_lineal, 'lineal fit', 'black')

            except:
                print 'fallo', interval[1] - interval[0]
            
        self.continuum_error = mean(std_arrays)   
                    
        return total_wavelength[self.clean_indeces], total_flux[self.clean_indeces]

    def model_parameters(self, obs_wave, obs_flux, err_obs):
    
        y_plus      = pymc.Uniform('He_abud', 0.050, 0.15)
#         Te          = pymc.Uniform('Te', self.ObjectData_dict['TOIII'] -2000, self.ObjectData_dict['TOIII']+2000)
        Te          = pymc.TruncatedNormal('Te', self.ObjectData_dict['TOIII'], self.ObjectData_dict['TOIII_error']**-2, a=self.ObjectData_dict['TOIII']- 6 * self.ObjectData_dict['TOIII_error'], b=self.ObjectData_dict['TOIII'] + 6 * self.ObjectData_dict['TOIII_error'])
#         Flux_Recomb = pymc.Uniform('Flux_Recomb', self.ObjectData_dict['Flux_Hbeta_Normalize'] - 4*self.ObjectData_dict['Error_Hbeta_Normalize'], self.ObjectData_dict['Flux_Hbeta_Normalize']+4 * self.ObjectData_dict['Error_Hbeta_Normalize'])
        Flux_Recomb = pymc.Normal('Flux_Recomb', self.ObjectData_dict['Flux_Hbeta_Normalize'], self.ObjectData_dict['Error_Hbeta_Normalize']**-2)
        
        #Calculate nebular continuum
        @pymc.deterministic
        def Calculate_Continuum(y_plus=y_plus, Te=Te, Flux_Recomb=Flux_Recomb, Wavelength_Range = obs_wave):

            return self.Calculate_Nebular_gamma(Te, Flux_Recomb, y_plus, HeIII_HII = 0.0, Wavelength_Range = Wavelength_Range)
         
        #Likelihood
        @pymc.stochastic(observed=True)
        def Likelihood_model(value=obs_flux, nebInt_cont = Calculate_Continuum, sigmaLines = err_obs):
            chi_F = np_sum(square(nebInt_cont - value) / square(sigmaLines))
            
            return - chi_F / 2
    
        return locals()
          
    def Calculate_Nebular_gamma(self, Te, Flux_Recomb, HeII_HII, HeIII_HII, Wavelength_Range):
        
        HII_HI          = 1.0
        
        Frac            = 1 + HeII_HII*4 + HeIII_HII*4
        
        Gamma_2q        = self.TwoPhotonContinuum(Te, Wavelength_Range)
             
        Gamma_FF_HI     = self.FreeFreeContinuum('HI', Te, Wavelength_Range)
        Gamma_FF        = Frac * Gamma_FF_HI
            
        Gamma_FB_HI     = self.FreeBoundContinuum_EP("HI", Te, Wavelength_Range)
        Gamma_FB_HeI    = self.FreeBoundContinuum_EP("HeI", Te, Wavelength_Range)
        Gamma_FB_HeII   = self.FreeBoundContinuum_EP("HeII", Te, Wavelength_Range)
            
        Gamma_FB        = HII_HI * Gamma_FB_HI + HeII_HII * Gamma_FB_HeI + HeIII_HII * Gamma_FB_HeII
    
        Gamma_Total     = Gamma_FB + Gamma_FF + Gamma_2q
        
        Gamma_lambda    = Gamma_Total * (self.c_AperS / power(Wavelength_Range,2))
        
#         NebularFlux_lambda    = self.Zanstra_Calibration_Hbeta(Te, Flux_Recomb, Gamma_Total, Wavelength_Range)
        NebularFlux_lambda    = self.Zanstra_Calibration_Halpha(Te, Flux_Recomb, Gamma_Total, Wavelength_Range)
        
        return NebularFlux_lambda
    
    def TwoPhotonContinuum(self, Te, Wavelength_Range):
        
        Gamma_2q = []
        
        #Coefficients for calculating A_2q The total radiative probability 2s -> 1s (s^-1)
        alpha_A = 0.88                                                                          
        beta_A  = 1.53
        gamma_A = 0.8
        C_A = 202.0                                                     # (s^-1)
        c_CGS = 2.99792458e10                                           #cm / s  
        h = 6.626068e-27                                                #erg s
        H0_Ionization_Energy = 13.6                                     #eV
        eV2erg = 1.602177e-12
        nu_0 = H0_Ionization_Energy * eV2erg / h                        #Hz
        
        A2q = 8.2249                                                    #(s^-1) Transition probability at lambda = 1215.7                                                              
        alpha_eff_2q = 6.5346e-11 * Te**-0.72315                   #(cm^3 s^-1) Effective Recombination coefficient (Este es mi ajuste. El de Molla es: 0.647e-10 * Te**-0.722)
        q2 = 5.92e-4 - 6.1e-9 * Te                                 #(cm^3 s^-1) Collisional transition rate coefficient for protons and electrons
        Lambda_0 =  1215.7
        
        for i in range(len(Wavelength_Range)):
            Lambda = float(Wavelength_Range[i])
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
            else:
                Gamma_2q.append(0.0)
    
        Np_Gamma_2q = array(Gamma_2q)    
        Np_Gamma_2q[isnan(Np_Gamma_2q)] = 0.0    
        
        return Np_Gamma_2q   

    def FreeFreeContinuum(self, Ion, Te, Wavelength_Range):
        
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
         
        CteTotal = (32 * (Z**2) * (e_proton**4)*h) / (3*(masseCGS**2)*(c_CGS**3)) * ((pi * H_RydbergEnergy / (3 * k * Te))**0.5) 
               
        for i in range(len(Wavelength_Range)):
            Wave = float(Wavelength_Range[i])
            nu = 2.99792458e18 / Wave  
            
            Flux_Comp1 = exp(((-1*h*nu)/(k*Te)))                                                                    #No unit    
            Flux_Comp2 = h*nu/(Z**2*e_proton*13.6057)
            Flux_Comp3 = k*Te/(h*nu)
    
            gff_mio = 1 + 0.1728*(Flux_Comp2)**0.33333*(1+2*Flux_Comp3) - 0.0496*(Flux_Comp2)**0.66667*(1+0.66667*Flux_Comp3+1.33333*(Flux_Comp3)**2)
            
            FluxValue = CteTotal * Flux_Comp1 * gff_mio                 # This is gamma_nu
            Flux_List.append(FluxValue)
             
        return array(Flux_List)

    def FreeBoundContinuum_EP(self, Ion, Te, Wavelength_Range):
        
        planckCGS = 6.626068e-27 # cm2 g / s  ==  erg s
        cspeedA = 3.e18 # A / s
        rydbergerg = 2.1798741e-11 # Rydberg to erg
        t4 = Te/1.e4
    
        PixelNumber = len(Wavelength_Range)
        
        nu          = zeros(PixelNumber)
        eryd        = zeros(PixelNumber)
        eryd_low    = zeros(PixelNumber)
        
        for i in range(len(Wavelength_Range)): 
            wave = Wavelength_Range[i]                             # wavelength en Angstroms
            nu[i] = cspeedA/wave                                        # frecuencia en Hz
            eryd[i] = planckCGS*cspeedA/(rydbergerg*wave)               # energia en Rydbergs
        
        if Ion == "HI":
            a = self.a_HI
        if Ion == "HeI":
            a = self.a_HeI
        if Ion == "HeII":
            a = self.a_HeII
    

        nTe = int(split(a[0])[0])                                       # number of Te columns
        nener = int(split(a[0])[1])                                     # number of energy points rows 
        skip = int(1+ceil(nTe/8.))                                      # 8 es el numero de valores de Te por fila.
        
        temp = zeros(nTe)
        for i in range(1,skip) :
            tt = split(a[i])
            for j in range(0,len(tt)) :
                temp[8*(i-1)+j] = tt[j]
        
        rte2 = 0

        if bisect_right(temp,log10(Te)) - bisect_left(temp,log10(Te)) == 1 : 
            rte1 = bisect_left(temp,log10(Te))  # Te existe
        elif bisect_right(temp,log10(Te)) - bisect_left(temp,log10(Te)) == 0 :
            rte1 = bisect_left(temp,log10(Te))-1 # interpola Te
            rte2 = bisect_right(temp,log10(Te))  # 
        else :
            print 'ERROR IN Te COLUMN DETECTION FOR Gamma_fb', Te
    
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
                gamma[i-skip] = float(split(a[i])[2+rte1]) + (float(split(a[i])[2+rte2])-float(split(a[i])[2+rte1]))/(10**temp[rte2]-10**temp[rte1]) * (Te-10**temp[rte1])
            
            ener_low = thresh*ener
            atmp=thresh.nonzero()[0]
            for i in range(1,len(atmp),2) : ener_low[atmp[i]] = ener_low[atmp[i-1]]
            for i in range(0,len(ener_low)) : 
                if ener_low[i] == 0 : ener_low[i] = ener_low[i-1]
            
            eryd_low = interp(eryd, ener, ener_low)
            Gamma_fb_f = interp(eryd, ener, gamma)  # linear interpolation at grid of interest
            Gamma_fb_f = Gamma_fb_f * 1e-40*t4**-1.5*exp(-15.7887*(eryd-eryd_low)/t4)
    
        return Gamma_fb_f      

    def Zanstra_Calibration_Hbeta(self, Te, Flux_EmLine, Gamma_nu, Wavelength_Range):
        
        # Pequignot et al. 1991
        t4 = Te / 10000
        
        alfa_eff = 0.668e-13 * t4**-0.507 / (1 + 1.221*t4**0.653)
        lambda_EmLine = 4862.683

        NebularFlux_lambda = Gamma_nu * lambda_EmLine * Flux_EmLine / (alfa_eff * self.planckCGS * power(Wavelength_Range,2)) # erg s-1 [cm-2] A-1 
        
        return NebularFlux_lambda

    def Zanstra_Calibration_Halpha(self, Te, Flux_EmLine, Gamma_nu, Wavelength_Range):
        
        # Pequignot et al. 1991
        t4              = Te / 10000
        
        alfa_eff        = 2.708e-13 * t4**-0.648 / (1 + 1.315*t4**0.523)
        lambda_EmLine   = 6562.819

        NebularFlux_lambda = Gamma_nu * lambda_EmLine * Flux_EmLine / (alfa_eff * self.planckCGS * power(Wavelength_Range,2)) # erg s-1 [cm-2] A-1 
        
        return NebularFlux_lambda

