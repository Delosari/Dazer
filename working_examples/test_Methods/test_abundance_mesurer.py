import  pyneb as pn
from    dazer_methods import Dazer
from    numpy import random
from    uncertainties import ufloat


def Sinthetic_cloud():
    
    Temp_dict   = {}
    TSIII       = ufloat(12000, 500)
    TOIII       = (1.0807 * TSIII/10000 - 0.0846) * 10000
    nSII        = ufloat(150, 50)
    
    OII_HII     = ufloat(2.5e-5, 0.1e-6)
    OIII_HIII   = ufloat(2.0e-5, 0.2e-6)

    SII_HII     = ufloat(1.5e-6, 0.1e-7)
    SIII_HII    = ufloat(2.0e-6, 0.2e-7)

    ArIII_HII   = ufloat(1.5e-7, 0.1e-8)
    ArIV_HII    = ufloat(1.0e-7, 0.2e-8)
    
    NII_HII     = ufloat(3.2e-6, 0.5e-7)
    
    Temp_dict['TOIII']  = random.normal(flux, error, size = self.MC_array_len)
    Temp_dict['TSIII']  = random.normal(flux, error, size = self.MC_array_len)
    Temp_dict['nSII']   = random.normal(flux, error, size = self.MC_array_len)

    self.Properties_dict['TOIII_approxfrom_TSIII_pn'] = (1.0807 * self.Properties_dict['TSIII_pn']/10000 - 0.0846) * 10000

    return

#Generate dazer object
# dz = Dazer()

MC_array_len = 100

Hbeta_label = 'H1_4861A'

#Set atomic data 
pn.atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')
pn.atomicData.setDataFile('s_iii_coll_HRS12.dat')

#Declare hydrogen atom so it becomes readily available:
S2_atom     = pn.Atom('S', 2)
S3_atom     = pn.Atom('S', 3)
Ar3_atom    = pn.Atom('Ar', 3)
Ar4_atom    = pn.Atom('Ar', 4)
N2_atom     = pn.Atom('N', 2)
O2_atom     = pn.Atom('O', 2)
O3_atom     = pn.Atom('O', 3)
            
H1_atom     = pn.RecAtom('H', 1)
He1_atom    = pn.RecAtom('He', 1)
He2_atom    = pn.RecAtom('He', 2)

Te = 10000.0
ne = 100.0
abundance = 0.5

SIII_9069A_Emis = S3_atom.getEmissivity(Te, ne, wave = 9069)
SIII_9531A_Emis = S3_atom.getEmissivity(Te, ne, wave = 9531)

SIII_9069A_flux = abundance * SIII_9069A_Emis
SIII_9531A_flux = abundance * SIII_9531A_Emis
Hbeta_flux      = H1_atom.getEmissivity(Te, ne, wave = 4861.0)



random.normal(flux, error, size = MC_array_len)



# Lines_Sum   = SIII_9069A_flux + SIII_9531A_flux
# 
# print 'Variables', Lines_Sum, Te, ne, Hbeta_flux
# 
# SIII_HII_abundance = S3_atom.getIonAbundance(int_ratio=Lines_Sum, tem=Te, den=ne, to_eval = 'L(9069)+L(9531)', Hbeta = Hbeta_flux)
# 
# print 'Predicted abundance', SIII_HII_abundance 
# 
# print 'Simulation completed'