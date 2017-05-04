from numpy import linspace, zeros, exp, array, arange, histogram, random, concatenate
import pymc

from Astro_Libraries.he_bayesian            import He_Inference_Abundance
from CodeTools.PlottingManager                      import myPickle


#Generate dazer object
pv = myPickle()

#Define plot frame and colors
pv.FigFormat_One(ColorConf='Night1')

He_abun_tools = He_Inference_Abundance()

print 'the matrix'
for matrix in He_abun_tools.H_CollCoeff_Matrix:
    print He_abun_tools.H_CollCoeff_Matrix[0]

y_plus0, n_e0, a_He0, tau0, Te_0, cHbeta0, a_H0, xi0, EW_Hbeta0, err_Flux_H0, err_Flux_He0 = He_abun_tools.Example_SyntheticParameters2()

print 'Initial values'
print 'xi', xi0

H_Flux  =   He_abun_tools.H_Flux_theo(Te_0, n_e0, xi0, cHbeta0, a_H0, h_Hlambda=1, a_Hbeta=1, h_Hbeta=1, EW_Hbeta=EW_Hbeta0)
He_Flux =   He_abun_tools.He_Flux_theo(Te_0, n_e0, xi0, cHbeta0, a_He0, h_Helambda=1, a_Hbeta=1, h_Hbeta=1, EW_Hbeta=EW_Hbeta0, y_plus=y_plus0, tau=tau0)
He_Flux =   He_abun_tools.He_Flux_theo_nof(Te_0, n_e0, xi0, cHbeta0, a_He0, h_Helambda=1, a_Hbeta=1, h_Hbeta=1, EW_Hbeta=EW_Hbeta0, y_plus=y_plus0, tau=tau0)

print 'Hydrogen fluxes'
print H_Flux
print 'Helium Fluxes'
print He_Flux
print 'Helium Fluxes no helium colisional excitation'
print He_Flux
juntos = concatenate([H_Flux, He_Flux])
# print 'all together', len(juntos), type(juntos), juntos[0], juntos[-1]

# print 'Eo', He_abun_tools.He1.getEmissivity(Te_0, n_e0, label=He_abun_tools.Helium_Labels[0]) / He_abun_tools.H1.getEmissivity(Te_0, n_e0, label=He_abun_tools.Hydrogen_Labels[2]) 

# print 'Hydrogen theoretical fluxes'
# print array([He_abun_tools.H1.getEmissivity(Te_0, n_e0, label=He_abun_tools.Hydrogen_Labels[0]), He_abun_tools.H1.getEmissivity(Te_0, n_e0, label=He_abun_tools.Hydrogen_Labels[1]), He_abun_tools.H1.getEmissivity(Te_0, n_e0, label=He_abun_tools.Hydrogen_Labels[3])]) / He_abun_tools.H1.getEmissivity(Te_0, n_e0, label=He_abun_tools.Hydrogen_Labels[2]) 

