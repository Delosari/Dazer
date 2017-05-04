#!/usr/bin/python

from numpy                                              import square, concatenate, sum
from numpy.random                                       import normal
import pymc

from Astro_Libraries.Abundances_InferenceModel_Helium   import He_Inference_Abundance
from CodeTools.PlottingManager                          import myPickle


#Run an example with a smaller interactions
#Run an example with initial values
#Run an example with MAP
#Run an example with greater intensities
#Graph with all the traces evolution
def Bayesian_HeliumAbudnance_Analysis(H_He_ObsFlux, H_He_ObsError, TOIII, EW_Hbeta, EW_Hbeta_err, Hbeta_error):
    
    y_plus  =   pymc.Uniform(           'He_abud',  0.065,  0.095)
    ne      =   pymc.TruncatedNormal(   'n_e',      100,    250**-2,    a = 0.0 ,   b = 1000.0)
    a_He    =   pymc.TruncatedNormal(   'abs_He',   1.0,    0.5**-2,    a = 0.0,    b = 5.0)
    tau     =   pymc.TruncatedNormal(   'Tau',      0.75,   0.5**-2,    a = 0.0,    b = 7.0)
    Te      =   pymc.Normal(            'T_e',      17000,  2000**-2)
    cHbeta  =   pymc.TruncatedNormal(   'c_Hbeta',  0.15,   0.05**-2,   a = 0.0,    b = 3.0)
    a_H     =   pymc.TruncatedNormal(   'abs_H',    1.0,    1**-2,      a = 0.0,    b = 14.0)
    xi      =   pymc.TruncatedNormal(   'Xi',       10,     200**-2,    a = 0.0,    b = 1000.0)
        
    #Calculate Hydrogen theoretical flux
    @pymc.deterministic
    def HFlux_theo(Te = Te, ne = ne, xi=xi, cHbeta=cHbeta, a_H=a_H, EW_Hbeta=EW_Hbeta):
        return bm.H_Flux_theo(Te, ne, xi, cHbeta, a_H, h_Hlambda=1, a_Hbeta=1, h_Hbeta=1, EW_Hbeta=EW_Hbeta)
    
    #Calculate Helium theoretical flux    
    @pymc.deterministic
    def He_Flux_theo_nof(Te = Te, ne = ne, xi=xi, cHbeta=cHbeta, a_H=a_H, a_He=a_He, EW_Hbeta=EW_Hbeta, y_plus = y_plus, tau=tau):
        return bm.He_Flux_theo(Te, ne, xi, cHbeta, a_H, a_He, h_Helambda=1, a_Hbeta=1, h_Hbeta=1, EW_Hbeta=EW_Hbeta, y_plus=y_plus, tau=tau)

    #Combine theoretical fluxes into a single array
    @pymc.deterministic
    def H_He_TheoFlux(HFlux_theo=HFlux_theo, HeFlux_theo=He_Flux_theo_nof):
        return concatenate([HFlux_theo, HeFlux_theo])
    
    #Chi_EquivalentWidth
    @pymc.deterministic
    def ChiSq_EW(EW_Hbeta_obs=EW_Hbeta, EW_Hbeta_sig = EW_Hbeta_err):
        EW_HbetaSynth   = normal(EW_Hbeta_obs, EW_Hbeta_sig, 1)
        chi_ew          = square(EW_HbetaSynth - EW_Hbeta_obs) / square(EW_Hbeta_sig)
        L_chi_ew        = - chi_ew / 2
        return L_chi_ew
    
    #Chi_Temperature
    @pymc.deterministic
    def ChiSq_T(T_Synth=Te, T_Meas = TOIII, sigma = TOIII*0.2):
        chi_tem          = square(T_Synth - T_Meas) / square(sigma)
        L_chi_tem        = - chi_tem / 2
        return L_chi_tem

    #Likelihood
    @pymc.stochastic(observed=True)
    def Likelihood_model(value=H_He_ObsFlux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = H_He_ObsError, ChiSq_EW=ChiSq_EW, ChiSq_T=ChiSq_T):
        chi_F           = sum(square(H_He_TheoFlux - value) / square(sigmaLines))
        ChiSq_F         = - chi_F / 2
        return ChiSq_F + ChiSq_EW + ChiSq_T
    
    #Deterministic method to track the evolution of the chi:
    @pymc.deterministic()
    def ChiSq(H_He_ObsFlux=H_He_ObsFlux, H_He_TheoFlux = H_He_TheoFlux, sigmaLines = H_He_ObsError, ChiSq_EW=ChiSq_EW, ChiSq_T=ChiSq_T):
        chi_F           = sum(square(H_He_TheoFlux - H_He_ObsFlux) / square(sigmaLines))
        ChiSq_F         = - chi_F / 2
        return ChiSq_F + ChiSq_EW + ChiSq_T
   
    return locals()

#DeviationPrior

#Generate dazer object
pv = myPickle()
bm = He_Inference_Abundance()

#Define plot frame and colors
pv.FigFormat_One(ColorConf  =   'Night1')

#Import observed data:
H_Flux_Obs, He_Flux_Obs         = bm.Import_Synthetic_Fluxes(case = 'case2')
H_Flux_Obs_err, He_Flux_Obs_er  = H_Flux_Obs*0.01, He_Flux_Obs*0.02

Hydrogen_Helium_Fluxes          = concatenate([H_Flux_Obs, He_Flux_Obs])
Hydrogen_Helium_Error           = concatenate([H_Flux_Obs_err, He_Flux_Obs_er])

Hbeta_error                     = 0.05 * 250

#Declare the Bayesian dictionary and database storage location
MCMC_dict                       = Bayesian_HeliumAbudnance_Analysis(Hydrogen_Helium_Fluxes, Hydrogen_Helium_Error, TOIII = 17000, EW_Hbeta = 250, EW_Hbeta_err = 10, Hbeta_error = Hbeta_error)
StoringDataFolder               = '/home/vital/Workspace/X_ModelData/MCMC_databases/' 


#New physical values to calculate the fluxes
db_name                         = 'he_Abundance_10000_OldRemade'   
csv                             = 'he_Abundance_Global_10000_OldRemade'
M = pymc.MCMC(MCMC_dict, db = 'pickle', dbname = StoringDataFolder + db_name)
M.sample(iter=10000)
M.write_csv(StoringDataFolder + csv, variables=['ChiSq', 'He_abud', 'T_e', 'n_e','abs_H', 'abs_He', 'Tau', 'c_Hbeta', 'Xi'])
M.db.close()
#     y_plus  =   pymc.Uniform(           'He_abud',  0.065,  0.095)
#     ne      =   pymc.TruncatedNormal(   'n_e',      500,    250**-2,    a = 0.0 ,   b = 1000.0)
#     a_He    =   pymc.TruncatedNormal(   'abs_He',   0.5,    0.5**-2,    a = 0.0,    b = 5.0)
#     tau     =   pymc.TruncatedNormal(   'Tau',      1.0,    0.5**-2,    a = 0.0,    b = 7.0)
#     Te      =   pymc.Normal(            'T_e',      16000,  2000**-2)
#     cHbeta  =   pymc.TruncatedNormal(   'c_Hbeta',  0.1,    0.05**-2,   a = 0.0,    b = 3.0)
#     a_H     =   pymc.TruncatedNormal(   'abs_H',    1.0,    1**-2,      a = 0.0,    b = 14.0)
#     xi      =   pymc.TruncatedNormal(   'Xi',       1,      200**-2,    a = 0.0,    b = 1000.0)

#Run the MCMC for 10000 iterations all centered  gaussian
# db_name                         = 'hE_Abundance_7000_AllGaussianGraph3_nofinHe_valueobservedTrue'   
# csv                             = 'he_Abundance_Global_7000_AllGaussianGraph3_nofinHe_valueobservedTrue'
# M = pymc.MCMC(MCMC_dict, db = 'pickle', dbname = StoringDataFolder + db_name)
# M.sample(iter=10000, burn=500)
# M.write_csv(StoringDataFolder + csv, variables=['ChiSq', 'He_abud', 'T_e', 'n_e','abs_H', 'abs_He', 'Tau', 'c_Hbeta', 'Xi'])
# M.db.close()

#2373 segundos
#y_plus  =   pymc.Uniform(           'He_abud',  0.065,  0.095)
#ne      =   pymc.TruncatedNormal(   'n_e',      100,    250**-2,    a = 0.0 ,   b = 1000.0, value=100)
#a_He    =   pymc.TruncatedNormal(   'abs_He',   1.0,    0.5**-2,    a = 0.0,    b = 5.0)
#tau     =   pymc.TruncatedNormal(   'Tau',      0.2,    0.5**-2,    a = 0.0,    b = 7.0, value=0.2)
#Te      =   pymc.Normal(            'T_e',      18000,  2000**-2)
#cHbeta  =   pymc.TruncatedNormal(   'c_Hbeta',  0.1,    0.05**-2,   a = 0.0,    b = 3.0)
#a_H     =   pymc.TruncatedNormal(   'abs_H',    1.0,    1**-2,      a = 0.0,    b = 14.0)
#xi      =   pymc.TruncatedNormal(   'Xi',       1,      200**-2,    a = 0.0,    b = 1000.0)

#Run the MCMC for 1000 iterations Uniform scheme
# db_name             = 'hE_Abundance_1000'   
# csv                 = 'he_Abundance_Global_1000'
# M = pymc.MCMC(MCMC_dict, db = 'pickle', dbname = StoringDataFolder + db_name)
# M.sample(iter=1000)
# M.write_csv(StoringDataFolder + csv, variables=['ChiSq', 'He_abud', 'T_e', 'n_e','abs_H', 'abs_He', 'Tau', 'c_Hbeta', 'Xi'])
# M.db.close()
#40 minutos

#Run the MCMC for 10000 iterations Uniform scheme udating the xi definition
# db_name             = 'hE_Abundance_10000_2'   
# csv                 = 'he_Abundance_Global_10000_2'
# M = pymc.MCMC(MCMC_dict, db = 'pickle', dbname = StoringDataFolder + db_name)
# M.sample(iter=10000, burn = 100)
# M.write_csv(StoringDataFolder + csv, variables=['ChiSq', 'He_abud', 'T_e', 'n_e','abs_H', 'abs_He', 'Tau', 'c_Hbeta', 'Xi'])
# M.db.close()

#Using the priors (the same as before
# Te      =   pymc.Uniform(   'T_e',     10000    ,22000)
# ne      =   pymc.Uniform(   'n_e',      10      ,500)
# xi      =   pymc.Uniform(   'Xi',       1       ,2500)    
# tau     =   pymc.Uniform(   'Tau',      0.0     ,1.8)
# a_H     =   pymc.Uniform(   'abs_H',    1       ,3)
# a_He    =   pymc.Uniform(   'abs_He',   0.8     ,1.2)
# cHbeta  =   pymc.Uniform(   'c_Hbeta',  0       ,1.4)
# y_plus  =   pymc.Uniform(   'He_abud',  0.065   ,0.095)
# sobre 20-30 minutos

#Run the MCMC for 10000 iterations Some gaussian distributions scheme
# db_name             = 'hE_Abundance_10000_Gaussian3'   
# csv                 = 'he_Abundance_Global_10000_Gaussian3'
# M = pymc.MCMC(MCMC_dict, db = 'pickle', dbname = StoringDataFolder + db_name)
# M.sample(iter=10000, burn = 1000,)
# M.write_csv(StoringDataFolder + csv, variables=['ChiSq', 'He_abud', 'T_e', 'n_e','abs_H', 'abs_He', 'Tau', 'c_Hbeta', 'Xi'])
# M.db.close()
#1700 segundos

#Run the MCMC for 10000 iterations Some gaussian distributions scheme
# db_name             = 'hE_Abundance_10000_GaussianGraph'   
# csv                 = 'he_Abundance_Global_10000_GaussianGraph'
# M = pymc.MCMC(MCMC_dict, db = 'pickle', dbname = StoringDataFolder + db_name)
# 
# #Make a slightly pretty graph out of it
# graph = pymc.graph.graph(M)
# graph.write_png(StoringDataFolder + "Bayesian_Model.png")

#     Te      =   pymc.Uniform(           'T_e',     10000   ,22000 )
#     ne      =   pymc.TruncatedNormal(   'n_e',      100,    400**-2, a = 0.0 , b = 9e99)
#     xi      =   pymc.Uniform(           'Xi',       1       ,2500 )    
#     tau     =   pymc.TruncatedNormal(   'Tau',      1.0,    0.5**-2, a = 0.0 , b = 9e99)
#     a_H     =   pymc.TruncatedNormal(   'abs_H',    1.0,    0.5**-2, a = 0.0 , b = 9e99)
#     a_He    =   pymc.TruncatedNormal(   'abs_He',   1.0,    0.5**-2, a = 0.0 , b = 9e99)
#     cHbeta  =   pymc.Uniform(           'c_Hbeta',  0       ,1.4 )
#     y_plus  =   pymc.Uniform(           'He_abud',  0.065   ,0.095 )

#1700 segundos
print 'Analyse completed'



    #Uniform priors
#     Te      =   pymc.Uniform(   'T_e',     10000   ,22000 )
#     ne      =   pymc.Uniform(   'n_e',      10      ,500 )
#     xi      =   pymc.Uniform(   'Xi',       1       ,2500 )    
#     tau     =   pymc.Uniform(   'Tau',      0.0     ,1.8 )
#     a_H     =   pymc.Uniform(   'abs_H',    1       ,3 )
#     a_He    =   pymc.Uniform(   'abs_He',   0.8     ,1.2 )
#     cHbeta  =   pymc.Uniform(   'c_Hbeta',  0       ,1.4 )
#     y_plus  =   pymc.Uniform(   'He_abud',  0.065   ,0.095 )