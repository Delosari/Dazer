from bin.lib.Astro_Libraries.Abundances_InferenceModel_Helium_v15 import Run_MCMC

#MCMC configuration run
db_name_ext         = 'v2_CorrectingFit'
synth_model         = 'Model3'
iterat, burn, thin  = 10000, 1500, 2

#Generate dazer object
bm = Run_MCMC()
  
#Generate the synthetic data
bm.calculate_synthFluxes(synth_model)
     
#Select right model according to data
bm.select_inference_model(synth_model)     
            
#Variables to save
params_list = ['ChiSq','y_plus', 'Te', 'ne', 'tau', 'cHbeta', 'xi']             
                   
#Run sampler
bm.run_pymc2(iterat, burn, thin, db_name_ext, variables_list=params_list)

