from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v23 import Run_MCMC

default_db_folder   = '/home/vital/Astrodata/Inference_output/'
iterat, burn, thin  = 10000, 0, 1
sim_model           = 'Reappeting_tests'
sim_components      = '_He_S_O_neb_stellar'
sim_name            = sim_model + sim_components
params_list         = ['y_plus', 'T_He', 'T_low','ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund','sigma_star','Av_star'] #'ChiSq_Recomb','ChiSq_S','ChiSq_O']      
burning             = 5000.0
                        
#Generate dazer object
dz = Dazer()
bm = Run_MCMC()

#Generate the synthetic data
bm.calculate_synthFluxes(sim_components)

#Declare line to measure #We might need to multiply this frames by the waves
line                    = 'H1_6563A'
bm.Current_Label        = bm.lick_idcs_df.loc[line].name
bm.Current_Ion          = bm.lick_idcs_df.loc[line].Ion
bm.Current_TheoLoc      = bm.lick_idcs_df.loc[line].lambda_theo * (1 + bm.obj_data['z_star'])
selections              = bm.lick_idcs_df.loc[line][3:9].values

#Measure the line
line_data = bm.measure_line(bm.obj_data['obs_wave_resam'], bm.emission_SED, selections, None, Measuring_Method = 'lmfit', store_data = False) 
print 'Flux measured'
print line_data['flux_intg']
print line_data['flux_gauss0']
print 'self.Halpha_norm', bm.Halpha_norm
print (bm.Hbeta_Flux/bm.stellar_SED['normFlux_stellar'])
dz.FigConf()
dz.data_plot(bm.obj_data['obs_wave_resam'], bm.emission_SED, label = 'Emision hydrogen')
dz.data_plot(line_data.x_resample, line_data.y_resample, label = 'manual fitting')
dz.FigWording(xlabel = 'Wavelength', ylabel = 'Flux', title = '')
dz.display_fig()

#Select right model according to data
bm.select_inference_model(sim_components)     
             
#Variables to save
db_address  = '{}{}_it{}_burn{}_thin{}'.format(default_db_folder, sim_name, iterat, burn, thin) 
    
#Run sampler
bm.run_pymc2(db_address, iterat, burn, thin, variables_list=params_list)
      
#Load database
pymc2_db, stat_db_dict = bm.load_pymc_database_manual(db_address, burning=burning)
     
#Traces plot
dz.traces_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_tracesPlot', save_pickle = False)
       
#Posteriors plot
dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_posteriorPlot', save_pickle = False)
       
#Posteriors plot
dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=4, n_rows=4)
dz.save_manager(db_address + '_acorrPlot', save_pickle = False)
      
#Corner plot
dz.corner_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_cornerPlot', save_pickle = False)
     
print '\nData treated'