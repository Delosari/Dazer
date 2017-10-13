from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v28 import Run_MCMC

iterat, burn, thin  = 20000, 0, 1
sim_model           = 'Just_metals_RS3'
sim_components      = '_He_S_O_neb_stellar'
obs_metals          = ['H', 'He1', 'He2', 'S2', 'S3', 'O2', 'O3', 'N2', 'Ar3', 'Ar4']
sim_name            = sim_model + sim_components
params_list         = ['T_low', 'ne','cHbeta','S2_abund','S3_abund','O2_abund','O3_abund','N2_abund','Ar3_abund','Ar4_abund'] #['y_plus', 'T_He', 'T_low', 'ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund','N2_abund','Ar3_abund','Ar4_abund','sigma_star','Av_star'] 
burning             = 0
                        
#Generate dazer object
dz = Dazer()
bm = Run_MCMC()

#Generate the synthetic data
bm.calculate_synthFluxes(sim_components, obs_lines = obs_metals)

# bm.obj_data['T_low']
# bm.obj_data['T_High']
# bm.obj_data['ne']
# bm.obj_data['cHbeta']
# bm.obj_data['']
# bm.calculate_colExcit_flux(T_low, T_High, ne, cHbeta, bm.obj_data['colLine_waves'], bm.obj_data['colLine_ions'], abund_dict, bm.obj_data['colLine_flambda'], bm.obj_data['colLine_pynebCode'])

#Select right model according to data
bm.select_inference_model(sim_components)     

#Variables to save
db_address  = '{}{}_it{}_burn{}'.format(bm.paths_dict['inference_folder'], sim_name, iterat, burn) 
 
#Run sampler
bm.run_pymc2(db_address, iterat, variables_list=params_list, prefit=False)
        
#Load database
pymc2_db, stat_db_dict = bm.load_pymc_database_manual(db_address, burning, params_list)
        
#Traces plot
print '-Generating traces plot'
dz.traces_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_tracesPlot_Test', save_pickle = False)
          
#Posteriors plot
print '-Generating posteriors plot'
dz.posteriors_plot(params_list, pymc2_db, stat_db_dict)
dz.save_manager(db_address + '_posteriorPlot', save_pickle = False)
           
#Posteriors plot
print '-Generating acorrelation plot'
dz.acorr_plot(params_list, pymc2_db, stat_db_dict, n_columns=4, n_rows=4)
dz.save_manager(db_address + '_acorrPlot', save_pickle = False)
         
#Corner plot
print '-Generating corner plot'
dz.corner_plot(params_list, pymc2_db, stat_db_dict, plot_true_values=True)
dz.save_manager(db_address + '_cornerPlot', save_pickle = False)
        
print '\nData treated'