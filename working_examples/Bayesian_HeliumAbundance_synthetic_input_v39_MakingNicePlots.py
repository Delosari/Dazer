# import seaborn as sns
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# sns.set_style("white")
# iris = sns.load_dataset("iris")    
#  
# g = sns.PairGrid(iris)
# g.map_diag(plt.hist, bins=20)
#  
# # def pairgrid_heatmap(x, y, **kws):
# #     cmap = sns.light_palette(kws.pop("color"), as_cmap=True)
# #     plt.hist2d(x, y, cmap=cmap, cmin=1, **kws)
# # g.map_offdiag(pairgrid_heatmap, bins=20, norm=LogNorm())
#      
# # def hexbin(x, y, color, **kwargs):
# #     cmap = sns.light_palette(color, as_cmap=True)
# #     plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)
# # g = sns.FacetGrid(iris, size=4)
# # g.map(hexbin, extent=[0, 50, 0, 10])
#  
# # def hexbin(x, y, color, **kwargs):
# #     cmap = sns.light_palette(color, as_cmap=True)
# #     plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)
# # g.map_offdiag(hexbin, bins=20, norm=LogNorm())
# 
# def kdebin(x, y, color, n_levels, **kwargs):
#     cmap = sns.light_palette(color, as_cmap=True)
#     sns.kdeplot(x, y, cmap='Blues', shade=True, shade_lowest=False)
# g.map_offdiag(kdebin, n_levels=20)
# 
# plt.show()

import pandas as pd
import numpy as np
from collections import OrderedDict
from dazer_methods import Dazer
from lib.Astro_Libraries.Abundances_InferenceModel_Helium_v51 import Run_MCMC
from seaborn import pairplot, light_palette, PairGrid, cubehelix_palette, kdeplot
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
  
  
iterat, burn, thin  = 15000, 0, 1
sim_model           = 'PuttingAllTogether'
sim_components      = '_neb_stars_Abunds'
obs_metals          = ['H1', 'He1', 'S2', 'S3', 'O2', 'O3', 'N2', 'Ar3', 'Ar4']
sim_name            = sim_model + sim_components
params_list         = ['He1_abund', 'T_He', 'T_low', 'ne','tau','cHbeta','xi','S2_abund','S3_abund','O2_abund','O3_abund', 'N2_abund', 'Ar3_abund', 'Ar4_abund', 'sigma_star', 'Av_star'] 
burning             = 7500
                       
#Generate dazer object
dz = Dazer()
bm = Run_MCMC()
  
# #Generate the synthetic data
bm.calculate_simObservation(sim_components, obs_lines = obs_metals)
  
#Variables to save
db_address = '{}{}_it{}_burn{}'.format(bm.paths_dict['inference_folder'], sim_name, iterat, burn) 
              
# #Load database
pymc2_db, stat_db_dict = bm.load_pymc_database_manual(db_address, burning, params_list)
 
rc={'axes.labelsize': 25, 'font.size': 32, 'legend.fontsize': 32.0, 'axes.titlesize': 32}
# plt.rcParams.update(**rc)
sns.set(rc=rc)
 
#Making the df for the seaborn plot
traces_dict = OrderedDict()
for entry in ['He1_abund', 'T_He', 'ne','tau','cHbeta','xi','tau']:
    entry_latex_key = dz.labels_latex_dic[entry]
    traces_dict[entry_latex_key] = stat_db_dict[entry]['trace']
df_seaborn = pd.DataFrame.from_dict(traces_dict)
g = PairGrid(df_seaborn)
g.map_diag(plt.hist, bins=20)
 
# def hexbin(x, y, color, **kwargs):
#     cmap = light_palette(color, as_cmap=True)
#     plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)
# g.map_offdiag(hexbin, bins=10)
def kdebin(x, y, color, n_levels, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    sns.kdeplot(x, y, cmap='Blues', shade=True, shade_lowest=False)
g.map_offdiag(kdebin, n_levels=20)


dz.savefig('/home/vital/Dropbox/Astrophysics/Seminars/Stasinska conference/helium_box.png', reset_fig=True)
 
rc={'axes.labelsize': 35, 'font.size': 35, 'legend.fontsize': 35.0, 'axes.titlesize': 35}
# plt.rcParams.update(**rc)
sns.set(rc=rc)
 
#Making the df for the seaborn plot
traces_dict = OrderedDict()
for entry in ['T_low', 'ne', 'S2_abund','S3_abund','O2_abund','O3_abund', 'N2_abund', 'Ar3_abund', 'Ar4_abund']:
    entry_latex_key = dz.labels_latex_dic[entry]
    traces_dict[entry_latex_key] = stat_db_dict[entry]['trace']
df_seaborn = pd.DataFrame.from_dict(traces_dict)
g = PairGrid(df_seaborn)
g.map_diag(plt.hist, bins=20)
 
# def hexbin(x, y, color, **kwargs):
#     cmap = light_palette(color, as_cmap=True)
#     plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)
# g.map_offdiag(hexbin, bins=20)
def kdebin(x, y, color, n_levels, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    sns.kdeplot(x, y, cmap='Blues', shade=True, shade_lowest=False)
g.map_offdiag(kdebin, n_levels=20)

 
dz.savefig('/home/vital/Dropbox/Astrophysics/Seminars/Stasinska conference/metals_box.png', resolution=150.0)
 
# print 'DONE'
# pairplot(df_seaborn)
# 
# g = sns.FacetGrid(tips, hue="time", col="time", size=4)
# g.map(hexbin, "total_bill", "tip", extent=[0, 50, 0, 10])
#   
# plt.show()
