import numpy as np
import pyneb as pn
import pyCloudy as pc
from dazer_methods import Dazer

import matplotlib.pyplot as plt
dz = Dazer()
dz.FigConf()

#Remove these pregenerated files
components_remove = [#['radius', '.rad'],
                     #['continuum', '.cont'],
                     #['physical conditions', '.phy'],
                     ['overview', '.ovr'],
                     ['heating', '.heat'],
                     ['cooling', '.cool'],
                     ['optical depth', '.opd']]

for item in components_remove:
    pc.config.SAVE_LIST.remove(item)

elements_remove = [['hydrogen','.ele_H'],
                   ['helium','.ele_He'],
                   ['carbon', '.ele_C'],
                   ['nitrogen', '.ele_N'],
                   ['oxygen', '.ele_O'],
                   ['argon', '.ele_Ar'],
                   ['neon', '.ele_Ne'],
                   ['sulphur', '.ele_S'],
                   ['chlorin', '.ele_Cl'],
                   ['iron', '.ele_Fe'],
                   ['silicon', '.ele_Si']]

for item in elements_remove:
    pc.config.SAVE_LIST_ELEMS.remove(item)

#These are the commands common to all the models (here only one ...)
emis_tab = ['H  1 4861.33A',
            'H  1 6562.81A',
            'HE 1 3888.63A',
            'HE 1 4026.20A',
            'HE 1 4471.49A',
            'HE 1 5875.64A',
            'HE 1 6678.15A',
            'HE 1 7065.22A']

# Tell pyCloudy where your cloudy executable is:
dir_ = '/home/vital/Desktop/tests_cloudy/'
pc.config.cloudy_exe = '/home/vital/cloudy/source/cloudy.exe'

# Define some parameters of the model:
model_name      = 'emissivity_grid'
full_model_name = '{0}{1}'.format(dir_, model_name)

# Physical parameters
dens    = 2.        #log cm-3
Teff    = 45000.    #K
qH      = 47.0      #s-1
dist    = 16.0      #kpc

def make_vary_Te_ne(Te, ne, simu_address):
    logTe = np.log10(Te)
    cloudy_model = pc.CloudyInput('{}M_i_Te{}_ne{}'.format(simu_address, Te, ne))
    cloudy_model.set_other(('init file "hheonly.ini"'))
    cloudy_model.set_other(('element helium abundance -3'))
    cloudy_model.set_BB(Teff=Teff, lumi_unit='q(H)', lumi_value=qH)
    cloudy_model.set_radius(r_in=dist)
    cloudy_model.set_other(('constant temperature {}'.format(logTe)))
    cloudy_model.set_cste_density(dens)
    cloudy_model.set_other(('database H-like levels large element hydrogen'))
    cloudy_model.set_other(('database H-like levels large element helium'))
    cloudy_model.set_other(('set dr 0'))
    cloudy_model.set_stop(('zone 1'))
    cloudy_model.set_emis_tab(emis_tab)  # better use read_emis_file(file) for long list of lines, where file is an external file
    cloudy_model.print_input()
    return cloudy_model

# #Run cloudy
Te_interval = [9000.0, 9500.0, 10000.0, 10500.0]
# for Te in Te_interval:
#     model_i = make_vary_Te_ne(Te, 100.0, dir_)
#     pc.log_.timer('Starting Cloudy', quiet = True, calling = 'test1')
#     model_i.run_cloudy()
#     pc.log_.timer('Cloudy ended after seconds:', calling = 'test1')

#Loading simulation results
Ms = pc.load_models(dir_)

#Reading the physical parameters in vector form
Te_vector = np.empty(len(Ms))
ne_vector = np.empty(len(Ms))
emis_dict = dict.fromkeys(emis_tab, np.empty(len(Ms)))

print 'somos ', len(Ms), ' y deberiamos ser', len(Te_interval)

for i in range(len(Ms)):
    ne_vector[i] = Ms[i].ne
    Te_vector[i] = Ms[i].T0
    for j in range(len(emis_dict)):
        line_label = emis_tab[j].replace(' ','_').replace('.','')
        emis_dict[emis_tab[j]][i] = Ms[i].get_emis(line_label)

# idx_sort = np.argsort(Te_vector)
# print Te_vector[idx_sort]
# print emis_dict['HE 1 7065.22A'][idx_sort]
#
# He1 = pn.RecAtom('He',1)
#
# HeI_7064A_emis = He1.getEmissivity(Te_vector[idx_sort], 100.2, label='7065.0')
#
# print HeI_7064A_emis.shape
# print emis_dict['HE 1 7065.22A'][idx_sort].shape
#
# dz.data_plot(Te_vector[idx_sort], emis_dict['HE 1 7065.22A'][idx_sort], 'Cloudy emissivities')
# dz.data_plot(Te_vector[idx_sort], HeI_7064A_emis, 'Pyneb emissivities')
#
#
# dz.FigWording(r'Temperature', 'Emissivity', 'Emissivities comparison')
#
# dz.display_fig()

idx_sort = np.argsort(Te_vector)
print Te_vector[idx_sort]

He1 = pn.RecAtom('He', 1)

HeI_5875_emis = He1.getEmissivity(Te_vector[idx_sort], 100.2, label='5876.0')

print HeI_5875_emis
print emis_dict['HE 1 5875.64A'][idx_sort]


dz.data_plot(Te_vector[idx_sort], emis_dict['HE 1 5875.64A'][idx_sort], 'Cloudy emissivities')
dz.data_plot(Te_vector[idx_sort], HeI_5875_emis, 'Pyneb emissivities')

dz.FigWording(r'Temperature', 'Emissivity', 'Emissivities comparison')

dz.display_fig()


# # Defining the object that will manage the input file for Cloudy
# c_input = pc.CloudyInput(full_model_name)
# c_input.set_other(('init file "hheonly.ini"'))
# c_input.set_other(('element helium abundance -3'))
# c_input.set_BB(Teff=Teff, lumi_unit='q(H)', lumi_value=qH)
# c_input.set_other(('constant temperature 4'))
# c_input.set_cste_density(dens)
# c_input.set_radius(r_in=dist)
# c_input.set_other(('set dr 0'))
# c_input.set_stop('zone 1')
# c_input.set_other(('stop zone 1'))
# c_input.set_emis_tab(emis_tab)  # better use read_emis_file(file) for long list of lines, where file is an external file
# c_input.print_input()
#
#Run cloudy simulation
# c_input.run_cloudy(dir_=dir_)

# title HeI_emissivity_grid
# Blackbody 45000.0
# init file "hheonly.ini"
# element helium abundance -3
# q(H) = 47.0
# radius = 16.0
# set dr 0
# constant temperature 4 vary
# grid range from 10000 to 20000 step 5000 linear
# hden = 2.0 vary
# grid range from 1.0 to 3.0 in 0.5 dex steps
# stop zone 1
# iterate
# save last lines emissivity "lines.emis"
# He 1  3888.63A
# He 1  4026.20A
# He 1  4471.49A
# He 1  5875.64A
# He 1  6678.15A
# He 1  7065.22A
# end of lines

# import numpy as np
# import matplotlib.pyplot as plt
# import pyCloudy as pc
#
# dir_ = '/home/vital/Desktop/tests_cloudy/'
# pc.config.cloudy_exe = '/home/vital/cloudy/source/cloudy.exe'
#
# def make_model(name, radius):
#     Min = pc.CloudyInput('/home/vital/Desktop/tests_cloudy/{}_{}'.format(name, radius))
#     Min.set_BB(Teff=43600, lumi_unit='Q(H)', lumi_value=49.34)
#     Min.set_cste_density(4)
#     Min.set_radius(radius)
#     Min.set_abund(predef='ism', nograins=False)
#     Min.set_other(('Cosmic Rays Background'))
#     Min.print_input()
#
# name = 'M1'
# for radius in np.linspace(13, 23,6):
#     make_model(name, radius)
#
# pc.print_make_file(dir_)
#
# pc.run_cloudy(dir_=dir_, n_proc=6, use_make=True)
#
# Ms = pc.load_models(dir_ + '/M1_', read_emis=False)
#
# for M in Ms:
#     print(M.model_name_s, M.out['Cloudy ends'])
#
#
# def make_varyZ(Z, folder_models):
#     Min = pc.CloudyInput('{}varyZ_{}'.format(folder_models, Z))
#     Min.set_BB(Teff=4e4, lumi_unit='Ionization parameter', lumi_value=-2)
#     Min.set_cste_density(0)
#     Min.set_stop(('zone = 1'))
#     Min.set_emis_tab(('O  3 5006.84A', 'H  1 4861.36A'))
#     Min.set_other(('metals {}'.format(Z),
#                    'set dr 0',
#                    'Cosmic Rays Background'))
#     Min.print_input()
#
# for Z in np.arange(-2, 1.1, 0.25):
#     make_varyZ(Z, dir_)
#
# pc.run_cloudy(dir_=dir_)
#
# Ms = pc.load_models(dir_)
# Ms = sorted(Ms, key=lambda x:x.abund['O'])
#
# print Ms[0].emis_labels
#
# O3Hb = [M.get_emis_vol('O__3_500684A') / M.get_emis_vol('H__1_486136A') for M in Ms]
# OH = [M.abund['O'] for M in Ms]
#
# f, ax = plt.subplots(figsize=(8,6))
# ax.plot(OH, O3Hb)
# ax.set_xlabel('log O/H')
# ax.set_ylabel(r'[OIII]/H$\beta$');
#
# T0 = [M.T0 for M in Ms]
# f, ax = plt.subplots(figsize=(10,8))
# ax.plot(OH, np.log(T0))
# ax.set_xlabel('log O/H')
# ax.set_ylabel('log Te')
#
# plt.show()