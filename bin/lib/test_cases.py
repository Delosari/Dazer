import pyneb as pn
import numpy as np
import pandas as pd
from inferenceModel import SpectraSynthesizer
from Astro_Libraries.spectrum_fitting.gasEmission_functions import TOIII_TSIII_relation

# Declare synthesizer object
specS = SpectraSynthesizer()

# Declare algorithm variables
# #TODO we need a default simulated observation to run of these test
obj_properties_file = 'C:\\Users\\Vital\\PycharmProjects\\thesis_pipeline\\article2_material\\synth_TestObjProperties.txt'
obj_prop_df = pd.read_csv(obj_properties_file, delim_whitespace=True, header=0, index_col=0)

obj_lines_file = 'C:\\Users\\Vital\\PycharmProjects\\thesis_pipeline\\article2_material\\synth_TestObjLines.txt'
obj_lines_df = pd.read_csv(obj_lines_file, delim_whitespace=True, header=0, index_col=0)

specS.obj_data = {}

specS.import_atomic_data()

specS.import_emission_line_data(obj_lines_df, obj_lines_df.index.values)

specS.lineFlambda = specS.gasExtincParams(specS.obj_data['lineWaves'], specS.config['R_v'], specS.config['reddenig_curve'])

specS.gasSamplerVariables(specS.obj_data['lineIons'], specS.config['high_temp_ions'])

# Dictionary with synthetic abundances
abund_df = obj_prop_df[obj_prop_df.index.str.contains('_abund')]
abund_keys = abund_df.index.str.rstrip('_abund').astype(str).values
abund_values = abund_df.variable_magnitude.values
specS.abund_dict = dict(zip(abund_keys, abund_values.T))

# Gas physical parameters
T_low = obj_prop_df.loc['T_low'][0]
n_e = obj_prop_df.loc['n_e'][0]
tau = obj_prop_df.loc['tau'][0]
cHbeta = obj_prop_df.loc['cHbeta'][0]

# Calculate T_high assuming we get T_low
T_high = TOIII_TSIII_relation(T_low)

# Prepare gas data
specS.lineLabels = specS.obj_data['lineLabels']
specS.lineIons = specS.obj_data['lineIons']
specS.linePnCode = specS.obj_data['linePynebCode']

lineFluxes = specS.calcEmFluxes(T_low, T_high, n_e, cHbeta, tau, specS.abund_dict, specS.emFluxTensors, np.zeros(specS.lineLabels.size))
lineFluxesPn = specS.calcEmFluxes_pyneb(T_low, T_high, n_e, cHbeta, tau, specS.abund_dict)

for i in range(lineFluxes.size):
    print specS.lineLabels[i], lineFluxes[i], lineFluxesPn[i]


# Te, ne, cHbeta, tau = 15000.0, 150.0, 0.125, 0.03
#
# specS.calcEmFluxes()


# flux = S3 * emis_R * np.power(10, -cHbeta * flambda) * ftau  + continuum
#
#
#
#
#
