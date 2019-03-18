from dazer_methods import Dazer
from lib.inferenceModel import SpectraSynthesizer

# Import functions
dz = Dazer()
specS = SpectraSynthesizer()

# Object data to prepare a synthetic observation
synth_data ={'obs_name'                 :'synthTestHIIgalaxy',
             'output_folder'            :'E:\\Research\\article_YpBayesian\\', #'/home/vital/article_YpBayesian/',#
             'obj_properties_file'      :'C:\\Users\\Vital\\PycharmProjects\\thesis_pipeline\\article2_material\\synth_TestObjProperties.txt',#'/home/vital/PycharmProjects/thesis_pipeline/article2_material/synth_TestObjProperties.txt',
             'obj_lines_file'           :'C:\\Users\\Vital\\PycharmProjects\\thesis_pipeline\\article2_material\\synth_TestObjLines.txt',#'/home/vital/PycharmProjects/thesis_pipeline/article2_material/synth_TestObjLines.txt',
             'wavelengh_limits'         :[4000, 6900],
             'resample_inc'             :1,
             'norm_interval'            :[5100, 5150],
             'ssp_lib_type'             :'starlight',
             'ssp_folder'               :'E:\\Research\\Starlight\\Bases\\',#'/home/vital/Starlight/Bases/', #
             'ssp_file'                 :'E:\\Research\\Starlight\\Bases\\Dani_Bases_Extra_short.txt',#'/home/vital/Starlight/Bases/Dani_Bases_Extra_short.txt',
             'obj_ssp_coeffs_file'      :'C:\\Users\\Vital\\PycharmProjects\\thesis_pipeline\\article2_material\\synth_StellarPop.txt',#'/home/vital/PycharmProjects/thesis_pipeline/article2_material/synth_StellarPop.txt', # ,
             'error_stellarContinuum'   :0.01,
             'error_lines'              :0.02,
             'atomic_data'              :None,
             'ftau_coeffs'              :None}

# Load atomic data
specS.import_atomic_data(atomic_data, ftau_coeffs, specS.config['temp_grid'], specS.config['den_grid'])

# Prepare data from emission line file (trick to put all the lines)
specS.import_emission_line_data(obj_lines_df, obj_lines_df.index.values)

# Reddening parameters
specS.obj_data['lineFlambda'] = specS.gasExtincParams(specS.obj_data['lineWaves'], specS.config['R_v'],
                                                    specS.config['reddenig_curve'])

# Create or load emissivity grids
emis_dict = specS.computeEmissivityDict(specS.obj_data['linePynebCode'], specS.obj_data['lineIons'],
                                            specS.obj_data['lineLabels'], output_folder)