import numpy as np
from pyneb import atomicData
from import_functions import make_folder
from specSynthesizer_tools import SpectrumFitter

class SpecSynthesizer(SpectrumFitter):

    def __init__(self, atomic_ref=None, nebular_data = None, ssp_type=None, ssp_folder=None, ssp_conf_file=None,
                 temp_grid=None, den_grid=None, high_temp_ions=None, R_v=3.4, reddenig_curve='G03_average', lowlimit_sspContribution=0.0):
