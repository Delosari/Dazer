import pyneb as pn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from inspect import getargspec
from collections import OrderedDict
from scipy.optimize import curve_fit
from inferenceModel import SpectraSynthesizer
from Astro_Libraries.spectrum_fitting.gasEmission_functions import EmissionComponents
from lib.Astro_Libraries.spectrum_fitting.plot_tools import MCMC_printer

# Tools for printing
mcPrinter = MCMC_printer()
outputFolder = 'E:\\Research\\article_YpBayesian\\'

# Generate the dictionary with pyneb ions
ionDict = pn.getAtomDict(atom_list=['H1r','He1r'])
pynebcode_list = np.array(['4471', '5876', '6678'])
ions_list = np.array(['He1r', 'He1r', 'He1r'])
labels_list = np.array(['He1_4471A', 'He1_5876A', 'He1_6678A'])
epm2017_emisCoeffs = {'He1_3889A': np.array([0.173, 0.00054, 0.904, 1e-5]),
                      'He1_4026A': np.array([-0.09, 0.0000063, 4.297, 1e-5]),
                      'He1_4471A': np.array([-0.1463, 0.0005, 2.0301, 1.5e-5]),
                      'He1_5876A': np.array([-0.226, 0.0011, 0.745, -5.1e-5]),
                      'He1_6678A': np.array([-0.2355, 0.0016, 2.612, 0.000146]),
                      'He1_7065A': np.array([0.368, 0.0017, 4.329, 0.0024]),
                      'He1_10830A': np.array([0.14, 0.00189, 0.337, -0.00027])}

porter2007_emisCoeffs = {'He1_4471A': np.array([-3.5209E+04, -4.5168E+02, 8.5367E+03, 9.1635E+03]),
                      'He1_5876A': np.array([2.0620E+05, 1.7479E+02, -1.3548E+04, -7.3492E+05]),
                      'He1_6678A': np.array([6.7315E+04, 7.5157E+01, -4.7101E+03, -2.3610E+05])}

# Temperature and density grid
temp_grid = np.array([9000, 20000, 251])
den_grid = np.array([1,600,101])
tempRange = np.linspace(temp_grid[0], temp_grid[1], temp_grid[2])
denRange = np.linspace(den_grid[0], den_grid[1], den_grid[2])
X, Y = np.meshgrid(tempRange, denRange)
tempGridFlatten, denGridFlatten = X.flatten(), Y.flatten()
te_ne_grid = (tempGridFlatten, denGridFlatten)

# Read lines file
obj_lines_file = 'C:\\Users\\Vital\\PycharmProjects\\thesis_pipeline\\article2_material\\synth_TestObjLines.txt'
obj_lines_df = pd.read_csv(obj_lines_file, delim_whitespace=True, header=0, index_col=0)
idcs_helium = obj_lines_df.index.isin(['He1_4471A', 'He1_5876A', 'He1_6678A'])
obj_lines_df = obj_lines_df[idcs_helium]

def emisEquation_HI(xy_space, a, b, c):
    temp_range, den_range = xy_space
    return a + b * np.log10(temp_range) + c * np.log10(temp_range) * np.log10(temp_range)

def emisEquation_HeI_fit(xy_space, a, b, c, d):
    temp_range, den_range = xy_space
    return (a + b * den_range) * np.log10(temp_range / 10000.0) - np.log10(c + d * den_range)

def emisEquation_HeII_fit(xy_space, a, b):
    temp_range, den_range = xy_space
    return a + b * np.log(temp_range / 10000)

def emisEquation_HeI_fitPor(xy_space, a, b, c, d):
    temp_range, den_range = xy_space
    return (a + b * np.log(tempRange) * np.log(tempRange) + c * np.log(tempRange) + d / np.log(tempRange)) * 1 / tempRange * 1.0e-25

def fitEmis(func_emis, xy_space, line_emis, p0=None):
    p1, p1_cov = curve_fit(func_emis, xy_space, line_emis, p0)
    return p1, p1_cov

def computeEmissivityDict(ionDict, pynebcode_list, ions_list, labels_list, tempGridFlatten, denGridFlatten):

    # Empty grid array
    emisDict = OrderedDict()
    emisDict['H1_4861A'] = ionDict['H1r'].getEmissivity(tempGridFlatten, denGridFlatten, wave=4861, product=False)

    for i in range(len(labels_list)):

        # Line emissivity references
        line_label = labels_list[i]

        # Line emissivity references
        emis_grid_i = ionDict[ions_list[i]].getEmissivity(tempGridFlatten, denGridFlatten, wave=float(pynebcode_list[i]), product=False)

        # Save along the number of points
        emisDict[line_label] = emis_grid_i

    return emisDict

def fitEmissivityPlane(labels_list, ionEmisEq_fit, emis_dict, tempGridFlatten, denGridFlatten, epm2017_emisCoeffs):

    # Dictionary to store the emissivity surface coeffients
    emisCoeffs = OrderedDict()
    for i in range(labels_list.size):

        lineLabel = labels_list[i]

        # Get equation type to fit the emissivity
        line_func  = ionEmisEq_fit[lineLabel]
        n_args = len(getargspec(line_func).args) - 2 #TODO Not working in python 2.7 https://stackoverflow.com/questions/847936/how-can-i-find-the-number-of-arguments-of-a-python-function

        # Compute emissivity functions coefficients
        emis_grid_i = emis_dict[lineLabel]
        p0 = epm2017_emisCoeffs[lineLabel] if lineLabel in epm2017_emisCoeffs else np.zeros(n_args)
        p1, cov1 = fitEmis(line_func, (tempGridFlatten, denGridFlatten), emis_grid_i, p0=p0)
        emisCoeffs[lineLabel] = p1

    return emisCoeffs

ionEmisEq_fit = {'H1_4102A': emisEquation_HI,
                  'H1_4341A': emisEquation_HI,
                  'H1_6563A': emisEquation_HI,
                  'He1_4471A': emisEquation_HeI_fit,
                  'He1_5876A': emisEquation_HeI_fit,
                  'He1_6678A': emisEquation_HeI_fit}

ionEmisEqPorter_fit = {'He1_4471A': emisEquation_HeI_fitPor,
                  'He1_5876A': emisEquation_HeI_fitPor,
                  'He1_6678A': emisEquation_HeI_fitPor}


# Compute emissivity grids
emisDict = computeEmissivityDict(ionDict, pynebcode_list, ions_list, labels_list, tempGridFlatten, denGridFlatten)

# Normalized log version
emisDict_log = {}
for label in labels_list:
    emisDict_log[label] = np.log10(emisDict[label] / emisDict['H1_4861A'])

# Fit equation original
emisCoeffs_log = fitEmissivityPlane(labels_list, ionEmisEq_fit, emisDict_log, tempGridFlatten, denGridFlatten, epm2017_emisCoeffs)

# Fit equation Porter
emisCoeffs_Porter = fitEmissivityPlane(labels_list, ionEmisEq_fit, emisDict_log, tempGridFlatten, denGridFlatten, porter2007_emisCoeffs)



# Plot the graphs
for label in labels_list:
    print 'Printing {}'.format(label)
    mcPrinter.emissivitySurfaceFit_2D(label, emisCoeffs_log[label], emisDict_log[label],
                             ionEmisEq_fit[label], te_ne_grid, denRange, tempRange)
    outputFile = '{}emissivityTeDe2DSecond{}_den{}-{}_temp{}-{}'.format(outputFolder, label, denRange[0], denRange[-1], tempRange[0], tempRange[-1])
    mcPrinter.savefig(outputFile, resolution=200)
    plt.clf()
