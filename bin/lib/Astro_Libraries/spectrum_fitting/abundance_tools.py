import numpy as np
from pymc3 import Deterministic
from collections import OrderedDict


class ChemicalModel():

    def __init__(self):

        # TODO should we use a dictionary?
        self.H1rCheck = None
        self.Ar3Check, self.Ar4Check= None, None
        self.He1rCheck, self.He2rCheck = None, None
        self.O2Check, self.O3Check = None, None
        self.N2Check = None
        self.S2Check, self.S3Check = None, None

        # Fernandez et al 2018 correction for the S^3+ fraction in the form:
        # log(Ar3/Ar4) = a * log(S3/S4) + b => logS4 =  (a * logS3 - logAr3 + logAr4 + b) / a
        # TODO this should be read from a text file
        self.a_S3corr, self.a_S3corrErr = 1.162, 0.00559
        self.b_S3corr, self.b_S3corrErr = 0.047, 0.0097

        # Ohrs 2016 relation for the OI_SI gradient
        # log(SI/OI) -1.53
        self.OI_SI = 0.029512
        self.OI_SIerr = 0.00339

    def checkIonObservance(self, ion, ionList):
        return True if ion in ionList else False

    def elementalChemicalModel(self, infParamDict, ionList, iterations):

        # Convert to natural scale
        tracesDict = {}
        for ion in ionList:
            if ion in ['He1r', 'He2r']:
                tracesDict[ion] = infParamDict[ion]
            else:
                tracesDict[ion] = np.power(10, infParamDict[ion] - 12)

        # Oxygen abundance
        if self.O2Check and self.O3Check:
           infParamDict['O_abund'] = self.oxygenAbundance(tracesDict)

        # Nitrogen abundance
        if self.N2Check and  self.O2Check and self.O3Check:
           infParamDict['N_abund'] = self.nitrogenAbundance(tracesDict)

        # Sulfur abundance
        if self.S2Check and self.S3Check:
           infParamDict['S_abund'] = self.sulfurAbundance(tracesDict, iterations)
           if 'S4' in tracesDict:
                infParamDict['ICF_SIV'] = infParamDict['S_abund'] / (tracesDict['S2'] + tracesDict['S3'])

        # Helium abundance
        if self.He1rCheck:
           infParamDict['He_abund'] = self.heliumAbundance(tracesDict)

        # Helium mass fraction by oxygen
        if self.He1rCheck and self.O2Check and self.O3Check:
            infParamDict['Ymass_O'] = self.heliumMassFractionOxygen(infParamDict['He_abund'], infParamDict['O_abund'])

        # Helium mass fraction by sulfur
        if self.He1rCheck and self.S2Check and self.S3Check:
            infParamDict['Ymass_S'] = self.heliumMassFractionSulfur(infParamDict['He_abund'], infParamDict['S_abund'], iterations)

        # Convert metal abundances to 12 + log(X^i+) notation
        for metal in ['O_abund', 'N_abund', 'S_abund']:
            if metal in infParamDict:
                infParamDict[metal] = 12 + np.log10(infParamDict[metal])

        return

    def oxygenAbundance(self, abundDict):

        O_abund = abundDict['O2'] + abundDict['O3']

        return O_abund

    def nitrogenAbundance(self, abundDict):

        NO_ratio = abundDict['N2'] / abundDict['O2']

        N_abund = NO_ratio * (abundDict['O2'] + abundDict['O3'])

        return N_abund

    def sulfurAbundance(self, abundDict, iterations):

        if self.Ar3Check and self.Ar4Check:
            aS3corrArray = np.random.normal(self.a_S3corr, self.a_S3corrErr, size=iterations)
            bS3corrArray = np.random.normal(self.b_S3corr, self.b_S3corrErr, size=iterations)

            S4_abund = (aS3corrArray * np.log10(abundDict['S3']) -  np.log10(abundDict['Ar3']) +
                        np.log10(abundDict['Ar4']) + bS3corrArray) / aS3corrArray

            abundDict['S4'] = np.power(10, S4_abund)

            S_abund = abundDict['S2'] + abundDict['S3'] + abundDict['S4']

        else:
            S_abund = abundDict['S2'] + abundDict['S3']

        return S_abund

    def heliumAbundance(self, abundDict):

        if self.He2rCheck:
            He_abund = abundDict['He1r'] + abundDict['He2r']
        else:
            He_abund = abundDict['He1r']

        return He_abund

    def heliumMassFractionOxygen(self, He_abund, O_abund):

        Y_fraction = (4 * He_abund * (1 - 20 * O_abund)) / (1 + 4 * He_abund)

        return Y_fraction

    def heliumMassFractionSulfur(self, He_abund, S_abund, iterations):

        OI_SIArray = np.random.normal(self.OI_SI, self.OI_SIerr, size=iterations)

        Y_fraction = (4 * He_abund * (1 - 20 * self.OI_SI * S_abund)) / (1 + 4 * He_abund)

        return Y_fraction