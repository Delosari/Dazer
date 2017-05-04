'''
Created on Sep 24, 2015

@author: vital
'''

import numpy as np
import pyneb as pn


O2 = pn.Atom('O',2)

233

10000

O2_3726A = 176.22 
O2_7319A = 2.33
O2_7330A = 1.94

t = 11697.45096
ne = 135.9073946
normalHbeta = 100

OII_HII = O2.getIonAbundance(int_ratio = O2_3726A, tem = t, den = ne,  to_eval='L(3726)+L(3729)', Hbeta = normalHbeta)
OII_HII_7319A = O2.getIonAbundance(int_ratio = O2_7319A, tem = t, den = ne, wave=7319, Hbeta = normalHbeta)
OII_HII_7330A = O2.getIonAbundance(int_ratio = O2_7330A, tem = t, den = ne, wave=7330, Hbeta = normalHbeta)
OII_HII_2_7330AA = O2.getIonAbundance(int_ratio = O2_7319A + O2_7330A, tem = t, den = ne, to_eval='L(7319)+L(7330)', Hbeta = normalHbeta)

print 'O2 abundance O2_3726A', np.log10(OII_HII)
print 'O2 abundance O2_7319A', np.log10(OII_HII_7319A)
print 'O2 abundance O2_7330A', np.log10(OII_HII_7330A)
print 'O2 abundance O2 7000A lines', np.log10(OII_HII_2_7330AA)