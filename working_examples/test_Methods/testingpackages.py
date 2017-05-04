import matplotlib.pyplot as plt
import pyneb as pn


pn.atomicData.setDataFile('he_i_rec_Pal12-Pal13.fits')

H1  = pn.RecAtom('H', 1)
He1 = pn.RecAtom('He', 1)
He2 = pn.RecAtom('He', 2)
S3 = pn.Atom('S', 3)


#Printing atomic data used in elements
print 'The H1 source file is', H1.printSources()
print 'The He1 source file is', He1.printSources()
print 'The He2 source file is', He2.printSources()
print 'The S3 source file is', S3.printSources()

#New Atom class with new atomic data set
pn.atomicData.setDataFile('s_iii_coll_HRS12.dat')
S3_new = pn.Atom('S', 3)

print 'The S3 source file is', S3_new.printSources()

tem = 10000
den = 100

print S3.printTransition(9069)
print S3.printTransition(9531)
print S3.getEmissivity(tem, den, wave = 9531) / S3.getEmissivity(tem, den, wave = 9069)


S3.plotGrotrian(tem, den)
plt.show()
#Warnings on missing labels. 3820 produces same emissivity for He1 and He2 atom

try:
    H13889A_Emis = H1.getEmissivity(tem, den, wave=3889)
except:
    H13889A_Emis = None

try:
    He13820A_Emis = He1.getEmissivity(tem, den, wave=3820.0)
except:
    He13820A_Emis = None

try:
    He23820A_Emis = He2.getEmissivity(tem, den, wave=3820.0)
except:
    He23820A_Emis = None


print 'H1 3889A emissivity',  H13889A_Emis
print 'He1 3820A emissivity',  He13820A_Emis
print 'He2 3820A emissivity',  He23820A_Emis



