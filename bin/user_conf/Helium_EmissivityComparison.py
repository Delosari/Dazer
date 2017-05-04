from numpy                          import array, linspace, zeros, pi

from CodeTools.PlottingManager      import myPickle
import pyneb                        as pn


S4             = pn.RecAtom('He',1)

He1_Emis_pyneb  = S4.getEmissivity(tem=10000.0, den=100.0, label = '3889.0')

He_Emis_article = 1.4005 * 1e-25 #Units: 4 Pi j / n_e / n_He+ (erg cm^3 s^-1).

#Emissivity coefficient
print 'Emissivity ratio', He1_Emis_pyneb / He_Emis_article

#Load script support classes
pv  = myPickle()
 
#Declare Figure format
pv.FigFormat_One(ColorConf='Night1')    
 
#Load pyneb atom
H1  = pn.RecAtom('H',1)
S4 = pn.RecAtom('He',1)
Helium_Labels = ['3889.0',        '4026.0',       '4471.0',       '5876.0',       '6678.0',       '7065.0']
 
#Declare initial conditions
Te_vector   = linspace(8000, 25000, 100)    #K
ne_vector   = linspace(10, 1000, 100)       #cm^-3
 
Te_0    = 10000
ne_0    = 100
Hbeta_0_Emis = H1.getEmissivity(tem=Te_0, den=ne_0, label='4_2') 
 
Normalizing_Constant = 1e-25 #erg cm^3 s^-1
 
#----Pyneb-----------------------------------------------
#Emissivity vector
HeEm_vector = zeros((len(Helium_Labels),len(ne_vector)))
 
#Calculate emissivities
for i in range(len(Helium_Labels)): 
    for j in range(len(Te_vector)):
        HeEm_vector[i][j] = S4.getEmissivity(tem=Te_vector[j], den=ne_0, label = Helium_Labels[i])
  
#Plot the data
for k in range(len(Helium_Labels)):
    Label = 'HeI ' + Helium_Labels[k] + r'$\AA$' + ' emissivity'
    pv.DataPloter_One(Te_vector, HeEm_vector[k]/1e-25, Label, pv.Color_Vector[2][k])
 
 
#-----PFSFD Emissivities from 2012 paper-----------------------------------------------
TableAddress = '/home/vital/git/Dazer_Local/Dazer/Astro_Libraries/PFSD_HeEmissivities_ne100'
Temps, He3889, He40226, He4471, He5876, He6678, He7065 = pv.get_TableColumn((0,2,3,4,5,6,7), TableAddress, HeaderSize=1, StringIndexes=False, unpack_check=True)
Emissivities_Table = [He3889, He40226, He4471, He5876, He6678, He7065]
 
print 'Helium emissivity ', 
 
#Plot the data
Conversionparameter = 1
 
#In the table the data is in units = 4 * pi * j / (ne * nHe+)
for k in range(len(Emissivities_Table)):
    Label = 'PFSD emissivity'
    Emissivity = Emissivities_Table[k] * Conversionparameter
    pv.DataPloter_One(Temps, Emissivity, Label, pv.Color_Vector[2][k], LineStyle=None)
 
# for k in range(len(Emissivities_Table)):
#     Label = 'PFSD emissivity'
#     pyneb_Counterpart = zeros(len(Temps))
#     for i in range(len(Temps)):
#         pyneb_Counterpart[i] = S4.getEmissivity(tem=Temps[i], den=ne_0, label = Helium_Labels[k])
#     print 'Linea', Helium_Labels[k]
#     print  Emissivities_Table[k]* 10e-25/ pyneb_Counterpart
#     Emissivity = Emissivities_Table[k]* 10e-25/ pyneb_Counterpart
#     pv.DataPloter_One(Temps, Emissivity, Label, pv.Color_Vector[2][k], LineStyle=None)
 
 
PlotTitle   = 'Helium lines emissivity evolution with Temperature'
x_label     = 'Temperature ' + r'$(K)$'
y_label     = r'$E_{He}\left(\lambda,\,T_{e},\,n_{e}=100\,cm^{-3}\right)_{He}\,/ 10^-25 ()$'
 
pv.Labels_Legends_One(PlotTitle, x_label, y_label)    
 
pv.DisplayFigure()
 
print 'Treatment completed'