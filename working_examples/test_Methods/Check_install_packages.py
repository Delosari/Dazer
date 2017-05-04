'''
Created on Aug 3, 2015

@author: vital
'''

print '--Starting'

#Approach 1
# help('modules')

#Approach 2
# import sys
# import os
# for p in sys.path:
#     print os.listdir( p )

#Approach 3
import pip


installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
for j in range(len(installed_packages_list)):
    print installed_packages_list[j]

#List of python libra ries



# #Install pip and main coding modules
# sudo apt-get install python-pip python-dev build-essential
# 
# #Install Easy install
# wget https://bootstrap.pypa.io/ez_setup.py -O - | sudo python
# 
# #Install scientific libraries
# #--Numpy
# sudo pip install numpy
# #--g fortran based libraries
# sudo apt-get install libatlas-base-dev gfortran
# #--scipy
# sudo pip install scipy
# #--matplotlib
# sudo pip install matplotlib
# #--pandas
# sudo pip install pandas
# #--ipython
# sudo pip install ipython 
# #--uncertaintines package
# sudo pip install uncertainties
# 
# #Install astronomical libraries
# #--pyfits
# sudo pip install pyfits
# #--pyraf
# sudo pip install pyraf
# #--pyneb
# sudo pip install pyneb
# #--astropy
# sudo pip install --no-deps astropy
# 
# #--Kmpc package: Download this file
# https://www.astro.rug.nl/software/kapteyn-beta/kapteyn-2.3b6.tar.gz
# #--Or the latest one which appears on this website:
# https://www.astro.rug.nl/software/kapteyn-beta/index.html
# #--You need to extract it and from a terminal go into this extracted folder and run:
# python setup.py install


print '-- Procedure tested'