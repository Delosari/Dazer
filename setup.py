from distutils.core import setup
setup(
  name = 'dazer',
  packages = ['dazer'], # this must be the same as the name above
  version = '0.1.3',
  description = 'Code for measuring emission lines in HII regions using a MC algorithm and tkinter interface',
  author = 'Vital Fernandez',
  author_email = 'vital.fernandez@gmail.com',
  url = 'https://github.com/Delosari/Dazer.git',
  download_url = 'https://github.com/Delosari/Dazer/tarball/0.1.3', 
  keywords = ['astrophysics', 'emission lines', 'MC'],
  classifiers = [], requires=['numpy', 'pandas', 'pymc3', 'theano']
)