from collections                    import OrderedDict, Sequence
from mimetypes                      import guess_type
from operator                       import itemgetter
from os                             import getcwd, walk, makedirs, mkdir, chdir, path
from os.path                        import isdir, isfile
from shutil                         import move as shutil_move, copyfile
from socket                         import gethostname
from subprocess                     import Popen, PIPE, STDOUT
from sys                            import argv, exit, stdout
from numpy                          import isnan, percentile, loadtxt, savetxt, where, zeros, searchsorted, ndarray, in1d, array, transpose, empty, nan as np_nan, float64 as np_float64, float32 as np_float32, full
from pandas                         import read_csv, DataFrame, read_pickle, ExcelFile, concat, ExcelWriter
#from pyfits                         import Column, ColDefs, open as pyfits_open, TableHDU, getdata
from astropy.io.fits                import Column, ColDefs, open as pyfits_open, TableHDU, getdata
from pylatex                        import Document, Figure, NewPage, NoEscape, Package, Tabular, Section, Tabu, Table, LongTable
from scipy                          import linspace
from scipy.interpolate              import interp1d
from uncertainties                  import UFloat, ufloat
from uncertainties.umath            import log10 as umath_log10, pow as unumath_pow
from uncertainties.unumpy           import uarray, nominal_values, std_devs, log10 as unum_log10, pow as unnumpy_pow
from astropy.io                     import fits
from string                         import ascii_uppercase
from pandas                         import notnull
from functools                      import partial
from sigfig                         import round_sig

class Images_Fits():
    
    def __init__(self):
        self.Opening_Procedure = None
        
    def Fits_to_Data(self, FolderName, FileName):
        
        #Get the fits file main data and header
        Data, Header_0 = getdata(FolderName + FileName, header=True)
        
        #In the case a fits file I made
        if ('WHTJOINW' in Header_0) or ("STALIGHT" in Header_0) or ("NEBUSPEC" in Header_0):
            x = Data['Wave']
            y = Data['Int']
    
#             FitsFile.close()
            
            return x, y, [Header_0]
        
        #In the case a dr10 file: Spectra redshifed and in absolute units 
        elif ("COEFF0" in Header_0 and "dr10" in FileName):
            FitsFile = pyfits_open(FolderName + FileName) 
            Spectra = FitsFile[1].data
            Header_2  = FitsFile[2].data
            Header_3 = FitsFile[3].data
            
            Int = Spectra['flux']
            Int_abs = Int / 1e17
            
            WavelenRange = 10.0**Spectra['loglam']
            SDSS_z = float(Header_2["z"][0] + 1)
            Wavelength_z = WavelenRange / SDSS_z
    
            Headers = (Header_0,Header_2,Header_3)
    
            FitsFile.close()
            
            return Wavelength_z, Int_abs, Headers
        
        #Any other fits file (This is a very old scheme)
        else:                        
            if Data.ndim == 1: 
                Int = Data
            else:
                Int = Data[0]
                
            if "COEFF0" in Header_0:
                dw              = 10.0**Header_0['COEFF1']                          # dw = 0.862936 INDEF (Wavelength interval per pixel)
                Wmin            = 10.0**Header_0['COEFF0']
                pixels          = Header_0['NAXIS1']                                # nw = 3801 number of output pixels
                Wmax            = Wmin + dw * pixels
                WavelenRange    = linspace(Wmin,Wmax,pixels,endpoint=False)
                     
                return WavelenRange, Int, [Header_0]
            
            elif "LTV1" in Header_0:
                StartingPix     = -1 * Header_0['LTV1']                   # LTV1 = -261. 
                Wmin_CCD        = Header_0['CRVAL1']
                dw              = Header_0['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
                pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
                Wmin            = Wmin_CCD + dw * StartingPix
                Wmax            = Wmin + dw * pixels
                WavelenRange    = linspace(Wmin,Wmax,pixels,endpoint=False)
                return WavelenRange, Int, [Header_0]
            
            else:
                Wmin            = Header_0['CRVAL1']
                dw              = Header_0['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
                pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
                Wmax            = Wmin + dw * pixels
                WavelenRange    = linspace(Wmin,Wmax,pixels,endpoint=False)
                return WavelenRange, Int, [Header_0]
     
    def get_spectra_data(self, file_address, ext=0, force_float64 = True):

        data_array, Header_0 = fits.getdata(file_address, header=True)    #This is the issue with the numpys created by myself
               
        #In the case a fits file I made
        if ('WHTJOINW' in Header_0) or ("STALIGHT" in Header_0) or ("NEBUSPEC" in Header_0):
            wavelength = data_array['Wave']
            Flux_array = data_array['Int']
                
        elif "COEFF0" in Header_0:
            dw              = 10.0**Header_0['COEFF1']                # dw = 0.862936 INDEF (Wavelength interval per pixel)
            Wmin            = 10.0**Header_0['COEFF0']
            pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
            Wmax            = Wmin + dw * pixels
            wavelength      = linspace(Wmin,Wmax,pixels,endpoint=False)
            Flux_array      = data_array
 
              
        elif "LTV1" in Header_0:
            StartingPix     = -1 * Header_0['LTV1']                   # LTV1 = -261. 
            Wmin_CCD        = Header_0['CRVAL1']
            dw              = Header_0['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
            pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
            Wmin            = Wmin_CCD + dw * StartingPix
            Wmax            = Wmin + dw * pixels
            wavelength      = linspace(Wmin,Wmax,pixels,endpoint=False)
            Flux_array      = data_array
        
        else:
            Wmin            = Header_0['CRVAL1']
            dw              = Header_0['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
            pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
            Wmax            = Wmin + dw * pixels
            wavelength      = linspace(Wmin,Wmax,pixels,endpoint=False)
            Flux_array      = data_array
  
        if force_float64:
            if isinstance(Flux_array[0], np_float32):
                Flux_array = Flux_array.astype(np_float64)
                wavelength = wavelength.astype(np_float64)
                
        return wavelength, Flux_array, Header_0    
    
    def getHeaderEntry(self, Entry, FitAddress, FitsType = None):
        
        FitsFile = pyfits_open(FitAddress) 
        Header_0 = FitsFile[0].header
        FitsFile.close()
           
        EntryValue = Header_0[Entry]
        
        return EntryValue
    
    def getArmColor(self, FitAddress, FitsType = 'WHT'):
        
        FitsFile = pyfits_open(FitAddress) 
        Header_0 = FitsFile[0].header
        FitsFile.close()
        
        Color = None
        
        if FitsType == 'WHT':
            Entry = 'ISIARM'
            EntryValue = Header_0[Entry]
        
            if EntryValue == 'Blue arm':
                Color = 'Blue'
            elif EntryValue == 'Red arm':
                Color = 'Red'
        
        return Color
    
    def Data_2_Fits(self, FileFolder, FitsName, Header, Wavelength, Intensity, NewKeyWord = None):
                
        Header[NewKeyWord[0]] = NewKeyWord[1]
        Column1     = fits.Column(name='Wave', format='E',array=Wavelength)
        Column2     = fits.Column(name='Int', format='E', array=Intensity)
        Columns     = fits.ColDefs([Column1, Column2])
        Table_HDU   = fits.TableHDU.from_columns(Columns, header=Header)
                
        Table_HDU.writeto(FileFolder + FitsName, overwrite = True)
        
        return

    def subSpectrum(self, Wave, Flux, Wlow, Whigh):
        
        indmin, indmax      = searchsorted(Wave, (Wlow, Whigh))
        #indmax              = min(len(Wave)-1, indmax) #WHAT IS THIS FOR???
        
        subWave   = Wave[indmin:indmax]
        subFlux   = Flux[indmin:indmax]
        
        return subWave, subFlux

class Tables_Txt():

    def __init__(self):
        
        self.Current_TableAddress           = None
        self.Current_HeaderSize             = None
        self.Current_HeaderKeys             = None
        
        self.Default_HeaderSize             = 1
        self.Default_ColumnDelimeter        = None

    def select_Table(self, TableAddress, HeaderSize = None, Delimeter = None, loadheaders_check = False):
        
        #In this method we define the table we are going to work with
        self.Current_TableAddress       = TableAddress
        
        #Defint header size... This is not as clear as it should
        if HeaderSize != None:
            self.Current_HeaderSize     = HeaderSize
                    
        #Define the separation criteria between columns
        if Delimeter == None:
            Delimeter = self.Default_ColumnDelimeter
        
        #Load the headers from the table
        if loadheaders_check == True:
            self.get_Headers_FromTable(self.Current_TableAddress, self.Current_HeaderSize, Delimeter)
        
        return
        
    def get_Headers_FromTable(self, TableAddress = None, TableHeaderSize = None, Delimeter = None):
        #WARNING: NOT VERY EFFICIENT AND NOT SURE IF I SHOULD DELETE THE DATA AFTER USING IT
        
        #Use default or declared delimiter for columns
        if Delimeter == None:
            Delimeter               = self.Default_ColumnDelimeter
        
        #Read the text file
        TextFile                    = open(TableAddress, 'r')
        TextFile_Lines              = TextFile.readlines()
        TextFile.close()
        
        #Import the headers (This assumes the header is the row just before the begining of the columns
        self.Current_HeaderKeys     = TextFile_Lines[TableHeaderSize - 1].split(Delimeter)
        
        return
                  
    def get_ColumnData(self, Columns, TableAddress = None, HeaderSize = None, StringIndexes = True, datatype = float, comments_icon = '#', unpack_check = True):
        #WARNING: Even if Columns is just one string or int, it must be in a row
        if TableAddress == None:
            TableAddress    = self.Current_TableAddress
        
        #In case the only rows to skip are the header
        if HeaderSize == None:
            HeaderSize      = self.Default_HeaderSize
                
        #This structure makes sure you can input either indexes or strings    
        if StringIndexes == True:
            self.get_Headers_FromTable(TableAddress, HeaderSize)
            List_Indexes = zeros(len(Columns)).tolist()
            for i in range(len(Columns)):
                List_Indexes[i] = int(self.Current_HeaderKeys.index(Columns[i]))
                
        else:
            List_Indexes = Columns
                
        #Import the data, just a single column
        if len(List_Indexes) == 1:
            return loadtxt(TableAddress, dtype = datatype, comments = comments_icon, skiprows = HeaderSize, usecols = List_Indexes)
        
        #Import data several columns with the same type
        else:
            return loadtxt(TableAddress, dtype = datatype, comments = comments_icon, skiprows = HeaderSize, usecols = (List_Indexes), unpack = unpack_check)
    
    def get_TableColumn(self, Columns, TableAddress, HeaderSize = None, StringIndexes = True, datatype = float, comments_icon = '#', Delimeter = None, unpack_check = True):
        
        if (StringIndexes) and (HeaderSize == None):
            HeaderSize = 1
        
        elif (HeaderSize== None) and (StringIndexes == False):
            exit('WARNING: Table does not have header: ' + TableAddress)

        #This structure makes sure you can input either indexes or strings    
        if StringIndexes == True:
            List_Indexes = zeros(len(Columns)).tolist()
            HeaderKeys = self.get_TableHeader(TableAddress, HeaderSize, Delimeter = Delimeter)
            
            for i in range(len(Columns)):
                List_Indexes[i] = int(HeaderKeys.index(Columns[i]))
                
        else:
            List_Indexes = Columns
            
        #Import the data, just a single column
        if len(List_Indexes) == 1:
            return loadtxt(TableAddress, dtype = datatype, comments = comments_icon, skiprows = HeaderSize, usecols = List_Indexes)
        
        #Import data several columns with the same type
        else:
            return loadtxt(TableAddress, dtype = datatype, comments = comments_icon, skiprows = HeaderSize, usecols = (List_Indexes), unpack = unpack_check)
           
        return
    
    def get_TableHeader(self, TableAddress, TableHeaderSize, Delimeter):
        #WARNING: NOT VERY EFFICIENT AND NOT SUER IF I SHOULD DELETE THE DATA AFTER USING IT
                
        #Read the text file
        TextFile                    = open(TableAddress, 'r')
        TextFile_Lines              = TextFile.readlines()
        TextFile.close()
        
        #Import the headers (This assumes the header is the row just before the begining of the columns
        Current_HeaderKeys     = TextFile_Lines[TableHeaderSize - 1].split(Delimeter)
        
        return Current_HeaderKeys

class table_formats():
    
    def __init__(self):
        
        self.header_dict = OrderedDict()
    
    def Catalogue_headers(self):
        
        #This dictionary stores the format of the headers in latex #It also serves as a list from the defined entry keys
        self.header_dict             = OrderedDict()
         
        #Object header
        self.header_dict['object']          = 'Object ID'
         
        #Emision lines ratios
        self.header_dict['O3_ratio']        = r'$\frac{\left[OIII\right]5007\AA}{\left[OIII\right]4959\AA}$'
        self.header_dict['N2_ratio']        = r'$\frac{\left[NII\right]6584\AA}{\left[NII\right]6548\AA}$'
        self.header_dict['S3_ratio']        = r'$\frac{\left[SIII\right]9531\AA}{\left[SIII\right]9069\AA}$'
     
        #Physical Parameters
        self.header_dict['TOIII_pn']           = r'$T_{\left[OIII\right]}$'
        self.header_dict['TSIII_pn']           = r'$T_{\left[SIII\right]}$'
        self.header_dict['TSII_pn']            = r'$T_{\left[SII\right]}$'
        self.header_dict['nSII_pn']            = r'$n_{e,\,\left[SII\right]}$'
        self.header_dict['cHBeta_red']          = r'$c\left(H\beta\right)$'
     
        #Abundances
        self.header_dict['OI_HI_pn']           = r'$12+log\left(\frac{O}{H}\right)$'
        self.header_dict['NI_HI_pn']           = r'$12+log\left(\frac{N}{H}\right)$'
#         self.header_dict['SI_HI_pn']           = r'$12+log\left(\frac{S}{H}\right)$'
        self.header_dict['SI_HI_ArCorr_pn']    = r'$12+log\left(\frac{S}{H}\right)_{Ar}$'
         
        self.header_dict['HeI_HI_pn']          = r'$\frac{He}{H}$'
        self.header_dict['HeII_HII_pn']        = r'$y^{+}$'
        self.header_dict['HeIII_HII_pn']       = r'$y^{++}$'
          
        self.header_dict['Y_Mass_O_pn']        = r'$Y_{\left(\frac{O}{H}\right)}$'
        self.header_dict['Y_Mass_S_pn']        = r'$Y_{\left(\frac{S}{H}\right)}$'
      
#         self.header_dict['Y_Inference_O_pn']   = r'$Y_{\left(\frac{O}{H}\right),\,inf}$'
#         self.header_dict['Y_Inference_S_pn']   = r'$Y_{\left(\frac{S}{H}\right),\,inf}$'
        
        #Physical Parameters
#         self.header_dict['y_plus_inf']      = r'$\left(\frac{HeI}{HI}\right)_{inf}$'
#         self.header_dict['Te_inf']          = r'$T_{e,\,inf}$'
#         self.header_dict['ne_inf']          = r'$n_{e,\,inf}$'
#         self.header_dict['cHbeta_inf']      = r'$c\left(H\beta\right)_{inf}$'
#         self.header_dict['ftau_inf']        = r'$\tau_{\inf}$'
#         self.header_dict['Xi_inf']          = r'$\xi_{inf}$'
    
        self.table_format = 'l' + ''.join([' c' for s in range(len(self.header_dict) - 1)])
    
        return
    
    def RedShifted_linesog_header(self):
        
        self.header_dict                    = OrderedDict()        
        self.header_dict['Object']          = 'Emission line'
        self.header_dict['Redshift']        = r'$z$'
        self.header_dict['Eqw']             = r'$Eqw(H\beta)$ $(\AA)$'
        self.header_dict['OIII4363']        = r'$[OIII]4363(\AA)$'
        self.header_dict['OIII5007']        = r'$[OIII]5007(\AA)$'
        self.header_dict['Halpha']          = r'$H\alpha6563(\AA)$'
        self.header_dict['HeI6678']         = r'$HeI6678\AA$'
        self.header_dict['SIIII9531']       = r'$[SIII]9531(\AA)$'        
        
        self.table_format = 'l' + ''.join([' c' for s in range(len(self.header_dict) - 1)])
        
        return
    
    def EmissionLinesLog_header(self):
        
        #This dictionary stores the format of the headers in latex #It also serves as a list from the defined entry keys
        self.header_dict                    = OrderedDict()        
        self.header_dict['Emission']        = 'Emission line'
        self.header_dict['f_lambda']         = r'$f(\lambda)$'
        self.header_dict['Eqw']             = r'$-EW(\AA)$'
        self.header_dict['Flux_undim']      = r'$F(\lambda)$'
        self.header_dict['Int_undim']       = r'$I(\lambda)$'
#         self.header_dict['Istellar']        = r'$I_{Stellar}(\lambda)$'
        
        self.table_format = 'l' + ''.join([' c' for s in range(len(self.header_dict) - 1)])
        
        return

    def galaxylog_v2_Total(self, thekey, Ratios_dict, RatiosColor_dict, variables_dict, object_column):

        #Treatment for the line ratios
        if thekey in ['O3_ratio', 'N2_ratio', 'S3_ratio']:
            value           = self.format_for_table(Ratios_dict[thekey])
            color           = RatiosColor_dict[thekey]
            cell            = r'\textcolor{{{color}}}{{{value}}}'.format(color = color, value = self.format_for_table(value, 3))
         
        #Treatment for the physical conditions
        elif thekey in ['TOIII_pn', 'TOII_pn', 'TSIII_pn', 'TSII_pn', 'nSII_pn', 'cHBeta_red']:

            value = object_column[thekey]
            
            if isinstance(value, UFloat):
                cell = self.format_for_table(value, 3)
            elif isnan(value):
                cell = None
            else:
                cell = None

        #Treatment for metallic abundances
        elif thekey in ['OI_HI_pn', 'NI_HI_pn', 'SI_HI_ArCorr_pn']:
            
            value = object_column[thekey]
            
            if isinstance(value, UFloat):
                abund_log   = 12 + umath_log10(value)
                cell        = self.format_for_table(abund_log)
            elif isnan(value):
                cell        = None
            else:
                cell        = None
                
        #Treatment for Helium abundances
        elif thekey in ['HeI_HI_pn', 'HeII_HII_pn', 'HeIII_HII_pn', 'Y_Mass_O_pn', 'Y_Mass_S_pn']:
            
            value = object_column[thekey]
            
            if isinstance(value, UFloat):
                cell = self.format_for_table(value, 3)
            elif isnan(value):
                cell = None
            else:
                cell = None


#         #Treatment for metallic abundances
#         elif thekey in ['OI_HI', 'NI_HI', 'SI_HI', 'SI_HI_ArCorr']:
#             value           = self.GetParameter_ObjLog(CodeName, FileFolder, thekey, 'float')
#             
#             if value != None:
#                 abund_log   = 12 + umath_log10(value)
#                 cell        = self.format_for_table(abund_log)
#             else:
#                 cell        = None
#             

#         
#         #Treatment for the inference parameters
#         elif thekey in ['y_plus_inf','Te_inf','ne_inf','cHbeta_inf','ftau_inf','Xi_inf']:
#             value    = self.GetParameter_ObjLog(CodeName, FileFolder, thekey, 'float')
#         
#             if value != None:
#                 error_type  = thekey[0:thekey.find('_inf')] + '_SD'
#                 error       = self.GetParameter_ObjLog(CodeName, FileFolder, error_type, 'float')
#                 value       = ufloat(value,error)
#                 color       = self.compare_variables(thekey, value, variables_dict, CodeName, FileFolder)
#                 cell        = r'\textcolor{{{color}}}{{{value}}}'.format(color = color, value = self.format_for_table(value, 3))
#             else:
#                 cell        = None
     
        return cell 

class Tables_Latex(table_formats):
    
    def __init__(self):
        
        table_formats.__init__(self)
        
    def color_evaluator(self, obs_value, theo_value, evaluation_pattern=[10.0, 25.0], evaluation_colors = ['ForestGreen', 'YellowOrange', 'Red']):
                
        if (theo_value * (1-evaluation_pattern[0]/100)) < obs_value < (theo_value * (1+evaluation_pattern[0]/100)):
            color = evaluation_colors[0]
        elif (theo_value * (1-evaluation_pattern[1]/100)) < obs_value < (theo_value * (1+evaluation_pattern[1]/100)):
            color = evaluation_colors[1]
        else:
            color = evaluation_colors[2]
             
        return color
    
    def color_evaluator_uplimits(self, obs_value, theo_value, evaluation_pattern, evaluation_colors = ['ForestGreen', 'YellowOrange', 'Red']):
        
        if (obs_value - theo_value) > evaluation_pattern[0]:
            color = evaluation_colors[0]
        elif (obs_value - theo_value) > evaluation_pattern[1]:
            color = evaluation_colors[1]            
        else:
            color = evaluation_colors[2]
        return color    
 
    def color_evaluator_lowlimits(self, obs_value, theo_value, evaluation_pattern, evaluation_colors = ['ForestGreen', 'YellowOrange', 'Red']):
       
        if (theo_value - obs_value) > evaluation_pattern[0]:
            color = evaluation_colors[0]
            
        elif (theo_value - obs_value) > evaluation_pattern[1]:
            color = evaluation_colors[1]            
        else:
            color = evaluation_colors[2]
           
        return color    
    
    def compare_variables(self, in_variable, in_variable_magnitude, variable_dict, CodeName, FileFolder):
                
        if in_variable in variable_dict.keys():
            variable_2_compare = variable_dict[in_variable]
        else:
            variable_2_compare = None
                
        if variable_2_compare != None:
            out_variable_magnitude = self.GetParameter_ObjLog(CodeName, FileFolder, variable_2_compare, 'float')
            
            if out_variable_magnitude != None:
                color = self.color_evaluator(in_variable_magnitude, out_variable_magnitude)
            else:
                color = 'black'

        else:
            color = 'black'
            
        return color
    
    def latex_header(self, table_address, table_type = 'standard', TitleColumn = None):

        if table_type == 'standard':

            #Generate object table
            self.doc = Document(table_address, documentclass='mn2e')             
            self.doc.packages.append(Package('preview', options=['active', 'tightpage',])) #Package to crop pdf to a figure
            self.doc.packages.append(Package('color', options=['usenames', 'dvipsnames',])) #Package to crop pdf to a figure
            self.doc.packages.append(Package('geometry', options=['pass', 'paperwidth=30in', 'paperheight=11in', ])) #Package to crop pdf to a figure
            self.doc.packages.append(Package('siunitx'))
            self.doc.packages.append(Package('booktabs'))
            self.doc.append(NoEscape(r'\sisetup{separate-uncertainty=true}'))
            
            #Table pre-commands
            self.doc.append(NoEscape(r'\begin{table*}[h]'))
            self.doc.append(NoEscape(r'\begin{preview}')) 
            self.doc.append(NoEscape(r'{\footnotesize'))
            self.doc.append(NoEscape(r'\centering'))
                 
            with self.doc.create(Tabular(self.table_format)) as self.table:
                if TitleColumn != None:
                    self.doc.append(NoEscape(r'\toprule'))
                    self.table.add_row(TitleColumn, escape=False)
                self.doc.append(NoEscape(r'\toprule'))
                #Declare the header
#                 self.table.add_hline()
                self.table.add_row(self.header_dict.values(), escape=False)
#                 self.table.add_hline()    
                self.doc.append(NoEscape(r'\midrule'))
                
                
        if table_type == 'lines log':
     
            #Generate object table
            self.doc = Document(table_address, documentclass='mn2e')             
            self.doc.packages.append(Package('preview', options=['active', 'tightpage',])) #Package to add new colors
            self.doc.packages.append(Package('color', options=['usenames', 'dvipsnames',])) #Package to crop pdf to a figure
            self.doc.packages.append(Package('geometry', options=['pass', 'paperwidth=30in', 'paperheight=11in', ])) #Package to crop pdf to a figure
            self.doc.packages.append(Package('siunitx'))
            self.doc.packages.append(Package('booktabs'))
            self.doc.append(NoEscape(r'\sisetup{separate-uncertainty=true}'))
            
            #Table pre-commands
            self.doc.append(NoEscape(r'\begin{table*}[h]'))
            self.doc.append(NoEscape(r'\begin{preview}')) 
            self.doc.append(NoEscape(r'{\footnotesize'))
            self.doc.append(NoEscape(r'\centering'))
                 
            with self.doc.create(Tabular(self.table_format)) as self.table:
                self.doc.append(NoEscape(r'\toprule'))
                self.table.add_row(['', '', 'HII Galaxy', CodeName, ''], escape=False)
                self.doc.append(NoEscape(r'\toprule'))
                #Declare the header
#                 self.table.add_hline()
                self.table.add_row(self.header_dict.values(), escape=False)
#                 self.table.add_hline()    
                self.doc.append(NoEscape(r'\midrule'))            
    
        return
    
    def table_header(self):
        
        with self.doc.create(Tabular(self.table_format)) as self.table:

            #Declare the header
            self.table.add_hline()
            self.table.add_row(self.header_dict.values(), escape=False)
            self.table.add_hline()
        
        return
    
    def table_footer(self, table_type = 'standard'):

        if table_type == 'standard':
            
            #Add a final line to the table 
            self.table.add_hline()
#             self.doc.append(NoEscape(r'\bottomrule'))

            #Close the preview
            self.doc.append(NoEscape('}'))
            self.doc.append(NoEscape(r'\end{preview}'))               
            self.doc.append(NoEscape(r'\end{table*}'))  
                 
            #Generate the document
#             self.doc.generate_tex()
            self.doc.generate_pdf(clean=True)

        return

class Pdf_printer():

    def __init__(self):
        
        self.pdf_type = None
        self.pdf_geometry_options = {'right'    : '1cm',
                                     'left'     : '1cm',
                                     'top'      : '1cm',
                                     'bottom'   : '2cm'}

    def create_pdfDoc(self, fname, pdf_type = 'graphs', geometry_options = None, document_class = u'article'):

        #TODO it would be nicer to create pdf object to do all these things

        self.pdf_type = pdf_type
        
        #Update the geometry if necessary (we coud define a dictionary distinction)
        if pdf_type == 'graphs':
            pdf_format = {'landscape':'true'}            
            self.pdf_geometry_options.update(pdf_format)
        
        elif pdf_type == 'table':
            pdf_format = {'landscape':'true',
                          'paperwidth':'30in',
                          'paperheight':'30in'}
            self.pdf_geometry_options.update(pdf_format)
            
        if geometry_options is not None:
            self.pdf_geometry_options.update(geometry_options)
        
        #Generate the doc
        self.pdfDoc = Document(fname, documentclass=document_class, geometry_options=self.pdf_geometry_options)

        if pdf_type == 'table':
            self.pdfDoc.packages.append(Package('preview', options=['active','tightpage',]))
            self.pdfDoc.packages.append(Package('hyperref', options=['unicode=true',]))
            self.pdfDoc.append(NoEscape(r'\pagenumbering{gobble}'))
            self.pdfDoc.packages.append(Package('nicefrac'))
            self.pdfDoc.packages.append(Package('color', options=['usenames', 'dvipsnames',])) #Package to crop pdf to a figure

        elif pdf_type == 'longtable':
            self.pdfDoc.append(NoEscape(r'\pagenumbering{gobble}'))

    def pdf_create_section(self, caption, add_page=False):
        
        with self.pdfDoc.create(Section(caption)):
        
            if add_page:
                self.pdfDoc.append(NewPage())
    
    def add_page(self):
        
        self.pdfDoc.append(NewPage())

        return
    
    def pdf_insert_image(self, image_address, fig_loc='htbp', width=r'1\textwidth'):
        
        with self.pdfDoc.create(Figure(position='h!')) as fig_pdf:
            fig_pdf.add_image(image_address, NoEscape(width))
        
        return

    def pdf_insert_table(self, column_headers=None, table_format = None, addfinalLine = True):

        #Set the table format
        if table_format is None:
            table_format = 'l' + 'c' * (len(column_headers) - 1)
        
        #Case we want to insert the table in a pdf
        if self.pdf_type != None:
        
            if self.pdf_type == 'table':
                self.pdfDoc.append(NoEscape(r'\begin{preview}')) 
            
                #Initiate the table
                with self.pdfDoc.create(Tabu(table_format)) as self.table:
                    if column_headers != None:    
                        self.table.add_hline()
                        self.table.add_row(map(str,column_headers), escape=False)
                        if addfinalLine:
                            self.table.add_hline()
                        
            elif self.pdf_type == 'longtable':
                
                #Initiate the table
                with self.pdfDoc.create(LongTable(table_format)) as self.table:
                    if column_headers != None:    
                        self.table.add_hline()
                        self.table.add_row(map(str,column_headers), escape=False)
                        if addfinalLine:
                            self.table.add_hline()

        #Table .tex without preamble
        else:
            self.table = Tabu(table_format)
            if column_headers != None:    
                self.table.add_hline()
                self.table.add_row(map(str,column_headers), escape=False)
                if addfinalLine:
                    self.table.add_hline()

    def pdf_insert_longtable(self, column_headers=None, table_format = None):

        #Set the table format
        if table_format is None:
            table_format = 'l' + 'c' * (len(column_headers) - 1)
        
        #Case we want to insert the table in a pdf
        if self.pdf_type != None:
        
            if self.pdf_type == 'table':
                self.pdfDoc.append(NoEscape(r'\begin{preview}')) 
            
            #Initiate the table
            with self.pdfDoc.create(Tabu(table_format)) as self.table:
                if column_headers != None:    
                    self.table.add_hline()
                    self.table.add_row(map(str,column_headers), escape=False)
                    self.table.add_hline()       

        #Table .tex without preamble
        else:
            self.table = LongTable(table_format)
            if column_headers != None:    
                self.table.add_hline()
                self.table.add_row(map(str,column_headers), escape=False)
                self.table.add_hline() 

    def addTableRow(self, input_row, row_format = 'auto', rounddig=4, rounddig_er=None, last_row=False):

        #Default formatting
        if row_format == 'auto':
            mapfunc = partial(self.format_for_table, rounddig=rounddig)
            output_row = map(mapfunc, input_row)
             
        #Append the row
        self.table.add_row(output_row, escape = False)
     
        #Case of the final row just add one line
        if last_row:
            self.table.add_hline()
                    
    def format_for_table(self, entry, rounddig = 4, rounddig_er=2, scientific_notation = False, nan_format = '-'):
                
        if rounddig_er == None:
            rounddig_er = rounddig
        
        #Check None entry
        if entry != None:
                
            #Check string entry
            if isinstance(entry, (str, unicode)): 
                formatted_entry = entry
               
            #Case of Numerical entry
            else:
                
                #Case of an array    
                scalarVariable = True
                if isinstance(entry, (Sequence, ndarray)):
                    
                    #Confirm is not a single value array
                    if len(entry) == 1:
                        entry           = entry[0]
                    #Case of an array
                    else:
                        scalarVariable  = False
                        formatted_entry = '_'.join(entry) # we just put all together in a "_' joined string    
                
                #Case single scalar        
                if scalarVariable:
                                  
                    #Case with error quantified
                    if isinstance(entry, UFloat):
                        formatted_entry = round_sig(nominal_values(entry), rounddig, scien_notation = scientific_notation) + r'$\pm$' +  round_sig(std_devs(entry), rounddig_er, scien_notation = scientific_notation)
                        
                    #Case single float
                    elif isnan(entry):
                        formatted_entry = nan_format
                        
                    #Case single float
                    else:
                        formatted_entry = round_sig(entry, rounddig, scien_notation = scientific_notation)
        else:
            #None entry is converted to None
            formatted_entry = 'None'
                
        return formatted_entry
        
    def fig_to_pdf(self, label=None, fig_loc='htbp', width=r'1\textwidth', add_page=False, *args, **kwargs):
        
        with self.pdfDoc.create(Figure(position=fig_loc)) as plot:
            plot.add_plot(width=NoEscape(width), placement='h', *args, **kwargs)
        
            if label is not None:
                plot.add_caption(label)
            
        if add_page:
            self.pdfDoc.append(NewPage())
        
    def generate_pdf(self, clean_tex = True, output_address=None):
        if output_address == None:
            if self.pdf_type == 'table':
                self.pdfDoc.append(NoEscape(r'\end{preview}')) 
            #self.pdfDoc.generate_pdf(clean_tex = clean_tex) # TODO this one does not work in windows
            self.pdfDoc.generate_pdf(clean_tex=clean_tex, compiler='pdflatex')
        else:
            self.table.generate_tex(output_address) 
            
        return

class Txt_Files_Manager(Images_Fits, Tables_Txt, Tables_Latex, Pdf_printer):

    def __init__(self):
        
        Images_Fits.__init__(self)
        Tables_Txt.__init__(self)
        Tables_Latex.__init__(self)
        Pdf_printer.__init__(self)
        
        self.Text_File_Type = None

    def Starlight_output_getdata(self, FileFolder, FileName):

        DataFile                    = open(FileFolder + FileName,"r")
        
        StarlightOutput             = DataFile.readlines()
        
        DataFile.close()
            
        # Synthesis Results - Best model #
        Chi2Line                    = self.LineFinder(StarlightOutput, "[chi2/Nl_eff]")
        AdevLine                    = self.LineFinder(StarlightOutput, "[adev (%)]")
        SumXdevLine                 = self.LineFinder(StarlightOutput, "[sum-of-x (%)]")
        v0_min_Line                 = self.LineFinder(StarlightOutput, "[v0_min  (km/s)]")
        vd_min_Line                 = self.LineFinder(StarlightOutput, "[vd_min  (km/s)]")
        Av_min_Line                 = self.LineFinder(StarlightOutput, "[AV_min  (mag)]")
    
        Nl_eff_line                 = self.LineFinder(StarlightOutput, "[Nl_eff]")
        
        
        SignalToNoise_Line          = self.LineFinder(StarlightOutput, "## S/N")
            
        l_norm_Line                 = self.LineFinder(StarlightOutput, "## Normalization info") + 1 
        llow_norm_Line              = self.LineFinder(StarlightOutput, "## Normalization info") + 2 
        lupp_norm_Line              = self.LineFinder(StarlightOutput, "## Normalization info") + 3     
        NormFlux_Line               = self.LineFinder(StarlightOutput, "## Normalization info") + 4 
        
        SpecLine                    = self.LineFinder(StarlightOutput, "## Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei")     #Location of my Spectrum in starlight output
       
        #Quality of fit                                         
        Chi2                        = float(StarlightOutput[Chi2Line].split()[0])                                        
        Adev                        = float(StarlightOutput[AdevLine].split()[0])
        SumXdev                     = float(StarlightOutput[SumXdevLine].split()[0])
        Nl_eff                      = float(StarlightOutput[Nl_eff_line].split()[0])
        v0_min                      = float(StarlightOutput[v0_min_Line].split()[0])
        vd_min                      = float(StarlightOutput[vd_min_Line].split()[0])
        Av_min                      = float(StarlightOutput[Av_min_Line].split()[0])
     
        #Signal to noise configuration                                           
        SignalToNoise_lowWave       = float(StarlightOutput[SignalToNoise_Line + 1].split()[0]) 
        SignalToNoise_upWave        = float(StarlightOutput[SignalToNoise_Line + 2].split()[0]) 
        SignalToNoise_magnitudeWave = float(StarlightOutput[SignalToNoise_Line + 3].split()[0]) 
        
        #Flux normailzation parameters                                          
        l_norm                      = float(StarlightOutput[l_norm_Line].split()[0])
        llow_norm                   = float(StarlightOutput[llow_norm_Line].split()[0])
        lupp_norm                   = float(StarlightOutput[lupp_norm_Line].split()[0])
        FluxNorm                    = float(StarlightOutput[NormFlux_Line].split()[0])
    
        Parameters                  = (Chi2, Adev, SumXdev, Nl_eff, v0_min, vd_min, Av_min, SignalToNoise_lowWave, SignalToNoise_upWave, SignalToNoise_magnitudeWave, l_norm, llow_norm, lupp_norm)
        
        #Spectra pixels location
        Pixels_Number               = int(StarlightOutput[SpecLine+1].split()[0])                     #Number of pixels in the spectra
        Ind_i                       = SpecLine+2                                                              #First pixel location      
        Ind_f                       = Ind_i + Pixels_Number                                                   #Final pixel location
        
        Input_Wavelength            = zeros(Pixels_Number)
        Input_Flux                  = zeros(Pixels_Number)
        Output_Flux                 = zeros(Pixels_Number)
        Output_Mask                 = zeros(Pixels_Number)
        
        for i in range(Ind_i, Ind_f): 
            Index                   = i - Ind_i
            Line                    = StarlightOutput[i].split()
            Input_Wavelength[Index] = float(Line[0])
            Input_Flux[Index]       = float(Line[1])*FluxNorm if Line[1] != '**********' else 0.0
            Output_Flux[Index]      = float(Line[2])*FluxNorm
            Output_Mask[Index]      = float(Line[3])
        
        MaskPixels                  = [[],[]]           #The 0 tag
        ClippedPixels               = [[],[]]           #The -1 tag
        FlagPixels                  = [[],[]]           #The -2 tag
        
        for j in range(len(Output_Mask)):
            PixelTag        = Output_Mask[j]
            Wave            = Input_Wavelength[j]
            if PixelTag == 0:
                MaskPixels[0].append(Wave)
                MaskPixels[1].append(Input_Flux[j])
            if PixelTag == -1:
                ClippedPixels[0].append(Wave)
                ClippedPixels[1].append(Input_Flux[j])            
            if PixelTag == -2:
                FlagPixels[0].append(Wave)
                FlagPixels[1].append(Input_Flux[j])            
      
        return Input_Wavelength, Input_Flux, Output_Flux, MaskPixels, ClippedPixels, FlagPixels, Parameters

    def GenerateStarlightFiles(self, FileFolder, FileName, CodeName, objData, X, Y, ExtraData = None, Velocity_Vector = ['FXK',  '0.0',  "10.0"], ComputerRoot = '/home/vital/', EmLinesFileExtension = "LickIndexes.txt", TwoArmsMode = False, ext_loop=''):
        #Velocity_Vector =                                                                               [Kinematics calculation type,    v0_start,    vd_start]
        
        Starlight_Folder    = ComputerRoot + 'Starlight/'
        Default_InputFolder = ComputerRoot + 'Starlight/Obs/'
        Default_MaskFolder  = ComputerRoot + 'Starlight/Masks/'
        Default_BasesFolder = ComputerRoot + 'Starlight/Bases/'
        Default_OutputFoler = ComputerRoot + 'Starlight/Output/'
    
        if TwoArmsMode == True:
            Color = None
            if 'Blue' in FileName:
                Color = 'Blue'
            elif 'Red' in FileName:
                Color = 'Red'
            
        #-----------------------     Generating Base File    ----------------------------------------------
#         BaseDataFile         = 'Dani_Bases_296.txt'
        BaseDataFile         = 'Dani_Bases_Extra.txt'

        #-----------------------     Generating Configuration File    -------------------------------------
        ConfigurationFile   = 'Sl_Config_v1.txt'
        
        #-----------------------Generating input spectra Textfile---------------------------------
        Interpolation       = interp1d(X, Y, kind = 'slinear')
        Wmin                = int(round(X[0],0))
        Wmax                = int(round(X[-1],0))
        
    #   Interpolate the new spectra to one angstrom per pixel resolution
        X_1Angs = range(Wmin+1,Wmax-1,1)
        Y_1Angs = Interpolation(X_1Angs)
        
        Sl_Input_Filename = FileName.replace(".fits", ".slInput")
        self.SaveSpectra_2_Starlight([X_1Angs, Y_1Angs], Default_InputFolder + Sl_Input_Filename)
        print '-- Starlight File:', Default_InputFolder + Sl_Input_Filename
        
        #-----------------------     Generating Mask File    ----------------------------------------------
            
        #Block Initial region
        Masks = []
        EdgeBoundary = 100
        Masks.append([Wmin, Wmin + EdgeBoundary, 'Lower_Edge'])
        Masks.append([Wmax - EdgeBoundary, Wmax , 'Upper_Edge'])
        
        #Import emision line location from lick indexes file
        LogName = CodeName + '_lick_indeces.txt'                             #Should change this to the Leak indexes plot  
        lick_idcs_df = read_csv(FileFolder + LogName, delim_whitespace = True, header = 0, index_col = 0, comment='L') #Dirty trick to avoid the Line_label row
        Labels_List, IniWave_List, FinWave_List  = lick_idcs_df.index.values, lick_idcs_df['Wave3'].values, lick_idcs_df['Wave4'].values
        
        #Loop through the lines and for the special cases increase the thickness
        for i in range(Labels_List.size):
            if (Wmin < IniWave_List[i])  and  (FinWave_List[i] < Wmax):
                Label = Labels_List[i]
                if (Label == 'H1_6563A') or (Label == 'O2_3726A') or (Label == 'H1_4340A') or (Label == 'N2_6584A') or (Label == 'O3_5007A') or (Label == 'S3_9531A'):
                    factor = 15
                else:
                    factor = 0
                    
                Masks.append([IniWave_List[i] - factor, FinWave_List[i] + factor, Label])
        
        if 'WHT' in FileName and TwoArmsMode == False:
            WHT_Wmatch_Wavelength   = objData.join_wavelength
            JoiningRegion_Begining  = WHT_Wmatch_Wavelength - 75
            JoiningRegion_End       = WHT_Wmatch_Wavelength + 75
            Not_MiddleRegionEncountered = False                #Boolean to check when we pass the match region... 
            Masks.append([JoiningRegion_Begining,JoiningRegion_End, 'WHT_Spectra_Joining'])
            
        else:
            Not_MiddleRegionEncountered = False                #Boolean to check when we pass the match region... 
            JoiningRegion_Begining = 0.0
            JoiningRegion_End = 0.0
                 
        if TwoArmsMode == True:
            MaskFileName = CodeName + "_" + Color + '_Mask.lineslog' + ext_loop
        else:
            MaskFileName = CodeName + '_Mask.lineslog' + ext_loop

        File = open(FileFolder + MaskFileName, "w")
        File.write(str(len(Masks)) + '\n')
        for k in range(len(Masks)):
            Line = str(Masks[k][0]) + '  ' + str(Masks[k][1]) + '  0.0  ' + str(Masks[k][2]) + '\n'
            File.write(Line)
         
        File.close()

        copyfile(FileFolder + MaskFileName, Default_MaskFolder + MaskFileName)
        print '-- Mask File:', Default_MaskFolder + MaskFileName
        
        #-----------------------     Generating output files    -------------------------------------------
    
        if ".fits" in FileName:
            Sl_Output_Filename = FileName.replace(".fits", ".slOutput")
        else:
            Sl_Output_Filename = FileName.replace(FileName[FileName.rfind("."):len(FileName)],'.slOutput')
    
        Sl_OutputFolder = Default_OutputFoler
        print '-- Output address:', Sl_OutputFolder + Sl_Output_Filename
    
        #-----------------------Generating Grid file---------------------------------
         
        GridLines = []
        GridLines.append("1")                               #"[Number of fits to run]"])
        GridLines.append(Default_BasesFolder)               #"[base_dir]"])
        GridLines.append(Default_InputFolder)               #"[obs_dir]"])
        GridLines.append(Default_MaskFolder)                #"[mask_dir]"])
        GridLines.append(Default_OutputFoler)               #"[out_dir]"])
        GridLines.append("-652338184")                      #"[your phone number]"])
        GridLines.append("4500.0 ")                         #"[llow_SN]   lower-lambda of S/N window"])
        GridLines.append("4550.0")                          #"[lupp_SN]   upper-lambda of S/N window"])
        GridLines.append("3400.0")                          #"[Olsyn_fin] upper-lambda for fit"])
        GridLines.append("12000.0")                         #"[Olsyn_fin] upper-lambda for fit"])
        GridLines.append("1.0")                             #"[Odlsyn]    delta-lambda for fit"])
        GridLines.append("1.0")                             #"[fscale_chi2] fudge-factor for chi2"])
        GridLines.append(Velocity_Vector[0])                #"[FIT/FXK] Fit or Fix kinematics"])
        GridLines.append("0")                               #"[IsErrSpecAvailable]  1/0 = Yes/No"])
        GridLines.append("0")                               #"[IsFlagSpecAvailable] 1/0 = Yes/No"])
    
        Redlaw = 'CCM'
        v0_start = Velocity_Vector[1]
        vd_start = Velocity_Vector[2]
    
        GridLines.append([Sl_Input_Filename, ConfigurationFile, BaseDataFile, MaskFileName, Redlaw, v0_start, vd_start, Sl_Output_Filename])     
        Grid_FileName = FileName.replace(FileName[FileName.rfind("."):len(FileName)],'.slGrid') + ext_loop
                
        File = open(Starlight_Folder + Grid_FileName,"w")        
        
        print '-- Grid File:', Starlight_Folder + Grid_FileName
        
        for i in range(len(GridLines) - 1):
            Parameter = GridLines[i]
            Element = str(Parameter) + "\n"
            File.write(Element)
        
        Element = "  ".join(GridLines[-1])+'\n'
        File.write(Element)
        File.close()
            
        return Grid_FileName, Sl_Output_Filename, Sl_OutputFolder, X_1Angs, Y_1Angs

    def GenerateStarlightFiles_ORIGINAL(self, FileFolder, FileName, CodeName, X, Y, ExtraData = None, Velocity_Vector = ['FIT',  '0.0',  "10.0"], ComputerRoot = '/home/vital/', EmLinesFileExtension = "LickIndexes.txt", TwoArmsMode = False):
        #Velocity_Vector =                                                                               [Kinematics calculation type,    v0_start,    vd_start]
        
        Starlight_Folder    = ComputerRoot + 'Starlight/'
        Default_InputFolder = ComputerRoot + 'Starlight/Obs/'
        Default_MaskFolder  = ComputerRoot + 'Starlight/Masks/'
        Default_BasesFolder = ComputerRoot + 'Starlight/Bases/'
        Default_OutputFoler = ComputerRoot + 'Starlight/Output/'
    
        if TwoArmsMode == True:
            Color = None
            if 'Blue' in FileName:
                Color = 'Blue'
            elif 'Red' in FileName:
                Color = 'Red'
            
        #-----------------------     Generating Base File    ----------------------------------------------
#         BaseDataFile         = 'Dani_Bases_296.txt'
        BaseDataFile         = 'Dani_Bases_Extra.txt'

        #-----------------------     Generating Configuration File    -------------------------------------
        ConfigurationFile   = 'Sl_Config_v1.txt'
        
        #-----------------------Generating input spectra Textfile---------------------------------
        Interpolation       = interp1d(X, Y, kind = 'slinear')
        Wmin                = int(round(X[0],0))
        Wmax                = int(round(X[-1],0))
        
    #   Interpolate the new spectra to one angstrom per pixel resolution
        X_1Angs = range(Wmin+1,Wmax-1,1)
        Y_1Angs = Interpolation(X_1Angs)
        
        Sl_Input_Filename = FileName.replace(".fits", ".slInput")
        self.SaveSpectra_2_Starlight([X_1Angs, Y_1Angs], Default_InputFolder + Sl_Input_Filename)
        print '-- Starlight File:', Default_InputFolder + Sl_Input_Filename
        
        #-----------------------     Generating Mask File    ----------------------------------------------
        
        if 'WHT' in FileName and TwoArmsMode == False:
            SpectraMeet             = self.GetParameter_ObjLog(CodeName, FileFolder, "Spectra_Meet")
            WHT_Wmatch_Wavelength   = self.GetParameter_ObjLog(CodeName, FileFolder, "WHT_Wmatch", Assumption='float')
            JoiningRegion_Begining  = WHT_Wmatch_Wavelength - 100
            JoiningRegion_End       = WHT_Wmatch_Wavelength + 100
            Not_MiddleRegionEncountered = False                #Boolean to check when we pass the match region... 
            
        else:
            Not_MiddleRegionEncountered = False                #Boolean to check when we pass the match region... 
            JoiningRegion_Begining = 0.0
            JoiningRegion_End = 0.0
    
    
        #Block Initial region
        Masks = []
        EdgeBoundary = 100
        Masks.append([Wmin, Wmin + EdgeBoundary, 'Lower_Edge'])
        
        #Import emision line location from lick indexes file
        LogName = CodeName + '_LickIndexes.txt'                             #Should change this to the Leak indexes plot  
        Labels_List = loadtxt(FileFolder + LogName, dtype=str, skiprows = 1, usecols = [0])
        IniWave_List, FinWave_List = loadtxt(FileFolder + LogName, dtype = float, skiprows = 1, usecols = (6,7), unpack = True)
        
        #In this first import we only load the masks within the spectra wavelength range
        for i in range(Labels_List.size):
            if (Wmin < IniWave_List[i])  and  (FinWave_List[i] < Wmax):
                Masks.append([IniWave_List[i],FinWave_List[i],Labels_List[i]])
        
        #Block final region
        Masks.append([Wmax - EdgeBoundary, Wmax , 'Upper_Edge'])
        MaskVector = [[],[],[]]
        
        #Generate the masks      
        for j in range(len(Masks)):
            Label               = Masks[j][2]
            Begining_point      = Masks[j][0]
            Ending_point        = Masks[j][1]
            
    
            if (Label == 'H1_6563A') or (Label == 'O2_3726A') or (Label == 'H1_4340A') or (Label == 'N2_6584A') or (Label == 'O3_5007A') or (Label == 'S3_9531A'):
                Increment = round((Masks[j][1] - Masks[j][0]) / 0.7, 0)
            elif (Label == 'Lower_Edge') or (Label == 'Upper_Edge'):  
                Increment = 0
            else:
                Increment = round((Masks[j][1] - Masks[j][0]) / 4, 0)
            
            IniWave = Masks[j][0] - Increment
            FinWave = Masks[j][1] + Increment
            
            if j > 0:
                PrevWave = MaskVector[1][-1]
                if IniWave <= PrevWave:
                    if FinWave <= PrevWave: 
                        MaskVector[2][-1] = MaskVector[2][-1] + ' ' +Label
                    else:    
                        MaskVector[1][-1] = FinWave
                        MaskVector[2][-1] = MaskVector[2][-1] + ' ' +Label
        
                else:
                    MaskVector[0].append(IniWave)
                    MaskVector[1].append(FinWave)
                    MaskVector[2].append(Label)
            else:
                MaskVector[0].append(IniWave)
                MaskVector[1].append(FinWave)
                MaskVector[2].append(Label)
                
            Case_Inside     = ((MaskVector[0][-1] >= JoiningRegion_Begining) and (MaskVector[1][-1] <= JoiningRegion_End))
            Case_WeMissedIt = ((MaskVector[0][-1] >= JoiningRegion_Begining) and (MaskVector[1][-1] >= JoiningRegion_End) and (Not_MiddleRegionEncountered == True))
            
            if Case_Inside:
                if Not_MiddleRegionEncountered == True:
                    MaskVector[0][-1] = JoiningRegion_Begining
                    MaskVector[1][-1] = JoiningRegion_End
                    MaskVector[2][-1] = 'Joining region'
                    Not_MiddleRegionEncountered = False
                else:
                    del MaskVector[0][-1]
                    del MaskVector[1][-1]
                    del MaskVector[2][-1]
            
            if Case_WeMissedIt:
                Ini_0   = MaskVector[0][-1]
                Fin_1   = MaskVector[1][-1]
                Lab     = MaskVector[2][-1]
                MaskVector[0][-1] = JoiningRegion_Begining
                MaskVector[1][-1] = JoiningRegion_End
                MaskVector[2][-1] = 'Joining region'
                MaskVector[0].append(Ini_0)
                MaskVector[1].append(Fin_1)
                MaskVector[2].append(Lab)            
                Not_MiddleRegionEncountered = False
        
        
        
        if TwoArmsMode == True:
            MaskFileName = CodeName + "_" + Color + '_Mask.lineslog'
        else:
            MaskFileName = CodeName + '_Mask.lineslog'
        
        #Esto como que jode el invento de antes....
        if SpectraMeet == 'True':
            MaskVector[0].append(JoiningRegion_Begining)
            MaskVector[1].append(JoiningRegion_End)
            MaskVector[2].append('Spectra meeting region')        
            
        
        File = open(FileFolder + MaskFileName, "w")
        File.write(str(len(MaskVector[0])) + '\n')
        for k in range(len(MaskVector[0])):
            Line = str(MaskVector[0][k]) + '  ' + str(MaskVector[1][k]) + '  0.0  ' + str(MaskVector[2][k]) + '\n'
            File.write(Line)
        
        File.close()
    
        copyfile(FileFolder + MaskFileName, Default_MaskFolder + MaskFileName)
        print '-- Mask File:', Default_MaskFolder + MaskFileName
        
        #-----------------------     Generating output files    -------------------------------------------
    
        if ".fits" in FileName:
            Sl_Output_Filename = FileName.replace(".fits", ".slOutput")
        else:
            Sl_Output_Filename = FileName.replace(FileName[FileName.rfind("."):len(FileName)],'.slOutput')
    
        Sl_OutputFolder = Default_OutputFoler
        print '-- Output address:', Sl_OutputFolder + Sl_Output_Filename
    
        #-----------------------Generating Grid file---------------------------------
         
        GridLines = []
        GridLines.append("1")                               #"[Number of fits to run]"])
        GridLines.append(Default_BasesFolder)               #"[base_dir]"])
        GridLines.append(Default_InputFolder)               #"[obs_dir]"])
        GridLines.append(Default_MaskFolder)                #"[mask_dir]"])
        GridLines.append(Default_OutputFoler)               #"[out_dir]"])
        GridLines.append("-652338184")                      #"[your phone number]"])
        GridLines.append("4500.0 ")                         #"[llow_SN]   lower-lambda of S/N window"])
        GridLines.append("4550.0")                          #"[lupp_SN]   upper-lambda of S/N window"])
        GridLines.append("3400.0")                          #"[Olsyn_fin] upper-lambda for fit"])
        GridLines.append("12000.0")                         #"[Olsyn_fin] upper-lambda for fit"])
        GridLines.append("1.0")                             #"[Odlsyn]    delta-lambda for fit"])
        GridLines.append("1.0")                             #"[fscale_chi2] fudge-factor for chi2"])
        GridLines.append(Velocity_Vector[0])                #"[FIT/FXK] Fit or Fix kinematics"])
        GridLines.append("0")                               #"[IsErrSpecAvailable]  1/0 = Yes/No"])
        GridLines.append("0")                               #"[IsFlagSpecAvailable] 1/0 = Yes/No"])
    
        Redlaw = 'CCM'
        v0_start = Velocity_Vector[1]
        vd_start = Velocity_Vector[2]
    
        GridLines.append([Sl_Input_Filename, ConfigurationFile, BaseDataFile, MaskFileName, Redlaw, v0_start, vd_start, Sl_Output_Filename])     
        Grid_FileName = FileName.replace(FileName[FileName.rfind("."):len(FileName)],'.slGrid')
                
        File = open(Starlight_Folder + Grid_FileName,"w")        
        
        print '-- Grid File:', Starlight_Folder + Grid_FileName
        
        for i in range(len(GridLines) - 1):
            Parameter = GridLines[i]
            Element = str(Parameter) + "\n"
            File.write(Element)
        
        Element = "  ".join(GridLines[-1])+'\n'
        File.write(Element)
        File.close()
            
        return Grid_FileName, Sl_Output_Filename, Sl_OutputFolder, X_1Angs, Y_1Angs
 
    def SaveSpectra_2_Starlight(self, TableOfValues, FileAddress):
        #WARNING: This is a very old approach
        
        File = open(FileAddress,"w")
        
        for i in range(len(TableOfValues[0])):
            Sentence = ''
            for j in range(len(TableOfValues)):
                Sentence = Sentence + ' ' + str(TableOfValues[j][i])
                
            File.write(Sentence + '\n')
        
        File.close()
    
        return    
    
    def LineFinder(self, myFile,myText):
        
        #THIS IS VERY INNEFFICIENT OPTION
        for i in range(len(myFile)):
            if myText in myFile[i]:
                return i

    def ImportDispersionVelocity(self, FileFolder, CodeName, LinesLogExtension):
  
        self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = "H1_4861A", Parameter_Header='sigma', LinesLog_suffix = LinesLogExtension)
  
#         HBeta_Sigma         = float(GetParameterFromDigitalLines(FileLines, "H1_4861A", "sigma", 2))
#         CaIII_8498_Sigma    = GetParameterFromDigitalLines(FileLines, "CaIII_8498", "sigma", 2)
#         CaIII_8542_Sigma    = GetParameterFromDigitalLines(FileLines, "CaIII_8542", "sigma", 2)
#         CaIII_8662_Sigma    = GetParameterFromDigitalLines(FileLines, "CaIII_8662", "sigma", 2)
        
        O3_5007             = self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = "O3_5007A",   Parameter_Header='sigma', LinesLog_suffix = LinesLogExtension)
        CaIII_8498_Sigma    = self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = "CaIII_8498", Parameter_Header='sigma', LinesLog_suffix = LinesLogExtension)
        CaIII_8542_Sigma    = self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = "CaIII_8542", Parameter_Header='sigma', LinesLog_suffix = LinesLogExtension)
        CaIII_8662_Sigma    = self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = "CaIII_8662", Parameter_Header='sigma', LinesLog_suffix = LinesLogExtension)
 
        Sigma               = None
        c_SI                = 300000.0
    
    #     if CaIII_8542_Sigma != None:
    #         Sigma =  float(CaIII_8542_Sigma) / 8542.0 * c_SI
    #     else:
        Sigma = float(O3_5007)  / 5007.0 * c_SI
        
        if (Sigma <= 0) or (Sigma == None):
            Sigma = 100.0
        
        # we always return HBeta_Sigma     
        return Sigma

    def GetParameter_LineLog(self, CodeName, FileFolder, LineLabel, Parameter_Header, LinesLog_suffix, typeParameter=float, LinesLogHeader_Address = None, verbose = False):
        
        #The lines log includes the datatype suffix... is this right????
        #LinesLog_suffix = _WHT_LinesLog_v3        
        #Not sure if I can read this file with 
#         if LinesLogHeader_Address == None:
#             LinesLogHeader_Address = self.Default_LinesLogHeaderAddress
        
#         print 'this is the file I am trying to open', LinesLogHeader_Address, 'hola'
#         print 'pero esta es buena', self.Default_LinesLogHeaderAddress
#         print 'aqui', self.Obj_LinesLog_Headers
        
#         Textfile_Headers                = loadtxt(LinesLogHeader_Address, dtype=str, skiprows = 1, usecols = [1], unpack = True)
        Header_Index                    = where(self.Obj_LinesLog_Headers==Parameter_Header)[0][0]
            
        Labels_Column                   = loadtxt(FileFolder + CodeName + LinesLog_suffix, dtype = str, skiprows = 2, usecols = [0] , unpack = True) 
        
        Parameter_Column                = loadtxt(FileFolder + CodeName + LinesLog_suffix, dtype = typeParameter, skiprows = 2, usecols = [Header_Index] , unpack = True) 
    
        Parameter = None   
                
        if len(where(Labels_Column==LineLabel)[0]) != 0:
            Label_Index                     = where(Labels_Column==LineLabel)[0][0]  
            Parameter                       = Parameter_Column[Label_Index]
            
        elif verbose == True:
            print '-- Parameter: ', Parameter_Header, 'not found for object', CodeName
                
        return Parameter

    def getFlux_LinesLog(self, FileAddress, Row, Column_identifier = None, flux_type = None, error_type = None, ufloat_check = True):
        
        #In case no Header is introduced to describe the recorded line the default is 'Line_Label' 
        if Column_identifier == None:
            Row_identifier = self.Labels_ColumnHeader
                
        #Lods the column which identify the parameter we want
        Labels_columnIndex  = where(self.Obj_LinesLog_Headers == Row_identifier)[0][0]
        
        Labels_column       = loadtxt(FileAddress, dtype=str, skiprows = self.LinesLog_HeaderLength, usecols = [Labels_columnIndex])
        Row_index           = where(Labels_column == Row)[0][0]
                
        #In case the flux type is not defined, the code will check if the line is blended, in that case it will use the gaussian fit, otherwise the integrated value
        if flux_type == None:
            Blended_columnIndex = where(self.Obj_LinesLog_Headers == self.BlendedLines_ColumnHeader)[0][0]
            Blended_column      = loadtxt(FileAddress, dtype=str, skiprows = self.LinesLog_HeaderLength, usecols = [Blended_columnIndex])
            Blended_value       = Blended_column[Row_index]
            
            if Blended_value == 'None':
                flux_type = 'FluxBrute'
            else:
                flux_type = 'FluxGauss'
        
        #If no error type is introduced we use ErrorEL_MCMC
        if error_type == None:
            error_type = self.GaussianError_ColumnHeader 
        
        #Load the right fluxes and error
        Flux_columnIndex            = where(self.Obj_LinesLog_Headers == flux_type)[0][0]
        Error_columnIndex           = where(self.Obj_LinesLog_Headers == error_type)[0][0]
        Flux_column, Error_column   = loadtxt(FileAddress, dtype=str, skiprows = self.LinesLog_HeaderLength, usecols = [Flux_columnIndex, Error_columnIndex], unpack=True)
        
        #Return the correct value
        if ufloat_check:
            return ufloat(Flux_column[Row_index], Error_column[Row_index])
        else:
            return Flux_column[Row_index], Error_column[Row_index]
        
    def getColumn_LinesLog(self, FileAddress, Column_Header, data_type = None, headersize = 0):
        
        #Decide the type of the data
        if data_type == None:
            if (Column_Header == self.Labels_ColumnHeader) or (Column_Header == self.Ion_ColumnHeader) or (Column_Header == self.BlendedLines_ColumnHeader):
                data_type = str
            else:
                data_type = float
        
        #Case single column: we return that column
        if isinstance(Column_Header, str):
            
            #Get the index of the header
            columnIndex = where(self.Obj_LinesLog_Headers == Column_Header)[0][0]            
            column      = loadtxt(FileAddress, dtype=data_type, skiprows = headersize, usecols = [columnIndex])
            
            return column
                
        #Case of several column we get a dictionary
        elif isinstance(Column_Header, (Sequence, ndarray)):
            if isinstance(Column_Header[0], str):          
                column_indicies = where(in1d(self.Obj_LinesLog_Headers, Column_Header, assume_unique=True))[0]
                columns         = transpose(loadtxt(FileAddress, dtype=data_type, skiprows = headersize, usecols = column_indicies))
                column_dict     = dict(zip(Column_Header, columns))
                
            else:
                columns         = transpose(loadtxt(FileAddress, dtype=data_type, skiprows = headersize, usecols = Column_Header))
                string_indexes  = array(map(str,Column_Header))
                column_dict     = dict(zip(string_indexes, columns))
                
            return column_dict
            
    def GetParameter_ObjLog(self, CodeName, FileFolder, Parameter, Assumption = None, sigfig = 5, strformat = '{:.5e}', logtype = 'Object'):
                
        #In this case we generate the address from the codename log
        if logtype == 'Object':
            ObjLog_Address          = FileFolder + CodeName + '_log.txt'
        
        #In this case we give directly the address
        if logtype == 'Catalogue':
            ObjLog_Address              = FileFolder + CodeName
                
        ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments    = loadtxt(ObjLog_Address, dtype=str, skiprows = 2, usecols = [0,1,2,3], unpack = True)
        
        Parameter_Index                                                         = where(ObjLog_Parameters == Parameter)[0][0]
        Parameter_Magnitude, Parameter_Error, Parameter_Comment                 = ObjLog_Magnitudes[Parameter_Index], ObjLog_Errors[Parameter_Index], ObjLog_Comments[Parameter_Index]
        
        #Basic test to check the quality of the analysis
        CheckPhysicality = ((Parameter_Magnitude != 'None') and (Parameter_Magnitude != 'nan') and (Parameter_Magnitude != '-'))
        
        #Special cases in which we want to convert de variable type or difine a default value
        if Assumption != None:
            
            # Basic float import
            if Assumption == 'float':
                if CheckPhysicality:
                    Converted_Parameter = float(Parameter_Magnitude)
                else:
                    Converted_Parameter = None
            
                if (Parameter_Error != '-') and (Converted_Parameter != None):
                    Converted_Parameter = ufloat(float(Parameter_Magnitude), float(Parameter_Error))
            
            if Assumption == 'string':                
                if CheckPhysicality:
                    if (Parameter_Error != '-') and (Parameter_Error != 'nan'):
                        Converted_Parameter = strformat.format(float(round_sig(float(Parameter_Magnitude),sigfig))) +  r'$\pm$' + strformat.format(float(round_sig(float(Parameter_Error),sigfig)))
                    else:
                        Converted_Parameter = strformat.format(float(round_sig(float(Parameter_Magnitude),sigfig)))   
                    
                else:
                    Converted_Parameter = '-'
            
            #Temperature needed case for HII region
            elif Assumption == 'Min_Temp':
                if CheckPhysicality:
                    Converted_Parameter = float(Parameter_Magnitude[0:Parameter_Magnitude.find('+/-')])
                else:
                    Converted_Parameter = 10000.0
    
            #Density needed case for HII region
            elif Assumption == 'Min_Den':
                if CheckPhysicality:
                    Converted_Parameter = float(Parameter_Magnitude[0:Parameter_Magnitude.find('+/-')])
                    if Converted_Parameter < 75.0:
                        Converted_Parameter = 50.0
                else:
                    Converted_Parameter = 100.0
                    
            elif Assumption == 'Min_HeII':
                #WARNING: Update for the new organizing error format                             
                if CheckPhysicality:
                    Converted_Parameter = float(Parameter_Magnitude[0:Parameter_Magnitude.find('+/-')])
                else:
                    Converted_Parameter = 0.1

            elif Assumption == 'Min_HeIII':
                #WARNING: Update for the new organizing error format                             
                if CheckPhysicality:
                    Converted_Parameter = float(Parameter_Magnitude[0:Parameter_Magnitude.find('+/-')])
                else:
                    Converted_Parameter = 0.0

            elif Assumption == 'cHbeta_min':
                #WARNING: Update for the new organizing error format                             
                if CheckPhysicality:
                    Converted_Parameter = float(Parameter_Magnitude)
                else:
                    Converted_Parameter = None
            
                if (Parameter_Error != '-') and (Converted_Parameter != None):
                    Converted_Parameter = ufloat(float(Parameter_Magnitude), float(Parameter_Error))
                    
                if Converted_Parameter != None:
                    if Converted_Parameter < 0.0:
                        Converted_Parameter = ufloat(0.0, 0.0)   
                    
            elif Assumption == "MatchingSpectra":
                
                if Parameter_Magnitude == 'True':
                    Wave_Index          = where(ObjLog_Parameters == 'WHT_Wmatch')[0][0]
                    Converted_Parameter = float(ObjLog_Magnitudes[Wave_Index])
        
                else:
                    Blue_lambdaF_ind    = ObjLog_Magnitudes[where(ObjLog_Parameters == 'WHT_Wmax_Blue')[0][0]]
                    Red_lambdaF_ind     = ObjLog_Magnitudes[where(ObjLog_Parameters == 'WHT_Wmin_Red')[0][0]]
                    
                    Converted_Parameter = (float(Blue_lambdaF_ind) + float(Red_lambdaF_ind)) / 2
            
            return Converted_Parameter
        
        else:
            
            if Parameter_Error != '-':
                Parameter_Magnitude = ufloat(float(Parameter_Magnitude), float(Parameter_Error))
            
            return Parameter_Magnitude
    
    def save_ChemicalData(self, FileAddress, Parameter, Magnitude, Error = None, Assumption = None, Log_extension = None):
        
        #HERE WE NEED TO ADD THE POSIBILITY OF COMMENTS (THE ERROR SHOULD BE DETECTED FROM A UNUMPY ARRAY)
        #Would be nice to make this a list updater
        #Loading the data from the text file
#         print 'Parameter', Parameter, isinstance(Magnitude, UFloat), type(Magnitude), 'goes to'
        
        #Open text file
        ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments    = loadtxt(FileAddress, dtype=str, usecols = [0,1,2,3], unpack = True)
    
        #Saving the new value for the given parameter
        Parameter_Index                                                         = where(ObjLog_Parameters == Parameter)[0][0]

        #Save the error in the case of a nparray
        if isinstance(Magnitude, UFloat) and (Error == None):
            ObjLog_Magnitudes[Parameter_Index]                                  = Magnitude.nominal_value
            ObjLog_Errors[Parameter_Index]                                      = Magnitude.std_dev
        elif Error != None:
            ObjLog_Magnitudes[Parameter_Index]                                  = str(Magnitude)
            ObjLog_Errors[Parameter_Index]                                      = str(Error)
        elif Magnitude != None:
            ObjLog_Magnitudes[Parameter_Index]                                  = str(Magnitude)
        elif Magnitude == None:
            ObjLog_Magnitudes[Parameter_Index]                                  = 'None'
            ObjLog_Errors[Parameter_Index]                                      = '-'
            
        #Saving the text file
        savetxt(FileAddress, transpose((ObjLog_Parameters, ObjLog_Magnitudes, ObjLog_Errors, ObjLog_Comments)), fmt='%s')
    
        return

    def prepareLog(self, LogFileName, LogFolder, LogFileFormat = None, ForceRemake = False):
        
        #Set the type of log file (THIS SHOULD BE CONFIGURE TO  USE STRINGS BINDED TO DEFAULT TEXT FILE-ADDRESSES AND FORMATS
        if LogFileFormat == None:
            LogFileFormat = self.Default_ObjectLogFormatAddress
        
        #We import the log format
        Format_Parameters= loadtxt(LogFileFormat, dtype=str, usecols = [0,1,2,3])
        
        #We import the data from previous log. If ForceRemake flag is on, we make new files
        log_address = LogFolder + LogFileName
        if isfile(log_address) and (ForceRemake == False):
            Log_Parameters = loadtxt(log_address, dtype=str, usecols = [0,1,2,3], unpack = True)

            #Loop through format file, if the same parameter is encountered in the log, it is stored
            CoincidenceArray = in1d(Format_Parameters[0],Log_Parameters[0],True)
#             for i in range(1, len(Format_Parameters)):
#                 if Format_Parameters[0]
        
        return
    
    def SetLogFile(self, LogFile, LogFolder, LogFileFormat = None, ForceRemake = False):
    
        #WARNING: THIS METHODOLOGY IS VERY OLD!! IT SHOULD BE UPDATED USING THE EXAMPLES IN SaveParameter_ObjLog AND GetParameter_ObjLog  
        if LogFileFormat == None:
            LogFileFormat = self.Default_ObjectLogFormatAddress
        
        Structure_Log = self.File2Table(LogFileFormat, "")
    
        #We import the data from previous log. If ForceRemake flag is on, we make new files
        if isfile(LogFolder + LogFile) and (ForceRemake == False):
            
            Obj_Log = self.File2Table(LogFolder, LogFile)
        
            for i in range(1,len(Structure_Log)):
                for j in range(len(Obj_Log)):
                    if Structure_Log[i][0] == Obj_Log[j][0]:
                        Structure_Log[i][1] = Obj_Log[j][1]
                        Structure_Log[i][2] = Obj_Log[j][2]
                        Structure_Log[i][3] = Obj_Log[j][3]
                        
        OutputFile = open(LogFolder + LogFile, 'w')
        ColumnSizes = []
        
        for i in range(len(Structure_Log[0])):
            Biggest_Length = 0
            for j in range(1,len(Structure_Log)):
                if len(Structure_Log[j][i]) > Biggest_Length:
                    Biggest_Length = len(Structure_Log[j][i])
            ColumnSizes.append(Biggest_Length + 2)
            
        NewFormatLine = "%" + str(ColumnSizes[0]) + "s" + "%" + str(ColumnSizes[1]) + "s" + "%" + str(ColumnSizes[2]) + "s" + "%" + str(ColumnSizes[3]) + "s"
            
        for z in range(1, len(Structure_Log)):  
            NewLine = NewFormatLine % (Structure_Log[z][0],Structure_Log[z][1],Structure_Log[z][2],Structure_Log[z][3])
            OutputFile.write(NewLine+"\n")
    
        OutputFile.close()
                    
        return

    def File2Table(self, Address,File):
        
        #THIS IS A VERY OLD APPROACH IT SHOULD BE UPDATED
        File = Address + File
        TextFile = open(File,"r")
        Lines = TextFile.readlines()
        TextFile.close()
        
        Table = []
        for line in Lines:
            Table.append(line.split())
                
        return Table

    def replace_line(self, file_name, line_num, text):
        
        Input_File = open(file_name, 'r')
        lines = Input_File.readlines()
        lines[line_num] = text
        
        Output_File = open(file_name, 'w')
        Output_File.writelines(lines)
        
        Input_File.close()
        Output_File.close()

    def Starlight_Launcher(self, Grid_FileName, ComputerRoot = '/home/vital/'):
        
        chdir(ComputerRoot + 'Starlight/')
        Command = './StarlightChains_v04.exe < ' + Grid_FileName
        print "Launch command:", Command
    #     Command = './Starlight_v04_Mac.exe < grid_example1.in'  
        
        p = Popen(Command, shell=True, stdout=PIPE, stderr=STDOUT)
             
        for line in p.stdout.readlines():
            print line,
            
        retval = p.wait()
        
        return

    def getEmLine_dict(self, Lines_List, Mode = 'Auto'):
        
        #The dictionary for the headers parameter should include the type
        Lines_dict = OrderedDict()
        Wavelength_Vector = zeros(len(Lines_List))
        
        #Decide the flux we must use: Integrated unless the line is blended, in that case we take the gaussian value
        for i in range(len(Lines_List)):
            
            Line = Lines_List[i]
 
            if Mode == 'Auto':
                
                Blended_Value = self.GetParameter_LineLog(self.CodeName, self.FileFolder, LineLabel = Line, Parameter_Header = 'Blended_Set',  typeParameter=str, LinesLog_suffix = self.AbundancesExtension)
                
                if Blended_Value == 'None':
                    FluxType = 'FluxBrute'
                else:
                    FluxType = 'FluxGauss'
                    
            elif Mode == 'Gauss':
                FluxType = 'FluxBrute'
                
            elif Mode == 'Integrated':
                FluxType = 'FluxBrute'
                        
            #Load the flux from the lines log. #WARNING: Not very efficient scheme 
            Flux                        = self.GetParameter_LineLog(self.CodeName, self.FileFolder, LineLabel = Line, Parameter_Header = FluxType,  LinesLog_suffix = self.AbundancesExtension)
            Flux_err                    = self.GetParameter_LineLog(self.CodeName, self.FileFolder, LineLabel = Line, Parameter_Header = 'ErrorEL_MCMC',  LinesLog_suffix = self.AbundancesExtension)
            Wavelength_Vector[i]        = self.GetParameter_LineLog(self.CodeName, self.FileFolder, LineLabel = Line, Parameter_Header = 'TheoWavelength',  LinesLog_suffix = self.AbundancesExtension)
                        
            #If the line was observed Line[i] = uarray(Flux, Flux_err). Otherwise Line = None
            Lines_dict[Line] = self.check_issues(magnitude = (Line, Flux, Flux_err), parameter_type = 'EmFlux')
            
        return Lines_dict, Wavelength_Vector

    def getEmLine_FluxDict(self, Lines_List, CodeName, FileFolder, AbundancesExtension, Mode = 'Auto'):
        
        #this is usefull should be universal with the previous one
        Lines_dict = OrderedDict()
        Wavelength_Vector = zeros(len(Lines_List))
        
        for i in range(len(Lines_List)):
            
            Line = Lines_List[i]
            #The dictionary for the headers parameter should include the type
            #Decide the flux we must use: Integrated unless the line is blended, in that case we take the gaussian value
            if Mode == 'Auto':
                
                Blended_Value = self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = Line, Parameter_Header = 'Blended_Set',  typeParameter=str, LinesLog_suffix = AbundancesExtension)
                
                if Blended_Value == 'None':
                    FluxType = 'FluxBrute'
                else:
                    FluxType = 'FluxGauss'
                    
            elif Mode == 'Gauss':
                FluxType = 'FluxBrute'
                
            elif Mode == 'Integrated':
                FluxType = 'FluxBrute'
                        
            #Load the flux from the lines log. #WARNING: Not very efficient scheme 
            Flux                        = self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = Line, Parameter_Header = FluxType,            LinesLog_suffix = AbundancesExtension)
            Flux_err                    = self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = Line, Parameter_Header = 'ErrorEL_MCMC',      LinesLog_suffix = AbundancesExtension)
            Wavelength_Vector[i]        = self.GetParameter_LineLog(CodeName, FileFolder, LineLabel = Line, Parameter_Header = 'TheoWavelength',    LinesLog_suffix = AbundancesExtension)
                        
            #Store the line (if it was observed
            if (Flux != None) and (Flux_err != None):
                Lines_dict[Line] = ufloat(Flux, Flux_err)
           
        return Lines_dict

class pd_Tools():
    
    def quick_indexing(self, df):

        df['quick_index'] = np_nan
    
        counter = 1
        for obj in df.index:
            
            if df.loc[obj,'Ignore_article'] != 'yes':
            
                if notnull(df.loc[obj,'Favoured_ref']):
                    df.loc[obj,'quick_index'] = df.loc[obj,'Favoured_ref']
                else:
                    df.loc[obj,'quick_index'] = "FTDTR-" + str(counter)
                    counter += 1
                    
        self.idx_include = notnull(df['quick_index'])
                 
class Dazer_Files(Txt_Files_Manager):

    def __init__(self):
           
        self.catalogue_properties_frame = None
        self.separators_rows            = [0, 17, 38, 40, 49, 58, 65, 68, 70, 76, 79, 84, 87, 91, 94, 97, 104, 109, 111, 118, 124, 127, 132, 135, 139, 142, 145, 152, 157, 160, 167, 175, 182, 189, 196, 203, 210, 213, 216, 226, 235, 244, 253, 262, 271, 280]
        
        #WARNING: This should be adapted so it can read the configuration from the installation folder
        #Dazer files structures
        Txt_Files_Manager.__init__(self)
        
        #All Headers
        self.ObjectLog_extension                = '_log.txt'
        self.Labels_ColumnHeader                = 'Line_Label'
        self.Ion_ColumnHeader                   = 'Ion'
        self.Wavelength_ColumnHeader            = 'TheoWavelength'
        self.BruteFlux_ColumnHeader             = 'Flux_Int'
        self.GaussianFlux_ColumnHeader          = 'Flux_Gauss'
        self.IntegratedError_ColumnHeader       = 'Error_FluxI'
        self.GaussianError_ColumnHeader         = 'Error_FluxG'
        self.EqW_ColumnHeader                   = 'Eqw'
        self.EqW_error_ColumnHeader             = 'Error_Eqw'        
        self.Line_Continuum_ColumnHeader        = 'Continuum_Median'
        self.Line_ContinuumSigma_Header         = 'Continuum_sigma'
        self.Helium_label_ColumnHeader          = 'HeI'
        self.BlendedLines_ColumnHeader          = 'group_label'
        
        #WARNING need a better way to find bin folder
        __location__ = path.realpath(path.join(getcwd(), path.dirname(__file__)))
        root_folder  = __location__[0:__location__.find('dazer')] + 'dazer/'
       
        #Dazer logs structure files location
        self.Default_ObjectLogFormatAddress     = root_folder + 'format/DZT_ObjectLog_Format.dz'
        self.Default_LinesLogHeaderAddress      = root_folder + 'format/DZT_LineLog_Headers.dz'
        self.list_lines_address                 = root_folder + 'format/DZT_EmLine_List_Emission.dz'

        self.LinesLog_HeaderFormatLength        = 1
        self.LinesLog_HeaderLength              = 2
        self.Obj_LinesLog_Headers               = loadtxt(self.Default_LinesLogHeaderAddress, dtype=str, skiprows = self.LinesLog_HeaderFormatLength, usecols = [1])
        self.LinesLogFormat_dict                = self.getColumn_LinesLog(self.list_lines_address, [0,1,2,3], data_type=str, headersize=0)
        
        #Files for the Helium inference models
        self.Hydrogen_CollCoeff_TableAddress    = root_folder + 'Astro_Libraries/Neutral_Hydrogen_Collisional_Correction_coef.txt'
        self.Helium_CollCoeff_TableAddress      = root_folder + 'Astro_Libraries/Neutral_Helium_Collisional_Correction_coef.txt'
        self.Helium_OpticalDepth_TableAddress   = root_folder + 'Astro_Libraries/Helium_OpticalDepthFunction_Coefficients.txt'
        
    def dazer_tableFormat(self, FileName, FileFolder, header_dict):
        
        #Declare document
        self.doc = Document(FileFolder + FileName, documentclass='mn2e')
        
        #Declare document packages
        self.doc.packages.append(Package('preview', options=['active', 'tightpage',])) #Package to crop pdf to a figure
        self.doc.packages.append(Package('color', options=['usenames', 'dvipsnames',])) #Package to crop pdf to a figure
        self.doc.packages.append(Package('geometry', options=['pass', 'paperwidth=30in', 'paperheight=11in', ])) #Package to crop pdf to a figure
         
        #Table pre-commands
        self.doc.append(NoEscape(r'\begin{table*}[h]'))
        self.doc.append(NoEscape(r'\begin{preview}')) 
        self.doc.append(NoEscape(r'{\footnotesize'))
        self.doc.append(NoEscape(r'\centering'))
        
        #Declare table format
        self.table_format = 'l' + ''.join([' c' for s in range(len(header_dict) - 1)])
        
        return
    
    def dazer_tableHeader(self, table, header_dict, mode = 'standard'):
        
        #Declare the header
        table.add_hline()
        table.add_row(header_dict.values(), escape=False)
        table.add_hline()        
        
        return

    def dazer_tableCloser(self, table_type='pdf', clean_check=False, mode='estandard'):
        
        #Close the preview
        self.doc.append(NoEscape('}'))
        self.doc.append(NoEscape(r'\end{preview}'))               
        self.doc.append(NoEscape(r'\end{table*}'))  
             
        #Generate the document
        # doc.generate_tex()
        
        if table_type == 'pdf':
            self.doc.generate_pdf(clean=clean_check)    
                
        elif table_type == 'tex':
            self.doc.generate_tex()
        
        return

    def load_Galaxy_Ratios(self, lineslog_frame, Atom_dict):
            
        Ratios_dict         = OrderedDict()
        RatiosColor_dict    = OrderedDict()
        
        lines_intensity = lineslog_frame['line_Int']
        
        #Calculate oxygen flux ratios
        if ('O3_4959A' in lineslog_frame.index) and ('O3_5007A' in lineslog_frame.index):
            Ratios_dict['O3_ratio']         = lines_intensity['O3_5007A'] / lines_intensity['O3_4959A']
            RatiosColor_dict['O3_ratio']    = self.color_evaluator(Ratios_dict['O3_ratio'], Atom_dict['O3_ratio'])
        else:
            Ratios_dict['O3_ratio']         = 'None'
            RatiosColor_dict['O3_ratio']    = 'Black'
             
        #Calculate nitrogen flux ratios
        if ('N2_6548A' in lineslog_frame.index) and ('N2_6584A' in lineslog_frame.index):
            Ratios_dict['N2_ratio']         = lines_intensity['N2_6584A'] / lines_intensity['N2_6548A']
            RatiosColor_dict['N2_ratio']    = self.color_evaluator(Ratios_dict['N2_ratio'], Atom_dict['N2_ratio'])
        else:
            Ratios_dict['N2_ratio']         = 'None'
            RatiosColor_dict['N2_ratio']    = 'Black'
             
        #Calculate sulfur flux ratios
        if ('S3_9069A' in lineslog_frame.index) and ('S3_9531A' in lineslog_frame.index):
            Ratios_dict['S3_ratio']         = lines_intensity['S3_9531A'] / lines_intensity['S3_9069A']
            RatiosColor_dict['S3_ratio']    = self.color_evaluator(Ratios_dict['S3_ratio'], Atom_dict['S3_ratio'])
        else:
            Ratios_dict['S3_ratio']         = 'None'
            RatiosColor_dict['S3_ratio']    = 'Black'
             
        return Ratios_dict, RatiosColor_dict

    def load_object_lines(self, FileFolder, CodeName, Log_extension, Mode = 'Auto', chbeta_coef = None):
        
        #Since the data comes in two formats it needs to be uploaded using two commands
        log_address         = FileFolder + CodeName + Log_extension
        String_columns      = [self.Labels_ColumnHeader, self.Ion_ColumnHeader, self.BlendedLines_ColumnHeader]
        Float_Columns       = ['TheoWavelength', 'FluxBrute', 'FluxGauss', 'Eqw', 'ErrorEqw', 'ErrorEL_MCMC']
        linelog_dict        = self.getColumn_LinesLog(log_address, String_columns, data_type = str, headersize = self.LinesLog_HeaderLength)
        float_Column_dict   = self.getColumn_LinesLog(log_address, Float_Columns, data_type = float, headersize = self.LinesLog_HeaderLength)
    
        #Update the first dictioary with the data from the second
        linelog_dict.update(float_Column_dict)
           
        #Empty array to store the right flux for each line
        Flux_array          = zeros(len(linelog_dict['FluxBrute']))
        
        for i in range(len(linelog_dict['FluxBrute'])):
            
            if Mode == 'Auto':
                
                Blended_Value = linelog_dict[self.BlendedLines_ColumnHeader][i]
                
                if Blended_Value == 'None':
                    Flux_mag= linelog_dict[self.BruteFlux_ColumnHeader][i]
                else:
                    Flux_mag = linelog_dict[self.GaussianFlux_ColumnHeader][i]
                    
            elif Mode == 'Gauss':
                Flux_mag = linelog_dict[self.GaussianFlux_ColumnHeader][i]
                
            elif Mode == 'Integrated':
                Flux_mag = linelog_dict[self.BruteFlux_ColumnHeader][i]        
        
            Flux_array[i] = Flux_mag        
        
        linelog_dict['Flux'] = uarray(Flux_array, linelog_dict[self.GaussianError_ColumnHeader])
        
        if chbeta_coef != None:
            if type(chbeta_coef) == str:
                cHbeta_mag = self.GetParameter_ObjLog(CodeName, FileFolder, Parameter = chbeta_coef, Assumption = 'cHbeta_min')
            else:
                cHbeta_mag = chbeta_coef

            f_lines                     = self.Reddening_curve(linelog_dict[self.Wavelength_ColumnHeader],  'Cardelli1989')
            Flux_derred                 = linelog_dict['Flux'] * unnumpy_pow(10,  f_lines * cHbeta_mag)
            
            linelog_dict['Flux']        = Flux_derred
            linelog_dict['f_lambda']    = f_lines
                
        return linelog_dict

    def load_lineslog_frame(self, lines_log_address, mode = 'Auto', chbeta_coef = None, key_check = 'group_label'):

        #Load a frame from the lines log       
        #lines_frame = read_csv(lines_log_address, skiprows = [0], delim_whitespace = True, header = 0, index_col = 0)
        lines_frame = read_csv(lines_log_address, delim_whitespace = True, header = 0, index_col = 0, comment='L')  #Dirty trick to avoid the Line_label row

        #Load the line flux           
        if mode == 'Auto':  #Gaussian flux for blended lines, integrated for the rest
    
            Int_indexes     = lines_frame[key_check] == 'None'
            Gauss_indexes   = lines_frame[key_check] != 'None'
                        
            F_Int_uarray    = uarray(lines_frame.loc[Int_indexes, 'flux_intg'].values,      lines_frame.loc[Int_indexes, 'flux_intg_er'].values) 
            F_Gauss_uarray  = uarray(lines_frame.loc[Gauss_indexes, 'flux_gauss'].values,   lines_frame.loc[Gauss_indexes, 'flux_gauss_er'].values)
            
            lines_frame.loc[Int_indexes, 'line_Flux']     = F_Int_uarray
            lines_frame.loc[Gauss_indexes, 'line_Flux']   = F_Gauss_uarray
 
        elif mode == 'Integrated':  #All integrated flux
            lines_frame['line_Flux'] = uarray(lines_frame['flux_intg'].values,  lines_frame['flux_intg_er'].values) 
        
        elif mode == 'Gauss': #All gaussian flux
            lines_frame['line_Flux'] = uarray(lines_frame['flux_gauss'].values,   lines_frame['flux_gauss_er'].values)
            
        #Load the line continuum         
        lines_frame['line_continuum'] = uarray(lines_frame['zerolev_mean'].values,  lines_frame['zerolev_std'].values) 
        
        #Load the line equivalent width
        lines_frame['line_Eqw'] = lines_frame['line_Flux'].values / lines_frame['line_continuum']
                    
        return lines_frame
    
    def load_catalogue_frame(self, Files_list):
        
        for i in range(len(Files_list)):
             
            CodeName, FileName, FileFolder = self.Analyze_Address(Files_list[i])

            #load object_frame
            obj_frame = read_csv(FileFolder + FileName, skiprows = self.separators_rows, delim_whitespace = True, names = ['mag', 'error', 'comments'])
            
            if self.catalogue_properties_frame is None:
                self.catalogue_properties_frame = DataFrame(index = obj_frame.index)

            #Filling the rows with nominal and error quantification
            index_Mag_and_error = (obj_frame['mag'] != '-') & (obj_frame['error'] != '-') & (obj_frame['mag'].notnull()) & (obj_frame['error'].notnull())
            nominal             = obj_frame.loc[index_Mag_and_error, 'mag'].values
            std_dev             = obj_frame.loc[index_Mag_and_error, 'error'].values
            self.catalogue_properties_frame.loc[index_Mag_and_error, CodeName] = uarray(nominal, std_dev)
            
            #Filling the rows with nominal quantification
            index_Mag           = (obj_frame['mag'] != '-') & (obj_frame['error'] == '-')
            nominal             = obj_frame.loc[index_Mag, 'mag'].values
            self.catalogue_properties_frame.loc[index_Mag, CodeName] = nominal

            #Filling the rows with None as np_nan
            index_Mag           = (obj_frame['mag'] == 'None') & (obj_frame['error'] == '-')
            nominal             = empty(index_Mag.sum())
            nominal.fill(np_nan)
            self.catalogue_properties_frame.loc[index_Mag, CodeName] = nominal

#             #Filling the rows with nan as np_nan
# #             index_Mag           = (obj_frame['mag'] == float('nan')) & (obj_frame['error'] == float('nan'))
#             index_Mag             = m_isnan(obj_frame['mag'])
# 
#             print 'indices', index_Mag
# #             print 'estos indices', index_Mag
# #             print 'Esta cosa', obj_frame.loc['OI_HI_pn']['mag'], type(obj_frame.loc['OI_HI_pn']['mag']), float('nan')
# #             print 'iguales', m_isnan(obj_frame.loc['OI_HI_pn']['mag'])
#             nominal             = empty(len(index_Mag))
#             nominal.fill(np_nan)
            
            self.catalogue_properties_frame.loc[index_Mag, CodeName] = nominal



           
            # #Code to read the number of separators
            # file = open(Folder + FileName)
            # lines = file.readlines()
            # values = []
            # 
            # for i in range(len(lines)):
            #     if lines[i][0] == '-':
            #         values.append(i)
            # 
            # print values
            
        return self.catalogue_properties_frame

    def getLineslog_frame(self, FileAddress):
        
        obj_frame = read_csv(Loglines, skiprows = [0], delim_whitespace = True, header =0, index_col=0)
        
        return obj_frame
    
    def get_line_value(self, linelog_dict, line_label, variable_in = 'Line_Label', variable_out = 'Flux'):
                
        line_Index  = where(linelog_dict[variable_in] == line_label)[0][0]
        magnitude   = linelog_dict[variable_out][line_Index]
        
        return magnitude

    def generate_catalogue_tree(self, catalogue_dict = None, obj = None, files_dict = None):

        if catalogue_dict != None:
             
            if isdir(catalogue_dict['Folder']) == False:
                mkdir(catalogue_dict['Folder'])
            
                if isdir(catalogue_dict['Obj_Folder']) == False:
                    mkdir(catalogue_dict['Obj_Folder'])
    
                if isdir(catalogue_dict['Data_Folder']) == False:
                    mkdir(catalogue_dict['Data_Folder'])
                    
            if obj != None:
                FileName = obj.replace('[', '').replace(']', '').replace('; ', '')

                if isdir(catalogue_dict['Obj_Folder'] + FileName + '/') == False:
                    mkdir(catalogue_dict['Obj_Folder'] + FileName + '/')
                
                return FileName
                    
        return

    def FamilyOfItemsInArray(self, Array):   
        
        d = {}
        for i in Array:
            if i in d:
                d[i] = d[i]+1
            else:
                d[i] = 1
        a = d.keys()
        a.sort()
        a[::-1]
        return list(a)    
            
    def File2Lines(self, Address,File):
        #THIS IS VERY OLD IT SHOULD BE UPDATED
        File = Address + File
        TextFile = open(File,"r")
        Lines = TextFile.readlines()
        TextFile.close()
        return Lines

    def SlOutput_2_StarData(self, FileFolder, FileName):
        
        #THIS IS VERY OLD IT SHOULD BE UPDATED
        Sl_Data = self.File2Lines(FileFolder, FileName)
        
        BasesLine = self.LineFinder(Sl_Data, "[N_base]")                                                          #Location of my normalization flux in starlight output
        Bases = int(Sl_Data[BasesLine].split()[0])
          
        Sl_DataHeader = self.LineFinder(Sl_Data, "# j     x_j(%)      Mini_j(%)     Mcor_j(%)     age_j(yr)")        #Location of my normalization flux in starlight output
        Ind_i = Sl_DataHeader + 1
        Ind_f = Sl_DataHeader + Bases                                    
        
        index = []
        x_j = []
        Mini_j = []
        Mcor_j = []
        age_j = []
        Z_j = []
        LbyM =[]
        Mstars = []
    
        for j in range(Ind_i, Ind_f + 1): 
            myDataLine = Sl_Data[j].split()
            index.append(float(myDataLine[0]))
            x_j.append(float(myDataLine[1]))
            Mini_j.append(float(myDataLine[2]))
            Mcor_j.append(float(myDataLine[3]))
            age_j.append(float(myDataLine[4]))
            Z_j.append(float(myDataLine[5]))
            LbyM.append(float(myDataLine[6]))
            Mstars.append(float(myDataLine[7]))
    
        return index, x_j, Mini_j, Mcor_j, age_j, Z_j, LbyM, Mstars

    def load_excel_DF(self, frame_address):
        
        #File which stores all the original file sheets #WARNING: this will not work if more than one excel DF loaded
        self.ipExcel_sheetColumns = OrderedDict() 
        
        #Load excel file:
        with ExcelFile(frame_address) as xlsx_file:
        
            #Load all sheets
            list_Df_sheet_i, sheets_names = [], xlsx_file.sheet_names
            for sheet in sheets_names:
                df_i = xlsx_file.parse(sheet, index_col=0)
                list_Df_sheet_i.append(df_i)
                self.ipExcel_sheetColumns[sheet] = list(df_i.columns.values)
        
        #Combine individual sheet-df into one
        df = concat(list_Df_sheet_i, axis=1)
        
        #Combine nominal and error columns for the same variables
        df_columns = df.columns.values
        for column in df_columns:
            if column + '_err' in df_columns:
                
                #Scheme only to combine rows which only contain value
                idcs_nan = df[column + '_err'].isnull() 
                
                #Empty error cells produce simple floats                             
                df.loc[idcs_nan, column] = df.loc[idcs_nan, column + '_err']
        
                #Otherwise they produce uarray                             
                df.loc[~idcs_nan, column] = uarray(df.loc[~idcs_nan, column], df.loc[~idcs_nan, column + '_err'])
                                      
                #df[column] = uarray(df[column].values, df[column + '_err'].values)
                
                #Remove error column from the dataframe
                df.drop(column + '_err', axis=1, inplace=True)
                                     
        return df

    def save_excel_DF(self, dataframe, frame_address, colsPerDict_dict = None, df_sheet_format = None, df_columns_format = []):
                
        if colsPerDict_dict == None:
            
            #Use stored format
            if df_sheet_format != None:
                colsPerDict_dict = OrderedDict()
                colsPerDict_dict['OBJ_ID']                      = ['objcode','obscode','SDSS_reference','Extra_IDs', 'Repeated_obs', 'Favoured_ref', 'SDSS_PLATE', 'SDSS_MJD', 'SDSS_FIBER', 'SDSS_Web', 'NED_Web']
                colsPerDict_dict['Data_location']               = ['Blue_file', 'Red_file',  'zBlue_file', 'zRed_file', 'tellRed_file', 'reduction_fits', 'emission_fits']
                colsPerDict_dict['OBJ_diagnostic']              = ['SIII_lines', 'T_low', 'T_high', 'O_valid', 'N_valid', 'S_valid', 'Ignore_article', '[OIII]5007A/[OIII]4959A','[NII]6548A/[NII]6584A','[SIII]9531A/[SIII]9069A','[OIII]5007A/[OIII]4959A_emis','[NII]6548A/[NII]6584A_emis','[SIII]9531A/[SIII]9069A_emis']
                colsPerDict_dict['Fits_properties']             = ['aperture','Blue_Grating','Red_Grating','Blue_CENWAVE','Red_CENWAVE','Dichroic','RA','DEC','UT_OBS','Wmin_Blue','Wmax_Blue','Wmin_Red','Wmax_Red']
                colsPerDict_dict['Reduction_data']              = ['obsfolder','calibration','calibration_star','telluric_star','Standard_stars','reduc_tag','join_wavelength','h_gamma_valid', 'z_SDSS', 'z_Blue', 'z_Blue_error', 'z_Red', 'z_Red_error']                
                colsPerDict_dict['Reddening']                   = ['E(B-V)_Galactic_dust', 'cHbeta_reduc', 'cHbeta_emis', 'cHbeta_G03_bar', 'cHbeta_G03_average', 'cHbeta_G03_supershell']
                colsPerDict_dict['Physical_Data']               = ['neSII','neOII','TeOII','TeSII','TeNII','TeOIII','TeSIII','TeOII_from_TeOIII','TeNII_from_TeOIII','TeSIII_from_TeOIII','TeOIII_from_TeSIII']
                colsPerDict_dict['Chemical_Abundances']         = ['SII_HII','SIII_HII','SIV_HII', 'ICF_SIV','OII_HII','OII_HII_3279A','OII_HII_7319A', 'OII_HII_ffO2', 'O_R3200', 'O_R3200_ffO2', 'O_R7300', 'O_R3', 'OIII_HII','NII_HII','ArIII_HII','ArIV_HII','HeII_HII_from_O','HeIII_HII_from_O','HeII_HII_from_S','HeIII_HII_from_S','SI_HI','OI_HI', 'OI_HI_ff02','NI_OI','NI_HI','HeI_HI_from_O','HeI_HI_from_S','Ymass_O','Ymass_S']
                
                colsPerDict_dict['Physical_Data_emis']          = map(lambda orig_string: orig_string + '_emis',    colsPerDict_dict['Physical_Data'])
                colsPerDict_dict['Chemical_Abundances_emis']    = map(lambda orig_string: orig_string + '_emis',    colsPerDict_dict['Chemical_Abundances'])
                colsPerDict_dict['Physical_Data_emis2nd']       = map(lambda orig_string: orig_string + '_emis2nd', colsPerDict_dict['Physical_Data'])
                colsPerDict_dict['Chemical_Abundances_emis2nd'] = map(lambda orig_string: orig_string + '_emis2nd', colsPerDict_dict['Chemical_Abundances'])
                
                colsPerDict_dict['Physical_Data_G03bar']   = map(lambda orig_string: orig_string + '_G03bar', colsPerDict_dict['Physical_Data'])
                colsPerDict_dict['Chemical_Abundances_G03bar'] = map(lambda orig_string: orig_string + '_G03bar', colsPerDict_dict['Chemical_Abundances'])
                
                colsPerDict_dict['Physical_Data_G03average']   = map(lambda orig_string: orig_string + '_G03average', colsPerDict_dict['Physical_Data'])
                colsPerDict_dict['Chemical_Abundances_G03average'] = map(lambda orig_string: orig_string + '_G03average', colsPerDict_dict['Chemical_Abundances'])

                colsPerDict_dict['Physical_Data_superS']   = map(lambda orig_string: orig_string + '_G03superS', colsPerDict_dict['Physical_Data'])
                colsPerDict_dict['Chemical_Abundances_superS'] = map(lambda orig_string: orig_string + '_G03superS', colsPerDict_dict['Chemical_Abundances'])


            #Everything into one sheet
            else:
                colsPerDict_dict = {'sheet' : list(dataframe.columns.values)}
                    
        #Check if dataframe includes unumpy arrays
        Ufloat_ocurrences   = dataframe.applymap(lambda x: isinstance(x, UFloat))        
        Ufloat_Columns      = Ufloat_ocurrences.apply(lambda x: (x == True).any())
        columns_WithError   = Ufloat_Columns[Ufloat_Columns].index.values
        
        if columns_WithError.shape[0] > 0:
            for column in columns_WithError:
                
                idcs_nan        = dataframe[column].isnull()
                original_data   = dataframe[column]
                
                #Add a new column
                dataframe.insert(dataframe.columns.get_loc(column) + 1, column + '_err', full(len(idcs_nan), np_nan)) #Place te err column after the its variable
                
                mag     = nominal_values(original_data[~idcs_nan].values)
                errors  = std_devs(original_data[~idcs_nan].values)
   
                dataframe.loc[~idcs_nan, column]           = mag
                dataframe.loc[~idcs_nan, column + '_err']  = errors
                   
                #Update the dictionary to add the error entry
                for entry in colsPerDict_dict:
                    if column in colsPerDict_dict[entry]:
                        colsPerDict_dict[entry].insert(colsPerDict_dict[entry].index(column) + 1, column + '_err')
               
        #List with the excel column
        columns_letters = list(ascii_uppercase)
        for i in range(3):
            for letter in ascii_uppercase:
                columns_letters.append(ascii_uppercase[i] + letter)
            
        #The final number of columns
        df_columns = list(dataframe.columns.values)
      
        with ExcelWriter(frame_address, engine='xlsxwriter') as writer:
             
            #Saving the sheets
            for sheet in colsPerDict_dict:
                 
                #Trick to remove not available variables from the dataframe while preserving sorting
                items = set(df_columns) & set(colsPerDict_dict[sheet])
                sheet_columns = sorted(items, key=lambda element: df_columns.index(element) + colsPerDict_dict[sheet].index(element))
                                 
                dataframe[sheet_columns].to_excel(writer, sheet_name=sheet)
                worksheet = writer.sheets[sheet]
                workbook  = writer.book
                 
                #Saving the columns
                for idx in range(len(sheet_columns) + 1): #We add one more to include the index column
                    if (idx != 0) or (sheet_columns[idx-1] in df_columns):
                        if sheet_columns[idx-1] not in df_columns_format:
                             
                            format = workbook.add_format()
                            format.set_align('right')
                            letter = '{columm_letter}:{columm_letter}'.format(columm_letter = columns_letters[idx])
                            
                            #For the index column
                            if letter == 'A:A':
                                header_maxlengh = len('Objects') + 2
                                data_maxlength  = dataframe.index.astype(str).map(len).max() + 2
                                column_width_set = header_maxlengh if header_maxlengh > data_maxlength else data_maxlength
                            #Rest of columns
                            else:
                                header_maxlengh = len(sheet_columns[idx-1]) + 2
                                data_maxlength  = dataframe[sheet_columns[idx-1]].astype(str).map(len).max() + 2
                                column_width_set = header_maxlengh if header_maxlengh > data_maxlength else data_maxlength
                            
                            #Set the format
                            worksheet.set_column(letter, column_width_set, format) 
                                     
            writer.save()

class File_Manager(Dazer_Files, pd_Tools):

    def __init__(self):
    
        #Run the classes for all the files we are able to deal with
        Dazer_Files.__init__(self)
    
        self.Arguments                  = []
        self.Arguments_Check            = None
        
        self.Tags                       = []        #MAYBE THIS COULD BE A DICTIONARY 
        
        self.Command_Launching_Files    = []
        self.Command_Launching_Folder   = None
        self.Flags_List                 = []
        
        self.RootFolder                 = None        
        self.verbosing                  = True
            
        self.ListFiles                  = []        #THIS ONE COULD BE A CONFLICT... CAREFULL IN THE PLOT MANAGER
        
        #Default extensions for the files we are treating. This should be loaded from text file
        self.Extensions_dictionary      = {
                                           'Starlight ouput'                : '.slOutput',
                                           'Reduction Instructions'         : '_Tasks_Configuration.txt',
                                           'Lineslog'                       : '_WHT_linesLog_reduc.txt'
                                           }
         
        #Launching methods
        self.Define_RootFolder()         
        
        self.ScriptCode                 = None
        self.ScriptName                 = None
        self.ScriptAddress              = None
        
        self.ErrorList                  = []
   
    def Arguments_Checker(self, LineaCommando):
        
        Num_Argum = len(LineaCommando)
                
        if Num_Argum == 1:
            self.Arguments_Check = False
            self.Arguments.append(LineaCommando[0])
            
        elif Num_Argum > 1:
            for i in range(len(LineaCommando)):
                self.Arguments_Check = True
                self.Arguments.append(LineaCommando[i])
       
    def Argument_Handler(self):   

        self.Command_Launching_Folder = getcwd()
        
        if self.Arguments_Check:
            for i in range(1, len(self.Arguments)):
                if "--" in self.Arguments[i]:
                    self.Flags_List.append(self.Arguments[i])
                else:
                    self.Command_Launching_Files.append(self.Command_Launching_Folder + '/' + self.Arguments[i])
            
            Num_Files = len(self.Command_Launching_Files)
            Num_Flags = len(self.Flags_List)
             
            if Num_Files > 0:
                print "-Files to treat:"
                print self.Command_Launching_Files
            if Num_Flags > 0:
                print "-FlagsList activated:"
                print self.Flags_List

    def Define_RootFolder(self):
        
        if self.RootFolder == None:
        
            MacName     = 'et'
            DellName    = 'foreshadowing'
            UbuntuName  = 'foreshadowing-G750JX'
    
            if gethostname() == MacName:
                self.RootFolder = '/Users/INAOE_Vital/'
            elif gethostname() == UbuntuName:
                self.RootFolder = '/home/delosari/'
            elif gethostname() == DellName:
                self.RootFolder = '/home/vital/'

    def File_Finder(self, Folder, myPattern):   
        
        #Define the list to store the files (should it be a self)
        FileList = []
            
        if type(myPattern) is not list:
            myPatternList = [myPattern]
            
        else:
            myPatternList = myPattern
                
        for Root, Dirs, Archives in walk(Folder):
            for Archive in Archives:
                Meets_one_Pattern = False
                for i in range(len(myPatternList)):
                    if (myPatternList[i] in Archive):
                        
                        #Security check to make sure we are not treating dummy files
                        if "~" not in Archive:                   
                            Meets_one_Pattern = True
                            
                if Meets_one_Pattern:        
                    if Root.endswith("/"):
                        FileList.append(Root + Archive)
                    else:
                        FileList.append(Root + "/" + Archive)
                                               
        return FileList

    def Analyze_Address(self, FileAddress, verbose = True):
        
        #Distinguish the three components from the address line
        FolderName      = FileAddress[0:FileAddress.rfind("/")+1]
        FileName        = FileAddress[FileAddress.rfind("/")+1:len(FileAddress)]
        CodeName        = FolderName[FolderName[0:-1].rfind("/")+1:len(FolderName)-1]
        
        #Special case for all spectra reduction nomenclature
        if FileName.startswith("obj") or FileName.startswith("std"):   
            CodeName    = FileName[3:FileName.find("_")]
            
        if verbose:
            print '--Treating file', CodeName, '(', FileName, ')', '\n'
        
        return CodeName, FileName, FolderName        

    def get_script_code(self):
        
        #Checking for arguments in terminal
        self.Arguments_Checker(argv)
        
        #Defining the script name, folder and order 
        self.ScriptName     = self.Arguments[0][self.Arguments[0].rfind("/")+1:len(self.Arguments[0])]
        self.ScriptFolder   = self.Arguments[0][0:self.Arguments[0].rfind("/")]        

        #WARNING THIS ORDER DEFINITION WILL NEED TO BE VARIABLE
        if self.ScriptName[2]   == '_':
            self.ScriptCode     =  self.ScriptName[0:2]
        else:
            self.ScriptCode     = ''

        return self.ScriptCode

    def Folder_Explorer(self, FilePattern, Containing_Folder, CheckComputer=False, Sort_Output = None, verbose = True):
        
        #Moving from an absolute to a relative structure
        if CheckComputer == True:
            Containing_Folder = self.RootFolder + Containing_Folder
        
        #Checking for arguments in terminal
        self.Arguments_Checker(argv)
        
        #Defining the script name, folder and order 
        self.ScriptName     = self.Arguments[0][self.Arguments[0].rfind("/")+1:len(self.Arguments[0])]
        self.ScriptFolder   = self.Arguments[0][0:self.Arguments[0].rfind("/")]
        
        #WARNING THIS ORDER DEFINITION WILL NEED TO BE VARIABLE
        if self.ScriptName[2]   == '_':
            self.ScriptCode     =  self.ScriptName[0:2]
                
        #Command from Terminal
        if self.Arguments_Check == True:
            self.Tags.append('Terminal')
            self.Argument_Handler()
            
            if len(self.Command_Launching_Files) == 1:
                self.Tags.append('SingleFile')
            
            self.ListFiles = self.Command_Launching_Files
                
        #Command from Eclipse
        else:
            self.Tags.append('Editor')
            
            FilesList = self.File_Finder(Containing_Folder, FilePattern)

            self.ListFiles = FilesList    
            
            if len(self.ListFiles) == 1:
                self.Tags.append('SingleFile')
                
        if Sort_Output == 'alphabetically':
            
            self.ListFiles = sorted(self.ListFiles)
        
        
        if verbose:
            print "Initiating "                         + self.Arguments[0][self.Arguments[0].rfind("/")+1:len(self.Arguments[0])] + " script\n"
            if self.ScriptCode != None: 
                print '- Code order:', self.ScriptCode, '\n'
              
            print "-Files found meeting the pattern: "   + str(FilePattern), "@", Containing_Folder, ':' 
            for File in self.ListFiles:
                print '\t-', File[File.rfind('/')+1:len(File)]
            print '\n'
        
        return self.ListFiles

    def File_to_data(self, FileFolder, FileName, InputProperties = None):

        #Fits files:
        if (".fit" in FileName) or ('fits' in FileName):
            Wave, Flux, ExtraData = self.Fits_to_Data(FileFolder, FileName)

            # I need this in case the Flux output is an array of arrays... need to deal with this in the lower level
            if type(Flux[0]) == type(Wave):                                      
                Y = Flux[0]
                X = Wave
            else:
                X, Y = Wave, Flux
            
            return X, Y, ExtraData
        
        #Text files
        elif self.Analyze_TextFile(FileFolder, FileName):
            
            if self.Extensions_dictionary['Starlight ouput'] in FileName:
                #                               0    1       2      3        4      5        6            7                    8                            9                10        11           12
                #       Parameters vector = (Chi2, Adev, SumXdev, Nl_eff, v0_min, vd_min, Av_min, SignalToNoise_lowWave, SignalToNoise_upWave, SignalToNoise_magnitudeWave, l_norm, llow_norm, lupp_norm)

                Input_Wavelength, Input_Flux, Output_Flux, MaskPixels, ClippedPixels, FlagPixels, Parameters = self.Starlight_output_getdata(FileFolder, FileName)
    
                return Input_Wavelength, Input_Flux, Output_Flux, MaskPixels, ClippedPixels, FlagPixels, Parameters
            
            elif self.Extensions_dictionary['Pyraf object reduction file'] in FileName:
                
                #In here the dictionary should look like : InputProperties = {'Columns' : ['Object', 'telluric_star']}
                InputProperties['HeaderSize']       = 2
                InputProperties['StringIndexes']    = True
                InputProperties['datatype']         = str
                InputProperties['unpack_check']     = True
                
                return self.get_ColumnData(**InputProperties)

    def Analyze_TextFile(self, FileFolder, FileName):
        
        FileCheck = False
            
        if (guess_type(FileFolder + FileName)[0] == 'text/plain') or (guess_type(FileFolder + FileName)[0] == 'application/x-ns-proxy-autoconfig') or (guess_type(FileFolder + FileName)[0] == None):
            FileCheck = True
    
        return FileCheck
    
    def GenerateFolder(self, FolderName, HostingAddress, FolderType = None):
        #WARNING: THERE CAN BE ISSUES IF THE FOLDER CANNOT BY CREATED
        
        if FolderType == 'Catalogue':
            if isdir(HostingAddress + FolderName) == False:
                makedirs(HostingAddress + FolderName)
            
        
        if isdir(HostingAddress + FolderName) == False:
            print 'WARNING: The folder could not be created'
            
        return
       
    def moveFile(self, FileName, HomeFolder, DestinationFolder, NewName = None, DeleteOriginal = True):
        
        #WARNING: THERE CAN BE ISSUES IF THE FOLDER CANNOT BY CREATED
        if isdir(DestinationFolder) == True:
            
            #Define initial address
            InitialAddress =  HomeFolder + FileName
            
            #Define destination address
            if NewName == None:
                FinalAddress = DestinationFolder + FileName
            else:
                FinalAddress = DestinationFolder + DestinationFolder
            
            #Move file
            shutil_move(InitialAddress, FinalAddress)
        
            #Check if file has been moved:
            if isfile(FinalAddress) == False:   
                print 'WARNING: File', FinalAddress, 'was not moved' 
        
        else:
            print 
            exit('WARNING: destination folder could not be found for:\n'+DestinationFolder+'\n'+FileName)
    
    
        return

    def copyFile(self, FileName, HomeFolder, DestinationFolder, NewName = None):
        
        #WARNING: THERE CAN BE ISSUES IF THE FOLDER CANNOT BY CREATED
        if isdir(DestinationFolder) == True:
            
            #Define initial address
            InitialAddress =  HomeFolder + FileName
            
            #Define destination address
            if NewName == None:
                FinalAddress = DestinationFolder + FileName
            else:
                FinalAddress = DestinationFolder + NewName
            
            #Move file
            copyfile(InitialAddress, FinalAddress)
        
            #Check if file has been moved:
            if isfile(FinalAddress) == False:   
                print 'WARNING: File', FinalAddress, 'was not moved' 
        
        else:
            print 
            exit('WARNING: destination folder could not be found for:\n'+DestinationFolder+'\n'+FileName)
    
    
        return    
    
 
#     def FindAndOrganize(self, FilePattern, Containing_Folder, unpack = False, CheckComputer=False, Sort_Output = 'Alpha'):
#         #THIS IS THE OLD CODE I USE I THE PLOTTINGMANAGER
#         
#         #Moving from an absolute to a relative structure
#         if CheckComputer == True:
#             FilesFolders = self.RootFolder + Containing_Folder
#         
#         #Checking for arguments in terminal
#         self.Argument_Check(argv)
#          
#         print "Initiating "                         + self.Arguments[0][self.Arguments[0].rfind("/")+1:len(self.Arguments[0])] + " task\n"
#            
#         print "Files found meeting the pattern: "   + str(FilePattern), "@", FilesFolders, ':\n' 
#            
#         #Command from Terminal
#         if self.ArgumentsCheck == True:
#             self.Tags.append('Terminal')
#             self.Argument_Handler()
#             
#             if len(self.Command_Launching_Files) == 1:
#                 self.Tags.append('SingleFile')
#             
#             self.ListFiles = (self.Command_Launching_Files,)
#             
#             return self.ListFiles
#     
#         #Command from Eclipse
#         else:
#             self.Tags.append('Editor')
#             
#             if unpack == False:     
#                 FilesList = self.File_Finder(FilesFolders, FilePattern)
# 
#                 Single_List = []
#                 Combined_List = []
#                 
#                 LeftPattern = '_@'
#                 RightPattern = '@_'
# 
#                 for FileAddress in FilesList:
#                     FileName = FileAddress[FileAddress.rfind("/")+1:len(FileAddress)]
#                     if (LeftPattern in FileName) and (RightPattern in FileName):
#                         FileLabel = FileName[FileName.find(LeftPattern)+2:FileName.find(RightPattern)]
#                         Matching_List = []
#                         Already_Seen = False
#                         
#                         for item in Combined_List:
#                             for family in item:
#                                 if family == FileAddress:
#                                     Already_Seen = True
#                         
#                         if Already_Seen == False:
#                             for FileLocation in FilesList:
#                                 FilNam = FileLocation[FileLocation.rfind("/")+1:len(FileLocation)]
#                                 if FileLabel in FilNam:
#                                     Matching_List.append(FileLocation)
#                                    
#                             Combined_List.append(tuple(Matching_List))
#                         
#                     else:
#                         Single_List.append((FileAddress,))
#                 
#                 self.ListFiles = tuple(Single_List + Combined_List)    
#                 
#                 if len(self.ListFiles) == 1:
#                     self.Tags.append('SingleFile')
#             else:
#                 FoldersList = []
#                 ArchivesList = []
#                 
#                 if type(FilePattern) is not list:
#                     myPatternList = [FilePattern]
#                      
#                 else:
#                     myPatternList = FilePattern
#                   
#                 for Root, Dirs, Archives in walk(FilesFolders):
#                     ValidArchives = []
#                     for Archive in Archives:
#                         Meets_one_Pattern = False
#                         for i in range(len(myPatternList)):
#                             if (myPatternList[i] in Archive):
#                                 if "~" not in Archive:                   
#                                     Meets_one_Pattern = True
#                                     print '--- File', Archive, '@', Dirs, Root               
#                 
#                         if Meets_one_Pattern:        
#                             if Root.endswith("/"):
#                                 FinalName = Root
#                             else:
#                                 FinalName = Root + "/"
#                 
#                             if FinalName in FoldersList:
#                                 ValidArchives.append(Archive)
#                             else:
#                                 FoldersList.append(FinalName)
#                                 ValidArchives.append(Archive)
#                      
#                     if len(ValidArchives) > 0:
#                         ValidArchives.sort()
#                         ArchivesList.append(ValidArchives)
#             
#                 if Sort_Output == 'Alpha':
#                     FoldersList, ArchivesList = zip(*sorted(zip(FoldersList, ArchivesList), key=itemgetter(0), reverse=False))
#                     
#                 self.ListFiles = (FoldersList, ArchivesList)    
# 
#         return self.ListFiles

    def load_dataframe(self, dataframe_address):
        
        df = read_pickle(dataframe_address)
    
        return df

    def save_dataframe(self, dataframe, dataframe_address):

        dataframe.to_pickle(dataframe_address)

        return

    def FindAndOrganize_dazer(self, FilePattern, FilesFolders, unpack = False, CheckComputer=False, Sort_Output = 'Alpha'):
        #THIS CODE IS NECESSARY SINCE DAZER DEPENDS ON THE [FOLDERLIST, ARCHIVELIST STRUCTURE] NEEDS TO BE UPDATED
        
        #Moving from an absolute to a relative structure
        if CheckComputer == True:
            FilesFolders = self.RootFolder + FilesFolders
        
        #Checking for arguments in terminal
        self.Arguments_Checker(argv)
                             
        #Command from Terminal
        if self.Arguments_Check == True:
            self.Tags.append('Terminal')
            self.Argument_Handler()
            
            if len(self.Command_Launching_Files) == 1:
                self.Tags.append('SingleFile')
            
            self.ListFiles = self.Command_Launching_Files
     
        #Command from eclipse
        if unpack == True:
            self.Tags.append('Editor')
            
            if unpack == True:     
                FoldersList = []
                ArchivesList = []
                 
                if type(FilePattern) is not list:
                    myPatternList = [FilePattern]
                      
                else:
                    myPatternList = FilePattern
                   
                for Root, Dirs, Archives in walk(FilesFolders):
                    ValidArchives = []
                    for Archive in Archives:
                        Meets_one_Pattern = False
                        for i in range(len(myPatternList)):
                            if (myPatternList[i] in Archive):
                                if "~" not in Archive:                   
                                    Meets_one_Pattern = True
                 
                        if Meets_one_Pattern:        
                            if Root.endswith("/"):
                                FinalName = Root
                            else:
                                FinalName = Root + "/"
                 
                            if FinalName in FoldersList:
                                ValidArchives.append(Archive)
                            else:
                                FoldersList.append(FinalName)
                                ValidArchives.append(Archive)
                      
                    if len(ValidArchives) > 0:
                        ValidArchives.sort()
                        ArchivesList.append(ValidArchives)
             
                if Sort_Output == 'Alpha':
                    FoldersList, ArchivesList = zip(*sorted(zip(FoldersList, ArchivesList), key=itemgetter(0), reverse=False))
                     
                self.ListFiles = (FoldersList, ArchivesList)    
 
        return self.ListFiles
        
    def extract_traces_statistics(self, traces_list = None):
        
        self.statistics_dict = OrderedDict()
                
        #If no list input we extract all the traces from the analysis        
        if traces_list == None:
            traces_list = self.traces_list
                
        for trace in traces_list:
            self.statistics_dict[trace] = OrderedDict()
            
            for stat in self.pymc_stats_keys:
                self.statistics_dict[trace][stat] = self.dbMCMC.trace(trace).stats()[stat]                
            
            Trace_array = self.pymc_database.trace(trace)[:] 
            self.statistics_dict[trace]['16th_p'] = percentile(Trace_array, 16)
            self.statistics_dict[trace]['84th_p'] = percentile(Trace_array, 84)
        
        return self.statistics_dict    
    
    def query_yes_no(self, question, default="yes"):
        """Ask a yes/no question via raw_input() and return their answer.
     
        "question" is a string that is presented to the user.
        "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).
     
        The "answer" return value is one of "yes" or "no".
        """
        valid = {"yes":True,   "y":True,  "ye":True,
                 "no":False,     "n":False}
        if default == None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)
     
        while True:
            stdout.write(question + prompt)
            choice = raw_input().lower()
            if default is not None and choice == '':
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                stdout.write("Please respond with 'yes' or 'no' "\
                                 "(or 'y' or 'n').\n")    
    
    
        
    