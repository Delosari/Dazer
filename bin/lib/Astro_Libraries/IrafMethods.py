from pyraf                      import iraf
from collections                import OrderedDict
from os                         import rename, remove
from os.path                    import isfile
from numpy                      import where
from lib.Astro_Libraries  import cosmics
from lib.Astro_Libraries  import f2n

class Iraf_Task_Configuration():
    
    def __init__(self):
                    
        #Attribute with the data type
        self.SpectraTreating = 'WHT'
        
        #Method with the employed fits header keys
        self.WHT_longslit_v0()
            
    def WHT_longslit_v0(self):
        
        self.Airmass_Key                = 'AIRMASS'
        self.Exptime_Key                = 'EXPTIME'
        
        self.GlobalSensitiviyCurve      = 'a_sen_wideCombAll_PYRAF'
        self.GlobalStandard_Ouput       = 'a_std_wideCombAll_PYRAF'
        self.Standard_Stars_Folder      = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/StandardStars_Calibration/'
        
        self.ArcCalibrationMap_Address  = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/Arc_Calibration/' + 'CuNe+CuAr.dat'
        self.ArcLampObject_Extension    = '_arc.fits'

        self.ReddeningFolder            = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/'
        self.ReddeningFile              = 'wht_extinction_file.dat'

    def Lacosmic_WHTData(self, InputFile, OutputFile, Fits_Folder, CodeName, gain, readnoise, sigclip, verbose):

        #Load the data
        (array, header) = cosmics.fromfits(Fits_Folder + InputFile, verbose = verbose)
        upsample = 1

        #Create the task object
        c = cosmics.cosmicsimage(array, pssl = 0.0, gain=1.16, readnoise=5, sigclip = 16, verbose = verbose)

        #Run the task
        c.run(maxiter = 4, verbose = verbose)
        
        #Proceed to generate the images
        
        #-- 1) The raw input array
        im = f2n.f2nimage(c.getrawarray(), verbose=False)
        im.setzscale()
        im.makepilimage("log")
        im.upsample(upsample)
        im.writetitle("Input image", colour = (0,255,0))
        im.tonet(Fits_Folder + CodeName + "_1_raw.png")        
        
        #-- 2) The raw input array
        #-- We output a list of the positions of detected cosmic ray hits.
        #-- This is made on purpose to be fed into f2n's drawstarslist :
        labeldict = c.labelmask()
        im = f2n.f2nimage(c.getrawarray(), verbose=False)
        im.setzscale()
        im.makepilimage("log")
        im.drawmask(c.getsatstars(), colour=(255, 0, 255))
        im.upsample(upsample)
        im.drawstarslist(labeldict, colour=(255,0,0))
        im.writetitle("Cosmic ray hits", colour = (0,255,0))
        im.tonet(Fits_Folder + CodeName + "_2_labels.png")       

        #-- 3) Tne png with the precise mask in green and a wider version in blue :
        im = f2n.f2nimage(c.getrawarray(), verbose=False)
        im.setzscale()
        im.makepilimage("log")
        im.drawmask(c.getdilatedmask(size=5), colour=(0, 0, 255))
        im.drawmask(c.getmask(), colour=(0, 255, 0))
        im.upsample(upsample)
        im.writetitle("Mask", colour = (0,255,0))
        im.tonet(Fits_Folder + CodeName + "_3_mask.png")

        #-- 4) And of course one png with the clean array :
        im = f2n.f2nimage(c.getcleanarray(), verbose=False)
        im.setzscale()
        im.makepilimage("log")
        im.upsample(upsample)
        im.writetitle("Cleaned image", colour = (0,255,0))
        im.tonet(Fits_Folder + CodeName + "_4_clean.png")
        
        #Save the FitsFile (CHECK THE HEADER):
        cosmics.tofits(Fits_Folder + OutputFile, c.cleanarray, header)
        
        return

    def IdentifyAttributes(self, ArcFile, Fits_Folder):
        
        identify_conf               = OrderedDict()
        
        if self.SpectraTreating == 'WHT':
                        
            identify_conf               = OrderedDict()
            identify_conf['images']     = Fits_Folder + ArcFile 
            identify_conf['section']    = 'middle column'
            identify_conf['niterat']    = 1
            identify_conf['fwidth']     = 7
            identify_conf['coordlist']  = self.ArcCalibrationMap_Address
            identify_conf['match']      = 10
            identify_conf['fwidth']     = 7.5
            identify_conf['cradius']    = 5
            identify_conf['threshold']  = 10
            identify_conf['minsep']     = 2

        return identify_conf

    def ReidentifyAttributes(self, ArcFile, Fits_Folder):
                
        reidentify_conf               = OrderedDict()
        
        if self.SpectraTreating == 'WHT':
                
            reidentify_conf['referenc']   = Fits_Folder + ArcFile 
            reidentify_conf['images']     = Fits_Folder + ArcFile
            
            reidentify_conf['section']    = 'middle column'
            reidentify_conf['coordlist']  = self.ArcCalibrationMap_Address
            
            reidentify_conf['interac']    = 'yes'
            reidentify_conf['overrid']    = 'yes'
            reidentify_conf['trace']      = 'no'
            reidentify_conf['nlost']      = 5
            reidentify_conf['threshold']  = 10
            reidentify_conf['match']      = 10
            reidentify_conf['verbose']    = 'yes'
            
        return reidentify_conf
    
    def FitcoordsAttributes(self, ArcFile, Fits_Folder):
        
        fitcoords_conf              = OrderedDict()
        
        if self.SpectraTreating == 'WHT':
            
            fitcoords_conf['images']        = Fits_Folder + ArcFile
            fitcoords_conf['interative']    = 'yes'
            fitcoords_conf['fitname']       = '""'
            fitcoords_conf['interative']    = 'yes'

        return fitcoords_conf

    def TransformAttributes(self, InputFile, OutputFile, Fits_Folder, ArcFile):
        
        transform_conf              = OrderedDict()
        
        if self.SpectraTreating == 'WHT':
            
            transform_conf['input']         = Fits_Folder + InputFile
            transform_conf['output']        = Fits_Folder + OutputFile
            transform_conf['fitname']       = Fits_Folder + ArcFile
            transform_conf['flux']          = 'no'

        return transform_conf

    def ApallAttributes(self, InputFile, OutputFile, Fits_Folder, line, nsum, lower_lim, upper_lim, extras):
        
        if self.SpectraTreating == 'WHT':
            
            apall_conf                   = OrderedDict()
            
            if lower_lim == None:
                lower_lim = 20

            if upper_lim == None:
                upper_lim = 20
                
            if nsum == None:
                nsum = 10

            if lower_lim == None:
                lower_lim = 20
            
            apall_conf['input']         = Fits_Folder + InputFile
            apall_conf['output']        = Fits_Folder + OutputFile
            apall_conf['width']         = 10
            apall_conf['nsum']          = nsum
            apall_conf['maxsep']        = 100
            apall_conf['t_niter']       = 1
            apall_conf['backgro']       = 'median'
            apall_conf['weights']       = 'variance'
            apall_conf['gain']          = 'gain'
            apall_conf['readnoi']       = 'readnoi'
            apall_conf['extras']        = extras
            apall_conf['interactive']   = 'yes'
            apall_conf['find']          = 'yes'
            apall_conf['nfind']         = 1  
            apall_conf['trace']         = 'yes'
            apall_conf['resize']        = 'no'  
            apall_conf['edit']          = 'yes'  
            
            return apall_conf
            
    def StandardAttributes(self, InputFile, OutputFile, starfile, Fits_Folder, airmass_value, exptime_value, Bandwidth, Bandsep ):
        
        if self.SpectraTreating == 'WHT':
            
            print '\n--- Calibration file used', self.Standard_Stars_Folder + starfile
            
            Observatory                     = 'lapalma'
            ExtinctionFileAddress           = self.ReddeningFolder + self.ReddeningFile     

            standard_conf                   = OrderedDict()
            standard_conf['input']          = Fits_Folder + InputFile
            standard_conf['output']         = Fits_Folder + OutputFile
            standard_conf['samestar']       = 'yes'
            standard_conf['beam_switch']    = 'no'
            standard_conf['apertures']      = '""' 
            standard_conf['bandwidth']      = Bandwidth
            standard_conf['bandsep']        = Bandsep
            standard_conf['fnuzero']        = 3.68000e-20
            standard_conf['extinction']     = ExtinctionFileAddress
            standard_conf['caldir']         = self.Standard_Stars_Folder
            standard_conf['observatory']    = Observatory
            standard_conf['interact']       = 'yes'
#           standard_conf['graphic']        = 'stdgraph'
#           standard_conf['cursor']         = '""'
            standard_conf['star_name']      = starfile              #MUST NOT INCLUDE EXTENSION IN THE STAR CALIBRATION FILE, NEITHER CAPITALS!!!!
            standard_conf['airmass']        = airmass_value
            standard_conf['exptime']        = exptime_value
            standard_conf['answer']         = 'yes'

        return standard_conf

    def SensfuncAttributes(self, standardfile_input, sensitivitycurve_output, Fits_Folder):
    
        if self.SpectraTreating == 'WHT':
                        
            Observatory                     = 'lapalma'
            ExtinctionFileAddress           = self.ReddeningFolder + self.ReddeningFile    
            
            sensfunc_conf                   = OrderedDict()
            sensfunc_conf['standards']      = Fits_Folder + standardfile_input
            sensfunc_conf['sensitivity']    = Fits_Folder + sensitivitycurve_output
            sensfunc_conf['apertures']      = '""'
            sensfunc_conf['ignoreaps']      = "yes"
#           sensfunc_conf['logfile']        = 00000000
            sensfunc_conf['extinction']     = ExtinctionFileAddress
#           sensfunc_conf['newextinction']  = 00000000
            sensfunc_conf['observatory']    = Observatory
            sensfunc_conf['function']       = 'spline3'
            sensfunc_conf['interactive']    = 'yes'
            sensfunc_conf['graphs']         = 'sri'
#           sensfunc_conf['marks']          = '"plus cross box"'
#           sensfunc_conf['colors']         = '"2 1 3 4"'
            sensfunc_conf['answer']         = 'yes'

        return sensfunc_conf

    def CalibrationAttributes(self,  InputFile, OutputFile, senstivityCurve, Fits_Folder, airmass_value, exptime_value, fnu_units):
        
        if self.SpectraTreating == 'WHT':
            
            Observatory                     = 'lapalma'
            ExtinctionFileAddress           = self.ReddeningFolder + self.ReddeningFile   
                           
            calibrate_conf                  = OrderedDict()
            calibrate_conf['input']         = Fits_Folder + InputFile
            calibrate_conf['output']        = Fits_Folder + OutputFile
            calibrate_conf['extinct']       = 'yes'
            calibrate_conf['flux']          = 'yes'
            calibrate_conf['extinction']    = ExtinctionFileAddress
            calibrate_conf['observatory']   = Observatory
            calibrate_conf['ignoreaps']     = 'yes'
            calibrate_conf['sensitivity']   = Fits_Folder + senstivityCurve
            calibrate_conf['fnu']           = fnu_units
            calibrate_conf['airmass']       = airmass_value
            calibrate_conf['exptime']       = exptime_value
            calibrate_conf['mode']          = 'ql'

            return calibrate_conf

    def ContinuumAttributes(self,  InputFile, OutputFile, Fits_Folder):
        
        continuum_conf      = OrderedDict()
        
        continuum_conf['input']= Fits_Folder + InputFile
        continuum_conf['output']= Fits_Folder + OutputFile
        
        return continuum_conf

    def SplotAttributes(self, InputFile, Fits_Folder, xmin = 'INDEF', xmax = 'INDEF', ymin = 'INDEF', ymax = 'INDEF'):
        
        splot_conf          = OrderedDict()
                
        splot_conf['images']= Fits_Folder + InputFile
        splot_conf['xmin']  = xmin
        splot_conf['xmax']  = xmax
        splot_conf['ymin']  = ymin
        splot_conf['ymax']  = ymax
        
        return splot_conf

    def SarithAttributes(self, InputFile1, InputFile2, OuputFile, Fits_Folder, operation, wmin = 'INDEF', wmax = 'INDEF', rebin='no', verbose = 'yes'):
        
        sarith_conf             = OrderedDict()
        
        sarith_conf['input1']   = Fits_Folder + InputFile1
        sarith_conf['op']       = operation
        sarith_conf['input2']   = Fits_Folder + InputFile2
        sarith_conf['output']   = Fits_Folder + OuputFile
        sarith_conf['w1']       = wmin
        sarith_conf['w2']       = wmax
        sarith_conf['rebin']    = rebin
        sarith_conf['verbose']  = verbose
        
        return sarith_conf

    def ScopyAttributes(self, InputFile, OuputFile, Fits_Folder, w1 = 'INDEF', w2 = 'INDEF', rebin='no', clobber = 'yes', verbose_output = 'yes'):
        
        scopy_conf              = OrderedDict()
        scopy_conf['input']     = Fits_Folder + InputFile
        scopy_conf['output']    = Fits_Folder + OuputFile
        scopy_conf['w1']        = w1
        scopy_conf['w2']        = w2
        scopy_conf['rebin']     = rebin
        scopy_conf['clobber']   = clobber
        scopy_conf['verbose']   = verbose_output

        return scopy_conf

    def WpspectextAttributes(self, InputFile, OuputFile, Fits_Folder, format='%0.1f'):
        # Export fits file into a text file
        # wspectext objXXblue_SlInput.fits objXXblue_SlInput.txt header- wformat="%0.1f" 
        
        wspectext_conf             = OrderedDict()
        
        return wspectext_conf

    def rspectextAttributes(self, InputFile, OutputFile, Fits_Folder):
        
        task_conf             = OrderedDict()
        task_conf['input']    = Fits_Folder + InputFile
        task_conf['output']   = Fits_Folder + OutputFile        
        task_conf['dtype']    = 'none'        

        return task_conf

    def DopcorAttributes(self, InputFile, OutputFile, Fits_Folder, z, FluxCorrection='yes', printoutput='yes'):
        
        dopcor_conf              = OrderedDict()
        dopcor_conf['input']     = Fits_Folder + InputFile
        dopcor_conf['output']    = Fits_Folder + OutputFile        
        dopcor_conf['redshift']  = z
        dopcor_conf['flux']      = FluxCorrection
        dopcor_conf['verbose']   = printoutput

        return dopcor_conf
    
    def DereddenAttributes(self, InputFile, OutputFile, Fits_Folder, Type, EBV, R):
        #"deredden "+ InputName+"_z " + InputName+"_z_EBV" + " " + "R=3.2 Type=E value="+myEBV

        deredden_conf               = OrderedDict()
        deredden_conf['input']      = Fits_Folder + InputFile
        deredden_conf['output']     = Fits_Folder + OutputFile        
        deredden_conf['type']       = Type
        deredden_conf['value']      = EBV
        deredden_conf['R']          = R
        
        return deredden_conf
           
class Iraf_Tools():
    
    def __init__(self):
        self.Reddening_Folder = 'None'
        self.IRAF_pypelineFolder = '/home/vital/Workspace/Thesis_Pipeline/SpectraReduction/'
    
    def outputNameGenerator(self, InputFile, Suffix):
        OutputFile =  InputFile[0:InputFile.find('.fits')] + '_' + Suffix + '.fits'
        return OutputFile
    
    def get_HeaderValue(self, Item, Fit_Address, Digitalize = False):
        Header = iraf.imhead(Fit_Address, long='yes', Stdout=1)
        for i in range(len(Header)):
            key = Header[i].replace(' ','').split('=')[0]
            if key == Item:
                value = Header[i].replace(' ','').split('=')[1]
                value = value.split('/')[0]
                if Digitalize == True:
                    return float(value)
                else:
                    return value 
    
    def printIrafCommand(self, Task, Attributes, printindividually = True):
        
        keys_attrib     = Attributes.keys()
        values_attrib   = Attributes.values()
    
        command_String  = Task
        
        if printindividually == False:
            for i in range(len(keys_attrib)):
                command_String = command_String + ' ' + keys_attrib[i] + '=' + str(values_attrib[i])
        else: 
            print '--Task atributtes:'
            for i in range(len(keys_attrib)):
                command_String = command_String + ' ' + keys_attrib[i] + '=' + str(values_attrib[i])
                print keys_attrib[i] + '=' + str(values_attrib[i])
                
        return command_String
    
    def StarnameFromFileName(self, FitsFile):
    
        #This function assumes a structure: 'stdBD17_wide_cr_f_a_bg_s'
        StarName = FitsFile[3:FitsFile.find('_')].lower()
        
        return StarName

    def getCalibrationFile(self, StarCode):
        
        if (StarCode == 'bd17') or (StarCode == 'sp2209+178') or (StarCode == 'sp2209178'):
            StarFile = 'bd17_d4708stisnic006_nobandpass'
            Bandwidth = '40'
            Bandsep  = '40'
 
        elif StarCode == 'wolf':
            StarFile = 'wolf_oke1974_40a'
            Bandwidth = 'INDEF'
            Bandsep  = 'INDEF'
                        
        elif StarCode == 'g191':
            StarFile = 'g191_b2bstisnic006_nobandpass'
            Bandwidth = '40'
            Bandsep  = '40'
          
        elif (StarCode == 'f34') or (StarCode == 'feige34'):
            StarFile = 'feige34_stis004_nobandpass'
            Bandwidth = '40'
            Bandsep  = '40'

        elif StarCode == 'hd84937':
            StarFile = 'hd84937_oke1983_40a'
            Bandwidth = 'INDEF'
            Bandsep  = 'INDEF'
            
        elif StarCode == 'g158':
            StarFile = 'g158_oke1990_1a'
            Bandwidth = '40'
            Bandsep  = '40' 

        elif (StarCode == 'hd19445') or (StarCode == 'sp0305+261') or (StarCode == 'sp0305261'):
            StarFile = 'hd19445_oke1983_40a'
            Bandwidth = 'INDEF'
            Bandsep  = 'INDEF'

        if (StarCode == 'bd28') or (StarCode == 'bd+284211') or (StarCode == 'sp2148286'):
            StarFile = 'bd28_d4211stis004_nobandpass'
            Bandwidth = '40'
            Bandsep  = '40'
            
        if (StarCode == 'bd33') or (StarCode == 'bd+332642') or (StarCode == 'sp1550330'):
            StarFile = 'bd33_d2642004_nobandpass'
            Bandwidth = '40'
            Bandsep  = '40'
                    
        return StarFile, Bandwidth, Bandsep
    
    def getTextHeaderRow(self, FileFolder, FileName):
    
        File        = open(FileFolder + FileName)
        FileData    = File.readlines()
        File.close()
            
        HeaderLine  = FileData[0].split() 
            
        return HeaderLine
    
    def movefiletofolder(self, FileFolder, FileName):
        
        if isfile(self.IRAF_pypelineFolder + FileName):
            rename(self.IRAF_pypelineFolder + FileName, FileFolder + FileName)
        
        return
      
    def getTaskConfiguration(self, Object, ObjList, ConfList, ConfAttributes):
       
        #This assumes a scheme with the type: Blue_5300_7000,Red_8400_10000  
        #WE SHOULD ADD A STEP TO MAKE SURE THE COLUMN IS READ IN SMALL LETTERS
        if ConfAttributes['Operation'] == 'Wavelength trimming':
            
            Object_index    = where(ObjList == Object.lower()) 
            Conf            = ConfList[Object_index]
            
            if 'Blue' in ConfAttributes['FileFolder']:
                Color_ind = 0
        
            elif 'Red' in ConfAttributes['FileFolder']:
                Color_ind = 1
            
            #WAR
            Color_Conf  = Conf[0].split(',')[Color_ind]
            
            Wmin = int(Color_Conf[Color_Conf.find('_') + 1:Color_Conf.rfind('_')])
            Wmax = int(Color_Conf[Color_Conf.rfind('_')+1:len(Color_Conf)])
            
            return Wmin, Wmax
        
        if ConfAttributes['Operation'] == 'Telluric correction':
            
            Object_index    = where(ObjList == Object.lower()) 
            Conf            = ConfList[Object_index][0]

            return Conf
  
        if ConfAttributes['Operation'] == 'Redshift correction':
            
            Object_index    = where(ObjList == Object.lower()) 
            Conf            = ConfList[Object_index]
            
            if 'Blue' in ConfAttributes['FileFolder']:
                Color_ind = 0
        
            elif 'Red' in ConfAttributes['FileFolder']:
                Color_ind = 1
                
            Color_Conf  = Conf[0].split(',')[Color_ind]
            
            z = float(Color_Conf[Color_Conf.find('_') + 1:len(Color_Conf)])
            
            return z 

    def getArmColor(self, FitAddress):
        
        #THIS COULD BE DONE WITH PYFITS  
        ArmColor = self.get_HeaderValue('Arm_value', Fit_Address = FitAddress, Digitalize = False)
        print 'the armcolor is', ArmColor
        Color = None
        if ArmColor == 'Blue_arm':
            Color = 'Blue'
        if ArmColor == 'Red_arm':
            Color = 'Red'
            
        return Color
   
    def deleteIrafFiles(self, IrafOutputAddress):
        
        if isfile(IrafOutputAddress):
            remove(IrafOutputAddress)    

class Pyraf_Workflow(Iraf_Task_Configuration, Iraf_Tools):
    
    def __init__(self, dataType):
        
        Iraf_Task_Configuration.__init__(self, dataType)
        Iraf_Tools.__init__(self)

    def Catalogue_Reduction_Format(self, Reduction_Operation, ObservationNight = '2011_12_Night2', Folder = None):
        
        self.FolderData                     = self.Catalogues_Folders(ObservationNight)
        
        self.PropertiesDictionary           = {}

        if Reduction_Operation == 'Cosmic_Ray_Removal':
            
            #This is because we might have combined images or not: '_c.fits' or not: '_nc.fits'
            self.PropertiesDictionary['Combined_Images'] = 'c'
    
            return self.FolderData, self.PropertiesDictionary
        
        if Reduction_Operation == 'Wavelength_Calibration':
            
            #This is because we might have applied the cosmic rays: '_r.fits' or not: '_nr.fits'
            self.PropertiesDictionary['Suffix_CosmicRaysRemoval'] = 'r'
    
            return self.FolderData, self.PropertiesDictionary

        if Reduction_Operation == 'Background_Extraction':

            self.PropertiesDictionary['Suffix_ArcCalibrated'] = '_a'
    
            return self.FolderData, self.PropertiesDictionary
    
        if Reduction_Operation == 'Std_FluxCalibration':
            
            if ObservationNight             == '2008_01_Night1':
                self.Stars_Extract_suffix   = 'wide_f_cr_a_bg.0001.fits'
            if ObservationNight             == '2008_01_Night2':
                self.Stars_Extract_suffix   = 'wide_f_cr_a_bg.0001.fits'            

            if ObservationNight             == '2009_07_Night1':
                self.Stars_Extract_suffix   = 'wide_f_cr_a_bg.0001.fits'
            if ObservationNight             == '2009_07_Night2':
                self.Stars_Extract_suffix   = 'wide_f_cr_a_bg.0001.fits'            

            if ObservationNight             == '2011_09_Night1':
                self.Stars_Extract_suffix   = 'wide_f_cr_a_bg.0001.fits'
            if ObservationNight             == '2011_09_Night2':
                self.Stars_Extract_suffix   = 'wide_f_cr_a_bg.0001.fits'
    
            if ObservationNight             == '2012_09_Night1':
                self.Stars_Extract_suffix   = 'wide_cr_f_a_bg.0001.fits'
            if ObservationNight             == '2012_09_Night2':
                self.Stars_Extract_suffix   = 'wide_f_cr_a_bg.0001.fits'   

            if ObservationNight             == '2011_12_Night1':
                self.Stars_Extract_suffix   = '_wide_cr_f_a_bg.0001.fits'
            if ObservationNight             == '2011_12_Night2':
                self.Stars_Extract_suffix   = '_wide_f_cr_a_bg.0001.fits'   
    
            return self.Stars_Extract_suffix, self.FolderData

        if Reduction_Operation == 'Obs_FluxCalibration':

            if ObservationNight             == '2008_01_Night1':
                self.Stars_Extract_suffix   = '_c_f_cr_a_bg.0001.fits'
            if ObservationNight             == '2008_01_Night2':
                self.Stars_Extract_suffix   = '_c_f_cr_a_bg.0001.fits'            

            if ObservationNight             == '2009_07_Night1':
                self.Stars_Extract_suffix   = '_c_f_cr_a_bg.0001.fits'
            if ObservationNight             == '2009_07_Night2':
                self.Stars_Extract_suffix   = '_c_f_cr_a_bg.0001.fits'            

            if ObservationNight             == '2011_09_Night1':
                self.Stars_Extract_suffix   = '_c_f_cr_a_bg.0001.fits'
            if ObservationNight             == '2011_09_Night2':
                self.Stars_Extract_suffix   = '_s_c_f_cr_a_bg.0001.fits'
    
            if ObservationNight             == '2011_12_Night1':
                self.Stars_Extract_suffix   = '_s_c_cr_f_a_bg.0001.fits'
            if ObservationNight             == '2011_12_Night2':
                self.Stars_Extract_suffix   = '_c_f_cr_a_bg.0001.fits'   
                
            if ObservationNight             == '2011_12_Night1':
                self.Stars_Extract_suffix   = '_s_c_cr_f_a_bg.0001.fits'
            if ObservationNight             == '2011_12_Night2':
                self.Stars_Extract_suffix   = '_c_f_cr_a_bg.0001.fits'                   
                
                
            return self.Stars_Extract_suffix, self.FolderData

        if Reduction_Operation == 'Obs_TrimCalibrated':
                
            self.PropertiesDictionary['Suffix_FluxCalibrated']      = '_fxW.fits'

            return self.FolderData, self.PropertiesDictionary

        if Reduction_Operation == 'Std_Clean':
            #WARNING: telluric correction with the non-flux calibrated spectrum
            if ObservationNight             == '2008_01_Night1':
                self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']      = '_narrow_f_cr_a_bg.0001.fits'
            if ObservationNight             == '2008_01_Night2':
                self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']      = '_narrow_f_cr_a_bg.0001.fits'            

            if ObservationNight             == '2009_07_Night1':
                self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']      = 'wide_f_cr_a_bg.0001_fxW_t.fits'
            if ObservationNight             == '2009_07_Night2':
                self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']      = 'wide_f_cr_a_bg.0001_fxW_t.fits'            

            if ObservationNight             == '2011_09_Night1':
                self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']      = '_wide_f_cr_a_bg.0001.fits'
            if ObservationNight             == '2011_09_Night2':
                self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']      = '_narrow_f_cr_a_bg.0001.fits'
    
            if ObservationNight             == '2011_12_Night1':
                self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']      = '_narrow_cr_f_a_bg.0001.fits'
            if ObservationNight             == '2011_12_Night2':
                self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']      = '_narrow_cr_f_a_bg.0001.fits'   


            self.PropertiesDictionary['Suffix_CleanCalibration']                = 'clean'  
            
            return self.FolderData, self.PropertiesDictionary

        if Reduction_Operation == 'deredden':

            self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']          = '_fxW_t.fits'
            self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed_Telluric'] = '_fxW_t_tell.fits'            

        if Reduction_Operation == 'dopcor':
            
            self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed']          = '_fxW_t.fits'
            self.PropertiesDictionary['Suffix_FluxCalibrated_Trimmed_Telluric'] = '_fxW_t_tell.fits'

            return self.FolderData, self.PropertiesDictionary

    def Catalogues_Folders(self, Catalogue_Keyword):
                    
#         if Catalogue_Keyword            == 'Night1_2011_12':
#             FolderData             = '/home/vital/Desktop/Flux_Calibration_Steps/Night1_2011_12/'                    
#                                
#         if Catalogue_Keyword            == 'Flux_Night2_2012':
#             FolderData             = '/home/vital/Desktop/Flux_Calibration_Steps/Night2_2012_Blue_Pyraf/'
# 
#         elif Catalogue_Keyword           == 'Flux_Night1_2012':
#             FolderData             = '/home/vital/Desktop/Flux_Calibration_Steps/PyRaf_Testing/'
# 
#         elif Catalogue_Keyword           == 'Night2_2009_07':
#             FolderData             = '/home/vital/Desktop/Flux_Calibration_Steps/Night2_2009_07/'
# 
#         elif Catalogue_Keyword           == 'Night2_2008_01':
#             FolderData             = '/home/vital/Desktop/Flux_Calibration_Steps/Night2_2008_01/'

        CatalogueRoot   = '/home/vital/Astrodata/WHT_HIIGalaxies_FluxCalibration/'
        FolderData      = CatalogueRoot + Catalogue_Keyword + '/'
               
        return FolderData

    def LacosmicTask(self, InputFile, OutputFile, Fits_Folder, CodeName, gain, readnoise, sigclip, verbose= True, Suffix = 'c'):
 
        #Incase no output name is given, we generate one with the provided suffix
        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile, Suffix)  
 
        if self.SpectraTreating == 'WHT':
        
            self.Lacosmic_WHTData(InputFile, OutputFile, Fits_Folder, CodeName, gain, readnoise, sigclip, verbose)
        
        return OutputFile

    def IdentifyTask(self, ArcFile, Fits_Folder):

        iraf.noao(_doprint=0)
        iraf.twodspec(_doprint=0)         
        iraf.longslit(_doprint=0)
        
        IdentifyConf = self.IdentifyAttributes(ArcFile, Fits_Folder)

        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('identify', IdentifyConf)
        print '--- Using the command'
        print Command
        
        iraf.twodspec.longslit.identify(**IdentifyConf)
        
        return

    def ReidentifyTask(self, ArcFile, Fits_Folder, ):

        iraf.noao(_doprint=0)
        iraf.twodspec(_doprint=0)         
        iraf.longslit(_doprint=0)
        
        ReidentifyConf = self.ReidentifyAttributes(ArcFile, Fits_Folder)

        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('reidentify', ReidentifyConf)
        print '--- Using the command'
        print Command
        
        iraf.twodspec.longslit.reidentify(**ReidentifyConf)

        return
        
    def FitcoordsTask(self, ArcFile, Fits_Folder):

        iraf.noao(_doprint=0)
        iraf.twodspec(_doprint=0)         
        iraf.longslit(_doprint=0)
        
        FitcoordsConf = self.FitcoordsAttributes(ArcFile, Fits_Folder)

        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('fitcords', FitcoordsConf)
        print '--- Using the command'
        print Command
        
        iraf.twodspec.longslit.fitcoords(**FitcoordsConf)

        return
    
    def TransformTask(self, InputFile, OutputFile, Fits_Folder, ArcFile, Suffix = 'a'):

        iraf.noao(_doprint=0)
        iraf.twodspec(_doprint=0)         
        iraf.longslit(_doprint=0)

        #Incase no output name is given, we generate one with the provided "preffix" (The defaul format is a_std_wolf.dat)
        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile, Suffix)  

        TransConf = self.TransformAttributes(InputFile, OutputFile, Fits_Folder, ArcFile)

        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('transform', TransConf)
        print '--- Using the command'
        print Command

        iraf.twodspec.longslit.transform(**TransConf)

        return OutputFile

    def ApallTask(self, InputFile, OutputFile, Fits_Folder, line, nsum = None, lower_lim = None, upper_lim = None, extras = 'yes', Suffix = 'bg'):
        
        iraf.noao(_doprint=0)
        iraf.twodspec(_doprint=0)         
        iraf.apextract(_doprint=0)
        
        #Incase no output name is given, we generate one with the provided "preffix" (The defaul format is a_std_wolf.dat)
        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile + '[1]', Suffix)  
            
        ApallConf = self.ApallAttributes(InputFile, OutputFile, Fits_Folder, line, nsum, lower_lim, upper_lim, extras)  
      
        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('apall', ApallConf)
        print '--- Using the command'
        print Command
      
#         p = Popen(iraf.twodspec.apextract.apall(**ApallConf), shell=True, stdout=PIPE, stderr=STDOUT)
        
        #Run the task        
        iraf.twodspec.apextract.apall(**ApallConf)
      
        return OutputFile
      
    def StandardTask(self, InputFile, OutputFile, FitsFolder, airmass_value, exptime_value):
   
        iraf.noao(_doprint=0)
        iraf.onedspec(_doprint=0)               
        
        #From the fits file determine which is the star being treated
        StarName = self.StarnameFromFileName(InputFile)
        
        #Get the corresponding Calibration file
        CalibrationFile, Bandwidth, Bandsep   = self.getCalibrationFile(StarName)
        
        #Incase no output name is given, we generate one with the provided "preffix" (The defaul format is a_std_wolf.dat)
        if OutputFile == None:
            OutputFile    = 'a_std_' + StarName
        
        #Prepare dictionary with the configuration for the tasks           
        Standard_conf_Comb  = self.StandardAttributes(InputFile, OutputFile, CalibrationFile, FitsFolder, airmass_value, exptime_value, Bandwidth, Bandsep)

        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('standard', Standard_conf_Comb)
        print '--- Using the command'
        print Command
        
        #Run the task        
        iraf.onedspec.standard(**Standard_conf_Comb)
        
        return OutputFile 
 
    def SensfuncTask_IndividualStarCurve(self, Input_Fits, Output_SensitivityCurve, FitsFolder):

#         iraf.noao(_doprint=0)
#         iraf.onedspec(_doprint=0)     

        #From the fits file determine which is the star being treated
        StarName = self.StarnameFromFileName(Input_Fits)
        
        #Get file generated by Standard task. It is assumed to have a format 
        Input_StandardFile = 'a_std_' + StarName
            
        #Incase no output name is given, we generate one with the provided "preffix" (The defaul format is a_std_wolf.dat)
        if Output_SensitivityCurve == None:
            Output_SensitivityCurve    = 'a_sen_' + StarName       

        Sensfunc_conf = self.SensfuncAttributes(Input_StandardFile, Output_SensitivityCurve, FitsFolder)

        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('sensfunc', Sensfunc_conf)
        print '--Using command'
        print Command

        #Run task command
        iraf.onedspec.sensfunc(**Sensfunc_conf)

        return

    def SensfuncTask_AllStellarCurves(self, FitsFolder, Input_StandardFile = None, Output_SensitivityCurve = None):

#         iraf.noao(_doprint=0)
#         iraf.onedspec(_doprint=0)     
            
        #Incase no output name is given, we generate one with the provided "preffix" (The defaul format is a_std_wolf.dat)
        if Input_StandardFile == None:  
            Input_StandardFile = self.GlobalStandard_Ouput

        if Output_SensitivityCurve == None:
            Output_SensitivityCurve = self.GlobalSensitiviyCurve
        
        Sensfunc_conf = self.SensfuncAttributes(Input_StandardFile, Output_SensitivityCurve, FitsFolder)

        #Run task command
        iraf.onedspec.sensfunc(**Sensfunc_conf)
        
        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('sensfunc', Sensfunc_conf, printindividually=False)
        print '--Using command'
        print Command

        return
       
    def CalibrateTask(self, InputFile, OutputFile, senstivityCurve, Fits_Folder, airmass_value, exptime, fnu_units = 'no', Suffix = 'FxW'):

#         iraf.noao(_doprint=0)
#         iraf.onedspec(_doprint=0)     
        
        #Incase no output name is given, we generate one with the provided "Suffix" (The defaul one is Fx)
        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile, Suffix)
        
        #Prepare dictionary with the configuration for the tasks           
        calibrate_conf = self.CalibrationAttributes(InputFile, OutputFile, senstivityCurve, Fits_Folder, airmass_value, exptime, fnu_units)
            
        #Run task command
        iraf.onedspec.calibrate(**calibrate_conf)
        
        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('calibrate', calibrate_conf, printindividually=False)
        print '--Using command'
        print Command
        
        return OutputFile    
    
    def ContinuumTask(self, InputFile, OutputFile, Fits_Folder, Suffix = 'N'):
        
        #Incase no output name is given, we generate one with the provided "Suffix" (The defaul one is N)
        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile, Suffix)
            
        continuum_conf = self.ContinuumAttributes(InputFile, OutputFile, Fits_Folder)

        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('continuum', continuum_conf, printindividually=False)
        print '--Using command'
        print Command

        
        iraf.continuum(**continuum_conf)
              
    def ScopyTask(self, InputFile, OutputFile, Fits_Folder, Wmin = 'INDEF', Wmax = 'INDEF', Suffix = 't'):

        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile, Suffix)

        scopy_conf = self.ScopyAttributes(InputFile, OutputFile, Fits_Folder, w1 = Wmin, w2 = Wmax)

        iraf.scopy(**scopy_conf)

        return OutputFile

    def SarithTask(self, Input1, Input2, operation, OutputFile, FitsFolder, Suffix = 'tell'):
        
        if OutputFile == None:
            OutputFile = self.outputNameGenerator(Input1, Suffix)        
        
        self.deleteIrafFiles(FitsFolder + OutputFile)

        sarith_conf = self.SarithAttributes(InputFile1 = Input1, InputFile2 = Input2, OuputFile = OutputFile, Fits_Folder = FitsFolder, operation = operation)
        
        iraf.sarith(**sarith_conf)
        
        return OutputFile

    def SplotTask(self, InputFile, Fits_Folder, xmin = 'INDEF', xmax = 'INDEF', ymin = 'INDEF', ymax='INDEF'):

        splot_conf = self.SplotAttributes(InputFile, Fits_Folder,xmin, xmax, ymin, ymax)
        
        iraf.splot(**splot_conf)

        return
          
    def DereddenTask(self, InputFile, OutputFile, Fits_Folder, R, Type, EBV = None, Suffix = 'EBVdust'):
        
        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile, Suffix)    
        
        self.deleteIrafFiles(Fits_Folder + OutputFile)
        
        deredden_conf = self.DereddenAttributes(InputFile, OutputFile, Fits_Folder, Type, EBV, R)

        iraf.deredden(**deredden_conf)
        
        return OutputFile               

    def Dopcortask(self, InputFile, OutputFile, Fits_Folder, z, Suffix = 'z'):

        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile, Suffix)    
        
        self.deleteIrafFiles(Fits_Folder + OutputFile)
                 
        dopcor_conf = self.DopcorAttributes(InputFile, OutputFile, Fits_Folder, z)
        
        iraf.dopcor(**dopcor_conf)
                        
        return OutputFile

    def readspectext(self, InputFile, OutputFile, Fits_Folder, conf_dict = None, Suffix = ''):
        
        if OutputFile == None:
            OutputFile = self.outputNameGenerator(InputFile, Suffix)    

        self.deleteIrafFiles(Fits_Folder + OutputFile)
        
        conf_dict = self.rspectextAttributes(InputFile, OutputFile, Fits_Folder)

        #Display the equivalent command in IRAF
        Command = self.printIrafCommand('rspectext', conf_dict)
        print '--Using command'
        print Command

        iraf.rspectext(**conf_dict)
        
        return
        
# def IRAFCommands_Redshift_Dust_Trim(Pattern,FilesList):
#     
#     if Pattern == "fxW.fits":
#         
#         print "----------------------------------------------------------------------------------------"
#         print "                                REDSHIFT,DUST AND POST-TRIM                             "
#         print "----------------------------------------------------------------------------------------"
#         
#         myObj_z_EBV = DATA_HIIGalaxiesProperties()
#         
#         CurrentFolder = "None"
#         
#         for i in range(len(FilesList)):
#         
#             SpectrumAddress = FilesList[i]
# 
#             ObjName, SpecName, SpecFolder = AddressAnalyzer(SpectrumAddress)
#             
#             InputName = SpecName[0:SpecName.find(".fits")]
#             
#             if CurrentFolder != SpecFolder:                         
#                     print "------ cd "+ SpecFolder
#                     print ""
#                     CurrentFolder = SpecFolder
#             
#             CodeSpec = SpecName[SpecName.find("obj")+3:SpecName.find("_")]            
#             
#             for i in range(len(myObj_z_EBV)):
#                 element = myObj_z_EBV[i]
#                 CodeObj = element[element.find("obj")+3:element.find("_")]
#                 if CodeObj == CodeSpec:
#                     myZ = element[element.find("_z")+2:element.find("_EBV")]
#                     myEBV = element[element.find("_EBV")+4:len(element)]
#                     print "--" + ObjName + ")" 
#                     print
#                     print "dopcor "+ InputName + " " + InputName+"_z" + " redshift=" + myZ
#                     print "deredden "+ InputName+"_z " + InputName+"_z_EBV" + " " + "R=3.2 Type=E value="+myEBV
#                     print "splot " + InputName+"_z_EBV " + "ymin=0 ymax=1.00e-15"
#                     print "scopy " + InputName+"_z_EBV " + InputName+"_z_EBV_t " + "w1=XXXXXX w2=XXXXXX rebin=no"
#                     print "splot " + InputName+"_z_EBV_t "
#                     prin            
                        