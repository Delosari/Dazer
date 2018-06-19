from pandas                     import DataFrame, read_pickle
from os                         import makedirs, walk, chdir, getcwd, remove, close
from os.path                    import isfile, isdir
from numpy                      import arange, loadtxt, savetxt, transpose, logical_not, in1d, linspace, where, datetime64, empty, array_str, mgrid, tile, ma, hstack, copy
from numpy.core                 import defchararray as np_f
from matplotlib                 import pyplot as plt
from pylatex                    import Document, Package, Figure, NoEscape
from collections                import OrderedDict
from astropy.io.fits            import getheader, getdata
from astropy.wcs                import WCS
from astropy.io                 import fits
from astropy.visualization      import ZScaleInterval
from itertools                  import cycle, izip_longest
import subprocess
import matplotlib.patches       as patches
from shutil                     import move
from numpy.polynomial.legendre  import legval
from datetime                   import timedelta
from tempfile                   import mkstemp

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def commands_to_log(task, task_conf, commands_log_address):

    command_String = task
    for key in task_conf:
        command_String = command_String + ' ' + key + '=' + task_conf[key]
    
    tasks_logged = loadtxt(commands_log_address, dtype = 'str', comments='--', delimiter = ' ', skiprows=1, usecols = [0], ndmin=1)
        
    with open(commands_log_address, 'a') as file_log:
        if (tasks_logged.size == 0):
            file_log.write('\n--' + task + '\n')
            file_log.write(command_String + '\n')
        elif (tasks_logged[-1] != task):
            file_log.write('\n--' + task + '\n')
            file_log.write(command_String + '\n')
        else:
            file_log.write(command_String + '\n')

def equivalent_Iraf_comand(Task, Attributes, verbose_output):
    
    #Extrack configuration dictionary
    keys_attrib     = Attributes.keys()
    values_attrib   = Attributes.values()
    
    if verbose_output:
        
        print '--Task atributtes:'
        for i in range(len(keys_attrib)):
            print keys_attrib[i] + '=' + str(values_attrib[i])
        
        print '\n'
    
    return

def set_Iraf_package(task, task_conf, commands_log_address):
    
    import pyraf
    
    if task == 'zerocombine':
        pyraf.iraf.noao.imred()
        pyraf.iraf.noao.imred.ccdred()    
        pyraf.iraf.noao.imred.ccdred.zerocombine(**task_conf)
        
    elif task == 'ccdproc':
        #With the imred package loaded type this to avoid the error no instrument loaded
        #ccdred.instrument = "ccddb$kpno/camera.dat"
        pyraf.iraf.noao.imred()
        pyraf.iraf.noao.imred.ccdred()
        pyraf.iraf.noao.imred.ccdred.ccdproc(**task_conf)
    
    elif task == 'flatcombine':
        pyraf.iraf.noao.imred()
        pyraf.iraf.noao.imred.ccdred()
        pyraf.iraf.noao.imred.ccdred.flatcombine(**task_conf)
    
    elif task == 'response':        
        pyraf.iraf.noao.twodspec()
        pyraf.iraf.noao.twodspec.longslit()
        pyraf.iraf.noao.twodspec.longslit.response(**task_conf)
    
    elif task == 'imcombine':
        pyraf.iraf.imcombine(**task_conf)

    elif task == 'imcopy':
        pyraf.iraf.imutil()
        pyraf.iraf.imutil.imcopy(**task_conf)

    elif task == 'imshift':
        pyraf.iraf.immatch()
        pyraf.iraf.immatch.imshift(**task_conf)

    elif task == 'skyflat':
        pyraf.iraf.noao.imred()
        pyraf.iraf.noao.imred.ccdred()
        pyraf.iraf.noao.imred.ccdred.combine(**task_conf)
    
    elif task == 'illumination':
        pyraf.iraf.noao.twodspec()
        pyraf.iraf.noao.twodspec.longslit()
        pyraf.iraf.noao.twodspec.longslit.illumination(**task_conf)
    
    elif task == 'identify':
        input_address       = task_conf['images']
        input_file          = input_address[input_address.rfind("/")+1:len(input_address)]
        folder_input        = input_address[0:input_address.rfind('/')]
        task_conf['images'] = input_file
        chdir(folder_input)
        pyraf.iraf.noao.twodspec()
        pyraf.iraf.noao.twodspec.longslit()
        pyraf.iraf.noao.twodspec.longslit.identify(**task_conf)
    
    elif task == 'reidentify':
        #Special case for identify set the working folder to the location
        input_address                   = task_conf['images']
        input_file                      = input_address[input_address.rfind("/")+1:len(input_address)]
        folder_input                    = input_address[0:input_address.rfind('/')]
        task_conf['images']             = input_file
        task_conf['referenc']           = input_file
        chdir(folder_input)
        pyraf.iraf.noao.twodspec()
        pyraf.iraf.noao.twodspec.longslit()
        pyraf.iraf.noao.twodspec.longslit.reidentify(**task_conf)
    
    elif task == 'fitcoords':
        #Special case for identify set the working folder to the location
        input_address                   = task_conf['images']
        input_file                      = input_address[input_address.rfind("/")+1:len(input_address)]
        folder_input                    = input_address[0:input_address.rfind('/')]
        task_conf['images']             = input_file
        chdir(folder_input)
        pyraf.iraf.noao.twodspec()
        pyraf.iraf.noao.twodspec.longslit()
        pyraf.iraf.noao.twodspec.longslit.fitcoords(**task_conf)
    
    elif task == 'transform':
        #Special case for identify set the working folder to the location
        input_address                   = task_conf['input']
        output_address                  = task_conf['output']
        input_file                      = input_address[input_address.rfind("/")+1:len(input_address)]
        output_file                     = output_address[output_address.rfind("/")+1:len(output_address)]
        folder_input                    = input_address[0:input_address.rfind('/')]
        task_conf['input']              = input_file
        task_conf['output']             = output_file
        chdir(folder_input)
        pyraf.iraf.noao.twodspec()
        pyraf.iraf.noao.twodspec.longslit()        
        pyraf.iraf.noao.twodspec.longslit.transform(**task_conf)
    
    elif task == 'dispcor':
        pyraf.iraf.noao.onedspec()
        pyraf.iraf.noao.onedspec.dispcor(**task_conf)       
        
    elif task == 'background':
        pyraf.iraf.noao.twodspec()
        pyraf.iraf.noao.twodspec.longslit()
        pyraf.iraf.noao.twodspec.longslit.background(**task_conf)
    
    elif task == 'apall':
        pyraf.iraf.noao.twodspec()
        pyraf.iraf.noao.twodspec.apextract()
        pyraf.iraf.noao.twodspec.apextract.apall(**task_conf)
    
    elif task == 'standard':
        pyraf.iraf.noao.onedspec()
        pyraf.iraf.noao.onedspec.standard(**task_conf)
        
    elif task == 'sensfunc':
        pyraf.iraf.noao.onedspec()
        pyraf.iraf.noao.onedspec.sensfunc(**task_conf)
        
    elif task == 'calibrate':
        pyraf.iraf.noao.onedspec()
        pyraf.iraf.noao.onedspec.calibrate(**task_conf)     

    elif task == 'continuum':
        pyraf.iraf.noao.onedspec()
        pyraf.iraf.noao.onedspec.continuum(**task_conf)    

    elif task == 'splot':
        pyraf.iraf.noao.onedspec()
        pyraf.iraf.noao.onedspec.splot(**task_conf)    

    elif task == 'sarith':
        pyraf.iraf.noao.onedspec()
        pyraf.iraf.noao.onedspec.sarith(**task_conf) 

    elif task == 'imarith':
        pyraf.iraf.imutil()
        pyraf.iraf.imutil.imarith(**task_conf)
        
    elif task == 'imstat':
        pyraf.iraf.imutil()
        pyraf.iraf.imutil.imstat(**task_conf)

    elif task == 'dopcor':
        pyraf.iraf.noao.onedspec()
        pyraf.iraf.noao.onedspec.dopcor(**task_conf)

    #Store the final command to a text file
    commands_to_log(task, task_conf, commands_log_address)
    
    return

def execute_Iraf_task(launching_command, commands_log, verbose_output=False):
    
    #Decompose comand to files
    FileAddress = launching_command.split(' ')[2]
    FileName    = FileAddress[FileAddress.rfind("/")+1:len(FileAddress)]
    task        = FileName[0:FileName.find('_')]
    
    print '-- IRAF task:', task
        
    #load the task configuration
    dict_keys, dict_values  = loadtxt(FileAddress, dtype='str', delimiter = ';', usecols = (0,1), unpack = True)
    task_conf               = OrderedDict(zip(dict_keys, dict_values))
    
    #Reproduce task name
    equivalent_Iraf_comand(task, task_conf, verbose_output)
    
    #Run the task
    set_Iraf_package(task, task_conf, commands_log)

    return

class fits_plots():
    
    def __init__(self):
        
        self.fig_type           = None
        self.page_grid_width    = None
        self.page_grid_height   = None
        
        self.Blue_Arc_lines = [4426.2047,
                         4657.8615,
                         4726.87,
                         4735.91,
                         5400.5615,
                         5852.5224,
                         6402.25,
                         6717.0216,
                         7723.8                     
                         ]

        self.Red_Arc_lines = [7245.17,
                         7383.98,
                         7635.1,
                         7948.18,
                         8006.16,
                         8014.79,
                         8115.31,
                         9122.9588,
                         9657.706,
                         9784.3944,
                         10052.1,
                         10470.05                         
                         ]
        
        linestyles = ["-","-","--","--","-.","-.",":",":"]
        self.linecycler = cycle(linestyles)
        
        self.plots_dict     = {'Blue' : 0, 'Red' : 1}
        self.frames_colors  = {'Blue arm':'bone', 'Red arm':'gist_heat'}

    def fits_to_frame(self, frame_address, axis_plot, color = None, ext = 0, valid = True, section = None):
        
        if color not in self.frames_colors:
            cmap_frame = 'gray'
        else:
            cmap_frame = self.frames_colors[color]
        
        #Clear the axis from previous plot
        axis_plot.clear()
        axis_plot.axes.get_yaxis().set_visible(False)
        axis_plot.axes.get_xaxis().set_visible(False)

        #Open the image
        with fits.open(frame_address) as hdu_list:
            image_data = hdu_list[ext].data
   
        #Complete image plotted
        if section is None:
                                                             
            #Get zscale limits for plotting the image
            IntensityLimits     = ZScaleInterval()
            int_min, int_max    = IntensityLimits.get_limits(image_data)[0], IntensityLimits.get_limits(image_data)[1]
            
            #Plot the data
            axis_plot.imshow(image_data, cmap=cmap_frame, origin='lower', vmin = int_min, vmax = int_max, interpolation='nearest', aspect='auto')
            axis_plot.set_xlim(0, image_data.shape[1])
            axis_plot.set_ylim(0, image_data.shape[0])

        #Only a section
        else:
            
            x_max, y_max        = int(section[0]), int(section[1])
            max_region          = image_data[y_max-5:y_max+5:, x_max-5:x_max+5]

            #Get zscale limits for plotting the image
            IntensityLimits     = ZScaleInterval()
            int_min, int_max    = IntensityLimits.get_limits(max_region)[0], IntensityLimits.get_limits(max_region)[1]
            
            #Plot the data
            axis_plot.imshow(max_region, cmap=cmap_frame, origin='lower', vmin = int_min, vmax = int_max, interpolation='nearest', aspect='auto')
            axis_plot.set_xlim(0, max_region.shape[1])
            axis_plot.set_ylim(0, max_region.shape[0])
        
        return
    
    def plot_spectral_axis(self, frame_address, axis_plot, label, ext = 0, valid = True):

        with fits.open(frame_address) as hdu_list:
            image_data = hdu_list[ext].data
        
        line_Style = '-' if valid else ':' 
        #colorTitle = 'black' if validity_check[index_i] else 'red'                      

        y_values = image_data.mean(axis=1) 
        x_values = range(len(y_values))
        axis_plot.plot(x_values, y_values, label = label, linestyle=line_Style)
              
        return

    def compute_background_median(self, fits_address, trace_address):

        #Read the text file
        file_trace  = open(trace_address)
        file_lines  = file_trace.readlines()        
        file_trace.close()
                
        #Open the image
        with fits.open(fits_address) as hdu_list:
            image_data = hdu_list[0].data
        
        #Reference indeces
        idx_start   = file_lines.index('\taxis\t1\n') + 1 #That '\t1' is the orientation of the frame
        idx_match   = [i for i in range(len(file_lines)) if '\tcenter\t' in file_lines[i]]
        
        #Aperture data
        center_line = file_lines[idx_match[0]]
        lower_line  = file_lines[idx_match[0]+1]
        upper_line  = file_lines[idx_match[0]+2]
        aper_center = map(float, center_line.split()[1:])
        aper_low    = map(float, lower_line.split()[1:])
        aper_high   = map(float, upper_line.split()[1:])
        number_pixels = aper_high - aper_low
        print number_pixels, int(number_pixels)
        
        #Fitting coefficients
        coef_n      = file_lines[idx_start].split()[1]
        fit_type    = float(file_lines[idx_start + 1].split()[0])
        order       = int(float(file_lines[idx_start + 2].split()[0]))
        #xmin        = int(float(file_lines[idx_start + 3].split()[0]))
        xmin        = 0
        xmax        = int(float(file_lines[idx_start + 4].split()[0]))
        coefs       = empty(order)
        for i in range(len(coefs)):
            coefs[i] = float(file_lines[idx_start + 5 + i].split()[0]) 

        #Plot the polynomial
        y_range     = arange(float(xmin), float(xmax))
        n           = (2 * y_range - (xmax + xmin)) / (xmax - xmin)
        poly_leg    = legval(n, coefs)
        trace_curve = poly_leg + aper_center[0]
          
        #Plot Background region
        idx_background  = [i for i in range(len(file_lines)) if '\t\tsample' in file_lines[i]]
        background_line = file_lines[idx_background[0]].split()[1:]
        
        y_grid, x_grid = mgrid[0:image_data.shape[0], 0:image_data.shape[1]]
        
        #GridAxis.imshow(image_data, origin='lower', cmap='gray', interpolation='none')
        
        Background_array = None
        for idx, region in enumerate(background_line):
            limits_region       = sorted(map(float, region.split(':')))
            low_limit, up_limit = limits_region[0], limits_region[1]
            trace_matrix        = tile(trace_curve, (image_data.shape[1], 1)).T
            bg_mask             = (x_grid > trace_matrix + low_limit) & (x_grid < trace_matrix  + up_limit)
            masked_image        = ma.masked_where(~bg_mask, image_data)
            if Background_array is None:
                Background_array = copy(masked_image)
            else:
                Background_array = hstack((Background_array,masked_image))
            
        background_flux = ma.sum(Background_array, axis=1)
            
        return background_flux      
        
    def trace_to_frame(self, frame_address, axis_plot, trace_folder, ext = 0, title_reference = ''):
        
        original_file   = frame_address[frame_address.rfind('/')+1:frame_address.rfind('.fits')]        
        trace_address   = trace_folder + 'ap_' + trace_folder[1:].replace('/','_').replace('_database_','_') + original_file.replace('/','_')
        #trace_address   = trace_address.replace('_standard_stars_','_objects_')

        #Read the text file
        file_trace  = open(trace_address)
        file_lines  = file_trace.readlines()        
        file_trace.close()
        
        #Reference indeces
        idx_start   = file_lines.index('\taxis\t1\n') + 1 #That '\t1' is the orientation of the frame
        idx_match   = [i for i in range(len(file_lines)) if '\tcenter\t' in file_lines[i]]
        
        #Aperture data
        center_line = file_lines[idx_match[0]]
        lower_line  = file_lines[idx_match[0]+1]
        upper_line  = file_lines[idx_match[0]+2]
        aper_center = map(float, center_line.split()[1:])
        aper_low    = map(float, lower_line.split()[1:])
        aper_high   = map(float, upper_line.split()[1:])

        #Fitting coefficients
        coef_n      = file_lines[idx_start].split()[1]
        fit_type    = float(file_lines[idx_start + 1].split()[0])
        order       = int(float(file_lines[idx_start + 2].split()[0]))
        xmin        = int(float(file_lines[idx_start + 3].split()[0]))
        xmax        = int(float(file_lines[idx_start + 4].split()[0]))
        coefs       = empty(order)
        for i in range(len(coefs)):
            coefs[i] = float(file_lines[idx_start + 5 + i].split()[0]) 

        #Plot the polynomial
        y_range     = arange(float(xmin), float(xmax))
        n           = (2 * y_range - (xmax + xmin)) / (xmax - xmin)
        poly_leg    = legval(n, coefs)
        trace_curve = poly_leg + aper_center[0]
        low_limit   = trace_curve + aper_low[0]
        high_limit  = trace_curve + aper_high[0]
        axis_plot.plot(trace_curve, y_range, color='red', linestyle=':')
        axis_plot.fill_betweenx(y_range, low_limit, high_limit, alpha=0.3, facecolor='green', edgecolor='green', linewidth=3.0)
          
        #Plot Background region
        idx_background  = [i for i in range(len(file_lines)) if '\t\tsample' in file_lines[i]]
        background_line = file_lines[idx_background[0]].split()[1:]
        print 'The regions', background_line
        for region in background_line:
            limits_region       = map(float, region.split(':'))
            low_limit_region    = trace_curve + limits_region[0]
            high_limit_region   = trace_curve + limits_region[1]    
            axis_plot.fill_betweenx(y_range, low_limit_region, high_limit_region, alpha=0.2, facecolor='yellow')
        
        #Wording for the plot
        coefficients_line   = r'$c_{{i}}$ = [{:s}]'.format(array_str(coefs, precision=3).translate(None, "[]"))
        title_lines         = '{:s} {:s} (order {:d})\n({:s})\nApperture (pixels) {:.2f}'.format(title_reference, 'legendre', order, coefficients_line, aper_high[0] - aper_low[0])
        axis_plot.set_title(title_lines, fontsize = 12) 

        return
 
    def fits_compare(self, pdf_address, indeces_frame, ext = 0, columns_mean = True):
        
        #Get data
        sorting_pattern     = ['frame_tag','ISIARM']
        files_name          = self.reducDf[indeces_frame].sort_values(sorting_pattern).file_name.values
        files_address       = self.reducDf[indeces_frame].sort_values(sorting_pattern).file_location.values
        frames_color        = self.reducDf[indeces_frame].sort_values(sorting_pattern).ISIARM.values
        frames_object       = self.reducDf[indeces_frame].sort_values(sorting_pattern).OBJECT.values
        validity_check      = self.reducDf[indeces_frame].sort_values(sorting_pattern).valid_file.values
                
        #Generate an artificial list with the range of files
        Number_frames   = len(files_address)
        indices_list    = range(Number_frames)

        #Figure configuration
        page_grid_width = 2
        self.Pdf_Fig, self.GridAxis = plt.subplots(1, page_grid_width, figsize=(10, 14)) 
        self.GridAxis_list = self.GridAxis.ravel()

        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=0cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}'))
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for sublist in grouper(page_grid_width, indices_list):
            for i in range(len(sublist)):
                
                #Get the index corresponding to the image to plot
                index_i = sublist[i]
            
                if index_i is not None:
                    colorTitle = 'black' if validity_check[index_i] else 'red'                      
                    self.fits_to_frame(files_address[index_i] + files_name[index_i], self.GridAxis_list[i], color=frames_color[index_i], ext=ext, valid = validity_check[index_i])
                    self.GridAxis_list[i].set_title('{filename}\n{object}'.format(filename = files_name[index_i], object = frames_object[index_i]), fontsize = 15, color=colorTitle)   

            plt.tight_layout()

            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla()
        
        #Plot mean row value      
        if columns_mean:

            #Figure for the plot            
            self.Pdf_Fig, self.GridAxis = plt.subplots(2, 1, figsize=(10, 12))   #figsize=(20, 24)
            self.GridAxis_list          = self.GridAxis.ravel()
            
            #Loop through the files
            for j in range(len(files_name)):
                
                #Get label
                CodeName = files_name[j][0:files_name[j].find('.')]
                
                #Get axis according to the color
                axis_index = self.plots_dict[frames_color[j].replace(' arm','')]

                #Plot the data
                self.plot_spectral_axis(files_address[j] + files_name[j], self.GridAxis_list[axis_index], CodeName, ext=ext, valid = validity_check[j])
            
            #Plot layout
            for idx, val in enumerate(['Blue arm spectra', 'Red arm spectra']):
                self.GridAxis[idx].set_xlabel('Pixel value',          fontsize = 15)
                self.GridAxis[idx].set_ylabel('Mean spatial count',   fontsize = 15)
                self.GridAxis[idx].set_title(val, fontsize = 20) 
                self.GridAxis[idx].tick_params(axis='both', labelsize=14)
                self.GridAxis[idx].legend(loc='upper right', prop={'size':14}, ncol=2)
            
            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla() 
 
        self.doc.generate_pdf(clean_tex=True)
        
        return

    def masked_pixels(self, pdf_address, indeces_frame, ext = 0, columns_mean = True):
        
        fill=False
        
        #Get data
        sorting_pattern     = ['RUN', 'frame_tag']
        files_name          = self.reducDf[indeces_frame].sort_values(sorting_pattern).file_name.values
        files_address       = self.reducDf[indeces_frame].sort_values(sorting_pattern).file_location.values
        frames_color        = self.reducDf[indeces_frame].sort_values(sorting_pattern).ISIARM.values
        frames_object       = self.reducDf[indeces_frame].sort_values(sorting_pattern).OBJECT.values
        validity_check      = self.reducDf[indeces_frame].sort_values(sorting_pattern).valid_file.values
                
        #Generate an artificial list with the range of files
        Number_frames   = len(files_address)
        indices_list    = range(Number_frames)

        #Figure configuration
        page_grid_width = 2
        self.Pdf_Fig, self.GridAxis = plt.subplots(1, page_grid_width, figsize=(10, 14)) 
        self.GridAxis_list = self.GridAxis.ravel()

        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=0cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}'))
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for sublist in grouper(page_grid_width, indices_list):
            for i in range(len(sublist)):
                
                #Get the index corresponding to the image to plot
                index_i = sublist[i]
            
                if index_i is not None:
                    colorTitle = 'black' if validity_check[index_i] else 'red'                      
                    self.fits_to_frame(files_address[index_i] + files_name[index_i], self.GridAxis_list[i], color=frames_color[index_i], ext=ext, valid = validity_check[index_i])
                    self.GridAxis_list[i].set_title('{filename}\n{object}'.format(filename = files_name[index_i], object = frames_object[index_i]), fontsize = 15, color=colorTitle)   

                    #load the mask
                    mask_file = '{rootfolder}badpix_{color}mask'.format(rootfolder=self.Catalogue_folder, color=frames_color[index_i].replace(' arm',''))
                    
                    c1, c2, r1, r2 = loadtxt(mask_file, dtype=int, comments='#', delimiter=' ', usecols=(0,1,2,3), unpack=True, ndmin=2)
                    
                    for j in range(len(c1)):
                        self.GridAxis_list[i].add_patch(patches.Rectangle((c1[j], r1[j]), c2[j] - c1[j], r2[j] - r1[j], linewidth = 0.5, color='yellow', fill=False))      # remove background

            plt.tight_layout()

            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla()
        
        #Plot mean row value      
        if columns_mean:

            #Figure for the plot            
            self.Pdf_Fig, self.GridAxis = plt.subplots(2, 1, figsize=(10, 12))   #figsize=(20, 24)
            self.GridAxis_list          = self.GridAxis.ravel()
            
            #Loop through the files
            for j in range(len(files_name)):
                
                #Get label
                CodeName = files_name[j][0:files_name[j].find('.')]
                
                #Get axis according to the color
                axis_index = self.plots_dict[frames_color[j].replace(' arm','')]

                #Plot the data
                self.plot_spectral_axis(files_address[j] + files_name[j], self.GridAxis_list[axis_index], CodeName, ext=ext, valid = validity_check[j])
            
            #Plot layout
            for idx, val in enumerate(['Blue arm spectra', 'Red arm spectra']):
                self.GridAxis[idx].set_xlabel('Pixel value',          fontsize = 15)
                self.GridAxis[idx].set_ylabel('Mean spatial count',   fontsize = 15)
                self.GridAxis[idx].set_title(val, fontsize = 20) 
                self.GridAxis[idx].tick_params(axis='both', labelsize=14)
                self.GridAxis[idx].legend(loc='upper right', prop={'size':14}, ncol=2)
            
            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla() 
 
        self.doc.generate_pdf(clean_tex=True)
        
        return

    def frame_combine(self, pdf_address, indx_frame_original, ext = 0, sorting_pattern = ['reduc_tag']):
        
        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}')) 
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #Through through the objects to generate to create all the plots
        frames_objects  = self.reducDf[indx_frame_original].frame_tag.unique()
        frames_colors   = self.reducDf[indx_frame_original].ISIARM.unique()
    
        for frame_object in frames_objects:
            
            for arm_color in frames_colors:
            
                sub_indices     = indx_frame_original & (self.reducDf.frame_tag == frame_object) & (self.reducDf.ISIARM == arm_color)
                files_name      = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).file_name.values
                files_address   = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).file_location.values
                validity_check  = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).valid_file.values
                
                Number_frames       = len(files_address)
                indices_list        = range(Number_frames)
                center_line_cords   = self.observation_dict['{CodeName}_refline_{Color}'.format(CodeName=frame_object, Color=arm_color.replace(' arm',''))]                   
                
                #---------------Figure complete frames
                page_grid_width     = 3
                Pdf_Fig, GridAxis   = plt.subplots(1, page_grid_width, figsize=(10, 14))  
                GridAxis_list       = GridAxis.ravel()
                                
                #Start looping through the indeces in groups with the same lengths as the columns * rows
                for sublist in grouper(page_grid_width, indices_list):
                    for i in range(len(sublist)):
                        
                        #Get the index corresponding to the image to plot
                        index_i = sublist[i]
                        
                        if index_i is not None:
                            colorTitle = 'black' if validity_check[index_i] else 'red'                      
                            self.fits_to_frame(files_address[index_i] + files_name[index_i], GridAxis_list[i], color=arm_color, ext=ext, valid = validity_check[index_i])
                            GridAxis_list[i].set_title('{filename}\n{object}'.format(filename = files_name[index_i], object = frame_object), fontsize = 15, color=colorTitle)
                        
                        else:
                            GridAxis_list[i].clear() #This seems abit too much but the tight layout crashes otherwise
 
                    plt.tight_layout()

                    #Add the plot                               
                    with self.doc.create(Figure(position='htbp')) as plot:       
                        plot.add_plot(height=NoEscape(r'1\textheight'))
                        self.doc.append(NoEscape(r'\newpage'))
                        plt.cla()
                        
                #---------------Figure reference lines maps
                grid_width, grid_height = 3, 2
                Pdf_Fig, GridAxis       = plt.subplots(grid_height, grid_width, figsize=(10, 12))  
                GridAxis_list           = GridAxis.ravel()
     
                #Start looping through the indeces in groups with the same lengths as the columns * rows
                for sublist in grouper(grid_width * grid_height, indices_list):
                    for i in range(len(sublist)):
                         
                        #Get the index corresponding to the image to plot
                        index_i = sublist[i]
     
                        if index_i is not None:
                            colorTitle = 'black' if validity_check[index_i] else 'red'
                            self.fits_to_frame(files_address[index_i] + files_name[index_i], GridAxis_list[i], color=arm_color, ext=ext, valid = validity_check[index_i], section=center_line_cords)
                            GridAxis_list[i].set_title('{filename}\n{object}'.format(filename = files_name[index_i], object = frame_object), fontsize = 15, color=colorTitle)   
                        else:
                            GridAxis_list[i].clear() #This seems a bit too much but the tight layout crashes otherwise
     
                    #Add the plot                               
                    with self.doc.create(Figure(position='htbp')) as plot:       
                        plot.add_plot(height=NoEscape(r'1\textheight'))
                        self.doc.append(NoEscape(r'\newpage'))
                        plt.cla()
     
                #---------------Reference line cut 
                grid_width, grid_height = 1, 2
                Pdf_Fig, GridAxis       = plt.subplots(grid_height, grid_width, figsize=(10, 12))  
                         
                for m in range(len(files_name)):
                     
                    with fits.open(files_address[m] + files_name[m]) as hdu_list:
                        image_data = hdu_list[ext].data
                         
                    x_max, y_max        = int(center_line_cords[0]), int(center_line_cords[1])
                    spatial_mean        = image_data[y_max-1:y_max+1, :].mean(axis=0)
                    spatial_length      = range(len(spatial_mean))

                    spectral_mean       = image_data[:, x_max-2:x_max+2].mean(axis=1)
                    spectral_length     = range(len(spectral_mean))

                    #Check if file is rejected:
                    if validity_check[m]:
                        label = files_name[m]
                    else:
                        label = 'REJECTED ' + files_name[m]
                     
                    GridAxis[0].step(spatial_length, spatial_mean, label = label)
                    GridAxis[1].step(spectral_length, spectral_mean, label = label)
                                 
                #Plot layout
                GridAxis[0].set_xlabel('Pixel value',          fontsize = 15)
                GridAxis[0].set_ylabel('Mean spatial count',   fontsize = 15)
                GridAxis[0].set_title(frame_object + ' ' + arm_color, fontsize = 15) 
                GridAxis[0].tick_params(axis='both', labelsize=15)
                GridAxis[0].legend(loc='upper right', prop={'size':12})
                GridAxis[0].set_xlim(x_max-40, x_max+40)

                GridAxis[1].set_xlabel('Pixel value',          fontsize = 15)
                GridAxis[1].set_ylabel('Mean spectral count',   fontsize = 15)
                GridAxis[1].set_title(frame_object + ' ' + arm_color, fontsize = 15) 
                GridAxis[1].tick_params(axis='both', labelsize=15)
                GridAxis[1].legend(loc='upper right', prop={'size':12})
                GridAxis[1].set_yscale('log')
     
                #Add the plot                               
                with self.doc.create(Figure(position='htbp')) as plot:       
                    plot.add_plot(height=NoEscape(r'1\textheight'))
                    self.doc.append(NoEscape(r'\newpage'))
                    plt.cla() 
          
        self.doc.generate_pdf(clean_tex=True)

    def frame_combine_shifted(self, pdf_address, indx_frame_original, ext = 0, sorting_pattern = ['reduc_tag']):
        
        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}')) 
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #Through through the objects to generate to create all the plots
        frames_objects  = self.reducDf[indx_frame_original].frame_tag.unique()
        frames_colors   = self.reducDf[indx_frame_original].ISIARM.unique()
                
        for frame_object in frames_objects:
            
            for arm_color in frames_colors:
            
                sub_indices     = indx_frame_original & (self.reducDf.frame_tag == frame_object) & (self.reducDf.ISIARM == arm_color)
                files_name      = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).file_name.values
                files_address   = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).file_location.values
                validity_check  = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).valid_file.values
                
                Number_frames       = len(files_address)
                indices_list        = range(Number_frames)
                center_line_cords   = self.observation_dict['{CodeName}_refline_{Color}'.format(CodeName=frame_object, Color=arm_color.replace(' arm',''))]                   
                                
                #---------------Figure complete frames
                page_grid_width     = 3
                Pdf_Fig, GridAxis   = plt.subplots(1, page_grid_width, figsize=(10, 14))  
                GridAxis_list       = GridAxis.ravel()
                                
                #Start looping through the indeces in groups with the same lengths as the columns * rows
                for sublist in grouper(page_grid_width, indices_list):
                    for i in range(len(sublist)):
                        
                        #Get the index corresponding to the image to plot
                        index_i = sublist[i]
                        
                        if index_i is not None:
                            colorTitle = 'black' if validity_check[index_i] else 'red'                      
                            self.fits_to_frame(files_address[index_i] + files_name[index_i], GridAxis_list[i], color=arm_color, ext=ext, valid = validity_check[index_i])
                            GridAxis_list[i].set_title('{filename}\n{object}'.format(filename = files_name[index_i], object = frame_object), fontsize = 15, color=colorTitle)
                        
                        else:
                            GridAxis_list[i].clear() #This seems abit too much but the tight layout crashes otherwise
 
                    plt.tight_layout()

                    #Add the plot                               
                    with self.doc.create(Figure(position='htbp')) as plot:       
                        plot.add_plot(height=NoEscape(r'1\textheight'))
                        self.doc.append(NoEscape(r'\newpage'))
                        plt.cla()
                        
                #---------------Figure reference lines maps
                grid_width, grid_height = 3, 2
                Pdf_Fig, GridAxis       = plt.subplots(grid_height, grid_width, figsize=(10, 12))  
                GridAxis_list           = GridAxis.ravel()
     
                #Start looping through the indeces in groups with the same lengths as the columns * rows
                for sublist in grouper(grid_width * grid_height, indices_list):
                    for i in range(len(sublist)):
                         
                        #Get the index corresponding to the image to plot
                        index_i = sublist[i]
     
                        if index_i is not None:
                            colorTitle = 'black' if validity_check[index_i] else 'red'
                            self.fits_to_frame(files_address[index_i] + files_name[index_i], GridAxis_list[i], color=arm_color, ext=ext, valid = validity_check[index_i], section=center_line_cords)
                            GridAxis_list[i].set_title('{filename}\n{object}'.format(filename = files_name[index_i], object = frame_object), fontsize = 15, color=colorTitle)   
                        else:
                            GridAxis_list[i].clear() #This seems a bit too much but the tight layout crashes otherwise
     
                    #Add the plot                               
                    with self.doc.create(Figure(position='htbp')) as plot:       
                        plot.add_plot(height=NoEscape(r'1\textheight'))
                        self.doc.append(NoEscape(r'\newpage'))
                        plt.cla()
     
                #---------------Reference line cut 
                grid_width, grid_height = 1, 2
                Pdf_Fig, GridAxis       = plt.subplots(grid_height, grid_width, figsize=(10, 12))  
                         
                for m in range(len(files_name)):
                     
                    with fits.open(files_address[m] + files_name[m]) as hdu_list:
                        image_data = hdu_list[ext].data
                         
                    x_max, y_max        = int(center_line_cords[0]), int(center_line_cords[1])
                    spatial_mean        = image_data[y_max-1:y_max+1, :].mean(axis=0)
                    spatial_length      = range(len(spatial_mean))

                    spectral_mean       = image_data[:, x_max-2:x_max+2].mean(axis=1)
                    spectral_length     = range(len(spectral_mean))

                    #Check if file is rejected:
                    if validity_check[m]:
                        label = files_name[m]
                    else:
                        label = 'REJECTED ' + files_name[m]
                     
                    GridAxis[0].step(spatial_length, spatial_mean, label = label)
                    GridAxis[1].step(spectral_length, spectral_mean, label = label)
                                 
                #Plot layout
                GridAxis[0].set_xlabel('Pixel value',          fontsize = 15)
                GridAxis[0].set_ylabel('Mean spatial count',   fontsize = 15)
                GridAxis[0].set_title(frame_object + ' ' + arm_color, fontsize = 15) 
                GridAxis[0].tick_params(axis='both', labelsize=15)
                GridAxis[0].legend(loc='upper right', prop={'size':12})
                GridAxis[0].set_xlim(x_max-40, x_max+40)

                GridAxis[1].set_xlabel('Pixel value',          fontsize = 15)
                GridAxis[1].set_ylabel('Mean spectral count',   fontsize = 15)
                GridAxis[1].set_title(frame_object + ' ' + arm_color, fontsize = 15) 
                GridAxis[1].tick_params(axis='both', labelsize=15)
                GridAxis[1].legend(loc='upper right', prop={'size':12})
                GridAxis[1].set_yscale('log')
     
                #Add the plot                               
                with self.doc.create(Figure(position='htbp')) as plot:       
                    plot.add_plot(height=NoEscape(r'1\textheight'))
                    self.doc.append(NoEscape(r'\newpage'))
                    plt.cla() 
          
        self.doc.generate_pdf(clean_tex=True)

    def objects_focus(self, pdf_address, indx_frame_original, ext = 0, sorting_pattern = ['ISIARM']):
        
        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}')) 
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #Through through the objects to generate to create all the plots
        frames_objects  = self.reducDf[indx_frame_original].frame_tag.unique()
        frames_colors   = self.reducDf[indx_frame_original].ISIARM.unique()
                
        for frame_object in frames_objects:

            for j in range(len(frames_colors)):

                grid_width, grid_height = 1, 2
                Pdf_Fig, GridAxis       = plt.subplots(grid_height, grid_width, figsize=(10, 12))  
          
                sub_indices     = indx_frame_original & (self.reducDf.frame_tag == frame_object) & (self.reducDf.ISIARM == frames_colors[j])
                files_name      = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).file_name.values
                files_address   = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).file_location.values
                validity_check  = self.reducDf[sub_indices].sort_values(sorting_pattern, ascending = False).valid_file.values
                                
                for perspective in [0, 1]:
                        
                    for m in range(len(files_name)):
                         
                        with fits.open(files_address[m] + files_name[m]) as hdu_list:
                            image_data = hdu_list[ext].data

                        #Check if file is rejected:
                        if validity_check[m]:
                            label = files_name[m]
                        else:
                            label = 'REJECTED ' + files_name[m]
                 
                        if  perspective == 0: 
                            #x_max, y_max       = int(center_line_cords[0]), int(center_line_cords[1])
                            #spatial_mean       = image_data[y_max-1:y_max+1, :].mean(axis=0)
                            spatial_mean        = image_data.mean(axis=0)
                            spatial_length      = range(len(spatial_mean))
                            GridAxis[perspective].step(spatial_length, spatial_mean, label = label)
                        
                        else:
                            #spectral_mean       = image_data[:, x_max-2:x_max+2].mean(axis=1)
                            spectral_mean       = image_data.mean(axis=1)   
                            spectral_length     = range(len(spectral_mean))
                            GridAxis[perspective].step(spectral_length, spectral_mean, label = label)
                                     
                    #Plot layout
                    GridAxis[perspective].set_xlabel('Pixel',          fontsize = 15)
                    GridAxis[perspective].set_ylabel('Mean spatial count',   fontsize = 15)
                    GridAxis[perspective].set_title(frame_object + ' ' + frames_colors[j], fontsize = 15) 
                    GridAxis[perspective].tick_params(axis='both', labelsize=15)
                    GridAxis[perspective].legend(loc='upper right', prop={'size':12})
     
                plt.tight_layout()
         
                #Add the plot                               
                with self.doc.create(Figure(position='htbp')) as plot:       
                    plot.add_plot(height=NoEscape(r'1\textheight'))
                    self.doc.append(NoEscape(r'\newpage'))
                    plt.cla() 
          
        self.doc.generate_pdf(clean_tex=True)
 
    def cosmic_plot(self, pdf_address, indeces_frame, ext = 0, colors_dict = None, columns_mean = True):

        #Get data
        sorting_pattern = ['ISIARM', 'frame_tag']
        files_name      = self.reducDf[indeces_frame].sort_values(sorting_pattern).file_name.values
        files_address   = self.reducDf[indeces_frame].sort_values(sorting_pattern).file_location.values
        frames_color    = self.reducDf[indeces_frame].sort_values(sorting_pattern).ISIARM.values
        frames_object   = self.reducDf[indeces_frame].sort_values(sorting_pattern).OBJECT.values
        
        #Generate an artificial list with the range of files
        Number_frames   = len(files_address)

        #Figure configuration
        titles          = ['Original', 'Mask', 'cleaned']
        page_grid_width = 3
        self.Pdf_Fig, self.GridAxis = plt.subplots(1, page_grid_width, figsize=(10, 14)) 
        self.GridAxis_list = self.GridAxis.ravel()

        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=0cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}'))
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))
                  
        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for i in range(Number_frames):
            
            clean_name  = files_address[i] + files_name[i]
            input_name  = clean_name.replace('_cr.fits', '.fits')
            mask_name   = clean_name.replace('_cr.fits', '_mask.fits')
            plot_frames = [input_name, input_name, clean_name]
                        
            for j in range(3):
                
                #Plot the frame           
                self.fits_to_frame(plot_frames[j], self.GridAxis_list[j], color=frames_color[i], ext=ext)
                self.GridAxis_list[j].set_title('{filename}\n{object}\n{type_f}'.format(filename = files_name[i], object = frames_object[i] ,type_f = titles[j]), fontsize = 15)   

                #Overplot cosmic rays location
                if j == 1:
                    
                    with fits.open(mask_name) as hdu_list:
                        image_mask = hdu_list[ext].data
             
                    mask_pixels = where(image_mask == 1)
                    xValues, yValues = mask_pixels[1], mask_pixels[0]
 
                    self.GridAxis_list[j].scatter(xValues, yValues, s=15, edgecolor='yellow', facecolor='none')

            plt.tight_layout()

            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla()
                      
        self.doc.generate_pdf(clean_tex=True)

    def cosmic_plot2(self, pdf_address, indeces_frame, ext = 0, colors_dict = None, columns_mean = True):
        
        files_name      = self.reducDf[indeces_frame].sort_values(['frame_tag']).file_name.values
        files_address   = self.reducDf[indeces_frame].sort_values(['frame_tag']).file_location.values
        frames_color    = self.reducDf[indeces_frame].sort_values(['frame_tag']).ISIARM.values
        slit_widths     = self.reducDf[indeces_frame].sort_values(['frame_tag']).ISISLITW.values
        tags_frame      = self.reducDf[indeces_frame].sort_values(['frame_tag']).frame_tag.values
        titles          = ['Original', 'cleaned', 'masked']
        Number_frames   = len(files_address)

        #Figure configuration
        Pdf_Fig, GridAxis   = plt.subplots(1, self.page_grid_width, figsize=(20, 24))  
        GridAxis_list       = GridAxis.ravel()

        #Create the doc
        doc = Document(pdf_address)
        doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        doc.append(NoEscape(r'\extrafloats{400}')) 
        doc.append(NoEscape(r'\pagenumbering{gobble}'))
                  
        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for i in range(Number_frames):
            
            clean_name  = files_address[i] + files_name[i]
            input_name  = clean_name.replace('_cr.fits', '.fits')
            mask_name   = clean_name.replace('_cr.fits', '_mask.fits')
            plot_frames = [input_name, input_name, clean_name]
            
            for j in range(3):
                                
                #Clear the axis from previous plot
                GridAxis_list[j].clear()
                GridAxis_list[j].axes.get_yaxis().set_visible(False)
                GridAxis_list[j].axes.get_xaxis().set_visible(False)
             
                with fits.open(plot_frames[j]) as hdu_list:
                    image_data = hdu_list[ext].data
                                                                  
                #Get zscale limits for plotting the image
                IntensityLimits     = ZScaleInterval()
                int_min, int_max    = IntensityLimits.get_limits(image_data)[0], IntensityLimits.get_limits(image_data)[1]
                
                #Plot the data
                cmap_frame = colors_dict[frames_color[i]]
                GridAxis_list[j].imshow(image_data, cmap=cmap_frame, origin='lower', vmin = int_min, vmax = int_max, interpolation='nearest', aspect='auto')
                GridAxis_list[j].set_title(files_name[i] + '\n' + titles[j] + '\n' + tags_frame[i] + ' ' + str(slit_widths[i]), fontsize = 12)                    
                plt.axis('tight')
                
                #Overplot cosmic rays location
                if j == 1:
                    with fits.open(mask_name) as hdu_list:
                        image_mask = hdu_list[ext].data
             
                    mask_pixels = where(image_mask == 1)
                    xValues     = mask_pixels[1]
                    yValues     = mask_pixels[0]
 
                    GridAxis_list[j].scatter(xValues, yValues, s=15, edgecolor='yellow', facecolor='none')
                    plt.axis('tight')
                plt.axis('tight')

            #Add the plot                               
            with doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                doc.append(NoEscape(r'\newpage'))
                plt.cla()
                     
        if columns_mean:

            #Figure for the plot            
            Pdf_Fig, GridAxis   = plt.subplots(2, 1, figsize=(20, 24))  
            GridAxis_list       = GridAxis.ravel()
            
            linestyles          = ["-","-","--","--","-.","-.",":",":"]
            linecycler          = cycle(linestyles)
            
            for j in range(len(files_name)):
                
                CodeName = files_name[j][0:files_name[j].find('.')]
                
                with fits.open(files_address[j] + files_name[j]) as hdu_list:
                    image_data = hdu_list[ext].data
                
                y_values = image_data.mean(axis=1) 
                x_values = range(len(y_values))
                                
                if frames_color[j] == 'Blue arm':
                    GridAxis_list[0].plot(x_values, y_values, label = CodeName, linestyle=next(linecycler))
                else:
                    GridAxis_list[1].plot(x_values, y_values, label = CodeName, linestyle=next(linecycler))
            
            #Plot layout
            plt.axis('tight')
            for idx, val in enumerate(['Blue arm spectra', 'Red arm spectra']):
                GridAxis[idx].set_xlabel('Pixel value',          fontsize = 30)
                GridAxis[idx].set_ylabel('Mean spatial count',   fontsize = 30)
                GridAxis[idx].set_title(val, fontsize = 30) 
                GridAxis[idx].tick_params(axis='both', labelsize=20)
                GridAxis[idx].legend(loc='upper right', prop={'size':12})
            
            #Add the plot                               
            with doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                doc.append(NoEscape(r'\newpage'))
                plt.cla() 
 
        doc.generate_pdf(clean_tex=True)

    def extracted_frames(self, pdf_address, indeces_frame, ext = 0):

        files_name          = self.reducDf[indeces_frame].sort_values(['frame_tag','ISIARM']).file_name.values
        files_address       = self.reducDf[indeces_frame].sort_values(['frame_tag','ISIARM']).file_location.values
        frames_color        = self.reducDf[indeces_frame].sort_values(['frame_tag','ISIARM']).ISIARM.values
        frames_object       = self.reducDf[indeces_frame].sort_values(['frame_tag','ISIARM']).OBJECT.values

        #Generate an artificial list with the range of files
        Number_frames       = len(files_address)
        indices_list        = range(Number_frames)

        #Figure configuration
        page_grid_width = 2
        self.Pdf_Fig, self.GridAxis = plt.subplots(1, page_grid_width, figsize=(8, 12))  
        self.GridAxis_list = self.GridAxis.ravel()

        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}'))
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #------Plotting trace curve observations
        for sublist in grouper(page_grid_width, indices_list):
            for i in range(len(sublist)):
                
                #Get the index corresponding to the image to plot
                index_i = sublist[i]
                if index_i is not None:
                
                    #Plot the fits frame
                    file_name_2d = files_name[index_i].replace('_e.fits', '.fits')
                    self.fits_to_frame(files_address[index_i] + file_name_2d, self.GridAxis_list[i], color=frames_color[index_i])
                
                    #Plot the trace
                    self.trace_to_frame(files_address[index_i] + file_name_2d, self.GridAxis_list[i], files_address[index_i] + 'database/', title_reference = '{:s}\n'.format(files_name[index_i]))
            
            plt.tight_layout()
                            
            #Add the plot
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla()
                
        #------Plotting extracted observations
        self.Pdf_Fig    = plt.figure(figsize = (16,10))
        self.GridAxis   = self.Pdf_Fig.add_subplot(111)        

        #Loop through the fits data:
        for i in range(Number_frames):
            
            wavelength, Flux_array, Header_0 = self.get_spectra_data(files_address[i] + files_name[i])
 
            if ('NAXIS2' in Header_0) and ('NAXIS3' in Header_0):
                appertures_number   = Header_0['NAXIS2']
                spectra_number      = 1
                for apper_n in range(appertures_number):
                    self.GridAxis.plot(wavelength, Flux_array[0][apper_n], label = 'extracted spectrum apperture: {apper_n}'.format(apper_n=apper_n+1))
            
            elif ('NAXIS2' in Header_0) and ('NAXIS1' in Header_0):
                appertures_number   = Header_0['NAXIS2']
                spectra_number      = 1
                for apper_n in range(appertures_number):
                    self.GridAxis.plot(wavelength, Flux_array[apper_n], label = 'extracted spectrum apperture: {apper_n}'.format(apper_n=apper_n+1))                                        
            
            else:
                spectra_number      = 1
                appertures_number   = 1
                self.GridAxis.plot(wavelength, Flux_array, label = 'extracted spectrum')
            
            #Plot wording
            self.GridAxis.set_xlabel(r'Wavelength $(\AA)$',  fontsize = 20)
            self.GridAxis.set_ylabel('Flux' + r'$(erg\,cm^{-2} s^{-1} \AA^{-1})$',  fontsize = 20)
            self.GridAxis.set_title(files_name[i] + ' extracted spectrum', fontsize = 25)   
            self.GridAxis.legend()
           
            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textwidth'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla()
        
        #Generate pdf
        self.doc.generate_pdf(clean_tex=True)               

        return

    def flats_plotting(self, pdf_address, indeces_frame, ext = 0, columns_mean = False):
        
        files_name          = self.reducDf[indeces_frame].sort_values(['frame_tag','ISIARM']).file_name.values
        files_address       = self.reducDf[indeces_frame].sort_values(['frame_tag','ISIARM']).file_location.values
        frames_color        = self.reducDf[indeces_frame].sort_values(['frame_tag','ISIARM']).ISIARM.values
        frames_object       = self.reducDf[indeces_frame].sort_values(['frame_tag','ISIARM']).OBJECT.values
        
        #Generate an artificial list with the range of files
        Number_frames       = len(files_address)
        indices_list        = range(Number_frames)

        #Figure configuration
        page_grid_width = 4
        self.Pdf_Fig, self.GridAxis = plt.subplots(1, page_grid_width, figsize=(30, 40))  
        self.GridAxis_list = self.GridAxis.ravel()
        
        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}'))
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for sublist in grouper(page_grid_width, indices_list):
            for i in range(len(sublist)):
                
                #Get the index corresponding to the image to plot
                index_i = sublist[i]
                
                if index_i is not None:
                    self.fits_to_frame(files_address[index_i] + files_name[index_i], self.GridAxis_list[i], frames_color[index_i], ext = ext)
                    self.GridAxis_list[i].set_title('{filename}\n{object}'.format(filename = files_name[index_i], object = frames_object[index_i]), fontsize = 25)                    
                    
                #plot spectral limits for matching high an low orders:
                if 'Flat_limits' in self.observation_dict:
                    x_limits = self.GridAxis_list[i].get_xlim()
                    limits_flat = map(int, self.observation_dict['Flat_limits'])    
                    if 'Blue' in frames_color[index_i]:
                        self.GridAxis_list[i].plot(x_limits, (limits_flat[0], limits_flat[0]), color='purple', linewidth=3, linestyle='-.')
                        self.GridAxis_list[i].plot(x_limits, (limits_flat[1], limits_flat[1]), color='purple', linewidth=3, linestyle='-.')
                    else:
                        self.GridAxis_list[i].plot(x_limits, (limits_flat[2], limits_flat[2]), color='purple', linewidth=3, linestyle='-.')
                        self.GridAxis_list[i].plot(x_limits, (limits_flat[3], limits_flat[3]), color='purple', linewidth=3, linestyle='-.')                        
            
            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla()
                
        #Plotting spectral axis conf with narrow limits
        if columns_mean:

            #Figure configuration
            self.Pdf_Fig, self.GridAxis = plt.subplots(2, 1, figsize=(10, 12))   #figsize=(20, 24)
            self.GridAxis_list          = self.GridAxis.ravel()
            
            #Loop through files
            for j in range(len(files_name)):
                
                #Get label
                CodeName = files_name[j][0:files_name[j].find('.')]
                
                #Get axis according to the color
                axis_index = self.plots_dict[frames_color[j].replace(' arm','')]

                #Plot the data
                self.plot_spectral_axis(files_address[j] + files_name[j], self.GridAxis_list[axis_index], CodeName)
            
            #plot spectral limits for matching high an low orders:
            if 'Flat_limits' in self.observation_dict:
                limits_flat = map(int, self.observation_dict['Flat_limits'])
                self.GridAxis[0].axvline(limits_flat[0], color = 'red', linewidth = 1, linestyle='--')
                self.GridAxis[0].axvline(limits_flat[1], color = 'red', linewidth = 2, linestyle='--')                      
                self.GridAxis[1].axvline(limits_flat[2], color = 'red', linewidth = 1, linestyle='--')
                self.GridAxis[1].axvline(limits_flat[3], color = 'red', linewidth = 2, linestyle='--')                    
            
            #Plot cropping limits
            for idx, val in enumerate(['Blue', 'Red']):
                cropping = map(int, self.observation_dict[val + '_cropping'])
                self.GridAxis[idx].axvline(cropping[2], color = 'red', linewidth = 1)
                self.GridAxis[idx].axvline(cropping[3], color = 'red', linewidth = 2)
                self.GridAxis[idx].set_xlim(cropping[2] - 10, cropping[3] + 10)    
                self.GridAxis[idx].set_ylim(0.3, 1.3)
                
            #Plot layout
            for idx, val in enumerate(['Blue arm spectra', 'Red arm spectra']):
                self.GridAxis[idx].set_xlabel('Pixel value',          fontsize = 15)
                self.GridAxis[idx].set_ylabel('Mean spatial count',   fontsize = 15)
                self.GridAxis[idx].set_title(val, fontsize = 20) 
                self.GridAxis[idx].tick_params(axis='both', labelsize=14)
                self.GridAxis[idx].legend(loc='upper right', prop={'size':14}, ncol=2)
            
            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla() 
           
        #Plotting spectral axis conf with narrow limits
        if columns_mean:

            #Figure configuration
            self.Pdf_Fig, self.GridAxis = plt.subplots(2, 1, figsize=(10, 12))   #figsize=(20, 24)
            self.GridAxis_list          = self.GridAxis.ravel()
            
            #Loop through files
            for j in range(len(files_name)):
                
                #Get label
                CodeName = files_name[j][0:files_name[j].find('.')]
                
                #Get axis according to the color
                axis_index = self.plots_dict[frames_color[j].replace(' arm','')]

                #Plot the data
                self.plot_spectral_axis(files_address[j] + files_name[j], self.GridAxis_list[axis_index], CodeName)
            
            #plot spectral limits for matching high an low orders:
            if 'Flat_limits' in self.observation_dict:
                limits_flat = map(int, self.observation_dict['Flat_limits'])
                self.GridAxis[0].axvline(limits_flat[0], color = 'purple', linewidth = 1)
                self.GridAxis[0].axvline(limits_flat[1], color = 'purple', linewidth = 2)                      
                self.GridAxis[1].axvline(limits_flat[2], color = 'purple', linewidth = 1)
                self.GridAxis[1].axvline(limits_flat[3], color = 'purple', linewidth = 2)                    
            
            #Plot cropping limits
            for idx, val in enumerate(['Blue', 'Red']):
                cropping = map(int, self.observation_dict[val + '_cropping'])
                self.GridAxis[idx].axvline(cropping[2], color = 'red', linewidth = 1)
                self.GridAxis[idx].axvline(cropping[3], color = 'red', linewidth = 2)
                
            #Plot layout
            for idx, val in enumerate(['Blue arm spectra', 'Red arm spectra']):
                self.GridAxis[idx].set_xlabel('Pixel value',          fontsize = 15)
                self.GridAxis[idx].set_ylabel('Mean spatial count',   fontsize = 15)
                self.GridAxis[idx].set_title(val, fontsize = 20) 
                self.GridAxis[idx].tick_params(axis='both', labelsize=14)
                self.GridAxis[idx].legend(loc='upper right', prop={'size':14}, ncol=2)
            
            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla() 
           
        self.doc.generate_pdf(clean_tex=True)

    def fits_catalogue(self, pdf_address, indeces_frame, ext = 0, columns_mean = False, sorting_pattern = ['frame_tag','ISIARM']):
        
        files_name          = self.reducDf[indeces_frame].sort_values(sorting_pattern).file_name.values
        files_address       = self.reducDf[indeces_frame].sort_values(sorting_pattern).file_location.values
        frames_color        = self.reducDf[indeces_frame].sort_values(sorting_pattern).ISIARM.values
        frames_object       = self.reducDf[indeces_frame].sort_values(sorting_pattern).OBJECT.values
        validity_check      = self.reducDf[indeces_frame].sort_values(sorting_pattern).valid_file.values
        
        #Generate an artificial list with the range of files
        Number_frames   = len(files_address)
        indices_list    = range(Number_frames)

        #Figure configuration
        page_grid_width = 5
        self.Pdf_Fig, self.GridAxis = plt.subplots(1, page_grid_width, figsize=(10, 14)) 
        self.GridAxis_list = self.GridAxis.ravel()

        #Create the doc
        self.doc = Document(pdf_address)
        self.doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=0cm']))
        self.doc.append(NoEscape(r'\extrafloats{100}'))
        self.doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for sublist in grouper(page_grid_width, indices_list):
            for i in range(len(sublist)):
                
                #Get the index corresponding to the image to plot
                index_i = sublist[i]
            
                if index_i is not None:
                    colorTitle = 'black' if validity_check[index_i] else 'red'
                    self.fits_to_frame(files_address[index_i] + files_name[index_i], self.GridAxis_list[i], frames_color[index_i], ext=ext, valid = validity_check[index_i])
                    self.GridAxis_list[i].set_title('{filename}\n{object}'.format(filename = files_name[index_i], object = frames_object[index_i]), fontsize = 12, color=colorTitle)   

            plt.tight_layout()

            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla()
            
            for axis in self.GridAxis_list:
                axis.clear()

        if columns_mean:

            #Figure for the plot            
            self.Pdf_Fig, self.GridAxis = plt.subplots(2, 1, figsize=(10, 12))   #figsize=(20, 24)
            self.GridAxis_list          = self.GridAxis.ravel()
                        
            #Loop through the files
            for j in range(len(files_name)):
                
                #Get label
                CodeName = files_name[j][0:files_name[j].find('.')]
                
                #Get axis according to the color
                if frames_color[j] in self.frames_colors:
                    axis_index = self.plots_dict[frames_color[j].replace(' arm','')]
                else:
                    axis_index = 0

                #Plot the data
                self.plot_spectral_axis(files_address[j] + files_name[j], self.GridAxis_list[axis_index], CodeName, ext=ext, valid = validity_check[j])
            
            #Plot layout
            for idx, val in enumerate(['Blue arm spectra', 'Red arm spectra']):
                self.GridAxis[idx].set_xlabel('Pixel value',          fontsize = 15)
                self.GridAxis[idx].set_ylabel('Mean spatial count',   fontsize = 15)
                self.GridAxis[idx].set_title(val, fontsize = 20) 
                self.GridAxis[idx].tick_params(axis='both', labelsize=14)
                self.GridAxis[idx].legend(loc='upper right', prop={'size':14}, ncol=2)
            
            #Add the plot                               
            with self.doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                self.doc.append(NoEscape(r'\newpage'))
                plt.cla() 
             
        self.doc.generate_pdf(clean_tex=True)

    def arcs_compare(self, pdf_address, files_address, files_name, frames_color, ext = 0, colors_dict = None):
        
        #Generate an artificial list with the range of files
        Number_frames   = len(files_address)

        #Create the doc
        doc = Document(pdf_address)
        doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        doc.append(NoEscape(r'\pagenumbering{gobble}'))
            
        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for i in range(Number_frames):
                                             
            #Load the image #need to remove this check.. simply add the files to rejected images
            hdu_list    = fits.open(files_address[i] + files_name[i])
            image_head  = hdu_list[ext].header
            image_data  = hdu_list[ext].data  #This might be different for non frames
            wcs         = WCS(image_head)                            
            hdu_list.close()

            #Figure configuration
            #Pdf_Fig, GridAxis   = plt.subplots(figsize=(20, 24), projection=wcs)  
            Pdf_Fig, GridAxis   = plt.subplots(figsize=(20, 28), ncols=2, subplot_kw={'projection': wcs})
            GridAxis_list       = GridAxis.ravel()
            
            #Get zscale limits for plotting the image
            IntensityLimits     = ZScaleInterval()
            int_min, int_max    = IntensityLimits.get_limits(image_data)[0], IntensityLimits.get_limits(image_data)[1]
            
            #Plot the data
            cmap_frame = colors_dict[frames_color[i]]
            GridAxis_list[0].imshow(image_data, cmap=cmap_frame, origin='lower', vmin = int_min, vmax = int_max, interpolation='nearest', aspect='auto')
            GridAxis_list[1].imshow(image_data, cmap=cmap_frame, origin='lower', vmin = int_min, vmax = int_max, interpolation='nearest', aspect='auto')

            GridAxis_list[0].set_title('Combined arc', fontsize = 28)
            GridAxis_list[1].set_title('Combined arc\nwith reference wavelengths', fontsize = 28)
            
            #Plot the lines
            x_limits = GridAxis_list[1].get_xlim()
            if frames_color[i] == 'Blue arm':
                for line_wavelength in self.Blue_Arc_lines:                
                    GridAxis_list[1].plot(x_limits, (line_wavelength, line_wavelength), color='purple', linewidth=3, linestyle='-.', transform=GridAxis_list[1].get_transform('world'))
            else:
                for line_wavelength in self.Red_Arc_lines:                
                    GridAxis_list[1].plot(x_limits, (line_wavelength, line_wavelength), color='purple', linewidth=3, linestyle='-.', transform=GridAxis_list[1].get_transform('world'))

            plt.axis('tight')
                        
            #Add the plot                               
            with doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                doc.append(NoEscape(r'\newpage'))
                plt.cla()
   
        doc.generate_pdf(clean_tex=True)

    def fast_combined(self, pdf_address, indeces_frame, observation_dict, reducDf, ext = 0, colors_dict = None):
        
        files_name      = self.reducDf[indeces_frame].sort_values(['frame_tag']).file_name.values
        files_address   = self.reducDf[indeces_frame].sort_values(['frame_tag']).file_location.values
        frames_color    = self.reducDf[indeces_frame].sort_values(['frame_tag']).ISIARM.values
        frames_object   = self.reducDf[indeces_frame].sort_values(['frame_tag']).frame_tag.values
        
        #Generate an artificial list with the range of files
        Number_frames   = len(files_address)
        indices_list    = range(Number_frames)

        #Figure configuration
        Pdf_Fig, GridAxis   = plt.subplots(1, 2, figsize=(20, 24))  
        GridAxis_list       = GridAxis.ravel()

        #Create the doc
        doc = Document(pdf_address)
        doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        doc.append(NoEscape(r'\extrafloats{100}')) 
        doc.append(NoEscape(r'\pagenumbering{gobble}'))
         
        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for sublist in grouper(2, indices_list):
            for i in range(len(sublist)):
                
                #Get the index corresponding to the image to plot
                index_i = sublist[i]
                
                #Clear the axis from previous plot
                GridAxis_list[i].clear()
                GridAxis_list[i].axes.get_yaxis().set_visible(False)
                GridAxis_list[i].axes.get_xaxis().set_visible(False)

                if index_i is not None:

                    with fits.open(files_address[index_i] + files_name[index_i]) as hdu_list:
                        image_data = hdu_list[ext].data
                                                                      
                    #Get zscale limits for plotting the image
                    IntensityLimits     = ZScaleInterval()
                    int_min, int_max    = IntensityLimits.get_limits(image_data)[0], IntensityLimits.get_limits(image_data)[1]
                    
                    #Plot the data
                    GridAxis_list[i].clear()
                    GridAxis_list[i].imshow(image_data, cmap=colors_dict[frames_color[index_i]], origin='lower', vmin = int_min, vmax = int_max, interpolation='nearest', aspect='auto')
                    GridAxis_list[i].set_title(files_name[index_i], fontsize = 20)                    
               
                    #Loading the reference line location
                    arm_color = frames_color[index_i].split()[0]
                    obj_refLine_key = '{CodeName}_refline_{Color}'.format(CodeName=frames_object[index_i], Color=arm_color)
                    refLine_cords   = observation_dict[obj_refLine_key]
                    x_peak, y_peak  = int(refLine_cords[0]), int(refLine_cords[1])
                    GridAxis_list[i].scatter(x_peak, y_peak, s=30, facecolor='green')
                    
                    #Loading the cropping area
                    obj_crop_key    = '{Color}_cropping'.format(Color=arm_color)
                    crop_cords      = map(int, observation_dict[obj_crop_key])
                    GridAxis_list[i].add_patch(patches.Rectangle((crop_cords[0], crop_cords[2]), crop_cords[1] - crop_cords[0], crop_cords[3] - crop_cords[2], linewidth = 2, color='red', fill=False))      # remove background

                    #Loading the scaling area area
                    obj_scaling_key = '{Color}_scale_region'.format(Color=arm_color)
                    scale_cords     = map(int, observation_dict[obj_scaling_key])
                    GridAxis_list[i].add_patch(patches.Rectangle((scale_cords[0], scale_cords[2]), scale_cords[1] - scale_cords[0], scale_cords[3] - scale_cords[2], linewidth = 2, color='black', fill=False))      # remove background

#                     try:
#                     #Print the maxima
#                     if frames_color[index_i] == 'Blue arm':
#                         region_zone =   self.default_OIII5007_region
#                         scale_zone  =   self.default_blue_Scaleregion
#                     elif frames_color[index_i] == 'Red arm':
#                         region_zone =   self.default_SIII9531_region
#                         scale_zone  =   self.default_red_Scaleregion
#
#                         section         = image_data[region_zone[0]:region_zone[1],region_zone[2]:region_zone[3]]
#                         max_value_sec   = npmax(section)
#                         max_indeces_sec = where(image_data == max_value_sec)  
#                                             
#                         new_limit           = max_indeces_sec[0][0] - 50                            
#                         section_2           = image_data[region_zone[0]:new_limit, region_zone[2]:region_zone[3]]
#                         max_value_sec2      = npmax(section_2)
#                         max_indeces_sec2    = where(image_data == max_value_sec2)
#                          
#                         #plotting the points
#                         GridAxis_list[i].clear()
#                          
#                         #Get zscale limits for plotting the image
#                         IntensityLimits     = ZScaleInterval()
#                         int_min, int_max    = IntensityLimits.get_limits(image_data)[0], IntensityLimits.get_limits(image_data)[1]                          
#                         cmap_frame          = self.frames_colors[frames_color[index_i]]
#                         GridAxis_list[i].imshow(image_data, cmap=cmap_frame, origin='lower', vmin = int_min, vmax = int_max)
#                          
#                         GridAxis_list[i].scatter(max_indeces_sec[1], max_indeces_sec[0], s=40, edgecolor='black', facecolor='none')
#                         GridAxis_list[i].scatter(max_indeces_sec2[1], max_indeces_sec2[0], s=40, edgecolor='yellow', facecolor='none')
#                         GridAxis_list[i].text(max_indeces_sec[1] + 10, max_indeces_sec[0], '{x} {y}'.format(x =max_indeces_sec[1], y=max_indeces_sec[0]), fontsize=16)
#                         GridAxis_list[i].text(max_indeces_sec2[1] + 10, max_indeces_sec2[0], '{x} {y}'.format(x =max_indeces_sec2[1], y=max_indeces_sec2[0]), fontsize=16)
#                         GridAxis_list[i].add_patch(patches.Rectangle((scale_zone[2], scale_zone[0]), scale_zone[3] - scale_zone[2], scale_zone[1] - scale_zone[0], fill=False))      # remove background
#                          
#                         GridAxis_list[i].set_title(files_name[index_i], fontsize = 20)                    
#                         GridAxis_list[i].axes.get_yaxis().set_visible(False)
#                         GridAxis_list[i].axes.get_xaxis().set_visible(False)
#                      
#                     except:
#                         print 'Could not print'

            #Add the plot                               
            with doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                doc.append(NoEscape(r'\newpage'))
                plt.cla()
        
        doc.generate_pdf(clean_tex=True)
    
        return

    def spectra(self, pdf_address, indeces_frame,  ext = 0):

        files_name      = self.reducDf[indeces_frame].sort_values(['frame_tag']).file_name.values
        files_address   = self.reducDf[indeces_frame].sort_values(['frame_tag']).file_location.values
        frames_color    = self.reducDf[indeces_frame].sort_values(['frame_tag']).ISIARM.values
        Number_frames   = len(files_name)
                
        #Create the doc
        doc = Document(pdf_address)
        doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm', 'landscape']))
        doc.append(NoEscape(r'\extrafloats{100}'))
        doc.append(NoEscape(r'\pagenumbering{gobble}'))

        #Loop through the fits data:
        for i in range(Number_frames):
        
            wavelength, Flux_array, Header_0 = self.get_spectra_data(files_address[i] + files_name[i])
                   
            #Define image
            Fig     = plt.figure(figsize = (16,10))
            Axis    = Fig.add_subplot(111)
            
            Axis.set_xlabel(r'Wavelength $(\AA)$',  fontsize = 20)
            Axis.set_ylabel('Flux' + r'$(erg\,cm^{-2} s^{-1} \AA^{-1})$',  fontsize = 20)
            Axis.set_title(files_name[i] + ' extracted spectrum', fontsize = 25)   

            if ('NAXIS2' in Header_0) and ('NAXIS3' in Header_0):
                appertures_number   = Header_0['NAXIS2']
                spectra_number      = 1
                for apper_n in range(appertures_number):
                    Axis.plot(wavelength, Flux_array[0][apper_n], label = 'extracted spectrum apperture: {apper_n}'.format(apper_n=apper_n+1))
            elif ('NAXIS2' in Header_0) and ('NAXIS1' in Header_0):
                appertures_number   = Header_0['NAXIS2']
                spectra_number      = 1
                for apper_n in range(appertures_number):
                    Axis.plot(wavelength, Flux_array[apper_n], label = 'extracted spectrum apperture: {apper_n}'.format(apper_n=apper_n+1))                                        
            else:
                spectra_number      = 1
                appertures_number   = 1
                Axis.plot(wavelength, Flux_array, label = 'extracted spectrum')
   
            Axis.legend()
            
            with doc.create(Figure(position='htbp')) as plot:      
                          
#                     plot.add_plot(height=NoEscape(r'1\textheight'))
                plot.add_plot(width=NoEscape(r'1\textwidth'))

                #Reset the image
                doc.append(NoEscape(r'\newpage'))
                plt.clf()
                
        doc.generate_pdf(clean_tex=True)               
           
        return

class pyraf_task_configuration():
    
    def __init__(self):
        
        self.task_attributes        = OrderedDict()
        self.objects_configuration  = 'HII_galaxies'
        self.script_location        = '/home/vital/git/Thesis_Pipeline/Thesis_Pipeline/Spectra_Reduction/'
        self.Standard_Stars_Folder  = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/StandardStars_Calibration/'
        
    def standar_stars_calibrationFile(self, StdStar_code):
        
        calibration_dict = {}
        calibration_dict['calibration_folder'] = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/StandardStars_Calibration/'
    
        if StdStar_code in ['FEI34', 'F34', 'feige34', 'Feige']:
            calibration_dict['file name'] = 'feige34_stis004_nobandpass'
            calibration_dict['bandwidth'] = 40
            calibration_dict['bandsep']  = 40
            
        elif StdStar_code in ['BD29']:
            calibration_dict['file name'] = 'bd29d2091_stis_003_nobandpass'
            calibration_dict['bandwidth'] = 40
            calibration_dict['bandsep']  = 40

        elif StdStar_code in ['BD28', 'BD+28']:
            calibration_dict['file name'] = 'bd28_d4211stis004_nobandpass'
            calibration_dict['bandwidth'] = 40
            calibration_dict['bandsep']  = 40
      
        elif StdStar_code in ['BD33', 'BD+33']:
            calibration_dict['file name'] = 'bd33_d2642004_nobandpass'
            calibration_dict['bandwidth'] = 40
            calibration_dict['bandsep']  = 40
            
        elif StdStar_code in ['WOLF1346', 'wolf', 'Wolf_1346', 'WOLF1346A']:
            calibration_dict['file name'] = 'wolf_oke1974_40a'
            calibration_dict['bandwidth'] = 'INDEF'
            calibration_dict['bandsep']  = 'INDEF'

        elif StdStar_code in ['bd17', 'BD+17', 'sp2209+178', 'sp2209178', 'BD17', 'BD+17_4708', 'SP2209+178']:
            calibration_dict['file name'] = 'bd17_d4708stisnic006_nobandpass'
            calibration_dict['bandwidth'] = 40
            calibration_dict['bandsep']  = 40

        elif StdStar_code in ['g191', 'G191']:
            calibration_dict['file name'] = 'g191_b2bstisnic006_nobandpass'
            calibration_dict['bandwidth'] = 40
            calibration_dict['bandsep']  = 40

        elif StdStar_code in ['g158', 'G158']:
            calibration_dict['file name'] = 'g158_oke1990_1a'
            calibration_dict['bandwidth'] = 40
            calibration_dict['bandsep']  = 40

        elif StdStar_code in ['hd19445', 'sp0305+261', 'sp0305261', 'SP0305+261']:
            calibration_dict['file name'] = 'hd19445_oke1983_40a'
            calibration_dict['bandwidth'] = 'INDEF'
            calibration_dict['bandsep']  = 'INDEF'

        elif StdStar_code in ['hd84937', 'HD84937']:
            calibration_dict['file name']   = 'hd84937_oke1983_40a' #'hd84937_oke1983_40a'
            calibration_dict['bandwidth']   = 'INDEF'
            calibration_dict['bandsep']     = 'INDEF'

        return calibration_dict
    
    def load_task_configuration(self, task):

        #Load the task attributes
        if task == 'zerocombine':
            self.zerocombine_task_configuration()

        elif task == 'ccdproc':
            self.ccdproc_task_configuration()

        elif task == 'flatcombine':
            self.flatcombine_task_configuration()

        elif task == 'response':
            self.response_task_configuration()

        elif task == 'imcombine':
            self.imcombine_task_configuration()

        elif task == 'imshift':
            self.imshift_configuration()

        elif task == 'imcopy':
            self.imcopy_task_configuration()

        elif task == 'illumination':
            self.illumination_task_configuration()

        elif task == 'lacosmic':
            self.lacosmic_task_configuration()

        elif task == 'identify':
            self.indentify_task_configuration()

        elif task == 'reidentify':
            self.reidentify_task_configuration()
        
        elif task == 'fitcoords':
            self.fitcoords_task_configuration()

        elif task == 'refspectra':
            self.refspectra_task_configuration()

        elif task == 'dispcor':
            self.dispcor_task_configuration()

        elif task == 'transform':
            self.transform_task_configuration()
            
        elif task == 'background':
            self.background_task_configuration()

        elif task == 'apall':
            self.apall_task_configuration()

        elif task == 'standard':
            self.standard_task_configuration()

        elif task == 'sensfunc':
            self.sens_task_configuration()

        elif task == 'calibrate':
            self.calibrate_task_configuration()

        elif task == 'continuum':
            self.continuum_task_configuration()

        elif task == 'splot':
            self.splot_task_configuration()

        elif task == 'sarith':
            self.sarith_task_configuration()

        elif task == 'imarith':
            self.imarith_task_configuration()

        elif task == 'imstat':
            self.imstat_task_configuration()

        elif task == 'dopcor':
            self.dopcor_task_configuration()
   
        return
    
    def zerocombine_task_configuration(self):
                
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['input']      = self.task_attributes['input']
            self.task_data_dict['output']     = self.task_attributes['output']
            self.task_data_dict['combine']    = self.task_attributes['combine']
            self.task_data_dict['reject']     = 'minmax'
            self.task_data_dict['ccdtype']    = '""'
            self.task_data_dict['process']    = 'no'
            self.task_data_dict['delete']     = 'no'
        
        return

    def ccdproc_task_configuration(self):

        #Load default configuraion
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['images']   = self.task_attributes['input']
            self.task_data_dict['output']   = self.task_attributes['output'] #Need a list of files here
            self.task_data_dict['ccdtype']  = '""'
            self.task_data_dict['readaxi']  = 'column'

        #for ccdproc we set all the defaults to no
        for task_attrib in ['fixpix', 'oversca', 'trim', 'zerocor', 'darkcor', 'flatcor', 'illumco', 'fringec', 'readcor', 'scancor']:
            if task_attrib in self.task_attributes:
                self.task_data_dict[task_attrib] = self.task_attributes[task_attrib]
            else:
                self.task_data_dict[task_attrib] = 'no'

        #Load specific files for each treatment                
        for task_attrib in ['fixfile','biassec', 'trimsec', 'zero', 'dark', 'flat', 'illum', 'fringe']:
            if task_attrib in self.task_attributes:
                self.task_data_dict[task_attrib] =  self.task_attributes[task_attrib]
            else:
                self.task_data_dict[task_attrib] = '""'


        return

    def flatcombine_task_configuration(self):
 
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['input']    = self.task_attributes['input']
            self.task_data_dict['output']   = self.task_attributes['output'] #Flat_combined
            self.task_data_dict['combine']  = self.task_attributes['combine']
            self.task_data_dict['reject']   = self.task_attributes['reject']
            self.task_data_dict['ccdtype']  = self.task_attributes['ccdtype']
            self.task_data_dict['scale']    = self.task_attributes['scale']
            self.task_data_dict['gain']     = self.task_attributes['gain'] 
            self.task_data_dict['snoise']   = self.task_attributes['snoise']
        
        return

    def response_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['calibrat']  = self.task_attributes['input'] #Object_flat or combined flat       
            self.task_data_dict['normaliz']  = self.task_attributes['normalizing_flat'] #Combined_flat       
            self.task_data_dict['response']  = self.task_attributes['output'] #Object_n_flat or n_combined_flat
            self.task_data_dict['functio']   = 'spline3'
            self.task_data_dict['order']     = self.task_attributes['order']
            self.task_data_dict['low_rej']   = '3.0'
            self.task_data_dict['high_rej']  = '3.0'
            self.task_data_dict['niterate']  = '1'
            self.task_data_dict['interact']  = 'yes'
            self.task_data_dict['threshold'] = self.task_attributes['threshold']

            
        return

    def imcombine_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['input']    = self.task_attributes['input']
            self.task_data_dict['output']   = self.task_attributes['output']
            self.task_data_dict['combine']    = self.task_attributes['combine']
            self.task_data_dict['scale']      = self.task_attributes['scale']
            self.task_data_dict['statsec']    = self.task_attributes['statsec']
            self.task_data_dict['reject']     = self.task_attributes['reject']
            self.task_data_dict['weight']     = self.task_attributes['weight']
            self.task_data_dict['gain']       = self.task_attributes['gain']
            self.task_data_dict['snoise']     = self.task_attributes['snoise']
            self.task_data_dict['sigma']      = self.task_attributes['output'].replace('.fits', '_sigma.fits')

        return

    def illumination_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['images']       = self.task_attributes['input']
            self.task_data_dict['illumination'] = self.task_attributes['output'] #This is the illumFlat
            self.task_data_dict['interact']     = 'yes' 
            self.task_data_dict['nbins']        = 5
            self.task_data_dict['low_rej']      = '3.0'
            self.task_data_dict['high_rej']     = '3.0'        
        
        return

    def imcopy_task_configuration(self):

        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['input']    = self.task_attributes['input']
            self.task_data_dict['output']   = self.task_attributes['output']
            self.task_data_dict['verbose']  = 'yes'
    
        return

    def imshift_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['input']    = self.task_attributes['input']
            self.task_data_dict['output']   = self.task_attributes['output']
            self.task_data_dict['xshift']   = self.task_attributes['xshift']
            self.task_data_dict['yshift']   = self.task_attributes['yshift']
    
        return        

    def indentify_task_configuration(self):
        
        self.task_data_dict = OrderedDict()

        if self.objects_configuration == 'HII_galaxies':
            
            if self.task_attributes['color'] == 'Blue':
                coordlist_address = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/CuNeCuAr_ArcLampBlue.csv'
            else:
                coordlist_address = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/CuNeCuAr_ArcLampRed.csv'

            self.task_data_dict['images']     = self.task_attributes['input']
            self.task_data_dict['database']   = self.task_attributes['database']
            self.task_data_dict['section']    = 'middle column'
            self.task_data_dict['niterat']    = 1
            self.task_data_dict['fwidth']     = 7
            self.task_data_dict['coordlist']  = coordlist_address
            self.task_data_dict['match']      = 10
            self.task_data_dict['maxfeat']    = 75
            self.task_data_dict['fwidth']     = 7.5
            self.task_data_dict['cradius']    = 5
            self.task_data_dict['threshold']  = 10
            self.task_data_dict['minsep']     = 2
            
        return

    def reidentify_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            
            if self.task_attributes['color'] == 'Blue':
                coordlist_address = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/CuNeCuAr_ArcLampBlue.csv'
            else:
                coordlist_address = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/CuNeCuAr_ArcLampRed.csv'            
            
            self.task_data_dict['referenc']   = self.task_attributes['referenc']
            self.task_data_dict['images']     = self.task_attributes['input']
            self.task_data_dict['section']    = 'middle column'
            self.task_data_dict['coordlist']  = coordlist_address
            self.task_data_dict['interac']    = 'no'
            self.task_data_dict['overrid']    = 'yes'
            self.task_data_dict['trace']      = 'no'
            self.task_data_dict['nlost']      = 5
            self.task_data_dict['threshold']  = 10
            self.task_data_dict['match']      = 10
            self.task_data_dict['verbose']    = 'yes'

        return

    def fitcoords_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['images']   = self.task_attributes['input']
            self.task_data_dict['fitname']  = self.task_attributes['fitname']
            self.task_data_dict['interac']  = 'yes'

        return
    
    def refspectra_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['images']   = self.task_attributes['input']
            self.task_data_dict['fitname']  = self.task_attributes['fitname']
            self.task_data_dict['interac']  = 'yes'

        return

    def dispcor_task_configuration(self):

        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':    
            self.task_data_dict['input']    = self.task_attributes['input']
            self.task_data_dict['output']   = self.task_attributes['output']
            self.task_data_dict['lineari']  = self.task_attributes['lineari']
            self.task_data_dict['lineari']  = self.task_attributes['lineari']
            self.task_data_dict['flux']     = self.task_attributes['flux']
            self.task_data_dict['global']   = self.task_attributes['global']
    
        return    
    
    def transform_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['input']    = self.task_attributes['input']
            self.task_data_dict['output']   = self.task_attributes['output']
            self.task_data_dict['fitnames'] = self.task_attributes['fitnames']
            self.task_data_dict['flux']     = 'yes' #is this right?        
        
        return
    
    def background_task_configuration(self):

        self.task_data_dict = OrderedDict()
        if self.objects_configuration       == 'HII_galaxies':
            self.task_data_dict['input']    = self.task_attributes['input']
            self.task_data_dict['output']   = self.task_attributes['output']
            self.task_data_dict['function'] = 'chebyshev'
            self.task_data_dict['order']    = self.task_attributes['order']  
            self.task_data_dict['axis']     = self.task_attributes['axis']
  
        return
        
    def apall_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':
            self.task_data_dict['input']        = self.task_attributes['input']
            self.task_data_dict['output']       = self.task_attributes['output']
            self.task_data_dict['find']         = 'yes'
            self.task_data_dict['trace']        = self.task_attributes['trace']
            self.task_data_dict['resize']       = self.task_attributes['resize']
            self.task_data_dict['recente']      = 'yes'
            self.task_data_dict['referen']      = self.task_attributes['referen']
            self.task_data_dict['edit']         = self.task_attributes['edit']
            self.task_data_dict['interactive']  = 'yes'
            self.task_data_dict['extract']      = 'yes'
            self.task_data_dict['review']       = 'yes'
            self.task_data_dict['ylevel']       = self.task_attributes['ylevel']   
            
            
            self.task_data_dict['nsum']         = self.task_attributes['nsum']
            
            self.task_data_dict['b_sampl']      = self.task_attributes['b_sample'] 
#             self.task_data_dict['b_naver']       = '-100'
#             self.task_data_dict['b_order']       = 2
#             
#             self.task_data_dict['width']         = 10
            
            self.task_data_dict['maxsep']       = 1000
            self.task_data_dict['b_order']      = self.task_attributes['b_order']
            self.task_data_dict['t_niter']      = 1
            self.task_data_dict['backgro']      = self.task_attributes['backgro']
            self.task_data_dict['saturation']   = '32400'            
            self.task_data_dict['gain']         = self.task_attributes['gain_key']
            self.task_data_dict['readnoi']      = self.task_attributes['readnois_key']
            self.task_data_dict['extras']       = self.task_attributes['extras']
            self.task_data_dict['line']         = self.task_attributes['line']
            self.task_data_dict['order']        = self.task_attributes['order']
                  
    def standard_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':        
            Observatory             = 'lapalma'
            ExtinctionFileAddress   = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/a_ing_ext.dat'    

            self.task_data_dict['input']          = self.task_attributes['input']
            self.task_data_dict['output']         = self.task_attributes['output']
            self.task_data_dict['samestar']       = 'yes'
            self.task_data_dict['beam_switch']    = 'no'
            self.task_data_dict['apertures']      = '""' 
            self.task_data_dict['bandwidth']      = self.task_attributes['bandwidth']
            self.task_data_dict['bandsep']        = self.task_attributes['bandsep']
            self.task_data_dict['fnuzero']        = 3.68000e-20
            self.task_data_dict['extinction']     = ExtinctionFileAddress
            self.task_data_dict['caldir']         = self.Standard_Stars_Folder
            self.task_data_dict['observatory']    = Observatory
            self.task_data_dict['interact']       = 'yes'
            self.task_data_dict['star_name']      = self.task_attributes['star_name']              #MUST NOT INCLUDE EXTENSION IN THE STAR CALIBRATION FILE, NEITHER CAPITALS!!!!
            self.task_data_dict['airmass']       = self.task_attributes['airmass']
            self.task_data_dict['exptime']       = self.task_attributes['exptime']
            self.task_data_dict['answer']         = 'yes'
        
    def sens_task_configuration(self):

        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':     
      
            Observatory             = 'lapalma'
            ExtinctionFileAddress   = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/a_ing_ext.dat'    
            
            self.task_data_dict['standards']      = self.task_attributes['input']
            self.task_data_dict['sensitivity']    = self.task_attributes['output']
            self.task_data_dict['apertures']      = '""'
            self.task_data_dict['ignoreaps']      = "yes"
            self.task_data_dict['extinction']     = ExtinctionFileAddress
            self.task_data_dict['observatory']    = Observatory
            self.task_data_dict['functio']        = self.task_attributes['functio']
            self.task_data_dict['order']          = self.task_attributes['order']
            self.task_data_dict['interactive']    = 'yes'
            self.task_data_dict['graphs']         = self.task_attributes['graphs']
            self.task_data_dict['answer']         = 'yes'        
        
        return

    def calibrate_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        if self.objects_configuration == 'HII_galaxies':             
        
            Observatory             = 'lapalma'
            ExtinctionFileAddress   = '/home/vital/Dropbox/Astrophysics/Tools/PyRaf/Telescope_ReductionFiles/WHT/a_ing_ext.dat'    
                           
            self.task_data_dict['input']         = self.task_attributes['input']
            self.task_data_dict['output']        = self.task_attributes['output']
            self.task_data_dict['extinct']       = 'yes'
            self.task_data_dict['flux']          = 'yes'
            self.task_data_dict['extinction']    = ExtinctionFileAddress
            self.task_data_dict['observatory']   = Observatory
            self.task_data_dict['ignoreaps']     = 'yes'
            self.task_data_dict['sensitivity']   = self.task_attributes['senstivityCurve']
            self.task_data_dict['fnu']           = 'no'
            self.task_data_dict['airmass']       = self.task_attributes['airmass']
            self.task_data_dict['exptime']       = self.task_attributes['exptime']
            self.task_data_dict['mode']          = 'ql'

            return
       
    def continuum_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        
        self.task_data_dict['input']     = self.task_attributes['input']
        self.task_data_dict['output']    = self.task_attributes['output']    
        
        return

    def sarith_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        
        self.task_data_dict['input1']   = self.task_attributes['input1']
        self.task_data_dict['op']       = self.task_attributes['op']
        self.task_data_dict['input2']   = self.task_attributes['input2']
        self.task_data_dict['output']   = self.task_attributes['output']
        self.task_data_dict['w1']       = 'INDEF'
        self.task_data_dict['w2']       = 'INDEF'
        self.task_data_dict['rebin']    = 'no'
        self.task_data_dict['verbose']  = 'yes'        
        
        return

    def imarith_task_configuration(self):
        
        self.task_data_dict = OrderedDict()
        
        self.task_data_dict['operand1'] = self.task_attributes['operand1']
        self.task_data_dict['op']       = self.task_attributes['op']
        self.task_data_dict['operand2'] = self.task_attributes['operand2']
        self.task_data_dict['result']   = self.task_attributes['output']
        self.task_data_dict['verbose']  = 'yes'        
        
        return

    def splot_task_configuration(self):
    
        self.task_data_dict = OrderedDict()
        self.task_data_dict['images']   = self.task_attributes['input']
        self.task_data_dict['xmin']     = self.task_attributes['xmin'] if 'xmin' in self.task_attributes  else 'INDEF'
        self.task_data_dict['xmax']     = self.task_attributes['xmax'] if 'xmax' in self.task_attributes  else 'INDEF'
        self.task_data_dict['ymax']     = self.task_attributes['ymax'] if 'ymax' in self.task_attributes  else 'INDEF'

    def imstat_task_configuration(self):
        self.task_data_dict = OrderedDict()
        self.task_data_dict['images']   = self.task_attributes['input']
        self.task_data_dict['lower']    = 'INDEF'
        self.task_data_dict['upper']    = 'INDEF'

    def dopcor_task_configuration(self):
        self.task_data_dict = OrderedDict()
        self.task_data_dict['input']    = self.task_attributes['input']
        self.task_data_dict['output']   = self.task_attributes['output']
        self.task_data_dict['redshift'] = self.task_attributes['redshift']
        self.task_data_dict['flux']     = self.task_attributes['flux']
        self.task_data_dict['verbose']  = 'yes'
                                                       
class spectra_reduction(fits_plots, pyraf_task_configuration):
    
    def __init__(self):

        #Declare the support classes
        fits_plots.__init__(self)
        pyraf_task_configuration.__init__(self)

        #Declare the catalogue (This should be done from a graphical intervace)
        self.get_observation_folder()

        self.observation_properties_file_name   = 'observation_properties.txt'
        
        #Create text file to store the iraf commands
        self.commands_log = 'commands_log.txt'
        if isfile(self.Catalogue_folder + self.commands_log) == False:
            with open (self.Catalogue_folder + self.commands_log, 'a') as f: f.write ('Catalogue ({cataloge_code}) reduction commands\n'.format(cataloge_code = self.Catalogue_folder))
            
        #Data location
        self.red_root_folder            = None
        self.testing_extension          = '_TestingContinua_3slits' #'_mediumSlits' #'_narrowSlits' #'_test1' #'_HeliumSlit'
    
        #Frame with the reduced data
        self.reducDf                    = None
        self.frame_filename             = 'reduction_frame.dz'
        
        #Flags
        self.verbose_output             = True
        
        #Fits extensions
        self.Extensions_dict = {}
        self.Extensions_dict['bias']    = '_b'
        self.Extensions_dict['flat']    = '_f'
        self.Extensions_dict['arc']     = '_w'
        self.Extensions_dict['flux']    = '_fx'
        
        #Organized folders
        self.default_Folders_dict = {}
        self.default_Folders_dict['bias']       = 'bias'
        self.default_Folders_dict['flat lamp']  = 'flats_lamp'
        self.default_Folders_dict['flat sky']   = 'flats_sky'
        self.default_Folders_dict['arcs']       = 'arcs'
        self.default_Folders_dict['objects']    = 'objects'
        self.default_Folders_dict['raw data']   = 'raw_fits'
        self.default_Folders_dict['telescope']  = 'conf_frames'
        self.default_Folders_dict['reduc_data'] = 'reduc_data'
    
    def get_observation_folder(self):
        
        self.Catalogue_folder   ='/home/vital/Astrodata/WHT_2011_11_er_propagation/Night1/'
        self.Catalogue_folder   = '/home/vital/Astrodata/WHT_2016_04/Night1/'
        #self.Catalogue_folder  = '/home/vital/Astrodata/WHT_2008_01/Night2/'
        #self.Catalogue_folder  = '/home/vital/Astrodata/WHT_2009_07/Night1/'
        #self.Catalogue_folder                 = '/home/vital/Astrodata/WHT_2011_01/Night2/'
        #self.Catalogue_folder                  = '/home/vital/Astrodata/WHT_2011_09/Night1/'
        #self.Catalogue_folder                 = '/home/vital/Astrodata/WHT_2011_11/Night1/'
        #self.Catalogue_folder                 = '/home/vital/Astrodata/WHT_2016_10/'
      
        return 
             
    def target_validity_check(self):

        #Check if the targets are valid and they belong to the objects we want to treat
        if (self.Objects_to_treat is None) or (self.Objects_to_treat is 'All'):
            boolean_array = (self.reducDf.valid_file) & (self.reducDf.frame_tag.isin(self.Objects_to_treat))
        else:
            boolean_array = (self.reducDf.valid_file)
        
        return boolean_array
    
    def Folder_Explorer(self, myPattern, Folder, Sort_Output = None, verbose = True):
                    
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

        if Sort_Output == 'alphabetically':
            
            FileList = sorted(FileList)
                                               
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

    def get_spectra_data(self, file_address, ext=0):

        Header_0    = fits.getheader(file_address, ext=ext)
        Flux_array  = fits.getdata(file_address, ext=ext)

        if "COEFF0" in Header_0:
            dw              = 10.0**Header_0['COEFF1']                # dw = 0.862936 INDEF (Wavelength interval per pixel)
            Wmin            = 10.0**Header_0['COEFF0']
            pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
            Wmax            = Wmin + dw * pixels
            wavelength      = linspace(Wmin,Wmax,pixels,endpoint=False)
              
        elif "LTV1" in Header_0:
            StartingPix     = -1 * Header_0['LTV1']                   # LTV1 = -261. 
            Wmin_CCD        = Header_0['CRVAL1']
            dw              = Header_0['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
            pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
            Wmin            = Wmin_CCD + dw * StartingPix
            Wmax            = Wmin + dw * pixels
            wavelength      = linspace(Wmin,Wmax,pixels,endpoint=False)
        
        else:
            Wmin            = Header_0['CRVAL1']
            dw              = Header_0['CD1_1']                       # dw = 0.862936 INDEF (Wavelength interval per pixel)
            pixels          = Header_0['NAXIS1']                      # nw = 3801 number of output pixels
            Wmax            = Wmin + dw * pixels
            wavelength      = linspace(Wmin,Wmax,pixels,endpoint=False)
    
        return wavelength, Flux_array, Header_0

    def declare_folde_structre(self):

        self.reducFolders       = {}
        for key in self.default_Folders_dict.keys():
            self.reducFolders[key] = self.reduc_RootFolder + self.default_Folders_dict[key] + '/'
            if isdir(self.reducFolders[key]) == False:
                makedirs(self.reducFolders[key])

        return

    def save_observation_dict(self):
        
        #Separate the parameters
        parameters_keys, parameters_values = self.observation_dict.keys(), self.observation_dict.values()
        
        #Format the entries
        for i in range(len(parameters_values)):
            if isinstance(parameters_values[i], list):
                parameters_values[i] = ' '.join(map(str, parameters_values[i]))
        
        #Saving the file
        observation_properties_file = self.Catalogue_folder + self.observation_properties_file_name     
        savetxt(fname = observation_properties_file, X = transpose([parameters_keys, parameters_values]), fmt='%s', delimiter = '; ')

        return
    
    def load_observation_dict(self, observation_properties_address):
        
        observation_dict = OrderedDict()
        observation_properties_file = loadtxt(observation_properties_address, dtype=str, delimiter = ';', usecols = [0,1])     
        for i in range(len(observation_properties_file)):
            key         =  observation_properties_file[i][0]
            list_values = observation_properties_file[i][1].split()
            observation_dict[key] = list_values

        return observation_dict

    def load_reduction_dataframe(self, dataframe_folder, dataframe_name = 'reduction_frame.dz'):

        reduction_dataframe = None

        #Pandas frame where we store all the reduction data
        if isfile(dataframe_folder + dataframe_name):
            reduction_dataframe     = read_pickle(dataframe_folder + dataframe_name)
            
            #Special trick to update the location 
            self.columns_reducDf    = reduction_dataframe.index.values
            
            #Check if files have been added to the rejection list
            self.check_rejected_files(reduction_dataframe, dataframe_folder)
            
            #Check if files were deleted and remove then from the data frame
            to_delete   = []
            addresses   = reduction_dataframe.file_location.values.astype(str)
            names       = reduction_dataframe.file_name.values.astype(str)
            frame_address = np_f.add(addresses, names)
 
            for i in range(len(frame_address)):
                if isfile(frame_address[i]) == False:
                    to_delete.append(names[i]) 
                     
            mask_delete         = logical_not(in1d(reduction_dataframe['file_name'], to_delete)) 
            indeces_to_delete   = reduction_dataframe.index[mask_delete]
            reduction_dataframe = reduction_dataframe.loc[indeces_to_delete]        

        return reduction_dataframe

    def declare_catalogue(self, catalogue_address, data_origin = 'WHT', objects_type = 'HII galaxies', verbose = True):

        #If you want to treat only one file please add it here
        self.Objects_to_treat = None

        #Declare catalogue addresses
        if catalogue_address is not None:
            self.reduc_RootFolder = catalogue_address
            self.Catalogue_folder = catalogue_address
        else:
            self.reduc_RootFolder = self.Catalogue_folder

        #Declare the main dictionary addresses
        self.declare_folde_structre()

        #Load the observation characteristics for the reduction from
        self.observation_dict = self.load_observation_dict(self.Catalogue_folder + self.observation_properties_file_name)

        #Load objects to treat:
        if self.Objects_to_treat == None:
            self.Objects_to_treat = self.observation_dict['objects'] + self.observation_dict['Standard_stars']

        #Load reduction dataframe
        self.reducDf = self.load_reduction_dataframe(catalogue_address)

        #Dictionary with the keys which correspond to a given telescope header
        if data_origin == 'WHT':
            self.columns_reducDf = ['OBJECT', 'OBSTYPE', 'file_name', 'RUN', 'file_location', 'reduc_tag', 'frame_tag', 'ISIARM', 'RA', 'DEC', 'UT', 'EXPTIME', 'AIRMASS', 'ISISLITW']

            #Case this is the first time
            if self.reducDf is None:
                self.reducDf = DataFrame(columns = self.columns_reducDf)
        if verbose:
            print 'Reduction dataframe {catalogue_address} loaded'.format(catalogue_address = catalogue_address + 'reduction_frame.dz')

    def check_rejected_files(self, reduction_dataframe, dataframe_folder, reject_file_name = 'rejected_files.txt'):
        
        #Set all to True
        reduction_dataframe['valid_file'] = True
        
        #Load list of files
        rejected_file = dataframe_folder + 'rejected_files.txt'
        list_rejected = loadtxt(rejected_file, dtype = str, usecols = [0], ndmin=1)
        
        #Set the column to false
        reduction_dataframe.loc[reduction_dataframe.file_name.isin(list_rejected), 'valid_file'] = False

    def generate_step_pdf(self, indeces_frame, file_address, plots_type = 'fits', ext = 1, include_graph = False, verbose = True, limits = None, sorting_pattern = ['frame_tag','ISIARM']):
                    
        if verbose:
            print 'Printing these files: {number}'.format(number = indeces_frame.sum())
            print self.reducDf[indeces_frame].sort_values(['frame_tag']).file_name.values

#------------Plotting pairs of fits
        if plots_type == 'fits':
            self.fits_catalogue(file_address, indeces_frame, ext = ext, columns_mean=include_graph, sorting_pattern=sorting_pattern)

#------------Plotting flats within limits of fits
        if plots_type == 'flats':
            self.flats_plotting(file_address, indeces_frame, ext = ext, columns_mean=include_graph)

#------------Plotting flats within limits of fits
        if plots_type == 'flat_calibration':
            self.objects_focus(file_address, indeces_frame, ext = ext)

#------------Plotting extracted frames:
        if plots_type == 'extraction':
            self.extracted_frames(file_address, indeces_frame, ext = ext)

#------------Plotting the fast combined with notes        
        if plots_type == 'cosmic_removal':          
            self.cosmic_plot(file_address, indeces_frame, ext = ext, columns_mean=include_graph)

#------------Plotting pairs of fits
        if plots_type == 'fits_compare':        
            self.fits_compare(file_address, indeces_frame, ext = ext, columns_mean=include_graph)

#------------Plotting the spectra        
        if plots_type == 'spectra':
            self.spectra(file_address, indeces_frame, ext)

#------------Plotting the spectra        
        if plots_type == 'masked_pixels':
            self.masked_pixels(file_address, indeces_frame, ext = ext, columns_mean=include_graph)

#------------Fast combine
        if plots_type == 'fast_combine':             
            self.fast_combined(file_address, indeces_frame, ext = 0, observation_dict = self.observation_dict, reducDf = self.reducDf, colors_dict = self.frames_colors)

#------------Comparing the combined frames                
        if plots_type == 'frame_combine':              
            self.frame_combine(file_address, indeces_frame, ext = 0)

#------------Comparing the combined frames                
        if plots_type == 'frame_combine_shifted':              
            self.frame_combine_shifted(file_address, indeces_frame, ext = 0)
                               
    def save_task_parameter(self, task, parameter, entry):
        
        import pyraf
        
        if task == 'response':
#             from pyraf.iraf.noao.twodspec import response
#             response()
            
            p_list = pyraf.iraf.noao.twodspec.longslit.response.getParList()

        #Getting the parameter form the task
        par = p_list[5].get()
                
        #Clean the format
        if par[0] == ' ':
            entry_value = par[1:].split()
        else:
            entry_value = par[0:].split()        

        #Save to the observation properties file
        self.observation_dict[entry] = entry_value
        self.save_observation_dict()
        
    def run_iraf_task(self, task, overwrite = True, run_externally = True, verbose = True):
        
        #Establish configuration folder
        configuration_file_address  = self.task_attributes['run folder'] + '{task}_{color}_conf.txt'.format(task = task, color = self.task_attributes['color'])
        
        #Establish the task configuration
        self.load_task_configuration(task)
                
        #Save the input file list if necessary
        if 'input array' in self.task_attributes:
            savetxt(self.task_attributes['run folder'] + self.task_attributes['in_list_name'], self.task_attributes['input array'], fmt='%s')
        if 'output array' in self.task_attributes:
            savetxt(self.task_attributes['run folder'] + self.task_attributes['out_list_name'], self.task_attributes['output array'], fmt='%s')

        #Save configuration data
        parameter_list, value_list  = self.task_data_dict.keys(), self.task_data_dict.values()
        savetxt(fname = configuration_file_address, X = transpose([parameter_list, value_list]), fmt='%s', delimiter = ';')
        
        #Destroy file if it already exists
        if ('output' in self.task_attributes) and overwrite:
            if isfile(self.task_attributes['output']):
                remove(self.task_attributes['output'])
        
        #Destroy files if they are in a list
        if ('output array' in self.task_attributes) and overwrite:
            for output_file in self.task_attributes['output array']:
                clean_name = output_file.replace('[1]','')
                if isfile(clean_name):
                    remove(clean_name)
                    
        #Destroy files if the output is an IRAF list (, separated) #This is a bit dirty
        if ('output' in self.task_attributes) and overwrite:
            if ',' in self.task_attributes['output']:
                output_files = self.task_attributes['output'].split(',')
                for output_file in output_files:
                    clean_name = output_file.replace('[1]','')
                    if isfile(clean_name):
                        remove(clean_name)            
                    
        #Run the task
        launch_command = 'python ' + self.script_location + '0_RunIraf_task.py'  + ' ' + configuration_file_address 
        
#         if run_externally:
#             print '-- Launching command:', launch_command
#             p = subprocess.Popen(launch_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#             out = p.communicate()[0]
#             print out
#         else:

        #Always run internally
        execute_Iraf_task(launch_command, commands_log= self.Catalogue_folder + self.commands_log, verbose_output = verbose)
        
        return
            
    def save_reducDF(self):
        
        self.reducDf.to_pickle(self.reduc_RootFolder + self.frame_filename)
        
        return

    def create_skycor_fits(self, idx_file, output_folder = None, ext = 0):
                
        frame_name      = self.reducDf.loc[idx_file].file_name.values[0]
        frame_folder    = self.reducDf.loc[idx_file].file_location.values[0]
        target_code     = self.reducDf.loc[idx_file].frame_tag.values[0].replace('[','').replace(']','')
        color_arm       = self.reducDf.loc[idx_file].ISIARM.values[0]
        
        #Define folder where Skycorr is run
        if output_folder == None:
            output_folder = '/home/vital/Skycorr/WHT_fittings/'

        #Get fits data (WHT case)
        print 'file', frame_folder + frame_name
        wavelength, Flux_array, Header_0 = self.get_spectra_data(frame_folder + frame_name)
        UT_start    = Header_0['UTSTART'].split(':')
        UT_start_s  = timedelta(hours=int(UT_start[0]),minutes=int(UT_start[1]),seconds=float(UT_start[2])).total_seconds()
        
        #--Calculate the median sky
        parent_file = frame_name.replace('_eSkyCorr.fits', '.fits')
        apall_file  = frame_folder + 'database/ap' + frame_folder.replace('/', '_') + parent_file.replace('.fits', '')
        #background_flux = self.compute_background_median(frame_folder + parent_file, apall_file)
        Flux_Obj        = Flux_array[0][0]
        Flux_Sky        = Flux_array[1][0]
        Flux_Combined   = Flux_array[0][0] + Flux_array[1][0]
        
#         Fig, Axis = plt.subplots(1, 1, figsize=(10, 12))
#         Axis.plot(wavelength, Flux_Sky, label='sky')
#         Axis.plot(wavelength, Flux_Combined, label='target')
#         Axis.legend()
#         plt.show()
        
        #Generate the header
        prihdr                  = fits.Header()
        prihdr['MJD-OBS']       = float(Header_0['MJD-OBS'])
        prihdr['TM-START']      = float(UT_start_s)
        prihdr['ESO TEL ALT']   = float(Header_0['LATITUDE'])
        prihdu = fits.PrimaryHDU(header=prihdr)
                
        #Generate the sky fits
        colA_1      = fits.Column(name='lambda',format='F', array=wavelength)
        colA_2      = fits.Column(name='flux', format='F', array=Flux_Sky)
        colsA       = fits.ColDefs([colA_1, colA_2])
        tbhduA      = fits.BinTableHDU.from_columns(colsA, header=prihdr)
        thdulistA   = fits.HDUList([prihdu, tbhduA])
        skyfits_address = '/home/vital/Skycorr/WHT_fittings/Input_sky/' + target_code + '_background.fits'
        thdulistA.writeto(skyfits_address, clobber=True)
  
        #Generate the object + sky fits
        colB_1      = fits.Column(name='lambda', format='F', array=wavelength)
        colB_2      = fits.Column(name='flux', format='F', array=Flux_Combined)
        colsB       = fits.ColDefs([colB_1, colB_2])
        tbhduB      = fits.BinTableHDU.from_columns(colsB, header=prihdr)
        thdulistB   = fits.HDUList([prihdu, tbhduB])
        targetfits_address = '/home/vital/Skycorr/WHT_fittings/Input_objects/' + target_code + '_spectrum_and_background.fits'
        thdulistB.writeto(targetfits_address, clobber=True)

        #Generate configuration file
        template    = '/home/vital/Skycorr/WHT_fittings/template_parameter_file.par'
        conf_file   = '/home/vital/Skycorr/WHT_fittings/{codeTarget}_{armcolor}_skycorr.in'.format(codeTarget = target_code, armcolor=color_arm.replace(' ','_'))
        output_file = '{codeTarget}_{armcolor}_skycorfit'.format(codeTarget = target_code, armcolor=color_arm.replace(' ','_'))
        
        dict_modifications = {}
        dict_modifications['INPUT_OBJECT_SPECTRUM=']    = targetfits_address
        dict_modifications['INPUT_SKY_SPECTRUM=']       = skyfits_address
        dict_modifications['OUTPUT_NAME=']              = output_file
        dict_keys = dict_modifications.keys()
        
        #Temporary file
        fh, abs_path = mkstemp()
        with open(abs_path,'w') as new_file:
            with open(template) as old_file:
                for line in old_file.readlines():
                    if line.replace('\n', '') in dict_keys:
                        formatted_line = line.replace('\n', '')
                        new_line = formatted_line + dict_modifications[formatted_line] + '\n'
                    else:
                        new_line = line
                    new_file.write(line.replace(line, new_line))
                    
        close(fh)

        #Move new file
        move(abs_path, conf_file)
                
        #Run skycor
        launch_command = '/home/vital/Skycorr/bin/skycorr {skycorr_script}'.format(skycorr_script = conf_file)
        p = subprocess.Popen(launch_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out = p.communicate()[0]
        print out
     
        return

    def object_to_dataframe(self, list_new_file_address, data_dict, new_entry_type = False):
        
        if type(list_new_file_address) is str:
            list_new_file_address = [list_new_file_address]
        
        for i in range(len(list_new_file_address)):
            
            new_file_address = list_new_file_address[i]
            
            #Identify frame
            CodeName, FileName, FileFolder = self.Analyze_Address(new_file_address, verbose=False)
            CodeName = FileName[0:FileName.find('.')]
            
            #Get headers
            Header0 = getheader(new_file_address, ext=0)
                
            #Add variables
            for key in self.columns_reducDf:
                if key in Header0:
                    self.reducDf.loc[CodeName, key] = Header0[key]
            
            #Add file location to the columns
            self.reducDf.loc[CodeName, 'file_name']      = FileName
            self.reducDf.loc[CodeName, 'file_location']  = FileFolder
            self.reducDf.loc[CodeName, 'valid_file']     = True
            
            #Adjusting corresponding variables
            for key in data_dict:
                
                #If key not in dictionary you can create it and set it to false
                if key not in self.reducDf.columns:
                    self.reducDf[key] = new_entry_type
                    
                self.reducDf.loc[CodeName, key]  = data_dict[key]
    
            #Add the frame_tag to be the same as its corresponding number
            run_number      = self.reducDf.loc[CodeName, 'RUN']
            run_firstindex  = (self.reducDf.RUN == run_number)
            frame_tag_value = self.reducDf.loc[run_firstindex, 'frame_tag'].values[0]
            self.reducDf.loc[CodeName, 'frame_tag'] = frame_tag_value
     
            self.save_reducDF()
                    
        return

    def printIrafCommand(self, Task, Attributes, printindividually = True):
        
        keys_attrib     = Attributes.keys()
        values_attrib   = Attributes.values()
    
        command_String  = Task
                
        print '--Whole command'
        for i in range(len(keys_attrib)):
            command_String = command_String + ' ' + keys_attrib[i] + '=' + str(values_attrib[i])
        
        print '--Task atributtes:'
        for i in range(len(keys_attrib)):
            print keys_attrib[i] + '=' + str(values_attrib[i])

        return command_String

    def get_closest_time(self, df, idx, bool_cond, to_this):
        others       = df.loc[bool_cond, to_this].values
        target       = datetime64(df.loc[idx, to_this])
        idx_closest  =  (abs(others-target)).argmin()
        closet_value = others[idx_closest]
        return df.loc[bool_cond & (df[to_this] == closet_value)].file_name.values[0],  closet_value

    def get_closest(self, df, idx, bool_cond, to_this):
        others       = df.loc[bool_cond, to_this].values
        target       = df.loc[idx, to_this]
        idx_closest  = (abs(others-target)).argmin()
        closet_value = others[idx_closest]
        return df.loc[bool_cond & (df[to_this] == closet_value)].file_name.values[0],  closet_value   

    def reset_task_dict(self):

        self.task_attributes = OrderedDict()

        return

    def beep_alarmn(self):
    
        #subprocess.call(['speech-dispatcher'])
        subprocess.call(['spd-say', '"your process has finished"'])
    
        return