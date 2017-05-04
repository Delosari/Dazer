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
    def frame_combine(self, pdf_address, indeces_frame, observation_dict, reducDf, ext = 0, colors_dict = None):
         
        files_name      = self.reducDf[indeces_frame].sort_values(['frame_tag', 'ISIARM']).file_name.values
        files_address   = self.reducDf[indeces_frame].sort_values(['frame_tag', 'ISIARM']).file_location.values
        frames_color    = self.reducDf[indeces_frame].sort_values(['frame_tag', 'ISIARM']).ISIARM.values
         
        #Generate an artificial list with the range of files
        Number_frames   = len(files_address)
        indices_list    = range(Number_frames)
 
#-------Plotting frames
        #Figure configuration
        Pdf_Fig, GridAxis   = plt.subplots(1, self.page_grid_width, figsize=(20, 24))  
        GridAxis_list       = GridAxis.ravel()
 
        #Create the doc
        doc = Document(pdf_address)
        doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=1cm', 'bottom=1cm']))
        doc.append(NoEscape(r'\extrafloats{100}')) 
        doc.append(NoEscape(r'\pagenumbering{gobble}'))
         
        #Start looping through the indeces in groups with the same lengths as the columns * rows
        for sublist in zip(*[iter(indices_list)] * self.page_grid_width):
                         
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
                    cmap_frame = colors_dict[frames_color[index_i]]
                    GridAxis_list[i].imshow(image_data, cmap=cmap_frame, origin='lower', vmin = int_min, vmax = int_max, interpolation='nearest', aspect='auto')
                    GridAxis_list[i].set_title(files_name[index_i], fontsize = 24)                    
 
            #Add the plot                               
            with doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                doc.append(NoEscape(r'\newpage'))
                plt.cla()
 
#-------Plotting maxima maps
        #Loop through the scientific objects
        objects_list = observation_dict['objects']
        for j in range(len(objects_list)):
             
            obj_target = objects_list[j] 
             
            Pdf_Fig, GridAxis   = plt.subplots(2, 3, figsize=(20, 24))  
            GridAxis_list       = GridAxis.ravel()
             
            #Loop through the arms    
            colors = ['Blue', 'Red']
            for k in range(len(colors)):
                arm_color = colors[k]
                 
                idx_objectFrames = (reducDf.frame_tag == obj_target) & (reducDf.reduc_tag == 'cosmic_ray_removal') & (reducDf.ISIARM == '{color} arm'.format(color = arm_color))
                idx_objCombine   = (reducDf.frame_tag == obj_target) & (reducDf.reduc_tag == 'obj_combine') & (reducDf.ISIARM == '{color} arm'.format(color = arm_color))   
                 
                frames_name      = self.reducDf[idx_objectFrames].sort_values(['file_name']).file_name.values
                frames_folders   = self.reducDf[idx_objectFrames].sort_values(['file_name']).file_location.values                
                 
                combframe_name   = self.reducDf[idx_objCombine].sort_values(['file_name']).file_name.values
                combframe_folder = self.reducDf[idx_objCombine].sort_values(['file_name']).file_location.values
                 
                frames_name         = list(frames_name) + list(combframe_name)
                frames_folders      = list(frames_folders) +list(combframe_folder)
                                 
                for m in range(len(frames_name)):
                     
                    with fits.open(frames_folders[m] + frames_name[m]) as hdu_list:
                        image_data = hdu_list[ext].data
                            
                    center_line_cords = observation_dict['{CodeName}_refline_{Color}'.format(CodeName=obj_target, Color=arm_color)]
                     
                    x_max, y_max = int(center_line_cords[0]), int(center_line_cords[1])
                    max_region = image_data[y_max-5:y_max+5:, x_max-5:x_max+5]
                                         
                    #Get zscale limits for plotting the image
                    IntensityLimits     = ZScaleInterval()
                    int_min, int_max    = IntensityLimits.get_limits(max_region)[0], IntensityLimits.get_limits(max_region)[1]
                     
                    cmap_frame = colors_dict[arm_color + ' arm']
                    GridAxis_list[m].imshow(max_region, cmap=cmap_frame, origin='lower', vmin = int_min, vmax = int_max, interpolation='nearest', aspect='auto')
                       
                    #Plot layout
                    GridAxis_list[m].set_xlabel('x',            fontsize = 30)
                    GridAxis_list[m].set_ylabel('y',            fontsize = 30)
                    GridAxis_list[m].set_title(frames_name[m],  fontsize = 30) 
                    GridAxis_list[m].tick_params(axis='both',   labelsize=20)
                    GridAxis_list[m].legend(loc='upper right',  prop={'size':18})
                    plt.axis('tight')
              
                #Add the plot                               
                with doc.create(Figure(position='htbp')) as plot:       
                    plot.add_plot(height=NoEscape(r'1\textheight'))
                    doc.append(NoEscape(r'\newpage'))
                    plt.cla()
 
 
 
#-------Reference_line_cut  
        objects_list        = observation_dict['objects']
        for j in range(len(objects_list)):
             
            obj_target = objects_list[j]
             
            Pdf_Fig, GridAxis   = plt.subplots(2, 1, figsize=(20, 24))  
            GridAxis_list       = GridAxis.ravel()
 
            #Loop through the arms    
            colors = ['Blue', 'Red']
            for k in range(len(colors)):
                arm_color = colors[k]
 
                idx_objectFrames = (reducDf.frame_tag == obj_target) & (reducDf.reduc_tag == 'cosmic_ray_removal') & (reducDf.ISIARM == '{color} arm'.format(color = arm_color))
                idx_objCombine   = (reducDf.frame_tag == obj_target) & (reducDf.reduc_tag == 'obj_combine') & (reducDf.ISIARM == '{color} arm'.format(color = arm_color))   
                 
                frames_name      = self.reducDf[idx_objectFrames].sort_values(['file_name']).file_name.values
                frames_folders   = self.reducDf[idx_objectFrames].sort_values(['file_name']).file_location.values                
                frame_validity   = self.reducDf[idx_objectFrames].sort_values(['file_name']).valid_file.values
                 
                combframe_name   = self.reducDf[idx_objCombine].sort_values(['file_name']).file_name.values
                combframe_folder = self.reducDf[idx_objCombine].sort_values(['file_name']).file_location.values
                 
                center_line_cords   = observation_dict['{CodeName}_refline_{Color}'.format(CodeName=obj_target, Color=arm_color)]
                 
                #Plot object fromes              
                for m in range(len(frames_name)):
                    with fits.open(frames_folders[m] + frames_name[m]) as hdu_list:
                        image_data = hdu_list[ext].data
                    x_max, y_max        = int(center_line_cords[0]), int(center_line_cords[1])
                    row_max             = image_data[y_max-1:y_max+1, :].mean(axis=0)
                    pixel_idx           = range(len(row_max))
 
                    #Check if file is rejected:
                    if frame_validity[m]:
                        label = frames_name[m]
                    else:
                        label = 'REJECTED ' + frames_name[m]
                     
                    GridAxis_list[k].step(pixel_idx, row_max, label = label)
             
                #Plot combined frame
                with fits.open(combframe_folder[0] + combframe_name[0]) as hdu_list:
                    image_data = hdu_list[ext].data            
                 
                x_max, y_max        = int(center_line_cords[0]), int(center_line_cords[1])
                row_max             = image_data[y_max-1:y_max+1, :].mean(axis=0)
                pixel_idx           = range(len(row_max))
                GridAxis_list[k].step(pixel_idx, row_max, label = obj_target + ' combined frame', color = 'black', linestyle=':')
             
                #Plot layout
                GridAxis[k].set_xlabel('Pixel value',          fontsize = 30)
                GridAxis[k].set_ylabel('Mean spatial count',   fontsize = 30)
                GridAxis[k].set_title(obj_target + ' ' + arm_color, fontsize = 30) 
                GridAxis[k].tick_params(axis='both', labelsize=20)
                GridAxis[k].legend(loc='upper right', prop={'size':18})
                GridAxis[k].set_xlim(x_max-40, x_max+40)
 
            #Add the plot                               
            with doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                doc.append(NoEscape(r'\newpage'))
                plt.cla() 
 
#-------Plotting spectra region  
        #Loop through the scientific objects
        objects_list        = observation_dict['objects']
        for j in range(len(objects_list)):
             
            obj_target = objects_list[j]
             
            Pdf_Fig, GridAxis   = plt.subplots(2, 1, figsize=(20, 24))  
            GridAxis_list       = GridAxis.ravel()
 
            linestyles          = ["-","-","--","--","-.","-.",":",":"]
            linecycler          = cycle(linestyles)
             
            #Loop through the arms    
            colors = ['Blue', 'Red']
            for k in range(len(colors)):
                arm_color = colors[k]
                 
                idx_objectFrames = (reducDf.frame_tag == obj_target) & (reducDf.reduc_tag == 'cosmic_ray_removal') & (reducDf.ISIARM == '{color} arm'.format(color = arm_color))
                idx_objCombine   = (reducDf.frame_tag == obj_target) & (reducDf.reduc_tag == 'obj_combine') & (reducDf.ISIARM == '{color} arm'.format(color = arm_color))   
                 
                frames_name      = self.reducDf[idx_objectFrames].sort_values(['file_name']).file_name.values
                frames_folders   = self.reducDf[idx_objectFrames].sort_values(['file_name']).file_location.values                
                frame_validity   = self.reducDf[idx_objectFrames].sort_values(['file_name']).valid_file.values
                 
                center_line_cords = observation_dict['{CodeName}_refline_{Color}'.format(CodeName=obj_target, Color=arm_color)]
                 
                combframe_name   = self.reducDf[idx_objCombine].sort_values(['file_name']).file_name.values
                combframe_folder = self.reducDf[idx_objCombine].sort_values(['file_name']).file_location.values
                                 
                for m in range(len(frames_name)):
                     
                    with fits.open(frames_folders[m] + frames_name[m]) as hdu_list:
                        image_data = hdu_list[ext].data
                            
                    x_max = int(center_line_cords[0])
                    y_values = image_data[:, x_max-2:x_max+2].mean(axis=1) 
                    x_values = range(len(y_values))
 
                    #Check if file is rejected:
                    if frame_validity[m]:
                        label = frames_name[m]
                    else:
                        label = 'REJECTED ' + frames_name[m]
                                         
                    GridAxis_list[k].plot(x_values, y_values, label = label, linestyle=next(linecycler))
 
                #Plot combined frame
                with fits.open(combframe_folder[0] + combframe_name[0]) as hdu_list:
                    image_data = hdu_list[ext].data            
                 
                x_max    = int(center_line_cords[0])
                y_values = image_data[:, x_max-2:x_max+2].mean(axis=1) 
                x_values = range(len(y_values))
                GridAxis_list[k].step(x_values, y_values, label = obj_target + ' combined frame', color = 'black', linestyle=':')
                      
                #Plot layout
                GridAxis[k].set_xlabel('Pixel value',          fontsize = 30)
                GridAxis[k].set_ylabel('Mean spatial count',   fontsize = 30)
                GridAxis[k].set_title(obj_target + ' ' + arm_color, fontsize = 30) 
                GridAxis[k].tick_params(axis='both', labelsize=20)
                GridAxis[k].legend(loc='upper right', prop={'size':18})
                GridAxis[k].set_yscale('log')
                 
            #Add the plot                               
            with doc.create(Figure(position='htbp')) as plot:       
                plot.add_plot(height=NoEscape(r'1\textheight'))
                doc.append(NoEscape(r'\newpage'))
                plt.cla() 
  
        doc.generate_pdf(clean_tex=True)