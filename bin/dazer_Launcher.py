#!/usr/bin/python

import  matplotlib.widgets  as widgets
from    bin.DZ_GUI              import Tk_GUI
from    bin.DZ_DataExplorer     import Plots_Manager

def Span_Manager(Wlow, Whig):
    
    #Check to make sure we are at the measuring lines screen
    if (pv.CurrentScreen == 'LINEMESURER'):    

        #Check we are not just clicking on the plot 
        if (Wlow != Whig):
            
            #Clear the plot
            pv.reset_fig()

            #Case selecting 1/3 region
            if len(pv.Selections) == 0:
                pv.Selections.append(Wlow)
                pv.Selections.append(Whig)
                
            #Case selecting 2/3 region   
            elif len(pv.Selections) == 2:
                pv.Selections.append(Wlow)
                pv.Selections.append(Whig)
                pv.Selections.sort()
                
            #Case selecting 3/3 region
            elif len(pv.Selections) == 4:
                pv.Selections.append(Wlow)
                pv.Selections.append(Whig)
                pv.Selections.sort()
                
            elif len(pv.Selections) == 6:
                 
                #Caso que se corrija la region de la linea
                if (Wlow > pv.Selections[1] and Whig < pv.Selections[4]):     
                    pv.Selections[2] = Wlow
                    pv.Selections[3] = Whig
                
                #Caso que se corrija el continuum izquierdo    
                elif (Wlow < pv.Selections[2] and Whig < pv.Selections[2]):   
                    pv.Selections[0] = Wlow
                    pv.Selections[1] = Whig
                                
                #Caso que se corrija el continuum derecho     
                elif (Wlow > pv.Selections[3] and Whig > pv.Selections[3]):   
                    pv.Selections[4] = Wlow
                    pv.Selections[5] = Whig
                
                #Case we want to select the complete region
                elif (Wlow < pv.Selections[0] and Whig > pv.Selections[5]):
                    
                    #Remove line from dataframe and save it
                    pv.remove_lines_df(pv.current_df, pv.Current_Label)
                    
                    #Save lines log df
                    pv.save_lineslog_dataframe(pv.current_df, pv.lineslog_df_address)
                    
                    #Clear the selections
                    del pv.Selections[:]
                    
            elif len(pv.Selections) > 6:
                print "WARNING: Selections size over limits 6"
            
            pv.ManagePlots()

        pv.FigCanvas.draw()

    return

def Key_Manager(event): 
    
    #Clear the screen
    pv.Axis.clear()
    
    #DO WE NEED THIS ONE HERE?
    pv.Axis.autoscale(True)
    
    #Move along data types
    if (event.key == "down") or (event.key == "s"):                                                        
        pv.Manage_DepthIndexes('NextObject')    
    elif (event.key == "up") or  (event.key == "w"):
        pv.Manage_DepthIndexes('PreviousObject')
        
    #Move along data equals
    elif (event.key == "right") or (event.key == "d"):     
        pv.Manage_SubSectionIndexes('NextObject')        
    elif (event.key == "left") or (event.key == "a"):
        pv.Manage_SubSectionIndexes('PreviousObject')
        
    #Quit the code
    elif (event.key == '0'):
        
        GUI.quit()
        print "Cerrando el programa"

        return
    
    #Unknown Key combination
    else:
        print '--', event.key, '-- has not been assigned a command'

    #Plot the corresponding images        
    pv.ManagePlots()

    #Redraw the image
    pv.FigCanvas.draw()

pv = Plots_Manager()

#---- Find files within root folder according to log files ".log"
pv.FilesList = pv.FindAndOrganize_dazer(pv.Pattern_PlotFiles, pv.Catalogue_Folder, unpack=True, CheckComputer = False)

#---- Generate Initial Frame
pv.FigConf('night')
pv.Axis.text(0.95, 0.05, 'Initiating visualizer', verticalalignment='bottom', horizontalalignment='right',transform = pv.Axis.transAxes, fontsize=15)

GUI = Tk_GUI(PlottingVector = pv)

pv.FigCanvas.show()
pv.FigCanvas.mpl_connect('key_press_event', Key_Manager)
Span = widgets.SpanSelector(pv.Axis, Span_Manager, 'horizontal', useblit=False, rectprops=dict(alpha=1, facecolor='Blue'))

GUI.mainloop()

print "Dazer Closed"

