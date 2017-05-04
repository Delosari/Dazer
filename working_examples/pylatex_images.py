# import matplotlib.pyplot as plt
# import matplotlib 
# 
# N = 50
# x = np.random.rand(N)
# y = np.random.rand(N)
# colors = np.random.rand(N)
# area = np.pi * (15 * np.random.rand(N))**2  # 0 to 15 point radiuses
# 
# plt.scatter(x, y, s=area, c=colors, alpha=0.5)
# plt.show()

import numpy as np
import matplotlib.pyplot as plt
from pylatex import Document, Section, Figure, NoEscape, NewPage, Package

class pylatex_plotter():

    def __init__(self):
        
        self.pdf_geometry_options = {'right'    : '1cm',
                                     'left'     : '1cm',
                                     'top'      : '1cm',
                                     'bottom'   : '2cm'}

    def create_pdfDoc(self, fname, pdf_type = 'graphs', geometry_options = None):
        
        #Update the geometry if necessary (we coud define a dictionary distinction)
        if pdf_type == 'graphs':
            pdf_format = {'landscape':'true'}
            #Package('geometry', ['landscape'])
            
            self.pdf_geometry_options.update(pdf_format)
            
        if geometry_options is not None:
            self.pdf_geometry_options.update(geometry_options)

        self.pdfDoc = Document(fname, geometry_options=self.pdf_geometry_options)

    def fig_to_pdf(self, label=None, fig_loc='htbp', width=r'1\textwidth', add_page=False, *args, **kwargs):
        
        with self.pdfDoc.create(Figure(position=fig_loc)) as plot:
            plot.add_plot(width=NoEscape(width), *args, **kwargs)
        
            if label is not None:
                plot.add_caption(label)
            
        if add_page:
            self.pdfDoc.append(NewPage())
        
    def generate_pdf(self, clean_tex = True):
        
        self.pdfDoc.generate_pdf(clean_tex = clean_tex)
        
        return
                   
dz = pylatex_plotter()

dz.create_pdfDoc('/home/vital/Desktop/matplotlib_ex-dpi')

for i in range(3):

    N = 50
    x = np.random.rand(N)
    y = np.random.rand(N)
    colors = np.random.rand(N)
    area = np.pi * (15 * np.random.rand(N))**2  # 0 to 15 point radiuses
     
    plt.scatter(x, y, s=area, c=colors, alpha=0.5)

    dz.fig_to_pdf(add_page=True)

dz.generate_pdf(False)


