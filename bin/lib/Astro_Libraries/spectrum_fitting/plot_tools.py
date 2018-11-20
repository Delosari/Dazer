import numpy as np
import corner
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib import rcParams
from matplotlib.mlab import detrend_mean
from mpl_toolkits.mplot3d import Axes3D
from numpy import array, reshape, empty, ceil, percentile, median, nan, flatnonzero, core
from lib.Plotting_Libraries.dazer_plotter import Plot_Conf
from lib.Plotting_Libraries.sigfig import round_sig
from lib.CodeTools.File_Managing_Tools import Pdf_printer
from scipy import stats
from scipy.integrate    import simps

def label_formatting(line_label):

    label = line_label.replace('_', '\,\,')
    if label[-1] == 'A':
        label = label[0:-1] + '\AA'
    label = '$' + label + '$'

    return label

def _histplot_bins(column, bins=100):
    """Helper to get bins for histplot."""
    col_min = np.min(column)
    col_max = np.max(column)
    return range(col_min, col_max + 2, max((col_max - col_min) // bins, 1))

class Basic_plots(Plot_Conf):

    def __init__(self):

        # Class with plotting tools
        Plot_Conf.__init__(self)

    def prefit_input(self):

        size_dict = {'figure.figsize': (20, 14), 'axes.labelsize': 16, 'legend.fontsize': 18}
        self.FigConf(plotSize=size_dict)

        #Input continuum
        self.data_plot(self.inputWave, self.inputContinuum, 'Input object continuum')

        #Observed continuum
        self.data_plot(self.obj_data['wave_resam'], self.obj_data['flux_norm'], 'Observed spectrum', linestyle=':')

        # Nebular contribution removed if available
        self.data_plot(self.nebDefault['wave_neb'], self.nebDefault['synth_neb_flux'], 'Nebular contribution removed', linestyle = '--')

        # #In case of a synthetic observation:
        # if 'neb_SED' in self.obj_data:
        #     self.data_plot(self.input_wave, self.obj_data['neb_SED']['neb_int_norm'], 'Nebular continuum')
        #     title_label = 'Observed spectrum'
        # if 'stellar_flux' in self.obj_data:
        #     self.data_plot(self.obj_data['obs_wave'], self.obj_data['stellar_flux']/self.obj_data['normFlux_coeff'], 'Stellar continuum')
        #     self.data_plot(self.obj_data['obs_wave'], self.obj_data['stellar_flux_err']/ self.obj_data['normFlux_coeff'], 'Stellar continuum with uncertainty', linestyle=':')
        #     title_label = 'Synthetic spectrum'

        self.FigWording(xlabel='Wavelength $(\AA)$', ylabel = 'Observed flux', title = 'Observed spectrum and prefit input continuum')

        return

    def prefit_comparison(self, obj_ssp_fit_flux):

        size_dict = {'figure.figsize': (20, 14), 'axes.labelsize': 22, 'legend.fontsize': 22}
        self.FigConf(plotSize=size_dict)

        self.data_plot(self.inputWave, self.obj_data['flux_norm'], 'Object normed flux')
        self.data_plot(self.inputWave, self.nebDefault['synth_neb_flux'], 'Nebular continuum')
        self.data_plot(self.inputWave, obj_ssp_fit_flux + self.nebDefault['synth_neb_flux'], 'Prefit continuum output', linestyle='-')
        self.Axis.set_yscale('log')
        self.FigWording(xlabel='Wavelength $(\AA)$', ylabel='Normalised flux', title='')

        return

    def prefit_ssps(self, sspPrefitCoeffs):

        size_dict = {'figure.figsize': (20, 14), 'axes.labelsize': 22, 'legend.fontsize': 22}
        self.FigConf(plotSize=size_dict)

        #TODO This function should go to my collection
        ordinal_generator = lambda n: "%d%s" % (n, "tsnrhtdd"[(n / 10 % 10 != 1) * (n % 10 < 4) * n % 10::4])
        ordinal_bases = [ordinal_generator(n) for n in range(len(sspPrefitCoeffs))]

        counter = 0
        for i in range(len(self.onBasesFluxNorm)):
            if sspPrefitCoeffs[i] > self.lowlimit_sspContribution:
                counter += 1
                label_i = '{} base: flux coeff {}, norm coeff {:.2E}'.format(ordinal_bases[i], sspPrefitCoeffs[i], self.onBasesFluxNormCoeffs[i])
                label_i = '{} base'.format(ordinal_bases[i])
                self.data_plot(self.onBasesWave, self.onBasesFluxNorm[i], label_i)

        self.area_fill(self.ssp_lib['norm_interval'][0], self.ssp_lib['norm_interval'][1], 'Norm interval: {} - {}'.format(self.ssp_lib['norm_interval'][0], self.ssp_lib['norm_interval'][1]), alpha=0.5)

        title = 'SSP prefit contributing stellar populations {}/{}'.format((sspPrefitCoeffs > self.lowlimit_sspContribution).sum(), len(self.onBasesFluxNorm))
        self.FigWording(xlabel='Wavelength $(\AA)$', ylabel = 'Normalised flux', title = '')

        return

    def continuumFit(self, db_dict):

        outputSSPsCoefs, Av_star = np.median(db_dict['w_i']['trace'], axis=0), db_dict['Av_star']['median']
        outputSSPsCoefs_std, Av_star_std = np.std(db_dict['w_i']['trace'], axis=0), db_dict['Av_star']['standard deviation']

        stellarFit = outputSSPsCoefs.dot(self.onBasesFluxNorm) * np.power(10, -0.4 * Av_star * self.Xx_stellar)
        stellarPrefit = self.sspPrefitCoeffs.dot(self.onBasesFluxNorm) * np.power(10, -0.4 * self.stellarAv_prior[0] * self.Xx_stellar)
        nebularFit = self.nebDefault['synth_neb_flux']

        stellarFit_upper = (outputSSPsCoefs + outputSSPsCoefs_std).dot(self.onBasesFluxNorm) * np.power(10, -0.4 * (Av_star - Av_star_std) * self.Xx_stellar) + nebularFit
        stellarFit_lower = (outputSSPsCoefs - outputSSPsCoefs_std).dot(self.onBasesFluxNorm) * np.power(10, -0.4 * (Av_star + Av_star_std) * self.Xx_stellar) + nebularFit

        size_dict = {'figure.figsize': (20, 14), 'axes.labelsize': 16, 'legend.fontsize': 18}
        self.FigConf(plotSize=size_dict)

        self.data_plot(self.inputWave, self.inputContinuum, 'Object normed flux')
        self.data_plot(self.inputWave, nebularFit, 'Nebular continuum')
        self.data_plot(self.inputWave, stellarFit + nebularFit,'Continuum fit', color='tab:green')
        self.Axis.fill_between(self.inputWave, stellarFit_lower, stellarFit_upper, color='tab:green', alpha=0.5)
        self.data_plot(self.inputWave, stellarPrefit + nebularFit, 'prefit fit', color='tab:red', linestyle=':')

        self.Axis.set_yscale('log')
        self.FigWording(xlabel='Wavelength $(\AA)$', ylabel='Observed flux', title='SSPs continuum prefit')

        return

    def masked_observation(self):

        size_dict = {'figure.figsize': (20, 14), 'axes.labelsize': 16, 'legend.fontsize': 16}
        self.FigConf(plotSize=size_dict)

        nLineMasks = self.boolean_matrix.shape[0]
        inputContinuum = self.obj_data['flux_norm'] - self.nebDefault['synth_neb_flux']
        self.data_plot(self.inputWave, inputContinuum, 'Unmasked input continuum')
        self.data_plot(self.inputWave, inputContinuum*np.invert(self.object_mask), 'Object mask', linestyle=':')

        for i in range(nLineMasks):
            self.data_plot(self.inputWave, inputContinuum * self.boolean_matrix[i,:], 'Line mask ' + self.obj_data['lineLabels'][i], linestyle='--')

        self.Axis.set_xlabel('Wavelength $(\AA)$')
        self.Axis.set_ylabel('Observed flux')
        self.Axis.set_title('Spectrum masks')
        self.Axis.set_xscale('log')
        self.Axis.legend(bbox_to_anchor = (0.95, 1), loc = 2, borderaxespad = 0.)

        return

    def resampled_observation(self):

        size_dict = {'figure.figsize': (20, 14), 'axes.labelsize': 16, 'legend.fontsize': 18}
        self.FigConf(plotSize=size_dict)

        self.data_plot(self.obj_data['obs_wavelength'], self.obj_data['obs_flux'], 'Observed spectrum')
        self.data_plot(self.obj_data['wave_resam'], self.obj_data['flux_resam'], 'Resampled spectrum', linestyle='--')
        self.data_plot(self.obj_data['wave_resam'], self.obj_data['flux_norm'] * self.obj_data['normFlux_coeff'], r'Normalized spectrum $\cdot$ {:.2E}'.format(self.obj_data['normFlux_coeff']), linestyle=':')
        self.area_fill(self.obj_data['norm_interval'][0], self.obj_data['norm_interval'][1], 'Norm interval: {} - {}'.format(self.obj_data['norm_interval'][0], self.obj_data['norm_interval'][1]), alpha=0.5)

        self.FigWording(xlabel='Wavelength $(\AA)$', ylabel = 'Observed flux', title = 'Resampled observation')

        return

    def linesGrid(self, linesDf, wave, flux, plotAddress):
        # Get number of lines to generate the figure
        lineLabels = linesDf.index.values
        nLabels = lineLabels.size
        if nLabels <= 56:
            nRows, nColumns = 7, 8
        else:
            nRows, nColumns = 8, 8


        # Generate figure
        size_dict = {'figure.figsize': (30, 20), 'axes.titlesize': 12, 'axes.labelsize': 10, 'legend.fontsize': 10}
        self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=nColumns, n_rows=nRows)

        # Get line regions
        wavesMatrix = linesDf.loc[:,'w1':'w6'].values
        wavesIdcs = np.searchsorted(wave, wavesMatrix)
        idcsLines = (wave[wavesIdcs[:, 2]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 3]])
        idcsContinua = ((wave[wavesIdcs[:, 0]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 1]])) | ((wave[wavesIdcs[:, 4]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 5]]))
        idcsRedContinuum = ((wave[wavesIdcs[:, 0]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 1]]))
        idcsBlueContinuum = ((wave[wavesIdcs[:, 4]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 5]]))

        # Loop through the line and plot them
        for i in range(nLabels):

            lineLabel = lineLabels[i]
            lineType = linesDf.iloc[i].region_label

            # Get line_i wave and continua
            lineWave, lineFlux = wave[idcsLines[:, i]], flux[idcsLines[:, i]]
            continuaWave, continuaFlux = wave[idcsContinua[:, i]], flux[idcsContinua[:, i]]
            continuaRedWave, continuaRedFlux = wave[idcsRedContinuum[:, i]], flux[idcsRedContinuum[:, i]]
            continuaBlueWave, continuaBlueFlux = wave[idcsBlueContinuum[:, i]], flux[idcsBlueContinuum[:, i]]
            lineRes = lineWave[1] - lineWave[0] # TODO need to understand this better

            # Compute linear line continuum and get the standard deviation on the continuum
            slope, intercept, r_value, p_value, std_err = stats.linregress(continuaWave, continuaFlux)
            linearLineContinua = lineWave * slope + intercept

            # lineFlux_i = linesDf.loc[lineLabel, 'line_Flux'].nominal_value
            # lineFlux_iSimps = simps(lineFlux, lineWave) - simps(linearLineContinua, lineWave)
            # lineFlux_iSum = (lineFlux.sum() - linearLineContinua.sum()) * lineRes

            #Excluding Hbeta
            recombCheck = True if (('H1' in lineLabel) or ('He1' in lineLabel) or ('He2' in lineLabel)) and (lineLabel != 'H1_4861A') and ('_w' not in lineLabel) else False
            blendedCheck = True if (lineType is not 'continuum_mask') and (lineType != 'None') else False
            # print i, lineLabels[i], recombCheck, blendedCheck #linesDf.iloc[i].obs_flux, simps(lineFlux, lineWave), lineFlux.sum()

            #Plot the data
            self.Axis[i].plot(continuaRedWave, continuaRedFlux, color='tab:orange')
            self.Axis[i].plot(continuaBlueWave, continuaBlueFlux, color='tab:orange')
            self.Axis[i].plot(lineWave, lineFlux, color='tab:blue')

            # Add flux for isolated recombination lines
            if recombCheck and blendedCheck is False:

                # Assign new values
                fluxContinuum = linearLineContinua.sum() * lineRes
                linesDf.loc[lineLabel, 'obs_flux'] = linesDf.iloc[i].obs_flux + fluxContinuum

                self.Axis[i].plot(lineWave, linearLineContinua, color='tab:green')

            # Add flux for blended recombination lines
            if recombCheck and blendedCheck:
                muLine, sigmaLine = linesDf.iloc[i].mu, linesDf.iloc[i].sigma
                lineHalfWidth = 3 * sigmaLine
                w3, w4 = muLine - lineHalfWidth, muLine + lineHalfWidth
                idx3, idx4 = np.searchsorted(wave, [w3, w4])
                idcsLines_trim = (wave[idx3] <= wave) & (wave <= wave[idx4])
                lineWaveTrim, linearLineContinuaTrim = wave[idcsLines_trim], wave[idcsLines_trim] * slope + intercept
                fluxContinuum = linearLineContinuaTrim.sum() * lineRes

                # Assign new values
                linesDf.loc[lineLabel, 'obs_flux'] = linesDf.iloc[i].obs_flux + fluxContinuum
                linesDf.loc[lineLabel, 'w3'], linesDf.loc[lineLabel, 'w4'] = w3, w4

                # Generate the plot
                self.Axis[i].plot(lineWaveTrim, linearLineContinuaTrim, color='tab:purple')

            # Adjust the flux in N2_6548
            if lineLabel == 'N2_6548A':
                if linesDf.loc[lineLabel, 'region_label'] != 'None':
                    linesDf.loc[lineLabel, 'obs_fluxErr'] = linesDf.loc['N2_6584A', 'obs_fluxErr']
                    if (linesDf.loc['N2_6548A', 'obs_fluxErr'] / linesDf.loc['N2_6548A', 'obs_flux']) > 0.1:
                        linesDf.loc['N2_6548A', 'obs_fluxErr'] = linesDf.loc['N2_6548A', 'obs_flux'] * 0.1
                    if (linesDf.loc['N2_6584A', 'obs_fluxErr'] / linesDf.loc['N2_6584A', 'obs_flux']) > 0.1:
                        linesDf.loc['N2_6584A', 'obs_fluxErr'] = linesDf.loc['N2_6584A', 'obs_flux'] * 0.1

            # Format the plot
            self.Axis[i].get_yaxis().set_visible(False)
            self.Axis[i].set_yticks([])
            self.Axis[i].get_xaxis().set_visible(False)
            self.Axis[i].set_xticks([])
            self.Axis[i].set_yscale('log')

            #Wording plot
            self.Axis[i].set_title(lineLabels[i])

        # Plot the data
        plt.savefig(plotAddress, dpi=200, bbox_inches='tight')

        return

    def linesGrid_noContinuum(self, linesDf, wave, flux, plotAddress):
        # Get number of lines to generate the figure
        lineLabels = linesDf.index.values
        nLabels = lineLabels.size
        if nLabels <= 56:
            nRows, nColumns = 7, 8
        else:
            nRows, nColumns = 8, 8


        # Generate figure
        size_dict = {'figure.figsize': (30, 20), 'axes.titlesize': 12, 'axes.labelsize': 10, 'legend.fontsize': 10}
        self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=nColumns, n_rows=nRows)

        # Get line regions
        wavesMatrix = linesDf.loc[:,'w1':'w6'].values
        wavesIdcs = np.searchsorted(wave, wavesMatrix)
        idcsLines = (wave[wavesIdcs[:, 2]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 3]])
        idcsContinua = ((wave[wavesIdcs[:, 0]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 1]])) | ((wave[wavesIdcs[:, 4]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 5]]))
        idcsRedContinuum = ((wave[wavesIdcs[:, 0]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 1]]))
        idcsBlueContinuum = ((wave[wavesIdcs[:, 4]] <= wave[:, None]) & (wave[:, None] <= wave[wavesIdcs[:, 5]]))

        # Loop through the line and plot them
        for i in range(nLabels):

            lineLabel = lineLabels[i]
            lineType = linesDf.iloc[i].region_label

            # Get line_i wave and continua
            lineWave, lineFlux = wave[idcsLines[:, i]], flux[idcsLines[:, i]]
            continuaWave, continuaFlux = wave[idcsContinua[:, i]], flux[idcsContinua[:, i]]
            continuaRedWave, continuaRedFlux = wave[idcsRedContinuum[:, i]], flux[idcsRedContinuum[:, i]]
            continuaBlueWave, continuaBlueFlux = wave[idcsBlueContinuum[:, i]], flux[idcsBlueContinuum[:, i]]
            lineRes = lineWave[1] - lineWave[0] # TODO need to understand this better

            # Compute linear line continuum and get the standard deviation on the continuum
            slope, intercept, r_value, p_value, std_err = stats.linregress(continuaWave, continuaFlux)
            linearLineContinua = lineWave * slope + intercept

            # lineFlux_i = linesDf.loc[lineLabel, 'line_Flux'].nominal_value
            # lineFlux_iSimps = simps(lineFlux, lineWave) - simps(linearLineContinua, lineWave)
            # lineFlux_iSum = (lineFlux.sum() - linearLineContinua.sum()) * lineRes

            #Excluding Hbeta
            recombCheck = True if (('H1' in lineLabel) or ('He1' in lineLabel) or ('He2' in lineLabel)) and (lineLabel != 'H1_4861A') and ('_w' not in lineLabel) else False
            blendedCheck = True if (lineType is not 'continuum_mask') and (lineType != 'None') else False
            # print i, lineLabels[i], recombCheck, blendedCheck #linesDf.iloc[i].obs_flux, simps(lineFlux, lineWave), lineFlux.sum()

            #Plot the data
            self.Axis[i].plot(continuaRedWave, continuaRedFlux, color='tab:orange')
            self.Axis[i].plot(continuaBlueWave, continuaBlueFlux, color='tab:orange')
            self.Axis[i].plot(lineWave, lineFlux, color='tab:blue')

            # Add flux for isolated recombination lines
            if recombCheck and blendedCheck is False:

                # Assign new values
                fluxContinuum = 0.0
                linesDf.loc[lineLabel, 'obs_flux'] = linesDf.iloc[i].obs_flux + fluxContinuum

                self.Axis[i].plot(lineWave, linearLineContinua, color='tab:green')

            # Add flux for blended recombination lines
            if recombCheck and blendedCheck:
                muLine, sigmaLine = linesDf.iloc[i].mu, linesDf.iloc[i].sigma
                lineHalfWidth = 3 * sigmaLine
                w3, w4 = muLine - lineHalfWidth, muLine + lineHalfWidth
                idx3, idx4 = np.searchsorted(wave, [w3, w4])
                idcsLines_trim = (wave[idx3] <= wave) & (wave <= wave[idx4])
                lineWaveTrim, linearLineContinuaTrim = wave[idcsLines_trim], wave[idcsLines_trim] * slope + intercept
                fluxContinuum = 0.0

                # Assign new values
                linesDf.loc[lineLabel, 'obs_flux'] = linesDf.iloc[i].obs_flux + fluxContinuum
                linesDf.loc[lineLabel, 'w3'], linesDf.loc[lineLabel, 'w4'] = w3, w4

                # Generate the plot
                self.Axis[i].plot(lineWaveTrim, linearLineContinuaTrim, color='tab:purple')

            # Adjust the flux in N2_6548
            if lineLabel == 'N2_6548A':
                if linesDf.loc[lineLabel, 'region_label'] != 'None':
                    linesDf.loc[lineLabel, 'obs_fluxErr'] = linesDf.loc['N2_6584A', 'obs_fluxErr']
                    if (linesDf.loc['N2_6548A', 'obs_fluxErr'] / linesDf.loc['N2_6548A', 'obs_flux']) > 0.1:
                        linesDf.loc['N2_6548A', 'obs_fluxErr'] = linesDf.loc['N2_6548A', 'obs_flux'] * 0.1
                    if (linesDf.loc['N2_6584A', 'obs_fluxErr'] / linesDf.loc['N2_6584A', 'obs_flux']) > 0.1:
                        linesDf.loc['N2_6584A', 'obs_fluxErr'] = linesDf.loc['N2_6584A', 'obs_flux'] * 0.1

            # Adjust the flux in N2_6548
            if lineLabel == 'O2_7319A':
                if linesDf.loc[lineLabel, 'region_label'] != 'None':
                    if (linesDf.loc['O2_7319A', 'obs_fluxErr'] / linesDf.loc['O2_7319A', 'obs_flux']) > 0.1:
                        linesDf.loc['O2_7319A', 'obs_fluxErr'] = linesDf.loc['O2_7319A', 'obs_flux'] * 0.1
                    if (linesDf.loc['O2_7330A', 'obs_fluxErr'] / linesDf.loc['O2_7330A', 'obs_flux']) > 0.1:
                        linesDf.loc['O2_7330A', 'obs_fluxErr'] = linesDf.loc['O2_7330A', 'obs_flux'] * 0.1

            # Format the plot
            self.Axis[i].get_yaxis().set_visible(False)
            self.Axis[i].set_yticks([])
            self.Axis[i].get_xaxis().set_visible(False)
            self.Axis[i].set_xticks([])
            self.Axis[i].set_yscale('log')

            #Wording plot
            self.Axis[i].set_title(lineLabels[i])

        # Plot the data
        plt.savefig(plotAddress, dpi=200, bbox_inches='tight')

        return

    def emissivitySurfaceFit_2D(self, line_label, emisCoeffs, emisGrid, funcEmis, te_ne_grid):

        # Plot format
        size_dict = {'figure.figsize': (20, 14),  'axes.titlesize': 16, 'axes.labelsize': 16, 'legend.fontsize': 18}
        rcParams.update(size_dict)

        # Generate figure
        fig, ax = plt.subplots(1, 1)

        # Generate fitted surface points
        matrix_edge = int(np.sqrt(te_ne_grid[0].shape[0]))
        surface_points = funcEmis(te_ne_grid, *emisCoeffs)

        # Plot plane
        plt.imshow(surface_points.reshape((matrix_edge, matrix_edge)), aspect=0.03,
                   extent=(te_ne_grid[1].min(), te_ne_grid[1].max(), te_ne_grid[0].min(), te_ne_grid[0].max()))

        # Compare pyneb values with values from fitting
        percentage_difference = (1 - surface_points/emisGrid.flatten()) * 100

        # Points with error below 1.0 are transparent:
        idx_interest = percentage_difference < 1.0
        x_values, y_values = te_ne_grid[1][idx_interest], te_ne_grid[0][idx_interest]
        ax.scatter(x_values, y_values, c="None", edgecolors='black', linewidths=0.35, label='Error below 1%')

        if idx_interest.sum() < emisGrid.size:

            # Plot grid points
            plt.scatter(te_ne_grid[1][~idx_interest], te_ne_grid[0][~idx_interest], c=percentage_difference[~idx_interest],
                        edgecolors='black', linewidths=0.1, cmap=cm.OrRd , label='Error above 1%')

            # Color bar
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('% difference', rotation=270, fontsize=15)

        #Trim the axis
        ax.set_xlim(te_ne_grid[1].min(), te_ne_grid[1].max())
        ax.set_ylim(te_ne_grid[0].min(), te_ne_grid[0].max())

        # Add labels
        ax.update({'xlabel': 'Density ($cm^{-3}$)', 'ylabel': 'Temperature $(K)$', 'title': line_label})

        return

    def emissivitySurfaceFit_3D(self, line_label, emisCoeffs, emisGrid, funcEmis, te_ne_grid):

        # Plot format
        size_dict = {'figure.figsize': (20, 14), 'axes.titlesize': 16, 'axes.labelsize': 16, 'legend.fontsize': 18}
        rcParams.update(size_dict)

        # Plot the grid points
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Generate fitted surface points
        matrix_edge = int(np.sqrt(te_ne_grid[0].shape[0]))

        # Plotting pyneb emissivities
        x_values, y_values = te_ne_grid[0].reshape((matrix_edge, matrix_edge)), te_ne_grid[1].reshape((matrix_edge, matrix_edge))
        ax.plot_surface(x_values, y_values, emisGrid.reshape((matrix_edge, matrix_edge)), color='g', alpha=0.5)

        # Plotting emissivity parametrization
        fit_points = funcEmis(te_ne_grid, *emisCoeffs)
        ax.scatter(te_ne_grid[0], te_ne_grid[1], fit_points, color='r', alpha=0.5)

        # Add labels
        ax.update({'ylabel': 'Density ($cm^{-3}$)', 'xlabel': 'Temperature $(K)$', 'title': line_label})

        return

    def traces_plot(self, traces_list, stats_dic):

        #Remove operations from the parameters list
        traces = traces_list[[i for i, v in enumerate(traces_list) if ('_Op' not in v) and ('_log__' not in v) and ('w_i' not in v)]]

        # Number of traces to plot
        n_traces = len(traces)

        # Declare figure format
        size_dict = {'figure.figsize': (14, 20), 'axes.titlesize':26, 'axes.labelsize': 24, 'legend.fontsize': 18}
        self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=1, n_rows=n_traces)

        # Generate the color map
        self.gen_colorList(0, n_traces)

        # Plot individual traces
        for i in range(n_traces):

            # Current trace
            trace_code = traces[i]
            trace_array = stats_dic[trace_code]['trace']

            # Label for the plot
            mean_value = stats_dic[trace_code]['mean']
            std_dev = stats_dic[trace_code]['standard deviation']
            if mean_value > 0.001:
                label = r'{} = ${}$ $\pm${}'.format(self.labels_latex_dic[trace_code], round_sig(mean_value, 4),
                                                    round_sig(std_dev, 4))
            else:
                label = r'{} = ${:.3e}$ $\pm$ {:.3e}'.format(self.labels_latex_dic[trace_code], mean_value, std_dev)

            # Plot the data
            self.Axis[i].plot(trace_array, label=label, color=self.get_color(i))
            self.Axis[i].axhline(y=mean_value, color=self.get_color(i), linestyle='--')
            self.Axis[i].set_ylabel(self.labels_latex_dic[trace_code])

            if i < n_traces - 1:
                self.Axis[i].set_xticklabels([])

            # Add legend
            self.legend_conf(self.Axis[i], loc=2)

        return

    def tracesPosteriorPlot(self, traces_list, stats_dic):

        #Remove operations from the parameters list # TODO addapt this line to discremenate better
        traces = traces_list[[i for i, v in enumerate(traces_list) if ('_Op' not in v) and ('_log__' not in v) and ('w_i' not in v)]]

        # Number of traces to plot
        n_traces = len(traces)

        # Declare figure format
        size_dict = {'axes.titlesize':20, 'axes.labelsize': 20, 'legend.fontsize': 10}
        rcParams.update(size_dict)#self.FigConf(plotSize=size_dict, Figtype='tracePosterior', n_columns=2, n_rows=n_traces)
        fig = plt.figure(figsize=(8, n_traces))

        # # Generate the color map
        self.gen_colorList(0, n_traces)
        gs = gridspec.GridSpec(n_traces*2, 4)

        gs.update(wspace=0.2, hspace=1.8)
        for i in range(n_traces):

            #Creat figure axis
            axTrace = fig.add_subplot(gs[2*i:2*(1+i),:3])
            axPoterior = fig.add_subplot(gs[2*i:2*(1+i), 3])

            # Current trace
            trace_code = traces[i]
            trace_array = stats_dic[trace_code]['trace']

            # Label for the plot
            mean_value = stats_dic[trace_code]['mean']
            std_dev = stats_dic[trace_code]['standard deviation']
            if mean_value > 0.001:
                label = r'{} = ${}$ $\pm${}'.format(self.labels_latex_dic[trace_code], round_sig(mean_value, 4),
                                                    round_sig(std_dev, 4))
            else:
                label = r'{} = ${:.3e}$ $\pm$ {:.3e}'.format(self.labels_latex_dic[trace_code], mean_value, std_dev)

            # Plot the traces
            axTrace.plot(trace_array, label=label, color=self.get_color(i))
            axTrace.axhline(y=mean_value, color=self.get_color(i), linestyle='--')
            axTrace.set_ylabel(self.labels_latex_dic[trace_code])

            #Plot the histograms
            axPoterior.hist(trace_array, bins=50, histtype='step', color=self.get_color(i), align = 'left')

            #Add true value if available
            if trace_code + '_true' in self.obj_data:
                value_true = self.obj_data[trace_code + '_true']
                axPoterior.axvline(x=value_true, color=self.get_color(i), linestyle='solid')

            # Add legend
            axTrace.legend(loc=7)

            #Remove ticks and labels
            if i < n_traces - 1:
                axTrace.get_xaxis().set_visible(False)
                axTrace.set_xticks([])
            axPoterior.yaxis.set_major_formatter(plt.NullFormatter())
            axPoterior.set_yticks([])

        return

    def posteriors_plot(self, traces_list, stats_dic):

        #Remove operations from the parameters list
        traces = traces_list[[i for i, v in enumerate(traces_list) if ('_Op' not in v) and ('_log__' not in v) and ('w_i' not in v)]]

        # Number of traces to plot
        n_traces = len(traces)

        # Declare figure format
        size_dict = {'figure.figsize': (14, 20), 'axes.titlesize':22, 'axes.labelsize': 22, 'legend.fontsize': 14}
        self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=1, n_rows=n_traces)

        # Generate the color map
        self.gen_colorList(0, n_traces)

        # Plot individual traces
        for i in range(len(traces)):

            # Current trace
            trace_code = traces[i]
            mean_value = stats_dic[trace_code]['mean']
            trace_array = stats_dic[trace_code]['trace']

            # Plot HDP limits
            HDP_coords = stats_dic[trace_code]['95% HPD interval']
            for HDP in HDP_coords:

                if mean_value > 0.001:
                    label_limits = 'HPD interval: {} - {}'.format(round_sig(HDP_coords[0], 4),
                                                                  round_sig(HDP_coords[1], 4))
                    label_mean = 'Mean value: {}'.format(round_sig(mean_value, 4))
                else:
                    label_limits = 'HPD interval: {:.3e} - {:.3e}'.format(HDP_coords[0], HDP_coords[1])
                    label_mean = 'Mean value: {:.3e}'.format(mean_value)

                self.Axis[i].axvline(x=HDP, label=label_limits, color='grey', linestyle='dashed')

            self.Axis[i].axvline(x=mean_value, label=label_mean, color='grey', linestyle='solid')
            self.Axis[i].hist(trace_array, histtype='stepfilled', bins=35, alpha=.7, color=self.get_color(i), normed=False)

            #Add true value if available
            if 'true_value' in stats_dic[trace_code]:
                value_true = stats_dic[trace_code]['true_value']
                if value_true is not None:
                    label_true = 'True value {:.3e}'.format(value_true)
                    self.Axis[i].axvline(x=value_true, label=label_true, color='black', linestyle='solid')

            # Figure wording
            self.Axis[i].set_ylabel(self.labels_latex_dic[trace_code])
            self.legend_conf(self.Axis[i], loc=2)

    def fluxes_distribution(self, lines_list, function_key, db_dict, true_values = None):

        # Number of traces to plot
        n_lines = len(lines_list)

        # Declare figure format
        size_dict = {'figure.figsize': (14, 20), 'axes.titlesize':11, 'axes.labelsize': 8, 'legend.fontsize': 11}
        self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=1, n_rows=n_lines)

        # Generate the color map
        self.gen_colorList(0, n_lines)

        # Flux statistics
        traces_array    = db_dict[function_key]['trace']
        median_values   = median(db_dict[function_key]['trace'], axis=0)
        p16th_fluxes    = percentile(db_dict[function_key]['trace'], 16, axis=0)
        p84th_fluxes    = percentile(db_dict[function_key]['trace'], 84, axis=0)

        # Plot individual traces
        for i in range(n_lines):

            #Current line
            label       = label_formatting(lines_list[i])
            trace       = traces_array[:,i]
            median_flux = median_values[i]

            # Plot HDP limits
            for HDP in [p16th_fluxes[i], p84th_fluxes[i]]:

                if median_flux > 0.001:
                    label_limits = 'HPD interval: {} - {}'.format(round_sig(p16th_fluxes[i], 4), round_sig(p84th_fluxes[i], 4))
                    label_mean = 'Mean value: {}'.format(round_sig(median_flux, 4))
                else:
                    label_limits = 'HPD interval: {:.3e} - {:.3e}'.format(p16th_fluxes[i], p84th_fluxes[i])
                    label_mean = 'Mean value: {:.3e}'.format(median_flux)

                self.Axis[i].axvline(x=HDP, label=label_limits, color='grey', linestyle='dashed')

            self.Axis[i].axvline(x=median_flux, label=label_mean, color='grey', linestyle='solid')
            self.Axis[i].hist(trace, histtype='stepfilled', bins=35, alpha=.7, color=self.get_color(i), normed=False)

            if  true_values is not None:
                true_value = true_values[i]
                label_true = 'True value: {}'.format(round_sig(true_value, 3))
                self.Axis[i].axvline(x=true_value, label=label_true, color='black', linestyle='solid')

            # Plot wording
            self.Axis[i].set_ylabel(label)
            self.legend_conf(self.Axis[i], loc=2)

        return

    def acorr_plot(self, traces_list, stats_dic, n_columns=4, n_rows=2):

        #Remove operations from the parameters list
        traces = traces_list[[i for i, v in enumerate(traces_list) if ('_Op' not in v) and ('_log__' not in v) and ('w_i' not in v)]]

        # Number of traces to plot
        n_traces = len(traces)

        # Declare figure format
        size_dict = {'figure.figsize': (14, 14), 'axes.titlesize': 20, 'legend.fontsize': 10}
        self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=n_columns, n_rows=n_rows)

        # Generate the color map
        self.gen_colorList(0, n_traces)

        # Plot individual traces
        for i in range(n_traces):

            # Current trace
            trace_code = traces[i]

            label = self.labels_latex_dic[trace_code]

            trace_array = stats_dic[trace_code]['trace']

            if trace_code != 'ChiSq':
                maxlags = min(len(trace_array) - 1, 100)
                self.Axis[i].acorr(x=trace_array, color=self.get_color(i), detrend=detrend_mean, maxlags=maxlags)

            else:
                # Apano momentaneo
                chisq_adapted = reshape(trace_array, len(trace_array))
                maxlags = min(len(chisq_adapted) - 1, 100)
                self.Axis[i].acorr(x=chisq_adapted, color=self.get_color(i), detrend=detrend_mean, maxlags=maxlags)

            self.Axis[i].set_xlim(0, maxlags)
            self.Axis[i].set_title(label)

        return

    def corner_plot(self, traces_list, stats_dic, true_values=None):

        #Remove operations from the parameters list
        traces = traces_list[[i for i, v in enumerate(traces_list) if ('_Op' not in v) and ('_log__' not in v) and ('w_i' not in v)]]

        # Number of traces to plot
        n_traces = len(traces)

        # Set figure conf
        sizing_dict = {}
        sizing_dict['figure.figsize'] = (14, 14)
        sizing_dict['legend.fontsize'] = 30
        sizing_dict['axes.labelsize'] = 30
        sizing_dict['axes.titlesize'] = 15
        sizing_dict['xtick.labelsize'] = 12
        sizing_dict['ytick.labelsize'] = 12

        rcParams.update(sizing_dict)

        # Reshape plot data
        list_arrays = []
        labels_list = []
        for trace_code in traces:
            trace_array = stats_dic[trace_code]['trace']
            list_arrays.append(trace_array)
            labels_list.append(self.labels_latex_dic[trace_code])
        traces_array = array(list_arrays).T

        # true_values = empty(len(traces))
        # for i in range(len(traces)):
        #     if 'true_value' in stats_dic[traces[i]]:
        #         true_values[i] = stats_dic[traces[i]]['true_value']
        #     else:
        #         true_values[i] = 'None'

        #Generate the plot
        self.Fig = corner.corner(traces_array[:, :], fontsize=30, labels=labels_list, quantiles=[0.16, 0.5, 0.84],
                                 show_titles=True, title_args={"fontsize": 200}, truths=true_values, truth_color='#ae3135',
                                 title_fmt='0.3f')



        # # Make
        # if plot_true_values:
        #
        #     true_values = empty(len(traces))
        #     for i in range(len(traces)):
        #         if 'true_value' in stats_dic[traces[i]]:
        #             true_values[i] = stats_dic[traces[i]]['true_value']
        #         else:
        #             true_values[i] = None
        #
        #     self.Fig = corner.corner(traces_array[:, :], fontsize=30, labels=labels_list, quantiles=[0.16, 0.5, 0.84],
        #                              show_titles=True, title_args={"fontsize": 200}, truths=true_values,
        #                              title_fmt='0.3f')
        # else:
        #     self.Fig = corner.corner(traces_array[:, :], fontsize=30, labels=labels_list, quantiles=[0.16, 0.5, 0.84],
        #                              show_titles=True, title_args={"fontsize": 200}, title_fmt='0.3f')

        return

class Basic_tables(Pdf_printer):

    def __init__(self):

        # Class with pdf tools
        Pdf_printer.__init__(self)

    def table_mean_outputs(self, table_address, parameters_list, db_dict, true_values = None):

        # Table headers
        headers = ['Parameter', 'True value', 'Mean', 'Standard deviation', 'Number of points', 'Median', r'$16^{th}$ $percentil$', r'$84^{th}$ $percentil$', r'$Difference$']

        # Generate pdf
        self.create_pdfDoc(table_address, pdf_type='table')
        self.pdf_insert_table(headers)

        #Loop around the parameters
        for param in parameters_list:

            if ('_Op' not in param) and param not in ['w_i']:

                # Label for the plot
                label       = self.labels_latex_dic[param]
                mean_value  = db_dict[param]['mean']
                std         = db_dict[param]['standard deviation']
                n_traces    = db_dict[param]['n']
                median      = db_dict[param]['median']
                p_16th      = db_dict[param]['16th_p']
                p_84th      = db_dict[param]['84th_p']
                true_value  = 'None'
                perDif      = 'None'

                if param + '_true' in true_values:
                    true_value = true_values[param + '_true']
                    perDif = 1 - (true_value/median * 100)

                self.addTableRow([label, true_value, mean_value, std, n_traces, median, p_16th, p_84th, perDif], last_row=False if parameters_list[-1] != param else True)

        self.generate_pdf(clean_tex=True)

        return

    def table_line_fluxes(self, table_address, lines_list, function_key, db_dict, true_values = None):

        # Generate pdf
        self.create_pdfDoc(table_address, pdf_type='table')

        # Table headers
        headers = ['Line Label', 'Mean', 'Standard deviation', 'Median', r'$16^{th}$ $percentil$', r'$84^{th}$ $percentil$']

        if true_values is not None:
            headers.insert(1, 'True value')
            headers.insert(2, 'Diff percentage')

        self.pdf_insert_table(headers)

        mean_line_values    = db_dict[function_key]['trace'].mean(axis=0)
        std_line_values     = db_dict[function_key]['trace'].std(axis=0)
        median_line_values  = median(db_dict[function_key]['trace'], axis=0)
        p16th_line_values   = percentile(db_dict[function_key]['trace'],16, axis=0)
        p84th_line_values   = percentile(db_dict[function_key]['trace'],84, axis=0)
        diff_Percentage     = 1 - (median_line_values/true_values * 100)

        for i in range(len(lines_list)):

            #Operations to format the label name "H1_4861A"
            label = label_formatting(lines_list[i])

            row_i = [label, mean_line_values[i], std_line_values[i], median_line_values[i], p16th_line_values[i], p84th_line_values[i]]

            if true_values is not None:
                row_i.insert(1, true_values[i])
                row_i.insert(2, diff_Percentage[i])

            self.addTableRow(row_i, last_row=False if lines_list[-1] != lines_list[i] else True)

        self.generate_pdf(clean_tex=True)

class MCMC_printer(Basic_plots, Basic_tables):

    def __init__(self):

        # Supporting classes
        Basic_plots.__init__(self)
        Basic_tables.__init__(self)

    def plot_emisFits(self, linelabels, emisCoeffs_dict, emisGrid_array, output_folder):

        # Temperature and density meshgrids
        X, Y = np.meshgrid(self.tem_grid_range, self.den_grid_range)
        XX, YY = X.flatten(), Y.flatten()
        te_ne_grid = (XX, YY)

        for i in range(linelabels.size):

            lineLabel = linelabels[i]
            print '--Fitting surface', lineLabel

            # 2D Comparison between PyNeb values and the fitted equation
            self.emissivitySurfaceFit_2D(lineLabel, emisCoeffs_dict[lineLabel], emisGrid_array[:,i], self.ionEmisEq[lineLabel], te_ne_grid)

            output_address = '{}{}_{}'.format(output_folder,'emissivityTeDe2D',lineLabel)
            self.savefig(output_address, resolution=200)
            plt.clf()

            # 3D Comparison between PyNeb values and the fitted equation
            self.emissivitySurfaceFit_3D(lineLabel, emisCoeffs_dict[lineLabel], emisGrid_array[:,i], self.ionEmisEq[lineLabel], te_ne_grid)

            output_address = '{}{}_{}'.format(output_folder,'emissivityTeDe3D',lineLabel)
            self.savefig(output_address, resolution=200)
            plt.clf()

        return

    def plot_emisRatioFits(self, diagnoslabels, emisCoeffs_dict, emisGrid_array, output_folder):

        # Temperature and density meshgrids
        X, Y = np.meshgrid(self.tem_grid_range, self.den_grid_range)
        XX, YY = X.flatten(), Y.flatten()
        te_ne_grid = (XX, YY)

        for i in range(diagnoslabels.size):

            lineLabel = diagnoslabels[i]
            print '--Fitting surface', lineLabel

            # 2D Comparison between PyNeb values and the fitted equation
            self.emissivitySurfaceFit_2D(lineLabel, emisCoeffs_dict[lineLabel], emisGrid_array[:,i], self.EmisRatioEq_fit[lineLabel], te_ne_grid)

            output_address = '{}{}_{}'.format(output_folder,'emissivityTeDe2D',lineLabel)
            self.savefig(output_address, resolution=200)
            plt.clf()

            # 3D Comparison between PyNeb values and the fitted equation
            self.emissivitySurfaceFit_3D(lineLabel, emisCoeffs_dict[lineLabel], emisGrid_array[:,i], self.EmisRatioEq_fit[lineLabel], te_ne_grid)

            output_address = '{}{}_{}'.format(output_folder,'emissivityTeDe3D',lineLabel)
            self.savefig(output_address, resolution=200)
            plt.clf()

        return

    def plotInputSSPsynthesis(self):

        # Plot input data
        self.prefit_input()
        self.savefig(self.input_folder + self.objName + '_prefit_input', resolution=200)

        # Plot resampling
        self.resampled_observation()
        self.savefig(self.input_folder + self.objName + '_resampled_spectra', resolution=200)

        # Spectrum masking
        self.masked_observation()
        self.savefig(self.input_folder + self.objName + '_spectrum_mask', resolution=200)

        return

    def plotOutputSSPsynthesis(self, pymc2_dbPrefit, pymc2_db_dictPrefit, obj_ssp_fit_flux = None, sspPrefitCoeffs=None):

        # Continua components from ssp prefit
        if obj_ssp_fit_flux is not None:
            self.prefit_comparison(obj_ssp_fit_flux)
            self.savefig(self.input_folder + self.objName + '_prefit_Output', resolution=200)

            # Prefit SSPs norm plot
            if sspPrefitCoeffs is not None:
                self.prefit_ssps(sspPrefitCoeffs)
                self.savefig(self.input_folder + self.objName + '_prefit_NormSsps', resolution=200)

        # SSP prefit traces # TODO increase flexibility for a big number of parameter
        traces_names = np.array(['Av_star', 'sigma_star'])

        # Stellar prefit Traces
        self.traces_plot(traces_names, pymc2_db_dictPrefit)
        self.savefig(self.input_folder + self.objName + '_PrefitTracesPlot', resolution=200)

        # Stellar prefit Posteriors
        # self.posteriors_plot(traces_names, pymc2_db_dictPrefit)
        # self.savefig(self.input_folder + self.objName + '_PrefitPosteriorPlot', resolution=200)

        # Stellar prefit Posteriors
        self.tracesPosteriorPlot(traces_names, pymc2_db_dictPrefit)
        self.savefig(self.input_folder + self.objName + '_ParamsTracesPosterios_StellarPrefit', resolution=200)

        return

    def plotOuputData(self, database_address, db, db_dict, params_list):

        if self.stellarCheck:
            self.continuumFit(db_dict)
            self.savefig(database_address + '_ContinuumFit', resolution=200)


        if self.emissionCheck:

            # Posteriors
            print '--Model parameters posterior diagram'
            self.tracesPosteriorPlot(params_list, db_dict)
            self.savefig(database_address + '_ParamsTracesPosterios', resolution=200)

            # # Table mean values
            print '--Model parameters table'
            self.table_mean_outputs(database_address + '_meanOutput', params_list, db_dict, true_values = self.obj_data)

            # # Line fluxes values
            print '--Line fluxes table' # TODO divide in two
            self.table_line_fluxes(database_address + '_LineFluxes', self.lineLabels, 'calcFluxes_Op', db_dict, true_values=self.obsLineFluxes)
            self.fluxes_distribution(self.lineLabels, 'calcFluxes_Op', db_dict, true_values=self.obsLineFluxes)
            self.savefig(database_address + '_LineFluxesPosteriors', resolution=200)

            # Traces
            print '--Traces diagram'
            self.traces_plot(params_list, db_dict)
            self.savefig(database_address + '_TracesPlot', resolution=200)

            # Posteriors
            print '--Model parameters posterior diagram'
            self.posteriors_plot(params_list, db_dict)
            self.savefig(database_address + '_PosteriorPlot', resolution=200)

            # Acorrelation
            print '--Acorrelation plot'
            self.acorr_plot(params_list, db_dict, n_columns=4, n_rows=int(ceil(len(params_list)/4.0)))
            self.savefig(database_address + '_Acorrelation', resolution=200)

            # Corner plot
            print '--Scatter plot matrix'
            self.corner_plot(params_list, db_dict)
            self.savefig(database_address + '_CornerPlot', resolution=100)

        return
