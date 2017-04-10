'''
Designed to plot comparision figures.
'''

from . import style
import matplotlib.collections
import numpy
import logging

# from ForceMode import ForceMode
# from ForceData import ForceData
# from ForceAnalysis import ForceAnalysis


class ForceVisualizationComparison:
    '''
    Visualize comparison of analytical results for multi cases.

    Parameters
    ----------
    force_analysis_list : a list of ForceAnalysis instance
        All analytical results of force data prepared for visualization.

    '''
    MIN_NUM_FOR_LARGE_FIGURE = 6

    def __init__(self, force_analysis_list, legend, frequency_st_list):
        assert len(force_analysis_list) == len(legend) == len(
            frequency_st_list)

        self._force_analysis_list = force_analysis_list
        self._legend = legend
        self._frequency_st_list = frequency_st_list

    @property
    def force_analysis_list(self):
        '''
        The list of force analysis instances for all concerned cases.
        '''

        return self._force_analysis_list

    def _plot_along_span(self,
                         variables,
                         out_filenames,
                         xlabel,
                         xmin,
                         xmax,
                         return_type,
                         grid=True,
                         figsize=(style.SINGLE_COLUMN_WIDTH / 2,
                                  style.SINGLE_COLUMN_LONG_HEIGHT / 2)):
        '''
        General function for plots comparison along the span axis.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        variables : a list of rank-1 array
            Varialbes for x axis.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.
        return_type : string ("std" or "rms")
            Dict type for return values.
        grid : bool
            True if grid is necessary.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        matplotlib.pyplot.clf()
        matplotlib.pyplot.gcf().set_size_inches(figsize)
        matplotlib.pyplot.grid(grid)
        matplotlib.pyplot.xlim(xmin, xmax)
        axis = matplotlib.pyplot.gca()
        # axis.xaxis.get_offset_text().set_x(0.1)

        axis.locator_params(axis='x', nbins=style.SHORT_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)
        plots = []
        return_dict = {}
        if return_type == 'max':
            for force_analysis, variable, legend in zip(
                    self._force_analysis_list, variables, self._legend):
                z = numpy.linspace(0, 1, num=force_analysis.node_number)
                plots.append(matplotlib.pyplot.plot(variable, z)[0])
                index = numpy.argmax(variable)
                return_dict[legend] = {
                    'maximun value': variable[index],
                    'position': z[index]
                }
        elif return_type == 'std':
            for force_analysis, variable, legend in zip(
                    self._force_analysis_list, variables, self._legend):
                z = numpy.linspace(0, 1, num=force_analysis.node_number)
                plots.append(matplotlib.pyplot.plot(variable, z)[0])
                std = numpy.std(variable)
                rms = numpy.sqrt(numpy.mean(numpy.square(variable)))

                return_dict[legend] = {
                    'standard deviation with respect to the span': std,
                    'root mean square with respect to the span': rms
                }

        matplotlib.pyplot.savefig('')

        x_sci_notaion = axis.xaxis.get_offset_text()
        x_sci_notaion.set_visible(False)
        if x_sci_notaion.get_text():
            xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                          x_sci_notaion.get_text()[1:])
        matplotlib.pyplot.xlabel(xlabel)
        matplotlib.pyplot.ylabel(r'$z\cdot L^{-1}$')

        matplotlib.pyplot.legend(plots, self._legend, loc='best')
        matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            matplotlib.pyplot.savefig(out_filename)

        return return_dict

    def _subplot_along_span(self,
                            out_filenames,
                            variables_list,
                            xlabel_list,
                            xmin_list,
                            xmax_list,
                            grid=False,
                            figsize=(style.SINGLE_COLUMN_WIDTH,
                                     style.SINGLE_COLUMN_SHORT_HEIGHT)):
        '''
        General function for subplots along the span axis.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        variable_list : a list of rank-1 array
            One varialbe for x axis of each column.
        xlabel_list : string
            The labels shown on the x axes.
        xmin_list : a list of floats
            Lower limit for x axes.
        xmax_list : a list of floats
            Upper limit for x axes.
        grid : bool
            True if grid is necessary.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        assert len(variables_list) == len(xlabel_list) == len(
            xmin_list) == len(xmax_list)
        figure, axis_tuple = matplotlib.pyplot.subplots(
            1, len(variables_list), figsize=figsize, sharey=True)
        plots = []
        for axis, variables, xlabel, xmin, xmax in zip(
                axis_tuple, variables_list, xlabel_list, xmin_list, xmax_list):
            axis.set_xlabel(xlabel)
            axis.set_xlim(xmin, xmax)
            # axis.xaxis.set_major_formatter(nice_math_text_form)
            # axis.xaxis.get_offset_text().set_x(0.2)

            axis.locator_params(axis='x', nbins=style.SHORT_XTICK_MAX_LENGTH)
            axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

            axis.grid(grid)
            for force_analysis, variable in zip(
                    self._force_analysis_list, variables):
                z = numpy.linspace(0, 1, num=force_analysis.node_number)
                plots.append(axis.plot(variable, z)[0])

            matplotlib.pyplot.savefig('')

            x_sci_notaion = axis.xaxis.get_offset_text()
            x_sci_notaion.set_visible(False)
            if x_sci_notaion.get_text():
                xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                              x_sci_notaion.get_text()[1:])
            matplotlib.pyplot.xlabel(xlabel)

        axis_tuple[0].set_ylabel(r'$z\cdot L^{-1}$')
        ncol = len(self._legend)
        while ncol > (len(xlabel_list) + 1):
            ncol /= 2
            ncol = int(numpy.ceil(ncol))
        legend = matplotlib.pyplot.figlegend(
            plots,
            self._legend,
            'lower center',
            bbox_to_anchor=[0.5, 1.0],
            ncol=ncol,
            frameon=False)
        figure.tight_layout()
        for out_filename in out_filenames:
            matplotlib.pyplot.savefig(
                out_filename,
                bbox_extra_artists=(legend, ),
                bbox_inches='tight')

    def plot_force_max(self, out_filenames, xlabel, xmin, xmax):
        '''
        Plot force max value along the whole time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.

        '''
        print('plot_force_max')
        force_max = []
        for force_analysis in self.force_analysis_list:
            force_max.append(force_analysis.force_max)
        return self._plot_along_span(
            force_max,
            out_filenames,
            xlabel,
            xmin,
            xmax,
            return_type='max')

    def plot_force_max_abs(self, out_filenames, xlabel, xmin, xmax):
        '''
        Plot force max abs value along the whole time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.

        '''
        print('plot_force_max_abs')
        force_max_abs = []
        for force_analysis in self.force_analysis_list:
            force_max_abs.append(
                numpy.maximum(
                    numpy.absolute(force_analysis.force_max),
                    numpy.absolute(force_analysis.force_min)))
        return self._plot_along_span(
            force_max_abs,
            out_filenames,
            xlabel,
            xmin,
            xmax,
            return_type='max')

    def plot_force_mean(self, out_filenames, xlabel, xmin, xmax):
        '''
        Plot force mean value along the whole time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.

        '''
        print('plot_force_mean')
        force_mean = []
        for force_analysis in self.force_analysis_list:
            force_mean.append(force_analysis.force_mean)
        return self._plot_along_span(
            force_mean,
            out_filenames,
            xlabel,
            xmin,
            xmax,
            return_type='max')

    def plot_curvature_mean(self, out_filenames, xlabel, xmin, xmax):
        '''
        Plot curvature mean value along the whole time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.

        '''
        print('plot_curvature_mean')
        curvature_mean = []
        for force_analysis in self.force_analysis_list:
            curvature_mean.append(-force_analysis.curvature_mean)
            # here set the mean curvature to the opposite value
        return self._plot_along_span(
            curvature_mean,
            out_filenames,
            xlabel,
            xmin,
            xmax,
            return_type='max')

    def plot_force_std(self, out_filenames, xlabel, xmin, xmax):
        '''
        Plot force standard deviations along the whole time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.

        '''
        print('plot_std')
        std = []
        for force_analysis in self.force_analysis_list:
            std.append(force_analysis.force_std)
        return self._plot_along_span(
            std, out_filenames, xlabel, xmin, xmax, return_type='std')

    def plot_curvature_std(self, out_filenames, xlabel, xmin, xmax):
        '''
        Plot curvature standard deviations for all concerned cases along the whole time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.

        '''
        print('plot_curvature_std')
        std = []
        for force_analysis in self.force_analysis_list:
            std.append(force_analysis.curvature_std)
        return self._plot_along_span(
            std, out_filenames, xlabel, xmin, xmax, return_type='std')

    def subplot_force_curvature_std(self, out_filenames, xlabel_list,
                                           xmin_list, xmax_list):
        '''
        Subplot standard deviations of force
        and curvature along the whole time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        xlabel_list : a list of string
            The labels shown on the x axes.
        xmin_list : float
            Lower limit for x axis.
        xmax_list : float
            Upper limit for x axis.

        '''
        print('subplot_force_curvature_std')
        std = []
        curvature_std = []
        for force_analysis in self.force_analysis_list:
            std.append(force_analysis.force_std)
            curvature_std.append(force_analysis.curvature_std)
        self._subplot_along_span(
            out_filenames, [std, curvature_std],
            xlabel_list,
            xmin_list,
            xmax_list,
            grid=True,
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_LONG_HEIGHT / 2))

    def subplot_force_curvature_mean_std(
            self, out_filenames, xlabel_list, xmin_list, xmax_list):
        '''
        Subplot mean values and standard deviations of
        force and curvature along the whole time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        xlabel_list : a list of string
            The labels shown on the x axes.
        xmin_list : float
            Lower limit for x axis.
        xmax_list : float
            Upper limit for x axis.

        '''
        print('subplot_force_curvature_mean_std')
        mean = []
        curvature_mean = []
        std = []
        curvature_std = []
        for force_analysis in self.force_analysis_list:
            mean.append(force_analysis.force_mean)
            curvature_mean.append(-force_analysis.curvature_mean)
            std.append(force_analysis.force_std)
            curvature_std.append(force_analysis.curvature_std)
        self._subplot_along_span(
            out_filenames, [mean, curvature_mean, std, curvature_std],
            xlabel_list,
            xmin_list,
            xmax_list,
            grid=True,
            figsize=(style.ONE_AND_HALF_COLUMN_WIDTH,
                     style.ONE_AND_HALF_COLUMN_SHORT_HEIGHT))

    def plot_modal_weight_force_std(self, out_filenames, ylabel):
        '''
        Plot the standard deviations of modal weight of force.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        ylabel : string
            The label shown on the y axis for the lowest figure.

        '''
        print('plot_modal_weight_force_std')
        matplotlib.pyplot.clf()
        if len(self.force_analysis_list) > self.MIN_NUM_FOR_LARGE_FIGURE:
            matplotlib.pyplot.gcf().set_size_inches(
                style.ONE_AND_HALF_COLUMN_WIDTH,
                style.ONE_AND_HALF_COLUMN_SHORT_HEIGHT)
        else:
            matplotlib.pyplot.gcf().set_size_inches(
                style.SINGLE_COLUMN_WIDTH, style.SINGLE_COLUMN_SHORT_HEIGHT)
        matplotlib.pyplot.xlim(
            self.force_analysis_list[0].frequency_min /
            max(self._frequency_st_list),
            self.force_analysis_list[0].frequency_max /
            min(self._frequency_st_list), )
        axis = matplotlib.pyplot.gca()
        # axis.xaxis.get_offset_text().set_x(0.1)

        axis.locator_params(axis='x', nbins=style.LONG_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

        plots = []
        return_dict = {}
        #rms -- f_n
        for force_analysis, frequency_st, marker_style, legend in\
                zip(
                    self.force_analysis_list,
                    self._frequency_st_list,
                    style.MARKER_STYLES,
                    self._legend
                ):
            plots.append(
                matplotlib.pyplot.
                plot(force_analysis.modal_natural_frequencies[
                    force_analysis.mode_number_min - 1:force_analysis.
                    mode_number_max] / frequency_st, force_analysis.
                     modal_weight_force_std, **marker_style)[0])
            return_dict[legend] = dict(
                zip(
                    list(
                        range(force_analysis.mode_number_min,
                              force_analysis.mode_number_max + 1)),
                    force_analysis.modal_weight_force_std))

            max_index = numpy.argmax(
                force_analysis.modal_weight_force_std)
            return_dict[legend]['dominant_mode_number'] = int(
                force_analysis.mode_number_min + max_index)

            if max_index == 0:
                secondary_max_index = numpy.argmax(
                    force_analysis.modal_weight_force_std[1:]) + 1
                return_dict[legend]['secondary_dominant_mode_number'] = int(
                    force_analysis.mode_number_min + secondary_max_index)
                logging.warning(
                    '''The dominant mode number can be wrong! Maybe the lower limit of mode number should be increased.  And secondary dominant mode number has been stored as an alternative.'''
                )

        y_sci_notaion = axis.yaxis.get_offset_text()
        y_sci_notaion.set_visible(False)
        if y_sci_notaion.get_text():
            ylabel = r"{:s} / {:s}".format(ylabel[:-1],
                                           y_sci_notaion.get_text()[1:])

        matplotlib.pyplot.xlabel(r'$f_n^m\cdot f_s^{-1}$')
        matplotlib.pyplot.ylabel(ylabel)
        matplotlib.pyplot.title('Mode {:d} to {:d}'.format(
            self._force_analysis_list[0].mode_number_min,
            self._force_analysis_list[0].mode_number_max))
        matplotlib.pyplot.grid()
        matplotlib.pyplot.legend(plots, self._legend, loc='best')
        matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            matplotlib.pyplot.savefig(out_filename)

        return return_dict

    def plot_modal_oscillatory_frequency(
            self, out_filenames, ylabel,
            diagonal_line_legend=r'$f_o^m=f_n^m$'):
        '''
        Plot the dominant frequency of force for each mode of each cases.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        ylabel : string
            The label shown on the y axis for the lowest figure.
        diagonal_line_legend : string
            The legend lable shown for the diagonal line.

        '''
        print('plot_modal_oscillatory_frequency')
        matplotlib.pyplot.clf()
        if len(self.force_analysis_list) > self.MIN_NUM_FOR_LARGE_FIGURE:
            matplotlib.pyplot.gcf().set_size_inches(
                style.ONE_AND_HALF_COLUMN_WIDTH,
                style.ONE_AND_HALF_COLUMN_SHORT_HEIGHT)
        else:
            matplotlib.pyplot.gcf().set_size_inches(
                style.SINGLE_COLUMN_WIDTH, style.SINGLE_COLUMN_SHORT_HEIGHT)
        matplotlib.pyplot.xlim(
            self.force_analysis_list[0].frequency_min /
            max(self._frequency_st_list),
            self.force_analysis_list[0].frequency_max /
            min(self._frequency_st_list), )
        axis = matplotlib.pyplot.gca()
        # axis.xaxis.get_offset_text().set_x(0.1)

        axis.locator_params(axis='x', nbins=style.LONG_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

        plots = []
        return_dict = {}
        #rms -- f_n
        for force_analysis, frequency_st, marker_style, legend in\
                zip(
                    self.force_analysis_list,
                    self._frequency_st_list, style.MARKER_STYLES, self._legend):
            plots.append(
                matplotlib.pyplot.plot(
                    force_analysis.modal_natural_frequencies[
                        force_analysis.mode_number_min - 1:
                        force_analysis.mode_number_max] / frequency_st,
                    force_analysis.modal_oscillatory_frequencies /
                    frequency_st, **marker_style)[0])
            return_dict[legend] = dict(
                zip(
                    list(
                        range(force_analysis.mode_number_min,
                              force_analysis.mode_number_max + 1)),
                    force_analysis.modal_oscillatory_frequencies))

        plots.append(
            matplotlib.pyplot.plot(
                matplotlib.pyplot.xlim(),
                matplotlib.pyplot.xlim(),
                marker='',
                color=style.REFERENCE_LIGHT_COLOR,
                linestyle='--')[0])

        y_sci_notaion = axis.yaxis.get_offset_text()
        y_sci_notaion.set_visible(False)
        if y_sci_notaion.get_text():
            ylabel = r"{:s} / {:s}".format(ylabel[:-1],
                                           y_sci_notaion.get_text()[1:])

        matplotlib.pyplot.xlabel(r'$f_n^m\cdot f_s^{-1}$')
        matplotlib.pyplot.ylabel(ylabel)
        matplotlib.pyplot.title('Mode {:d} to {:d}'.format(
            self._force_analysis_list[0].mode_number_min,
            self._force_analysis_list[0].mode_number_max))
        matplotlib.pyplot.grid()
        #matplotlib.pyplot.xlim(0, 3)
        matplotlib.pyplot.legend(
            plots, self._legend + [diagonal_line_legend], loc='best')
        matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            matplotlib.pyplot.savefig(out_filename)

        return return_dict
