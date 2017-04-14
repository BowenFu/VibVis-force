'''
Designed to plot data.
'''
import numpy
import scipy.signal
import scipy.interpolate
import style
import matplotlib.collections
import matplotlib.animation


def _subplot_row(out_filenames,
                 x_variable,
                 y_variable_list,
                 xlabel,
                 ylabel,
                 text_list,
                 xmin,
                 xmax,
                 ymin,
                 ymax,
                 color=style.DARK_COLOR,
                 grid=True,
                 figsize=(style.ONE_AND_HALF_COLUMN_WIDTH,
                          style.ONE_AND_HALF_COLUMN_SHORT_HEIGHT)):
    '''
    General function for subplots of multi rows.

    Parameters
    ----------
    out_filenames : a list of string
        The filenames for saving figures.
    x_variable : rank-1 array
        Varialbe for x axis.
    y_variable_list : rank-2 array or a list of rank-1 array
        A list of varialbes for y axis. One variable for one row.
    xlabel : string
        The label shown on the x axis.
    ylabel : string
        The label shown on the y axis for the lowest figure.
    text_list : a list of string
        Texts to label each row or each y variable.
    xmin : float
        Lower limit for x axis.
    xmax : float
        Upper limit for x axis.
    ymin : float
        Lower limit for y axis.
    ymax : float
        Upper limit for y axis.
    color : string (must be acceptable for matplotlib )
        Line color used in the figure.
    grid : bool
        True if grid is necessary.
    figsize : float tuple (float1, float2)
        Figure size (width, height).

    '''
    assert len(y_variable_list) == len(text_list), print(
        len(y_variable_list), '!=', len(text_list))
    figure, axis_tuple = matplotlib.pyplot.subplots(
        len(y_variable_list), 1, figsize=figsize, sharex=True)

    for i, (axis, y_variable,
            text) in enumerate(zip(axis_tuple, y_variable_list, text_list)):
        axis.plot(x_variable, y_variable, color=color)
        axis.set_xlim(xmin, xmax)
        axis.set_ylim(ymin, ymax)
        # axis.yaxis.set_major_formatter(nice_math_text_form)
        axis.locator_params(axis='x', nbins=style.LONG_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.SHORT_YTICK_MAX_LENGTH)

        if i < len(text_list) - 1:
            matplotlib.pyplot.setp(axis.get_yticklabels(), visible=False)
            axis.yaxis.set_major_formatter(style.NO_POWER_FORM)
        else:
            figure.savefig('')
            x_sci_notaion = axis.xaxis.get_offset_text()
            x_sci_notaion.set_visible(False)
            if x_sci_notaion.get_text():
                xlabel = r"{:s} / {:s}".format(xlabel[:-1],
                                               x_sci_notaion.get_text()[1:])

            y_sci_notaion = axis.yaxis.get_offset_text()
            y_sci_notaion.set_visible(False)
            if y_sci_notaion.get_text():
                ylabel = r"{:s} / {:s}".format(ylabel[:-1],
                                               y_sci_notaion.get_text()[1:])

            axis.set_xlabel(xlabel)
            axis.set_ylabel(ylabel)
        axis.text(
            style.LEFT_CORNER,
            style.TOP_CORNER,
            text,
            horizontalalignment='left',
            verticalalignment='top',
            transform=axis.transAxes)
        axis.grid(grid)
    # matplotlib.pyplot.tight_layout()
    for out_filename in out_filenames:
        figure.savefig(out_filename)
    matplotlib.pyplot.close()


def _contour(
        out_filenames,
        variable,
        x,
        y,
        xlabel,
        ylabel,
        colorbar_min,
        colorbar_max,
        contour_num,
        color=style.DARK_COLOR,
        figsize=(style.SINGLE_COLUMN_WIDTH, style.SINGLE_COLUMN_SHORT_HEIGHT)):
    '''
    General function for contours.

    Parameters
    ----------
    out_filenames : a list of string
        The filenames for saving figures.
    variable : rank-2 array
        Concerned varialbe for contours.
    x : rank-1 array or a list of rank-1 array
        varialbe for x axis.
    y : rank-1 array or a list of rank-1 array
        varialbe for y axis.
    xlabel : string
        The label shown on the x axis.
    ylabel : string
        The label shown on the y axis.
    colorbar_min : float
        Lower limit for the contour.
    colorbar_max : float
        Upper limit for the contour.
    contour_num : int
        Line number for the contour.
    color : string (must be acceptable for matplotlib )
        Line color used in the figure.
    figsize : float tuple (float1, float2)
        Figure size (width, height).

    '''

    figure, axis = matplotlib.pyplot.subplots(figsize=figsize)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)

    axis.locator_params(
            axis='x', nbins=style.LONG_XTICK_MAX_LENGTH)
    axis.locator_params(
            axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

    X, Y = numpy.meshgrid(x, y, indexing='ij')
    contour_range = numpy.linspace(colorbar_min, colorbar_max, contour_num)
    # variable[variable > colorbar_max] = colorbar_max
    # variable[variable < colorbar_min] = colorbar_min
    axis.contour(X, Y, variable, contour_range, colors=color, extend='both')
    # matplotlib.pyplot.tight_layout()
    for out_filename in out_filenames:
        figure.savefig(out_filename)
    matplotlib.pyplot.close('all')


def _contourf(
        out_filenames,
        variable,
        x,
        y,
        xlabel,
        ylabel,
        colorbar_min,
        colorbar_max,
        contourf_num,
        colorbar=False,
        colorbar_zlabel='',
        cmap=style.CMAP_SINGLE,
        figsize=(style.SINGLE_COLUMN_WIDTH, style.SINGLE_COLUMN_SHORT_HEIGHT)):
    '''
    General function for contourfs.

    Parameters
    ----------
    out_filenames : a list of string
        The filenames for saving figures.
    variable : rank-2 array
        Concerned varialbe for contours.
    x : rank-1 array or a list of rank-1 array
        varialbe for x axis.
    y : rank-1 array or a list of rank-1 array
        varialbe for y axis.
    xlabel : string
        The label shown on the x axis.
    ylabel : string
        The label shown on the y axis.
    colorbar_min : float
        Lower limit for the contour.
    colorbar_max : float
        Upper limit for the contour.
    contourf_num : int
        Line number for the contourf.
    colorbar : bool
        True if colorbar is necessary.
    colorbar_zlabel : string
        Label shown on the colorbar.
    cmap : matplotlib cmap
        Color maps used for the contourfs.
    figsize : float tuple (float1, float2)
        Figure size (width, height).

    '''
    figure, axis = matplotlib.pyplot.subplots(figsize=figsize)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)

    axis.locator_params(
        axis='x', nbins=style.LONG_XTICK_MAX_LENGTH)
    axis.locator_params(
        axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

    X, Y = numpy.meshgrid(x, y, indexing='ij')
    contourf_range = numpy.linspace(colorbar_min, colorbar_max, contourf_num)
    # variable[variable > colorbar_max] = colorbar_max
    # variable[variable < colorbar_min] = colorbar_min
    contourf_ = axis.contourf(X, Y, variable, contourf_range, cmap=cmap, extend='both')

    if colorbar:
        colorbar_ax = figure.colorbar(contourf_).ax
        if colorbar_zlabel:
            figure.savefig('')
            colorbar_sci_notaion = colorbar_ax.yaxis.get_offset_text()
            colorbar_sci_notaion.set_visible(False)
            if colorbar_sci_notaion.get_text():
                colorbar_zlabel = "{:s} / {:s}".format(
                    colorbar_zlabel[:-1], colorbar_sci_notaion.get_text()[1:])
            colorbar_ax.set_title(colorbar_zlabel)

    # matplotlib.pyplot.tight_layout()
    for out_filename in out_filenames:
        figure.savefig(out_filename)
    matplotlib.pyplot.close()


def _contourf_contour(
        out_filenames,
        variable,
        x,
        y,
        xlabel,
        ylabel,
        colorbar_min,
        colorbar_max,
        contourf_num,
        contour_num,
        colorbar_zlabel='',
        cmap=style.CMAP_DOUBLE,
        color=style.LIGHT_COLOR,
        figsize=(style.SINGLE_COLUMN_WIDTH, style.SINGLE_COLUMN_SHORT_HEIGHT)):
    '''
    General function for contourfs.

    Parameters
    ----------
    out_filenames : a list of string
        The filenames for saving figures.
    variable : rank-2 array
        Concerned varialbe for contours.
    x : rank-1 array or a list of rank-1 array
        varialbe for x axis.
    y : rank-1 array or a list of rank-1 array
        varialbe for y axis.
    xlabel : string
        The label shown on the x axis.
    ylabel : string
        The label shown on the y axis.
    colorbar_min : float
        Lower limit for the contour.
    colorbar_max : float
        Upper limit for the contour.
    contourf_num : int
        Line number for the contourf.
    contour_num : int
        Line number for the contour.
    colorbar_zlabel : string
        Label shown on the colorbar.
    cmap : matplotlib cmap
        Color maps used for the contourfs.
    color : string (must be acceptable for matplotlib )
        Line color used in the figure.
    figsize : float tuple (float1, float2)
        Figure size (width, height).

    '''
    figure, axis = matplotlib.pyplot.subplots(figsize=figsize)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)

    axis.locator_params(
        axis='x', nbins=style.LONG_XTICK_MAX_LENGTH)
    axis.locator_params(
        axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

    X, Y = numpy.meshgrid(x, y, indexing='ij')
    contourf_range = numpy.linspace(colorbar_min, colorbar_max, contourf_num)
    contour_range = numpy.linspace(colorbar_min, colorbar_max, contour_num)
    # variable[variable > colorbar_max] = colorbar_max
    # variable[variable < colorbar_min] = colorbar_min
    axis.contour(X, Y, variable, contour_range, colors=color, extend='both')
    contourf_ = axis.contourf(X, Y, variable, contourf_range, cmap=cmap, extend='both')
    contour_range = numpy.round(
        contour_range[::2],
        decimals=-int(
            numpy.floor(numpy.log10(numpy.diff(contour_range[::2])[0]))))
    colorbar_ax = figure.colorbar(contourf_, ticks=contour_range).ax
    if colorbar_zlabel:
        figure.savefig('')
        colorbar_sci_notaion = colorbar_ax.yaxis.get_offset_text()
        colorbar_sci_notaion.set_visible(False)
        if colorbar_sci_notaion.get_text():
            colorbar_zlabel = "{:s} / {:s}".format(
                colorbar_zlabel[:-1], colorbar_sci_notaion.get_text()[1:])
        colorbar_ax.set_title(colorbar_zlabel)

    # colorbar_ax.yaxis.set_major_formatter(style.NO_POWER_FORM)

    matplotlib.pyplot.tight_layout()
    for out_filename in out_filenames:
        figure.savefig(out_filename)
    matplotlib.pyplot.close()


class ForceVisualization:
    '''
    Visualize all analytical results of force data obtained from ForceAnalysis.

    Parameters
    ----------
    force_analysis : ForceAnalysis instance
        All analytical results of force data prepared for visualization.

    '''

    def __init__(self, force_analysis):
        self._force_analysis = force_analysis

    @property
    def force_analysis(self):
        '''
        ForceAnalysis instance containing all analytical results of force data.
        '''

        return self._force_analysis

    def _plot_along_time_function(self,
            out_filenames,
            time_function,
            variable,
            xlabel,
            start_time,
            end_time,
            ylabel,
            ymin,
            ymax,
            color=style.DARK_COLOR,
            grid=True,
            figsize=(style.SINGLE_COLUMN_WIDTH,
                style.SINGLE_COLUMN_SHORT_HEIGHT)):
        '''
        General function for plots y along x.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        time_function : a ufunc of time
            Function results for x axis.
        variable : rank-1 array
            Varialbe for y axis.
        xlabel : string
            The label shown on the x axis for the lowest figure.
        start_time : float
            Lower limit for time function.
        end_time : float
            Upper limit for time function.
        ylabel : string
            The label shown on the y axis for the lowest figure.
        ymin : float
            Lower limit for y axis.
        ymax : float
            Upper limit for y axis.
        color : string (must be acceptable for matplotlib )
            Line color used in the figure.
        grid : bool
            True if grid is necessary.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        figure, axis = matplotlib.pyplot.subplots(figsize=figsize)

        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        selected_time = self.force_analysis.time[time_index]
        selected_variable = variable[time_index]
        matplotlib.pyplot.plot(
            selected_time,
            selected_variable,
            color=color)

        # matplotlib.pyplot.xlim(xmin, xmax)
        axis.set_ylim(ymin, ymax)
        matplotlib.pyplot.grid(grid)

        axis.locator_params(axis='x', nbins=style.LONG_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.SHORT_YTICK_MAX_LENGTH)

        figure.savefig('')

        x_sci_notaion = axis.xaxis.get_offset_text()
        x_sci_notaion.set_visible(False)
        if x_sci_notaion.get_text():
            xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                          x_sci_notaion.get_text()[1:])

        y_sci_notaion = axis.yaxis.get_offset_text()
        y_sci_notaion.set_visible(False)
        if y_sci_notaion.get_text():
            ylabel = r"{:s} / {:s}".format(ylabel[:-1],
                                           y_sci_notaion.get_text()[1:])

        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        # matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            figure.savefig(out_filename)
        matplotlib.pyplot.close()


    def _plot_along_time(self,
                         out_filenames,
                         variable,
                         start_time,
                         end_time,
                         ylabel,
                         ymin,
                         ymax,
                         xlabel=r'$t\mathrm{\ (s)}$',
                         color=style.DARK_COLOR,
                         reference_line_factor_tuple=(),
                         reference_line_share_maximum=False,
                         tol=1e-2,
                         input_time_referece_tuple=(),
                         grid=True,
                         figsize=(style.SINGLE_COLUMN_WIDTH,
                                  style.SINGLE_COLUMN_SHORT_HEIGHT)):
        '''
        General function for plots along the time axis.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        variable : rank-1 array
            Varialbe for y axis.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        ylabel : string
            The label shown on the y axis for the lowest figure.
        ymin : float
            Lower limit for y axis.
        ymax : float
            Upper limit for y axis.
        color : string (must be acceptable for matplotlib )
            Line color used in the figure.
        grid : bool
            True if grid is necessary.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        figure, axis = matplotlib.pyplot.subplots(figsize=figsize)

        xlabel = xlabel
        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        selected_time = self.force_analysis.time[time_index]
        selected_variable = variable[time_index]
        axis.plot(
            selected_time,
            selected_variable,
            # self.force_analysis.time[start_index: end_index],
            # variable[start_index: end_index],
            color=color)

        assert len(reference_line_factor_tuple) in (0, 2)
        if len(reference_line_factor_tuple) is 2:
            reference_max = numpy.amax(selected_variable)
            reference_min = numpy.amin(selected_variable)
            if reference_line_share_maximum:
                reference_max = max(reference_max, abs(reference_min))
                reference_min = -reference_max

            reference_factor_max, reference_factor_min = reference_line_factor_tuple

            upper_reference = reference_factor_max * reference_max
            lower_reference = reference_factor_min * reference_min

            upper_intersection_times = []
            lower_intersection_times = []
            if upper_reference >= 0:
                upper_intersection_times = selected_time[(
                    selected_variable > upper_reference * (1 - tol)) & (
                        selected_variable < upper_reference * (1 + tol))]
            else:
                upper_intersection_times = selected_time[(
                    selected_variable < upper_reference * (1 - tol)) & (
                        selected_variable > upper_reference * (1 + tol))]
            if lower_reference >= 0:
                lower_intersection_times = selected_time[(
                    selected_variable > lower_reference * (1 - tol)) & (
                        selected_variable < lower_reference * (1 + tol))]
            else:
                lower_intersection_times = selected_time[(
                    selected_variable < lower_reference * (1 - tol)) & (
                        selected_variable > lower_reference * (1 + tol))]
            intersection_time_min = numpy.amin(numpy.r_[
                upper_intersection_times, lower_intersection_times])
            intersection_time_max = numpy.amax(numpy.r_[
                upper_intersection_times, lower_intersection_times])

            axis.axhline(
                y=upper_reference, **style.REFERENCE_LINE_STYLE)
            axis.axhline(
                y=lower_reference, **style.REFERENCE_LINE_STYLE)
            axis.axvline(
                x=intersection_time_min, **style.REFERENCE_LINE_STYLE)
            axis.axvline(
                x=intersection_time_max, **style.REFERENCE_LINE_STYLE)
        # none reference
        elif input_time_referece_tuple:
            assert len(input_time_referece_tuple) is 2
            axis.axvline(
                x=input_time_referece_tuple[0], **style.REFERENCE_LINE_STYLE)
            axis.axvline(
                x=input_time_referece_tuple[1], **style.REFERENCE_LINE_STYLE)

        axis.set_xlim(start_time, end_time)
        axis.set_ylim(ymin, ymax)
        matplotlib.pyplot.grid(grid)

        axis.locator_params(axis='x', nbins=style.LONG_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.SHORT_YTICK_MAX_LENGTH)

        figure.savefig('')

        x_sci_notaion = axis.xaxis.get_offset_text()
        x_sci_notaion.set_visible(False)
        if x_sci_notaion.get_text():
            xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                          x_sci_notaion.get_text()[1:])

        y_sci_notaion = axis.yaxis.get_offset_text()
        y_sci_notaion.set_visible(False)
        if y_sci_notaion.get_text():
            ylabel = r"{:s} / {:s}".format(ylabel[:-1],
                                           y_sci_notaion.get_text()[1:])

        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        # matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            figure.savefig(out_filename)
        matplotlib.pyplot.close()
        if reference_line_factor_tuple:
            return (intersection_time_min, intersection_time_max)

    def _subplot_along_time(self,
                            out_filenames,
                            variable_list,
                            start_time,
                            end_time,
                            ylabel,
                            text_list,
                            ymin,
                            ymax,
                            xlabel=r'$t\mathrm{\ (s)}$',
                            grid=False,
                            figsize=(style.SINGLE_COLUMN_WIDTH,
                                     style.SINGLE_COLUMN_SHORT_HEIGHT)):
        '''
        General function for plots along the time axis.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        variable_list : a list of rank-1 array
            One varialbe for each row.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        ylabel : string
            The label shown on the y axis for the lowest figure.
        ymin : float
            Lower limit for y axis.
        ymax : float
            Upper limit for y axis.
        grid : bool
            True if grid is necessary.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        _subplot_row(
            out_filenames=out_filenames,
            x_variable=self.force_analysis.time[time_index],
            y_variable_list=variable_list[:, time_index],
            # x_variable=self.force_analysis.time[start_index: end_index],
            # y_variable_list=variable_list[:, start_index: end_index],
            xlabel=xlabel,
            ylabel=ylabel,
            text_list=text_list,
            xmin=start_time,
            xmax=end_time,
            ymin=ymin,
            ymax=ymax,
            grid=grid,
            figsize=figsize)

    def _plot_along_span(self,
                         out_filenames,
                         variable,
                         xlabel,
                         xmin,
                         xmax,
                         color=style.DARK_COLOR,
                         clf=True,
                         save=True,
                         grid=False,
                         figsize=(style.SINGLE_COLUMN_WIDTH / 2,
                                  style.SINGLE_COLUMN_LONG_HEIGHT / 2)):
        '''
        General function for plots along the span axis.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        variable : rank-1 array
            Varialbe for x axis.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.
        color : string (must be acceptable for matplotlib )
            Line color used in the figure.
        clf : bool
            True for starting a new figure.
        save : bool
            True for ending a figure.
        grid : bool
            True if grid is necessary.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        if clf:
            matplotlib.pyplot.clf()
            matplotlib.pyplot.gcf().set_size_inches(figsize)
        axis = matplotlib.pyplot.gca()
        figure = matplotlib.pyplot.gcf()
        axis.locator_params(axis='x', nbins=style.SHORT_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

        matplotlib.pyplot.grid(grid)
        matplotlib.pyplot.plot(
            variable, self.force_analysis.span, color=color)
        matplotlib.pyplot.xlim(xmin, xmax)

        figure.savefig('')

        x_sci_notaion = axis.xaxis.get_offset_text()
        x_sci_notaion.set_visible(False)
        if x_sci_notaion.get_text():
            xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                          x_sci_notaion.get_text()[1:])

        y_sci_notaion = axis.yaxis.get_offset_text()
        y_sci_notaion.set_visible(False)
        if y_sci_notaion.get_text():
            ylabel = r"{:s} / {:s}".format(ylabel[:-1],
                                           y_sci_notaion.get_text()[1:])

        matplotlib.pyplot.xlabel(xlabel)
        matplotlib.pyplot.ylabel(r'$z\cdot L^{-1}$')

        matplotlib.pyplot.tight_layout()
        if save:
            for out_filename in out_filenames:
                figure.savefig(out_filename)
            matplotlib.pyplot.close()

    def _subplot_along_span(self,
                            out_filenames,
                            variable_list,
                            xlabel_list,
                            xmin_list,
                            xmax_list,
                            color=style.DARK_COLOR,
                            grid=False,
                            figsize=(style.ONE_AND_HALF_COLUMN_WIDTH,
                                     style.ONE_AND_HALF_COLUMN_SHORT_HEIGHT)):
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
        color : string (must be acceptable for matplotlib )
            Line color used in the figure.
        grid : bool
            True if grid is necessary.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        assert len(variable_list) == len(xlabel_list) == len(xmin_list) == len(
            xmax_list)
        figure, axis_tuple = matplotlib.pyplot.subplots(
            1, len(variable_list), figsize=figsize, sharey=True)
        axis_tuple[0].set_ylabel(r'$z\cdot L^{-1}$')
        for axis, variable, xlabel, xmin, xmax in\
                zip(axis_tuple, variable_list, xlabel_list, xmin_list, xmax_list):
            axis.plot(variable, self.force_analysis.span, color=color)
            axis.set_xlim(xmin, xmax)
            axis.locator_params(axis='x', nbins=style.SHORT_XTICK_MAX_LENGTH)
            axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

            figure.savefig('')

            x_sci_notaion = axis.xaxis.get_offset_text()
            x_sci_notaion.set_visible(False)
            if x_sci_notaion.get_text():
                xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                              x_sci_notaion.get_text()[1:])

            axis.set_xlabel(xlabel)

            # axis.xaxis.set_major_formatter(nice_math_text_form)
            axis.grid(grid)
        # matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            figure.savefig(out_filename)
        matplotlib.pyplot.close()

    def plot_force_along_time_function(
            self,
            out_filenames,
            node_i,
            start_time,
            end_time,
            xlabel,
            force_label,
            force_min,
            force_max,
            time_function,
            grid=False,
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_SHORT_HEIGHT)):
        '''
        Plot the time history of force for specific node.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        node_i : int
            Node number of the force to be plotted.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        force_label : string
            The label shown on the force axis for the lowest figure.
        force_min : float
            Lower limit for force axis.
        force_max : float
            Upper limit for force axis.

        '''
        print('plot_force_along_time_function')
        return self._plot_along_time_function(
            out_filenames=out_filenames,
            time_function=time_function,
            variable=self.force_analysis.force[:, node_i - 1],
            xlabel=xlabel,
            start_time=start_time,
            end_time=end_time,
            ylabel=force_label,
            ymin=force_min,
            ymax=force_max,
            grid=grid,
            figsize=figsize)

    def plot_time_history_force(
            self,
            out_filenames,
            node_i,
            start_time,
            end_time,
            force_label,
            force_min,
            force_max,
            xlabel=r'$t\mathrm{\ (s)}$',
            reference_line_factor_tuple=(),
            grid=False,
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_SHORT_HEIGHT / 2)):
        '''
        Plot the time history of force for specific node.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        node_i : int
            Node number of the force to be plotted.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        force_label : string
            The label shown on the force axis for the lowest figure.
        force_min : float
            Lower limit for force axis.
        force_max : float
            Upper limit for force axis.

        '''
        print('plot_time_history_force')
        return self._plot_along_time(
            out_filenames=out_filenames,
            variable=self.force_analysis.force[:, node_i - 1],
            start_time=start_time,
            end_time=end_time,
            ylabel=force_label,
            ymin=force_min,
            ymax=force_max,
            reference_line_factor_tuple=reference_line_factor_tuple,
            grid=grid,
            figsize=figsize)

    def plot_time_history_force_deviation(
            self,
            out_filenames,
            node_i,
            start_time,
            end_time,
            force_deviation_label,
            force_deviation_min,
            force_deviation_max,
            xlabel=r'$t\mathrm{\ (s)}$',
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_SHORT_HEIGHT / 2)):
        '''
        Plot the time history of force deviation for specific node.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        node_i : int
            Node number of the force deviation to be plotted.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        force_deviation_label : string
            The label shown on the force deviation axis for the lowest figure.
        force_deviation_min : float
            Lower limit for force deviation axis.
        force_deviation_max : float
            Upper limit for force deviation axis.

        '''
        print('plot_time_history_force_deviation')
        self._plot_along_time(
            out_filenames=out_filenames,
            variable=self.force_analysis.force_deviation[:, node_i -
                                                                    1],
            start_time=start_time,
            end_time=end_time,
            xlabel=xlabel,
            ylabel=force_deviation_label,
            ymin=force_deviation_min,
            ymax=force_deviation_max,
            figsize=figsize)

    def plot_time_history_velocity(
            self,
            out_filenames,
            node_i,
            start_time,
            end_time,
            velocity_label,
            velocity_min,
            velocity_max,
            reduced_velocity=None,
            xlabel=r'$t\mathrm{\ (s)}$',
            reference_line_factor_tuple=(),
            input_time_referece_tuple=(),
            grid=False,
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_SHORT_HEIGHT / 2)):
        '''
        Plot the time history of velocity for specific node.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        node_i : int
            Node number of the velocity to be plotted.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        velocity_label : string
            The label shown on the velocity axis for the lowest figure.
        velocity_min : float
            Lower limit for velocity axis.
        velocity_max : float
            Upper limit for velocity axis.
        reduced_velocity : string
            Means to obtain non-dimensional velocity.

        '''
        print('plot_time_history_velocity')
        if reduced_velocity is None:
            velocity = self.force_analysis.velocity
        elif reduced_velocity is 'fundamental_natural_frequency':
            velocity = self.force_analysis.velocity / (
                self.force_analysis.modal_natural_frequencies[0])
        else:
            raise TypeError('Wrong reduced velocity type.')
        return self._plot_along_time(
            out_filenames=out_filenames,
            variable=velocity[:, node_i - 1],
            start_time=start_time,
            end_time=end_time,
            xlabel=xlabel,
            ylabel=velocity_label,
            ymin=velocity_min,
            ymax=velocity_max,
            reference_line_factor_tuple=reference_line_factor_tuple,
            input_time_referece_tuple=input_time_referece_tuple,
            grid=grid,
            figsize=figsize)

    def plot_modal_weight_force(self, out_filenames, mode_i, start_time,
                                       end_time, ylabel, ymin, ymax,
                                       xlabel=r'$t\mathrm{\ (s)}$',
                                       ):
        '''
        Plot the modal weight history of force for specific mode.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        mode_i : int
            Mode number of the modal weight history to be plotted.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        ylabel : string
            The label shown on the y axis for the lowest figure.
        ymin : float
            Lower limit for y axis.
        ymax : float
            Upper limit for y axis.

        '''
        if not self.force_analysis.modal_analysis:
            return
        print('plot_modal_weight_force')
        self._plot_along_time(
            out_filenames=out_filenames,
            variable=self.force_analysis.
            modal_weight_force[:, mode_i -
                                      self.force_analysis.mode_number_min],
            start_time=start_time,
            end_time=end_time,
            xlabel=xlabel,
            ylabel=ylabel,
            ymin=ymin,
            ymax=ymax)

    def plot_outline(self,
                     out_filenames,
                     start_time,
                     end_time,
                     line_number,
                     xlabel,
                     xmin,
                     xmax,
                     show_min_max=False,
                     show_y=True,
                     figsize=(style.SINGLE_COLUMN_WIDTH / 2,
                              style.SINGLE_COLUMN_LONG_HEIGHT / 2)):
        '''
        Plot instaneous force at multi-times.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        line_number : int
            Line number of the figure.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.
        show_min_max : bool
            True for showing minimum and maximum force along the whole time.
        show_y : bool
            True if ylabel is needed.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        print('plot_outline')
        figure, axis = matplotlib.pyplot.subplots(figsize=figsize)

        ylabel = r'$z\cdot L^{-1}$'

        if show_min_max:
            axis.plot(self.force_analysis.force_min,
                                   self.force_analysis.span,
                                   **style.LIGHT_LINE_STYLE)
            axis.plot(self.force_analysis.force_max,
                                   self.force_analysis.span,
                                   **style.LIGHT_LINE_STYLE)

        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)

        # for time_index in numpy.linspace(
        #         start_index, end_index, line_number, endpoint=False):
        step = int(sum(time_index) // line_number)
        for force in self.force_analysis.force[
                time_index, :][::step, :]:
            # start_index, end_index, line_number, endpoint=False):
            # time_index = self.force_analysis.time_index(time)
            # time_index = int(numpy.around(time_index))
            axis.plot(
                force,
                # self.force_analysis.force[time_index, :],
                self.force_analysis.span,
                **style.SINGLE_LINE_STYLE)
        axis.set_xlim(xmin, xmax)
        axis.locator_params(axis='x', nbins=style.SHORT_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

        figure.savefig('')

        x_sci_notaion = axis.xaxis.get_offset_text()
        x_sci_notaion.set_visible(False)
        if x_sci_notaion.get_text():
            xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                          x_sci_notaion.get_text()[1:])

        axis.set_xlabel(xlabel)

        axis.set_ylabel(ylabel)
        if not show_y:
            matplotlib.pyplot.setp(axis.get_yticklabels(), visible=False)
            axis.set_ylabel(r' ')

        matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            figure.savefig(out_filename)
        matplotlib.pyplot.close()

    def plot_deviation_outline(self,
                               out_filenames,
                               start_time,
                               end_time,
                               line_number,
                               xlabel,
                               xmin,
                               xmax,
                               show_min_max=False,
                               show_y=True,
                               figsize=(style.SINGLE_COLUMN_WIDTH / 2,
                                        style.SINGLE_COLUMN_LONG_HEIGHT / 2)):
        '''
        Plot instaneous force deviation at multi-times.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        line_number : int
            Line number of the figure.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.
        show_min_max : bool
            True for showing minimum and maximum force along the whole time.
        show_y : bool
            True if ylabel is needed.
        figsize : float tuple (float1, float2)
            Figure size (width, height).

        '''
        print('plot_deviation_outline')
        matplotlib.pyplot.clf()
        matplotlib.pyplot.gcf().set_size_inches(figsize)
        if show_min_max:
            matplotlib.pyplot.plot(self.force_analysis.force_min -
                                   self.force_analysis.force_mean,
                                   self.force_analysis.span,
                                   **style.LIGHT_LINE_STYLE)
            matplotlib.pyplot.plot(self.force_analysis.force_max -
                                   self.force_analysis.force_mean,
                                   self.force_analysis.span,
                                   **style.LIGHT_LINE_STYLE)

        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)

        #  for specific time period
        # for time in numpy.linspace(start_time, end_time, line_number,
        # endpoint=False):
        # for time_index in numpy.linspace(
        #        start_index, end_index, line_number, endpoint=False):
        step = int(sum(time_index) // line_number)
        for force_deviation in self.force_analysis.force_deviation[
                time_index, :][::step, :]:
            # time_index = self.force_analysis.time_index(time)
            # time_index = int(numpy.around(time_index))
            matplotlib.pyplot.plot(
                # self.force_analysis.force_deviation[time_index,
                # :],
                force_deviation,
                self.force_analysis.span,
                **style.SINGLE_LINE_STYLE)
        matplotlib.pyplot.xlim(xmin, xmax)
        axis = matplotlib.pyplot.gca()
        axis.locator_params(axis='x', nbins=style.SHORT_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)

        matplotlib.pyplot.savefig('')

        x_sci_notaion = axis.xaxis.get_offset_text()
        x_sci_notaion.set_visible(False)
        if x_sci_notaion.get_text():
            xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                          x_sci_notaion.get_text()[1:])

        matplotlib.pyplot.xlabel(xlabel)

        matplotlib.pyplot.ylabel(r'$z\cdot L^{-1}$')
        if not show_y:
            matplotlib.pyplot.setp(axis.get_yticklabels(), visible=False)
            matplotlib.pyplot.ylabel(r' ')

        # matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            matplotlib.pyplot.savefig(out_filename)
        matplotlib.pyplot.close()

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
        self._plot_along_span(
            out_filenames,
            self.force_analysis.force_mean,
            xlabel,
            xmin,
            xmax,
            color=style.DARK_COLOR,
            clf=True,
            save=False,
            grid=True)
        self._plot_along_span(
            out_filenames,
            self.force_analysis.force_min,
            xlabel,
            xmin,
            xmax,
            color=style.LIGHT_COLOR,
            clf=False,
            save=False,
            grid=True)
        self._plot_along_span(
            out_filenames,
            self.force_analysis.force_max,
            xlabel,
            xmin,
            xmax,
            color=style.LIGHT_COLOR,
            clf=False,
            save=True,
            grid=True)
        return {
            'mean': self.force_analysis.force_mean.tolist(),
            'min': self.force_analysis.force_min.tolist(),
            'max': self.force_analysis.force_max.tolist()
        }

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
        # here plot mean curvature to the opposite value
        self._plot_along_span(
            out_filenames,
            -self.force_analysis.curvature_mean,
            xlabel,
            xmin,
            xmax,
            color=style.DARK_COLOR,
            clf=True,
            save=False,
            grid=True)
        self._plot_along_span(
            out_filenames,
            -self.force_analysis.curvature_min,
            xlabel,
            xmin,
            xmax,
            color=style.LIGHT_COLOR,
            clf=False,
            save=False,
            grid=True)
        self._plot_along_span(
            out_filenames,
            -self.force_analysis.curvature_max,
            xlabel,
            xmin,
            xmax,
            color=style.LIGHT_COLOR,
            clf=False,
            save=True,
            grid=True)
        #  return actual value rather than opposite value.
        return {
            'curvature_mean': self.force_analysis.curvature_mean.tolist(),
            'min': self.force_analysis.curvature_min.tolist(),
            'max': self.force_analysis.curvature_max.tolist()
        }

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
        print('plot_force_std')
        self._plot_along_span(
            out_filenames,
            self.force_analysis.force_std,
            xlabel,
            xmin,
            xmax,
            grid=True)
        return self.force_analysis.force_std

    def plot_curvature_std(self, out_filenames, xlabel, xmin, xmax):
        '''
        Plot curvature standard deviations along the whole time.

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
        self._plot_along_span(
            out_filenames,
            self.force_analysis.curvature_std,
            xlabel,
            xmin,
            xmax,
            grid=True)
        return self.force_analysis.curvature_std

    def subplot_force_curvature_std(self, out_filenames, xlabel_list,
                                           xmin_list, xmax_list):
        '''
        Subplot standard deviations of force and curvature along the whole time.

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
        self._subplot_along_span(
            out_filenames, [
                self.force_analysis.force_std,
                self.force_analysis.curvature_std,
            ],
            xlabel_list,
            xmin_list,
            xmax_list,
            grid=True,
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_LONG_HEIGHT / 2))

    def subplot_force_curvature_mean_std(
            self, out_filenames, xlabel_list, xmin_list, xmax_list):
        '''
        Subplot mean values and standard deviations of force and curvature along the whole time.

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
        self._subplot_along_span(
            out_filenames, [
                self.force_analysis.force_mean,
                -self.force_analysis.curvature_mean,
                self.force_analysis.force_std,
                self.force_analysis.curvature_std
            ],
            xlabel_list,
            xmin_list,
            xmax_list,
            grid=True,
            figsize=(style.ONE_AND_HALF_COLUMN_WIDTH,
                     style.ONE_AND_HALF_COLUMN_SHORT_HEIGHT))

    def subplot_modal_weight_force(
            self,
            out_filenames,
            start_time,
            end_time,
            ylabel,
            ymin,
            ymax,
            xlabel=r'$t\mathrm{\ (s)}$',
            figsize=(style.FULL_WIDTH * 2 / 3, style.FULL_SHORT_HEIGHT)):
        '''
        Plot the time history of modal weight of force.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        ylabel : string
            The label shown on the y axis for the lowest figure.
        ymin : float
            Lower limit for y axis.
        ymax : float
            Upper limit for y axis.

        '''
        if not self.force_analysis.modal_analysis:
            return
        print('subplot_modal_weight_force')
        self._subplot_along_time(
            out_filenames=out_filenames,
            variable_list=self.force_analysis.modal_weight_force.T,
            start_time=start_time,
            end_time=end_time,
            xlabel=xlabel,
            ylabel=ylabel,
            text_list=[
                r'Mode {:d}'.format(mode_i)
                for mode_i in
                range(self.force_analysis.mode_number_min,
                      self.force_analysis.mode_number_max + 1)
            ],
            ymin=ymin,
            ymax=ymax,
            grid=False,
            figsize=figsize)

    def contourf_contour_spatio_temporal_force(
            self,
            out_filenames,
            start_time,
            end_time,
            colorbar_min,
            colorbar_max,
            contourf_num,
            contour_num,
            xlabel=r'$t\mathrm{\ (s)}$',
            colorbar_zlabel='', ):
        '''
        Spatio temporal contour over contourf of force.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        colorbar_min : float
            Lower limit for the contour.
        colorbar_max : float
            Upper limit for the contour.
        contour_num : int
            Line number for the contour.
        contourf_num : int
            Line number for the contourf.
        colorbar_zlabel : string
            Label shown on the colorbar.

        '''
        print('contourf_contour_spatio_temporal_force')
        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        time = self.force_analysis.time[time_index]
        _contourf_contour(
            out_filenames=out_filenames,
            variable=self.force_analysis.force_deviation[
                time_index],
            x=time,
            y=self.force_analysis.span,
            xlabel=r'$t\mathrm{\ (s)}$',
            ylabel=r'$z\cdot L^{-1}$',
            colorbar_min=colorbar_min,
            colorbar_max=colorbar_max,
            contourf_num=contourf_num,
            contour_num=contour_num,
            colorbar_zlabel=colorbar_zlabel)

    def contour_spatio_temporal_force(self, out_filenames, start_time,
                                             end_time, colorbar_min,
                                             colorbar_max, contour_num,
                                             xlabel=r'$t\mathrm{\ (s)}$',
                                             ):
        '''
        Spatio temporal contour for force.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        colorbar_min : float
            Lower limit for the contour.
        colorbar_max : float
            Upper limit for the contour.
        contour_num : int
            Line number for the contour.

        '''
        print('contour_spatio_temporal_force')
        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        time = self.force_analysis.time[time_index]
        _contour(
            out_filenames=out_filenames,
            variable=self.force_analysis.force_deviation[
                time_index],
            x=time,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=r'$z\cdot L^{-1}$',
            colorbar_min=colorbar_min,
            colorbar_max=colorbar_max,
            contour_num=contour_num)

    def contourf_contour_spatio_temporal_curvature(
            self,
            out_filenames,
            start_time,
            end_time,
            colorbar_min,
            colorbar_max,
            contourf_num,
            contour_num,
            xlabel=r'$t\mathrm{\ (s)}$',
            colorbar_zlabel='', ):
        '''
        Spatio temporal contour over contourf of curvature.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        colorbar_min : float
            Lower limit for the contour.
        colorbar_max : float
            Upper limit for the contour.
        contour_num : int
            Line number for the contour.
        contourf_num : int
            Line number for the contourf.
        colorbar_zlabel : string
            Label shown on the colorbar.

        '''
        print('contourf_contour_spatio_temporal_curvature')
        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        time = self.force_analysis.time[time_index]
        _contourf_contour(
            out_filenames=out_filenames,
            variable=self.force_analysis.curvature_deviation[time_index],
            x=time,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=r'$z\cdot L^{-1}$',
            colorbar_min=colorbar_min,
            colorbar_max=colorbar_max,
            contourf_num=contourf_num,
            contour_num=contour_num,
            colorbar_zlabel=colorbar_zlabel)

    def contour_spatio_temporal_curvature(self, out_filenames, start_time,
                                          end_time, colorbar_min, colorbar_max,
                                          contour_num,
                                          xlabel=r'$t\mathrm{\ (s)}$',
                                          ):
        '''
        Spatio temporal contour for curvature.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        colorbar_min : float
            Lower limit for the contour.
        colorbar_max : float
            Upper limit for the contour.
        contour_num : int
            Line number for the contour.

        '''
        print('contour_spatio_temporal_curvature')
        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        time = self.force_analysis.time[time_index]
        _contour(
            out_filenames=out_filenames,
            variable=self.force_analysis.curvature_deviation[time_index],
            x=time,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=r'$z\cdot L^{-1}$',
            colorbar_min=colorbar_min,
            colorbar_max=colorbar_max,
            contour_num=contour_num)

    def subplot_span_force(self,
                                  out_filenames,
                                  force_label,
                                  start_time,
                                  end_time,
                                  force_min,
                                  force_max,
                                  xlabel=r'$t\mathrm{\ (s)}$',
                                  num=9,
                                  figsize=(style.SINGLE_COLUMN_WIDTH,
                                           style.SINGLE_COLUMN_LONG_HEIGHT)):
        '''
        Subplots of force along time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        force_label : string
            Label for force.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        force_min : float
            Lower limit for force axis.
        force_max : float
            Upper limit for force axis.
        num : int
            Number of nodes for subplots.
        '''

        num += 1
        print('subplot_span_force')
        step = int(self.force_analysis.node_number // num)
        span_list = numpy.arange(0, 1,
                                 step / self.force_analysis.node_number)
        self._subplot_along_time(
            out_filenames=out_filenames,
            variable_list=self.force_analysis.force[:, ::step].T[
                -2:0:-1],
            start_time=start_time,
            end_time=end_time,
            xlabel=xlabel,
            ylabel=force_label,
            text_list=[
                r'$z\cdot L^{{-1}}={:.1f}$'.format(span)
                for span in span_list[-2:0:-1]#inverse
            ],
            ymin=force_min,
            ymax=force_max,
            grid=False,
            figsize=figsize)

    def subplot_span_force_deviation(
            self,
            out_filenames,
            force_deviation_label,
            start_time,
            end_time,
            force_deviation_min,
            force_deviation_max,
            xlabel=r'$t\mathrm{\ (s)}$',
            num=9,
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_LONG_HEIGHT)):
        '''
        Subplots of force deviation along time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        force_deviation_label : string
            Label for force.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        force_deviation_min : float
            Lower limit for force axis.
        force_deviation_max : float
            Upper limit for force axis.
        num : int
            Number of nodes for subplots.
        '''

        num += 1
        print('subplot_span_force_deviation')
        step = int(self.force_analysis.node_number // num)
        span_list = numpy.arange(0, 1,
                                 step / self.force_analysis.node_number)
        self._subplot_along_time(
            out_filenames=out_filenames,
            variable_list=self.force_analysis.
            force_deviation[:, ::step].T[-2:0:-1],
            start_time=start_time,
            end_time=end_time,
            xlabel=xlabel,
            ylabel=force_deviation_label,
            text_list=[
                r'$z\cdot L^{{-1}}={:.1f}$'.format(span)
                for span in span_list[-2:0:-1]#inverse
            ],
            ymin=force_deviation_min,
            ymax=force_deviation_max,
            grid=False,
            figsize=figsize)

    def subplot_span_velocity(self,
                              out_filenames,
                              velocity_label,
                              start_time,
                              end_time,
                              velocity_min,
                              velocity_max,
                              num=9,
                              figsize=(style.SINGLE_COLUMN_WIDTH,
                                       style.SINGLE_COLUMN_LONG_HEIGHT)):
        '''
        Subplots of velocity along time.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        velocity_label : string
            Label for velocity.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        velocity_min : float
            Lower limit for velocity axis.
        velocity_max : float
            Upper limit for velocity axis.
        num : int
            Number of nodes for subplots.
        '''

        num += 1
        print('subplot_span_velocity')
        step = int(self.force_analysis.node_number // num)
        span_list = numpy.arange(0, 1,
                                 step / self.force_analysis.node_number)
        self._subplot_along_time(
            out_filenames=out_filenames,
            variable_list=self.force_analysis.velocity[:, ::step].T[-1:0:-1],
            start_time=start_time,
            end_time=end_time,
            ylabel=velocity_label,
            text_list=[
                r'$z\cdot L^{{-1}}={:.1f}$'.format(span)
                for span in span_list[-2:0:-1]#inverse
            ],
            ymin=velocity_min,
            ymax=velocity_max,
            grid=False,
            figsize=figsize)

    def subplot_fft_amplitude(
            self,
            out_filenames,
            num=9,
            xlabel=r'$f_o\mathrm{\ (Hz)}$',
            ylabel=r'$\mathrm{A_{FFT}}$',
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_LONG_HEIGHT)):
        '''
        Subplots of FFT amplitude.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        num : int
            number of nodes for subplots.
        '''

        if not self.force_analysis.frequency_domain_analysis:
            return

        num += 1
        print('subplot_fft_amplitude')
        step = int(self.force_analysis.node_number // num)
        span_list = numpy.arange(0, 1,
                                 step / self.force_analysis.node_number)
        _subplot_row(
            out_filenames=out_filenames,
            x_variable=self.force_analysis.fft_frequency,
            y_variable_list=self.force_analysis.fft_amplitude[:, ::step].T[1:-1],
            xlabel=xlabel,
            ylabel=ylabel,
            text_list=[
                r'$z\cdot L^{{-1}}={:.1f}$'.format(span)
                for span in span_list[1:-1]
            ],
            xmin=self.force_analysis.frequency_min,
            xmax=self.force_analysis.frequency_max,
            ymin=0,
            ymax=numpy.nanmax(self.force_analysis.fft_amplitude),
            grid=False,
            figsize=figsize, )

    def subplot_power_spectral_density(
            self,
            out_filenames,
            num=9,
            xlabel=r'$f_o\mathrm{\ (Hz)}$',
            ylabel=r'$\mathrm{PSD}$',
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_LONG_HEIGHT)):
        '''
        Subplots of power spectral density.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        num : int
            number of nodes for subplots.
        '''

        if not self.force_analysis.frequency_domain_analysis:
            return

        num += 1
        print('subplot_power_spectral_density')
        step = int(self.force_analysis.node_number // num)
        span_list = numpy.arange(0, 1,
                                 step / self.force_analysis.node_number)
        _subplot_row(
            out_filenames=out_filenames,
            x_variable=self.force_analysis.fft_frequency,
            y_variable_list=self.force_analysis.
            power_spectral_density[:, ::step].T[-2:0:-1],#inverse
            xlabel=xlabel,
            ylabel=ylabel,
            text_list=[
                r'$z\cdot L^{{-1}}={:.1f}$'.format(span)
                for span in span_list[-2:0:-1]#inverse
            ],
            xmin=self.force_analysis.frequency_min,
            xmax=self.force_analysis.frequency_max,
            ymin=0,
            ymax=numpy.nanmax(self.force_analysis.power_spectral_density),
            grid=False,
            figsize=figsize, )

    def subplot_modal_fft_amplitude(
            self,
            out_filenames,
            xlabel=r'$f_o^m\mathrm{\ (Hz)}$',
            ylabel=r'$\mathrm{A_{FFT}}$',
            figsize=(style.FULL_WIDTH * 1 / 3, style.FULL_SHORT_HEIGHT)):
        '''
        Subplots modal fft amplitude for multi modes.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        figsize : float tuple (float1, float2)
            Figure size (width, height).
        '''

        if not self.force_analysis.modal_analysis:
            return
        print('subplot_modal_fft_amplitude')
        _subplot_row(
            out_filenames=out_filenames,
            x_variable=self.force_analysis.fft_frequency,
            y_variable_list=self.force_analysis.modal_fft_amplitude.T,
            xlabel=xlabel,
            ylabel=ylabel,
            text_list=[
                r'Mode {:d}'.format(mode_i)
                for mode_i in
                range(self.force_analysis.mode_number_min,
                      self.force_analysis.mode_number_max + 1)
            ],
            xmin=self.force_analysis.frequency_min,
            xmax=self.force_analysis.frequency_max,
            ymin=0,
            ymax=numpy.nanmax(
                self.force_analysis.modal_fft_amplitude),
            grid=False,
            figsize=figsize)


    def subplot_modal_power_spectral_density(
            self,
            out_filenames,
            xlabel=r'$f_o^m\mathrm{\ (Hz)}$',
            ylabel=r'$\mathrm{PSD}$',
            figsize=(style.FULL_WIDTH * 1 / 3, style.FULL_SHORT_HEIGHT)):
        '''
        Subplots modal power spectral density for multi modes.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        figsize : float tuple (float1, float2)
            Figure size (width, height).
        '''

        if not self.force_analysis.modal_analysis:
            return
        print('subplot_modal_power_spectral_density')
        _subplot_row(
            out_filenames=out_filenames,
            x_variable=self.force_analysis.fft_frequency,
            y_variable_list=self.force_analysis.
            modal_power_spectral_density.T,
            xlabel=xlabel,
            ylabel=ylabel,
            text_list=[
                r'Mode {:d}'.format(mode_i)
                for mode_i in
                range(self.force_analysis.mode_number_min,
                      self.force_analysis.mode_number_max + 1)
            ],
            xmin=self.force_analysis.frequency_min,
            xmax=self.force_analysis.frequency_max,
            ymin=0,
            ymax=numpy.nanmax(
                self.force_analysis.modal_power_spectral_density),
            grid=False,
            figsize=figsize)

    def plot3d_power_spectral_density(self,
                                      out_filenames,
                                      xlabel=r'$f_o\mathrm{\ (Hz)}$',
                                      ylabel=r'$z\cdot L^{-1}$',
                                      figsize=(style.SINGLE_COLUMN_WIDTH / 2,
                                               style.SINGLE_COLUMN_WIDTH / 2)):
        '''
        3d power spectral density plots.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        figsize : float tuple (float1, float2)
            Figure size (width, height).
        '''

        if not self.force_analysis.frequency_domain_analysis:
            return
        print('plot3d_power_spectral_densityl density')
        from mpl_toolkits.mplot3d import Axes3D
        matplotlib.pyplot.clf()
        matplotlib.pyplot.gcf().set_size_inches(figsize)
        ax = matplotlib.pyplot.gca(projection='3d')

        ax.set_xlabel(xlabel)
        ax.set_xlim3d(self.force_analysis.frequency_min,
                      self.force_analysis.frequency_max)
        ax.set_ylim3d(0, 1)
        ax.set_ylabel(ylabel)
        ax.set_zlim3d(
            numpy.nanmin(self.force_analysis.power_spectral_density),
            numpy.nanmax(self.force_analysis.power_spectral_density))
        ax.set_zlabel(r'Power spectral density', labelpad=30)
        ax.zaxis.set_major_formatter(
            matplotlib.ticker.StrMethodFormatter('{x:.0e}'))

        num = 10
        step = self.force_analysis.node_number // num
        span = numpy.linspace(0, 1, num=num)
        line = matplotlib.collections.LineCollection(
            [
                list(zip(self.force_analysis.fft_frequency, A.T))
                for A in self.force_analysis.power_spectral_density[:, ::
                                                                        step].T
            ],
            edgecolor=matplotlib.colors.ColorConverter().to_rgb(
                style.DARK_COLOR))
        ax.add_collection3d(line, zs=span, zdir='y')
        ax.view_init(elev=45, azim=-120)
        ax.grid()
        matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            matplotlib.pyplot.savefig(out_filename)
        matplotlib.pyplot.close()

    def contourf_fft_amplitude(self, out_filenames,
            contourf_num,
            xlabel=r'$f_o\mathrm{\ (Hz)}$',
            ylabel=r'$z\cdot L^{-1}$',
            ):
        '''
        Contourfs of FFT amplitude for all nodes.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        contourf_num : int
            Line number for the contourf.

        '''
        if not self.force_analysis.frequency_domain_analysis:
            return
        print('contourf_fft_amplitude')
        _contourf(
            out_filenames=out_filenames,
            variable=self.force_analysis.fft_amplitude,
            x=self.force_analysis.fft_frequency,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=ylabel,
            colorbar_min=numpy.nanmin(
                self.force_analysis.fft_amplitude),
            colorbar_max=numpy.nanmax(
                self.force_analysis.fft_amplitude),
            contourf_num=contourf_num, )

    def contourf_power_spectral_density(self, out_filenames,
            contourf_num,
            xlabel=r'$f_o\mathrm{\ (Hz)}$',
            ylabel=r'$z\cdot L^{-1}$',
            ):
        '''
        Contourfs of power spectral density for all nodes.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        contourf_num : int
            Line number for the contourf.

        '''
        if not self.force_analysis.frequency_domain_analysis:
            return
        print('contourf_power_spectral_density')
        _contourf(
            out_filenames=out_filenames,
            variable=self.force_analysis.power_spectral_density,
            x=self.force_analysis.fft_frequency,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=ylabel,
            colorbar_min=numpy.nanmin(
                self.force_analysis.power_spectral_density),
            colorbar_max=numpy.nanmax(
                self.force_analysis.power_spectral_density),
            contourf_num=contourf_num, )

    def contour_fft_amplitude(self, out_filenames,
            contour_num,
            xlabel=r'$f_o\mathrm{\ (Hz)}$',
            ylabel=r'$z\cdot L^{-1}$',
            ):
        '''
        Contours of fft amplitude for all nodes.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        contour_num : int
            Line number for the contourf.

        '''
        if not self.force_analysis.frequency_domain_analysis:
            return
        print('contour_fft_amplitude')
        _contour(
            out_filenames=out_filenames,
            variable=self.force_analysis.fft_amplitude,
            x=self.force_analysis.fft_frequency,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=ylabel,
            colorbar_min=numpy.nanmin(
                self.force_analysis.fft_amplitude),
            colorbar_max=numpy.nanmax(
                self.force_analysis.fft_amplitude),
            contour_num=contour_num)

    def contour_power_spectral_density(self, out_filenames,
            contour_num,
            xlabel=r'$f_o\mathrm{\ (Hz)}$',
            ylabel=r'$z\cdot L^{-1}$',
            ):
        '''
        Contours of power spectral density for all nodes.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        contour_num : int
            Line number for the contourf.

        '''
        if not self.force_analysis.frequency_domain_analysis:
            return
        print('contour_power_spectral_density')
        _contour(
            out_filenames=out_filenames,
            variable=self.force_analysis.power_spectral_density,
            x=self.force_analysis.fft_frequency,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=ylabel,
            colorbar_min=numpy.nanmin(
                self.force_analysis.power_spectral_density),
            colorbar_max=numpy.nanmax(
                self.force_analysis.power_spectral_density),
            contour_num=contour_num)

    def contourf_wavelet(self, node_i, out_filenames, start_time, end_time,
                         contourf_num,
                         xlabel=r'$t\mathrm{\ (s)}$',
                         ylabel=r'$f_o\mathrm{\ (Hz)}$',
                         ):
        '''
        Contourfs of wavelet for specific node.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        contourf_num : int
            Line number for the contourf.

        '''
        if not self.force_analysis.wavelet_analysis:
            return
        print('contourf_wavelet')
        # time_index = (
        #     self.force_analysis.wavelet_time >= start_time) & (
        #         self.force_analysis.wavelet_time <= end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        _contourf(
            out_filenames,
            variable=self.force_analysis.wavelet_power[node_i - 1].T[
                time_index],
            #x=self.force_analysis.wavelet_time[time_index],
            x=self.force_analysis.time[time_index],
            y=self.force_analysis.wavelet_frequency,
            xlabel=xlabel,
            ylabel=ylabel,
            colorbar_min=numpy.nanmin(self.force_analysis.wavelet_power[
                node_i - 1]),
            colorbar_max=numpy.nanmax(self.force_analysis.wavelet_power[
                node_i - 1]),
            contourf_num=contourf_num, )

    def contour_wavelet(self, node_i, out_filenames, start_time, end_time,
                        contour_num,
                        xlabel=r'$t\mathrm{\ (s)}$',
                        ylabel=r'$f_o\mathrm{\ (Hz)}$',
                        ):
        '''
        Contours of wavelet for specific node.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        contour_num : int
            Line number for the contourf.

        '''
        if not self.force_analysis.wavelet_analysis:
            return
        print('contour_wavelet')
        # time_index = (
        #     self.force_analysis.wavelet_time >= start_time) & (
        #         self.force_analysis.wavelet_time <= end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        _contour(
            out_filenames,
            variable=self.force_analysis.wavelet_power[node_i - 1].T[
                time_index],
            #x=self.force_analysis.wavelet_time[time_index],
            x=self.force_analysis.time[time_index],
            y=self.force_analysis.wavelet_frequency,
            xlabel=xlabel,
            ylabel=ylabel,
            colorbar_min=numpy.nanmin(self.force_analysis.wavelet_power[
                node_i - 1]),
            colorbar_max=numpy.nanmax(self.force_analysis.wavelet_power[
                node_i - 1]),
            contour_num=contour_num)

    def plot_wavelet_dominant_frequency(self, out_filenames, node_i,
                                        start_time, end_time,
                                        xlabel=r'$t\mathrm{\ (s)}$',
                                        ylabel=r'$f_o\mathrm{\ (Hz)}$',
                                        ):
        '''
        Dominant frequency along time axis for specific node.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        node_i : int
            Node number of the dominant frequency to be plotted.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.

        '''
        if not self.force_analysis.wavelet_analysis:
            return
        print('plot_wavelet_dominant_frequency')
        self._plot_along_time(
            out_filenames=out_filenames,
            variable=self.force_analysis.wavelet_dominant_frequencies[
                node_i - 1],
            start_time=start_time,
            end_time=end_time,
            xlabel=xlabel,
            ylabel=ylabel,
            ymin=self.force_analysis.frequency_min,
            ymax=self.force_analysis.frequency_max)

    def contourf_wavelet_dominant_frequency(self, out_filenames, start_time,
                                            end_time, contourf_num,
                                            xlabel=r'$t\mathrm{\ (s)}$',
                                            ):
        '''
        Contourf of dominant frequency along time axis for all nodes.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        contourf_num : int
            Line number for the contourf.

        '''
        if not self.force_analysis.wavelet_analysis:
            return
        print('contourf_wavelet_dominant_frequency')

        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time = self.force_analysis.time[time_index]
        _contourf(
            out_filenames=out_filenames,
            variable=self.force_analysis.
            wavelet_dominant_frequencies[:, time_index].T,
            x=time,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=r'$z\cdot L^{-1}$',
            colorbar_min=self.force_analysis.frequency_min,
            colorbar_max=self.force_analysis.frequency_max,
            contourf_num=contourf_num, )

    def contour_wavelet_dominant_frequency(self, out_filenames, start_time,
                                           end_time, contour_num,
                                           xlabel=r'$t\mathrm{\ (s)}$',
                                           ):
        '''
        Contour of dominant frequency along time axis for all nodes.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        contour_num : int
            Line number for the contourf.

        '''
        if not self.force_analysis.wavelet_analysis:
            return
        print('contour_wavelet_dominant_frequency')
        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.next_time_index(end_time)
        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        time = self.force_analysis.time[time_index]
        _contour(
            out_filenames=out_filenames,
            variable=self.force_analysis.
            wavelet_dominant_frequencies[:, time_index].T,
            x=time,
            y=self.force_analysis.span,
            xlabel=xlabel,
            ylabel=r'$z\cdot L^{-1}$',
            colorbar_min=self.force_analysis.frequency_min,
            colorbar_max=self.force_analysis.frequency_max,
            contour_num=contour_num)

    def subplot_wavelet_dominant_frequency(
            self,
            out_filenames,
            start_time,
            end_time,
            ymin,
            ymax,
            num=9,
            xlabel=r'$t\mathrm{\ (s)}$',
            ylabel=r'$f_o\mathrm{\ (Hz)}$',
            figsize=(style.SINGLE_COLUMN_WIDTH,
                     style.SINGLE_COLUMN_LONG_HEIGHT)):
        '''
        Plot the time history of dominant frequency of force.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        ylabel : string
            The label shown on the y axis for the lowest figure.
        ymin : float
            Lower limit for y axis.
        ymax : float
            Upper limit for y axis.
        num : int
            number of nodes for subplots.
        '''
        if not self.force_analysis.wavelet_analysis:
            return
        print('subplot_wavelet_dominant_frequency')

        num += 1
        print('subplot_power_spectral_density')
        step = int(self.force_analysis.node_number // num)
        span_list = numpy.arange(0, 1,
                                 step / self.force_analysis.node_number)
        self._subplot_along_time(
            out_filenames=out_filenames,
            variable_list=self.force_analysis.
            wavelet_dominant_frequencies[::step, :][-2:0:-1],#inverse
            start_time=start_time,
            end_time=end_time,
            xlabel=xlabel,
            ylabel=ylabel,
            text_list=[
                r'$z\cdot L^{{-1}}={:.1f}$'.format(span)
                for span in span_list[-2:0:-1]#inverse
            ],
            ymin=ymin,
            ymax=ymax,
            grid=False,
            figsize=figsize)

    def _add_and_move_line_along_time(
            self,
            variable_along_span,
            out_filenames,
            start_time,
            end_time,
            xlabel,
            xmin,
            xmax,
            grid,
            data_step,
            fps,
            dpi=100,
            line_number=10,
            color=style.DARK_COLOR,
            figsize=(style.SINGLE_COLUMN_WIDTH / 2,
                     style.SINGLE_COLUMN_LONG_HEIGHT / 2)):
        matplotlib.pyplot.clf()
        figure = matplotlib.pyplot.gcf()
        figure.set_size_inches(figsize)

        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        variable_along_span = variable_along_span[time_index]
        lines = []
        lines.append(
            matplotlib.pyplot.plot(
                variable_along_span[0, :],
                self.force_analysis.span,
                color=color)[0])
        matplotlib.pyplot.xlim(xmin, xmax)
        axis = matplotlib.pyplot.gca()
        axis.locator_params(axis='x', nbins=style.SHORT_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)
        matplotlib.pyplot.grid(grid)

        figure.savefig('')

        x_sci_notaion = axis.xaxis.get_offset_text()
        x_sci_notaion.set_visible(False)
        if x_sci_notaion.get_text():
            xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                          x_sci_notaion.get_text()[1:])
        matplotlib.pyplot.xlabel(xlabel)
        matplotlib.pyplot.ylabel(r'$z\cdot L^{-1}$')

        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.time_index(end_time)

        add_line_interval = int(sum(time_index) // (data_step * line_number))

        # animation function.  This is called sequentially

        def animate(i):
            lines[0].set_xdata(variable_along_span[data_step * i, :])
            if i % add_line_interval == 0:
                lines.append(
                    matplotlib.pyplot.plot(variable_along_span[
                        data_step * i, :], self.force_analysis.span, **
                                           style.LIGHT_LINE_STYLE)[0])
            return lines

        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)
        animator = matplotlib.animation.FuncAnimation(
            figure,
            animate,
            frames=int(sum(time_index) // data_step),
            blit=False)

        matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            animator.save(out_filename, writer='imagemagick', fps=fps, dpi=dpi)

    def _move_line_along_time(self,
                              variable_along_span,
                              out_filenames,
                              start_time,
                              end_time,
                              xlabel,
                              xmin,
                              xmax,
                              grid,
                              data_step,
                              fps,
                              dpi=100,
                              color=style.DARK_COLOR,
                              figsize=(style.SINGLE_COLUMN_WIDTH / 2,
                                       style.SINGLE_COLUMN_LONG_HEIGHT / 2)):
        matplotlib.pyplot.clf()
        figure = matplotlib.pyplot.gcf()
        figure.set_size_inches(figsize)

        time_index = (self.force_analysis.time >= start_time) & (
            self.force_analysis.time <= end_time)

        variable_along_span = variable_along_span[time_index]

        line, = matplotlib.pyplot.plot(
            variable_along_span[0, :],
            self.force_analysis.span,
            color=color)
        matplotlib.pyplot.xlim(xmin, xmax)
        axis = matplotlib.pyplot.gca()
        axis.locator_params(axis='x', nbins=style.SHORT_XTICK_MAX_LENGTH)
        axis.locator_params(axis='y', nbins=style.LONG_YTICK_MAX_LENGTH)
        matplotlib.pyplot.grid(grid)

        figure.savefig('')

        x_sci_notaion = axis.xaxis.get_offset_text()
        x_sci_notaion.set_visible(False)
        if x_sci_notaion.get_text():
            xlabel = "{:s} / {:s}".format(xlabel[:-1],
                                          x_sci_notaion.get_text()[1:])
        matplotlib.pyplot.xlabel(xlabel)
        matplotlib.pyplot.ylabel(r'$z\cdot L^{-1}$')

        # start_index = self.force_analysis.time_index(start_time)
        # end_index = self.force_analysis.time_index(end_time)
        # initialization function: plot the background of each frame
        def init():
            line.set_data(variable_along_span[0, :],
                          self.force_analysis.span)
            return line,

        # animation function.  This is called sequentially
        def animate(i):
            line.set_xdata(variable_along_span[data_step * i, :])
            return line,

        # call the animator.  blit=True means only re-draw the parts that have
        # changed.
        animator = matplotlib.animation.FuncAnimation(
            figure,
            animate,
            init_func=init,
            frames=int(sum(time_index) // data_step),
            blit=True)

        matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            animator.save(out_filename, writer='imagemagick', fps=fps, dpi=dpi)

    def make_force_animation_along_time(self,
                                               out_filenames,
                                               start_time,
                                               end_time,
                                               xlabel,
                                               xmin,
                                               xmax,
                                               grid=True,
                                               data_step=50,
                                               fps=24):
        '''
        Make force animation.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.
        grid : bool
            True if grid is necessary.
        data_step : int
            Plot every data step
        fps : int
            Frame per second

        '''
        print('make_force_animation_along_time')
        self._move_line_along_time(self.force_analysis.force,
                                   out_filenames, start_time, end_time, xlabel,
                                   xmin, xmax, grid, data_step, fps)

    def make_curvature_animation_along_time(self,
                                            out_filenames,
                                            start_time,
                                            end_time,
                                            xlabel,
                                            xmin,
                                            xmax,
                                            grid=True,
                                            data_step=50,
                                            fps=24):
        '''
        Make curvature animation.

        Parameters
        ----------
        out_filenames : a list of string
            The filenames for saving figures.
        start_time : float
            Lower limit for time axis.
        end_time : float
            Upper limit for time axis.
        xlabel : string
            The label shown on the x axis.
        xmin : float
            Lower limit for x axis.
        xmax : float
            Upper limit for x axis.
        grid : bool
            True if grid is necessary.
        data_step : int
            Plot every data step
        fps : int
            Frame per second

        '''
        print('make_curvature_animation_along_time')
        self._move_line_along_time(-self.force_analysis.curvature,
                                   out_filenames, start_time, end_time, xlabel,
                                   xmin, xmax, grid, data_step, fps)
