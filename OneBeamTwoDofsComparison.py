'''
Designed to plot comparison of several one-cylinder cases with two dofs.
Users just need to instantiate this class,
thus avoiding interactions with more abstarct classes including
ForceData, ForceMode, ForceAnalysis, and ForceVisualization.
'''

import json
import numpy
from .ForceData import ForceData
from .Beam import Beam
from .ForceAnalysis import ForceAnalysis
from .ForceVisualizationComparison import ForceVisualizationComparison


class OneBeamTwoDofsComparison:
    def __init__(
            self,
            # for force_data
            force_filename_list,
            legend_list,
            save_directory,
            filetypes,
            # for beam
            mass_ratio_list,
            top_tension_list,
            horizontal_list,
            tension_include_boyancy_list,
            length_list,
            diameter_list,
            bending_stiffness_list,
            frequency_st_list):
        '''
        Visualize comparison of analytical results for multi cases.

        Parameters
        ----------
        force_filename_list : ForceAnalysis instance
            All analytical results of force data prepared for visualization.

        '''
        self._save_directory = save_directory
        if self._save_directory[-1] != '/':
            self._save_directory += '/'

        self._length = max(
            len(force_filename_list),
            len(legend_list),
            len(mass_ratio_list),
            len(top_tension_list),
            len(horizontal_list),
            len(tension_include_boyancy_list),
            len(length_list), len(diameter_list), len(bending_stiffness_list))

        force_filename_list = self._fill_if_one_element(
            force_filename_list)
        mass_ratio_list = self._fill_if_one_element(mass_ratio_list)
        top_tension_list = self._fill_if_one_element(top_tension_list)
        horizontal_list = self._fill_if_one_element(horizontal_list)
        tension_include_boyancy_list = self._fill_if_one_element(
            tension_include_boyancy_list)
        diameter_list = self._fill_if_one_element(diameter_list)
        length_list = self._fill_if_one_element(length_list)
        bending_stiffness_list = self._fill_if_one_element(
            bending_stiffness_list)
        bending_stiffness_list = self._fill_if_one_element(
            bending_stiffness_list)
        frequency_st_list = self._fill_if_one_element(frequency_st_list)

        self._filetypes = filetypes
        self._legend_list = legend_list
        self._frequency_st_list = frequency_st_list
        self._inline_force_data_list = []
        self._crossflow_force_data_list = []
        self._beam_list = []
        # self._line_styles = line_styles
        # self._MARKER_STYLES = MARKER_STYLES

        assert len(force_filename_list) == len(legend_list) ==\
            len(mass_ratio_list) == len(top_tension_list) ==\
            len(horizontal_list) ==\
            len(tension_include_boyancy_list) == len(length_list) ==\
            len(diameter_list) == len(bending_stiffness_list) == self._length,\
            print(
                len(force_filename_list),
                len(legend_list),
                len(mass_ratio_list),
                len(top_tension_list),
                len(horizontal_list),
                len(tension_include_boyancy_list),
                len(length_list),
                len(diameter_list),
                len(bending_stiffness_list),
                self._length
        )

        for i, (force_filename, mass_ratio, top_tension,
                horizontal,
                tension_include_boyancy, length, diameter,
                bending_stiffness) in enumerate(
                    zip(force_filename_list, mass_ratio_list,
                        top_tension_list, horizontal_list, tension_include_boyancy_list,
                        length_list, diameter_list, bending_stiffness_list)):
            time, inline_force, crossflow_force = numpy.loadtxt(
                force_filename, unpack=True)

            self._inline_force_data_list.append(
                ForceData(
                    time=time, force=inline_force))
            self._crossflow_force_data_list.append(
                ForceData(
                    time=time, force=crossflow_force))
            self._beam_list.append(
                Beam(
                    node_number=self._crossflow_force_data_list[i]
                    .node_number,
                    mass_ratio=mass_ratio,
                    top_tension=top_tension,
                    horizontal=horizontal,
                    tension_include_boyancy=tension_include_boyancy,
                    length=length,
                    diameter=diameter,
                    bending_stiffness=bending_stiffness))

    def __append_filetypes(self, filename):
        return [
            self._save_directory + filename + '.' + filetype
            for filetype in self._filetypes
        ]

    def _fill_if_one_element(self, list_to_be_checked):
        if len(list_to_be_checked) == 1:
            list_to_be_checked = list_to_be_checked * self._length
        return list_to_be_checked

    # x
    def plot_crossflow_force(
            self,
            # for force_analysis
            start_time,
            end_time,
            force_min,
            force_max,
            force_std,
            curvature_min,
            curvature_max,
            curvature_std,
            frequency_min,
            frequency_max,
            mode_number_min,
            mode_number_max, ):
        crossflow_force_analysis_list = []
        for force_data, beam in\
                zip(self._crossflow_force_data_list, self._beam_list):
            crossflow_force_analysis_list.append(
                ForceAnalysis(
                    force_data=force_data,
                    beam=beam,
                    start_time=start_time,
                    end_time=end_time,
                    frequency_min=frequency_min,
                    frequency_max=frequency_max,
                    mode_number_min=mode_number_min,
                    mode_number_max=mode_number_max,
                    modal_analysis=True, ))
        crossflow_force_visulization_comparison = ForceVisualizationComparison(
            crossflow_force_analysis_list, self._legend_list,
            self._frequency_st_list)
        # crossflow_force_visulization_comparison.pre_plot()

        if end_time == -1:
            end_time = min([
                force_analysis.end_time
                for force_analysis in
                crossflow_force_visulization_comparison.
                force_analysis_list
            ])

        out_dict = {}
        force_mean_dict = crossflow_force_visulization_comparison.plot_force_mean(
            self.__append_filetypes(
                'crossflow_force_mean_{:.1f}_{:.1f}'.format(start_time,
                                                                   end_time)),
            xlabel=r'$\bar{y}\cdot D^{-1}$',
            xmin=force_min,
            xmax=force_max, )
        force_std_dict = crossflow_force_visulization_comparison.plot_force_std(
            self.__append_filetypes('crossflow_force_std_{:.1f}_{:.1f}'.
                                    format(start_time, end_time)),
            xlabel=r'$\sigma_y\cdot D^{-1}$',
            xmin=0,
            xmax=force_std, )
        curvature_mean_dict = crossflow_force_visulization_comparison.plot_curvature_mean(
            self.__append_filetypes('crossflow_curvature_mean_{:.1f}_{:.1f}'.
                                    format(start_time, end_time)),
            xlabel=r'$\bar{c}_y\cdot D$',
            xmin=curvature_min,
            xmax=curvature_max, )
        curvature_std_dict = crossflow_force_visulization_comparison.plot_curvature_std(
            self.__append_filetypes('crossflow_curvature_std_{:.1f}_{:.1f}'.
                                    format(start_time, end_time)),
            xlabel=r'$\sigma_{c_y}\cdot D$',
            xmin=0,
            xmax=curvature_std, )
        crossflow_force_visulization_comparison.\
            subplot_force_curvature_std(
                self.__append_filetypes(
                    'crossflow_subplot_force_curvature_std_{:.1f}_{:.1f}'.format(
                        start_time, end_time)
                ),
                xlabel_list=[r'$\sigma_y\cdot D^{-1}$',
                             r'$\sigma_{c_y}\cdot D$'],
                xmin_list=[0, 0],
                xmax_list=[force_std, curvature_std],
            )

        # Save to fjile
        time_domain_dict = {}
        time_domain_dict['force_mean'] = force_mean_dict
        time_domain_dict['curvature mean'] = curvature_mean_dict
        time_domain_dict['force_std'] = force_std_dict
        time_domain_dict['curvature_std'] = curvature_std_dict
        time_domain_dict['time_interval'] = '({:.1f}, {:.1f})'.format(
            start_time, end_time)
        out_dict['time_domain_analysis'] = time_domain_dict

        modal_weight_force_std_dict = crossflow_force_visulization_comparison.plot_modal_weight_force_std(
            self.__append_filetypes('crossflow_modal_weight_force_std'),
            # ylabel=r'$v^m_{\mathrm{rms}}\cdot D^{-1}$'
            ylabel=r'$\sigma_{v^m}\cdot D^{-1}$')

        modal_fo_dict = crossflow_force_visulization_comparison.plot_modal_oscillatory_frequency(
            self.__append_filetypes('crossflow_modal_oscillatory_frequency'),
            ylabel=r'$f^m_{oy}/f_s$',
            diagonal_line_legend=r'$f^m_{oy}=f^m_n$')

        mode_dict = {}
        mode_dict['std'] = modal_weight_force_std_dict
        mode_dict['fo'] = modal_fo_dict
        out_dict['modal_analysis'] = mode_dict
        with open(self._save_directory + 'crossflow.txt', 'w') as outfile:
            json.dump(out_dict, outfile, indent=2)
        print('Crossflow plotting finished.')

    # y
    # x

    def plot_inline_force(
            self,
            # for force_analysis
            start_time,
            end_time,
            force_min,
            force_max,
            force_std,
            curvature_min,
            curvature_max,
            curvature_std,
            frequency_min,
            frequency_max,
            mode_number_min,
            mode_number_max, ):
        inline_force_analysis_list = []
        for force_data, beam in\
                zip(self._inline_force_data_list, self._beam_list):
            inline_force_analysis_list.append(
                ForceAnalysis(
                    force_data=force_data,
                    beam=beam,
                    start_time=start_time,
                    end_time=end_time,
                    frequency_min=frequency_min,
                    frequency_max=frequency_max,
                    mode_number_min=mode_number_min,
                    mode_number_max=mode_number_max,
                    modal_analysis=True, ))
        inline_force_visulization_comparison = ForceVisualizationComparison(
            inline_force_analysis_list, self._legend_list,
            self._frequency_st_list)
        # inline_force_visulization_comparison.pre_plot()

        if end_time == -1:
            end_time = min([
                force_analysis.end_time
                for force_analysis in
                inline_force_visulization_comparison.
                force_analysis_list
            ])

        out_dict = {}

        force_mean_dict = inline_force_visulization_comparison.plot_force_mean(
            self.__append_filetypes('inline_force_mean_{:.1f}_{:.1f}'.
                                    format(start_time, end_time)),
            xlabel=r'$\bar{x}\cdot D^{-1}$',
            xmin=force_min,
            xmax=force_max, )
        force_std_dict = inline_force_visulization_comparison.plot_force_std(
            self.__append_filetypes('inline_force_std_{:.1f}_{:.1f}'.
                                    format(start_time, end_time)),
            xlabel=r'$\sigma_x\cdot D^{-1}$',
            xmin=0,
            xmax=force_std, )
        curvature_mean_dict = inline_force_visulization_comparison.plot_curvature_mean(
            self.__append_filetypes('inline_curvature_mean_{:.1f}_{:.1f}'.
                                    format(start_time, end_time)),
            xlabel=r'$\bar{c}_x\cdot D$',
            xmin=curvature_min,
            xmax=curvature_max, )
        curvature_std_dict = inline_force_visulization_comparison.plot_curvature_std(
            self.__append_filetypes('inline_curvature_std_{:.1f}_{:.1f}'.
                                    format(start_time, end_time)),
            xlabel=r'$\sigma_{c_x}\cdot D$',
            xmin=0,
            xmax=curvature_std, )
        inline_force_visulization_comparison.\
            subplot_force_curvature_mean_std(
                self.__append_filetypes(
                    'inline_subplot_force_curvature_mean_std_{:.1f}_{:.1f}'.format(start_time, end_time)),
                xlabel_list=[r'$\bar{x}\cdot D^{-1}$',
                             r'$\bar{c}_x\cdot D$',
                             r'$\sigma_x\cdot D^{-1}$',
                             r'$\sigma_{c_x}\cdot D$'],
                xmin_list=[force_min,
                           curvature_min,
                           0,
                           0],
                xmax_list=[force_max,
                           curvature_max,
                           force_std,
                           curvature_std],
            )

        # Save to file
        time_domain_dict = {}
        time_domain_dict['force_mean'] = force_mean_dict
        time_domain_dict['curvature_mean'] = curvature_mean_dict
        time_domain_dict['force_std'] = force_std_dict
        time_domain_dict['curvature_std'] = curvature_std_dict
        time_domain_dict['time_interval'] = '({:.1f}, {:.1f})'.format(
            start_time, end_time)
        out_dict['time_domain_analysis'] = time_domain_dict

        modal_weight_force_std_dict = inline_force_visulization_comparison.plot_modal_weight_force_std(
            self.__append_filetypes('inline_modal_weight_force_std'),
            # ylabel=r'$u^m_{\mathrm{rms}}\cdot D^{-1}$'
            ylabel=r'$\sigma_{u^m}\cdot D^{-1}$')

        modal_fo_dict = inline_force_visulization_comparison.plot_modal_oscillatory_frequency(
            self.__append_filetypes('inline_modal_oscillatory_frequency'),
            ylabel=r'$f^m_{ox}/f_s$',
            diagonal_line_legend=r'$f^m_{ox}=f^m_n$')

        mode_dict = {}
        mode_dict['std'] = modal_weight_force_std_dict
        mode_dict['fo'] = modal_fo_dict
        out_dict['modal_analysis'] = mode_dict
        with open(self._save_directory + 'inline.txt', 'w') as outfile:
            json.dump(out_dict, outfile, indent=2)
        print('Inline plotting finished.')
