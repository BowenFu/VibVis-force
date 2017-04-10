'''
Designed to plot data of two cylinders with one dof.
Users just need to instantiate this class,
thus avoiding interactions with more abstarct classes including
ForceData, ForceMode, ForceAnalysis, and ForceVisualization.
'''

import json
import logging
import numpy
from .ForceData import ForceData
from .Beam import Beam
from .ForceAnalysis import ForceAnalysis
from .ForceVisualization import ForceVisualization


class TwoBeamsOneDof:
    '''
    node_i from 1 to node_num
    '''

    def __init__(
            self,
            # for force_data
            force_filename,
            save_directory,
            filetypes,
            # for beam
            mass_ratio,
            top_tension,
            horizontal,
            tension_include_boyancy,
            length,
            diameter,
            bending_stiffness):
        time, first_force, second_force = numpy.loadtxt(
            force_filename, unpack=True)

        self._first_force_data = ForceData(
            time=time, force=first_force)
        self._second_force_data = ForceData(
            time=time, force=second_force)
        self._beam = Beam(
            node_number=self._second_force_data.node_number,
            mass_ratio=mass_ratio,
            top_tension=top_tension,
            horizontal=horizontal,
            tension_include_boyancy=tension_include_boyancy,
            length=length,
            diameter=diameter,
            bending_stiffness=bending_stiffness)
        self._save_directory = save_directory
        if self._save_directory[-1] != '/':
            self._save_directory += '/'

        self._filetypes = filetypes

        self._beam.plot_modal_shapes(
            self._append_filetypes('modal_shapes'), max_order=5)

        with open(self._save_directory + 'natural_frequencies.txt',
                  'wb') as file_:
            numpy.savetxt(
                file_, self._beam.natural_frequencies[:20], fmt='%.4e')

    def _append_filetypes(self, filename):
        return [
            self._save_directory + filename + '.' + filetype
            for filetype in self._filetypes
        ]

    # x
    def plot_first_force(
            self,
            # for force_analysis
            start_time,
            end_time,
            medium_step,
            short_step,
            force_min,
            force_max,
            force_delta,
            force_std,
            curvature_min,
            curvature_max,
            curvature_delta,
            curvature_std,
            frequency_min,
            frequency_max,
            mode_number_min,
            mode_number_max,
            modal_weight_force_delta,
            make_animation=False,
            contour_curvature=False,
            wavelet_analysis=False):
        first_force_analysis =\
                ForceAnalysis(
                        force_data=self._first_force_data,
                        beam=self._beam,
                        start_time=start_time,
                        end_time=end_time,
                        frequency_min=frequency_min,
                        frequency_max=frequency_max,
                        mode_number_min=mode_number_min,
                        mode_number_max=mode_number_max,
                        modal_analysis=True,
                        wavelet_analysis=wavelet_analysis,
                        frequency_domain_analysis=True)
        first_force_visulization = ForceVisualization(
            first_force_analysis)

        start_time = first_force_visulization.force_analysis.start_time
        end_time = first_force_visulization.force_analysis.end_time

        first_force_visulization.subplot_power_spectral_density(
            self._append_filetypes('first_power_spectral_density_subplot'))
        first_force_visulization.contourf_power_spectral_density(
            self._append_filetypes('first_power_spectral_density_contourf'),
            contourf_num=51)
        first_force_visulization.contour_power_spectral_density(
            self._append_filetypes('first_power_spectral_density_contour'),
            contour_num=11)
        first_force_visulization.subplot_modal_power_spectral_density(
            self._append_filetypes(
                'first_modal_power_spectral_density_subplot'))
        first_force_visulization.contourf_wavelet_dominant_frequency(
            self._append_filetypes(
                'first_wavelet_dominant_frequency_contourf'),
            start_time=start_time,
            end_time=end_time,
            contourf_num=51)
        first_force_visulization.contour_wavelet_dominant_frequency(
            self._append_filetypes('first_wavelet_dominant_frequency_contour'),
            start_time=start_time,
            end_time=end_time,
            contour_num=11)

        first_force_visulization.subplot_wavelet_dominant_frequency(
            out_filenames=self._append_filetypes('first_subplot_wavelet'),
            start_time=start_time,
            end_time=end_time,
            ymin=frequency_min,
            ymax=frequency_max,
            num=9)

        # animation
        if make_animation:
            first_force_visulization.make_force_animation_along_time(
                [self._save_directory + 'first_force_animation.gif'],
                start_time=start_time,
                end_time=end_time,
                xlabel=r'$y_L\cdot D^{-1}$',
                xmin=force_min,
                xmax=force_max)

        for node_i in range(
                10,
                first_force_visulization.force_analysis.node_number,
                10):
            first_force_visulization.contourf_wavelet(
                node_i,
                self._append_filetypes(
                    'first_wavelet_contourf_node_{:d}'.format(node_i)),
                start_time=start_time,
                end_time=end_time,
                contourf_num=51)
            first_force_visulization.contour_wavelet(
                node_i,
                self._append_filetypes(
                    'first_wavelet_contour_node_{:d}'.format(node_i)),
                start_time=start_time,
                end_time=end_time,
                contour_num=11)
            first_force_visulization.plot_wavelet_dominant_frequency(
                self._append_filetypes(
                    'first_wavelet_dominant_frequency_node_{:d}'.format(
                        node_i)),
                node_i,
                start_time=start_time,
                end_time=end_time)

        first_force_visulization.subplot_span_force_deviation(
            out_filenames=self._append_filetypes(
                'first_subplot_force_deviation'),
            start_time=start_time,
            end_time=end_time,
            force_deviation_label=r'$(y_L-\bar{y_L})\cdot D^{-1}$',
            force_deviation_min=-force_delta,
            force_deviation_max=force_delta)

        first_force_visulization.subplot_modal_weight_force(
            out_filenames=self._append_filetypes('first_mode'),
            start_time=start_time,
            end_time=end_time,
            ylabel=r'$v_L^m\cdot D^{-1}$',
            ymin=-modal_weight_force_delta,
            ymax=modal_weight_force_delta)

        first_mean = first_force_visulization.plot_force_mean(
            self._append_filetypes('first_force_mean_{:.1f}_{:.1f}'.
                                   format(start_time, end_time)),
            xlabel=r'$y_L\cdot D^{-1}$',
            xmin=force_min,
            xmax=force_max, )
        with open(self._save_directory +
                  'first_force_mean_{:.1f}_{:.1f}.txt'.format(
                      start_time, end_time), 'w') as outfile:
            json.dump(first_mean, outfile)
            first_curvature_mean = first_force_visulization.plot_curvature_mean(
                self._append_filetypes('first_curvature_mean_{:.1f}_{:.1f}'.
                                       format(start_time, end_time)),
                xlabel=r'$y_L\cdot D^{-1}$',
                xmin=curvature_min,
                xmax=curvature_max, )
            with open(self._save_directory +
                      'first_curvature_mean_{:.1f}_{:.1f}.txt'.format(
                          start_time, end_time), 'w') as outfile:
                json.dump(first_curvature_mean, outfile)
            first_std = first_force_visulization.plot_force_std(
                self._append_filetypes('first_force_std_{:.1f}_{:.1f}'.
                                       format(start_time, end_time)),
                xlabel=r'$\sigma_{y_L}\cdot D^{-1}$',
                xmin=0,
                xmax=force_std, )
            numpy.savetxt(self._save_directory +
                          'first_force_std_{:.1f}_{:.1f}'.format(
                              start_time, end_time), first_std)
            first_std_curvature = first_force_visulization.plot_curvature_std(
                self._append_filetypes('first_curvature_std_{:.1f}_{:.1f}'.
                                       format(start_time, end_time)),
                xlabel=r'$\sigma_{c_{y_L}}\cdot D$',
                xmin=0,
                xmax=curvature_std, )
            numpy.savetxt(self._save_directory +
                          'first_curvature_std_{:.1f}_{:.1f}'.format(
                              start_time, end_time), first_std_curvature)

            first_force_visulization.subplot_force_curvature_mean_std(
                self._append_filetypes(
                    'first_mean_std_force_curvature_{:.1f}_{:.1f}'.
                    format(start_time, end_time)),
                xlabel_list=[
                    r'$\bar{y}_L\cdot D^{-1}$', r'$\bar{c}_{y_L}\cdot D$',
                    r'$\sigma_{y_L}\cdot D^{-1}$', r'$\sigma_{c_{y_L}}\cdot D$'
                ],
                xmin_list=[force_min, curvature_min, 0, 0],
                xmax_list=[
                    force_max, curvature_max, force_std,
                    curvature_std
                ], )

        end_time += -medium_step

        for time in numpy.arange(start_time, end_time, medium_step):
            first_force_visulization.contourf_contour_spatio_temporal_force(
                self._append_filetypes(
                    'first_spatio_temporal_force_{:.1f}'.format(time)),
                start_time=time,
                end_time=(time + medium_step),
                colorbar_min=-force_delta,
                colorbar_max=force_delta,
                contourf_num=51,
                contour_num=11,
                colorbar_zlabel=r'$y_L\cdot D^{-1}$')
            first_force_visulization.contour_spatio_temporal_force(
                self._append_filetypes(
                    'first_spatio_temporal_force_contour_{:.1f}'.format(
                        time)),
                start_time=time,
                end_time=(time + medium_step),
                colorbar_min=-force_delta,
                colorbar_max=force_delta,
                contour_num=11)

        if contour_curvature:
            for time in numpy.arange(start_time, end_time, medium_step):
                first_force_visulization.contourf_contour_spatio_temporal_curvature(
                    self._append_filetypes(
                        'first_spatio_temporal_curvature_{:.1f}'.format(time)),
                    start_time=time,
                    end_time=(time + medium_step),
                    colorbar_min=-curvature_delta,
                    colorbar_max=curvature_delta,
                    contourf_num=51,
                    contour_num=11,
                    colorbar_zlabel=r'$c_{y_L}\cdot D^{-1}$')
                first_force_visulization.contour_spatio_temporal_curvature(
                    self._append_filetypes(
                        'first_spatio_temporal_curvature_contour_{:.1f}'.
                        format(time)),
                    start_time=time,
                    end_time=(time + medium_step),
                    colorbar_min=-curvature_delta,
                    colorbar_max=curvature_delta,
                    contour_num=11)

        end_time += (medium_step - short_step)

        for time in numpy.arange(start_time, end_time, short_step):
            first_force_visulization.plot_outline(
                self._append_filetypes('first_outline_{:.1f}'.format(time)),
                xlabel=r'$y_L\cdot D^{-1}$',
                xmin=force_min,
                xmax=force_max,
                start_time=time,
                end_time=(time + short_step),
                line_number=25)
            first_force_visulization.plot_deviation_outline(
                self._append_filetypes('first_deviation_outline_{:.1f}'.format(
                    time)),
                xlabel=r'$(y_L-\bar{y}_L)\cdot D^{-1}$',
                xmin=-force_delta,
                xmax=force_delta,
                start_time=time,
                end_time=(time + short_step),
                line_number=25)

        with open(self._save_directory + 'first_oscillatory_frequencies.txt',
                  'wb') as file_:
            numpy.savetxt(
                file_,
                first_force_analysis.oscillatory_frequencies,
                fmt='%.4e')
        print('First plotting completed.')

    # y
    def plot_second_force(
            self,
            # for force_analysis
            start_time,
            end_time,
            medium_step,
            short_step,
            force_min,
            force_max,
            force_delta,
            force_std,
            curvature_min,
            curvature_max,
            curvature_delta,
            curvature_std,
            frequency_min,
            frequency_max,
            mode_number_min,
            mode_number_max,
            modal_weight_force_delta,
            make_animation=False,
            contour_curvature=False):
        second_force_analysis = ForceAnalysis(
            force_data=self._second_force_data,
            beam=self._beam,
            start_time=start_time,
            end_time=end_time,
            frequency_min=frequency_min,
            frequency_max=frequency_max,
            mode_number_min=mode_number_min,
            mode_number_max=mode_number_max,
            modal_analysis=True,
            frequency_domain_analysis=True)
        second_force_visulization = ForceVisualization(
            second_force_analysis)

        start_time = second_force_visulization.force_analysis.start_time
        end_time = second_force_visulization.force_analysis.end_time
        '''
        second_force_visulization.plot3d_power_spectral_density(
            self._append_filetypes('second_power_spectral_density')
        )
        '''

        second_force_visulization.subplot_power_spectral_density(
            self._append_filetypes('second_power_spectral_density_subplot'))
        second_force_visulization.subplot_modal_power_spectral_density(
            self._append_filetypes(
                'second_modal_power_spectral_density_subplot'))
        second_force_visulization.contourf_power_spectral_density(
            self._append_filetypes('second_power_spectral_density_contourf'),
            contourf_num=51)
        second_force_visulization.contour_power_spectral_density(
            self._append_filetypes('second_power_spectral_density_contour'),
            contour_num=11)
        second_force_visulization.contourf_wavelet_dominant_frequency(
            self._append_filetypes(
                'second_wavelet_dominant_frequency_contourf'),
            start_time=start_time,
            end_time=end_time,
            contourf_num=51)
        second_force_visulization.contour_wavelet_dominant_frequency(
            self._append_filetypes(
                'second_wavelet_dominant_frequency_contour'),
            start_time=start_time,
            end_time=end_time,
            contour_num=11)

        # animation
        if make_animation:
            second_force_visulization.make_force_animation_along_time(
                [self._save_directory + 'second_force_animation.gif'],
                start_time=start_time,
                end_time=end_time,
                xlabel=r'$y_R\cdot D^{-1}$',
                xmin=force_min,
                xmax=force_max)
        '''
        for node_i in range(1, second_force_visulization.force_analysis.node_number, 10):
            second_force_visulization.plot_time_history(
                self._append_filetypes('second_node_{:d}'.format(node_i)),
                node_i,
                ylabel=r'$y_R\cdot D^{-1}$',
                ymin=force_min,
                ymax=force_max,
                start_time=start_time,
                end_time=end_time
            )
        '''

        for node_i in range(
                10,
                second_force_visulization.force_analysis.node_number,
                10):
            second_force_visulization.contourf_wavelet(
                node_i,
                self._append_filetypes(
                    'second_wavelet_contourf_node_{:d}'.format(node_i)),
                start_time=start_time,
                end_time=end_time,
                contourf_num=51)
            second_force_visulization.contour_wavelet(
                node_i,
                self._append_filetypes(
                    'second_wavelet_contour_node_{:d}'.format(node_i)),
                start_time=start_time,
                end_time=end_time,
                contour_num=11)
            second_force_visulization.plot_wavelet_dominant_frequency(
                self._append_filetypes(
                    'second_wavelet_dominant_frequency_node_{:d}'.format(
                        node_i)),
                node_i,
                start_time=start_time,
                end_time=end_time)

        second_force_visulization.subplot_span_force(
            out_filenames=self._append_filetypes(
                'second_subplot_force'),
            start_time=start_time,
            end_time=end_time,
            force_label=r'$y_R\cdot D^{-1}$',
            force_min=force_min,
            force_max=force_max)

        second_force_visulization.subplot_modal_weight_force(
            out_filenames=self._append_filetypes('second_mode'),
            start_time=start_time,
            end_time=end_time,
            ylabel=r'$v_R^m\cdot D^{-1}$',
            ymin=-modal_weight_force_delta,
            ymax=modal_weight_force_delta)

        second_force_visulization.plot_force_mean(
            self._append_filetypes('second_force_mean_{:.1f}_{:.1f}'.
                                   format(start_time, end_time)),
            xlabel=r'$y_R\cdot D^{-1}$',
            xmin=force_min,
            xmax=force_max, )
        second_force_visulization.plot_curvature_mean(
            self._append_filetypes('second_curvature_mean_{:.1f}_{:.1f}'.
                                   format(start_time, end_time)),
            xlabel=r'$y_R\cdot D^{-1}$',
            xmin=curvature_min,
            xmax=curvature_max, )
        second_force_visulization.plot_force_std(
            self._append_filetypes('second_force_std_{:.1f}_{:.1f}'.
                                   format(start_time, end_time)),
            xlabel=r'$\sigma_{y_R}\cdot D^{-1}$',
            xmin=0,
            xmax=force_std, )
        second_force_visulization.plot_curvature_std(
            self._append_filetypes('second_curvature_std_{:.1f}_{:.1f}'.format(
                start_time, end_time)),
            xlabel=r'$\sigma_{c_{y_R}}\cdot D$',
            xmin=0,
            xmax=curvature_std, )
        second_force_visulization.subplot_force_curvature_std(
            self._append_filetypes(
                'second_std_force_curvature_{:.1f}_{:.1f}'.format(
                    start_time, end_time)),
            xlabel_list=[
                r'$\sigma_{y_R}\cdot D^{-1}$', r'$\sigma_{c_{y_R}}\cdot D$'
            ],
            xmin_list=[0, 0],
            xmax_list=[force_std, curvature_std], )

        end_time += -medium_step

        for time in numpy.arange(start_time, end_time, medium_step):
            second_force_visulization.contourf_contour_spatio_temporal_force(
                self._append_filetypes(
                    'second_spatio_temporal_force_{:.1f}'.format(time)),
                start_time=time,
                end_time=(time + medium_step),
                colorbar_min=-force_delta,
                colorbar_max=force_delta,
                contourf_num=51,
                contour_num=7,
                colorbar_zlabel=r'$y_R\cdot D^{-1}$')
            second_force_visulization.contour_spatio_temporal_force(
                self._append_filetypes(
                    'second_spatio_temporal_force_contour_{:.1f}'.
                    format(time)),
                start_time=time,
                end_time=(time + medium_step),
                colorbar_min=-force_delta,
                colorbar_max=force_delta,
                contour_num=7)

        if contour_curvature:
            for time in numpy.arange(start_time, end_time, medium_step):
                second_force_visulization.contourf_contour_spatio_temporal_curvature(
                    self._append_filetypes(
                        'second_spatio_temporal_curvature_{:.1f}'.format(
                            time)),
                    start_time=time,
                    end_time=(time + medium_step),
                    colorbar_min=-curvature_delta,
                    colorbar_max=curvature_delta,
                    contourf_num=51,
                    contour_num=11,
                    colorbar_zlabel=r'$c_{y_R}\cdot D^{-1}$')
                second_force_visulization.contour_spatio_temporal_curvature(
                    self._append_filetypes(
                        'second_spatio_temporal_curvature_contour_{:.1f}'.
                        format(time)),
                    start_time=time,
                    end_time=(time + medium_step),
                    colorbar_min=-curvature_delta,
                    colorbar_max=curvature_delta,
                    contour_num=11)

        end_time += medium_step - short_step

        for time in numpy.arange(start_time, end_time, short_step):
            second_force_visulization.plot_outline(
                self._append_filetypes('second_outline_{:.1f}'.format(time)),
                xlabel=r'$y_R\cdot D^{-1}$',
                xmin=force_min,
                xmax=force_max,
                start_time=time,
                end_time=(time + short_step),
                line_number=25)
            second_force_visulization.plot_deviation_outline(
                self._append_filetypes(
                    'second_deviation_outline_{:.1f}'.format(time)),
                xlabel=r'$(y_R-\bar{y}_R)\cdot D^{-1}$',
                xmin=-force_delta,
                xmax=force_delta,
                start_time=time,
                end_time=(time + short_step),
                line_number=25)

        with open(self._save_directory + 'second_oscillatory_frequencies.txt',
                  'wb') as file_:
            numpy.savetxt(
                file_,
                second_force_analysis.oscillatory_frequencies,
                fmt='%.4e')
        print('Second plotting completed.')
