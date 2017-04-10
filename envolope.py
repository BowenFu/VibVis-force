'''
Designed to plot data of one cylinder with two dofs.
Users just need to instantiate this class,
thus avoiding interactions with more abstarct classes including
ForceData, ForceMode, ForceAnalysis, and ForceVisualization.
'''
import json
import numpy
import os.path
from .ForceData import ForceData
from .Beam import Beam
from .ForceAnalysis import ForceAnalysis
from .ForceVisualization import ForceVisualization
from . import style


class OneBeamTwoDofs:
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
            tension_include_boyancy,
            length,
            diameter,
            bending_stiffness,
            ):
        time, inline_force, crossflow_force =\
            numpy.loadtxt(force_filename, unpack=True)

        self._inline_force_data = ForceData(
            time=time, force=inline_force)
        self._crossflow_force_data = ForceData(
            time=time, force=crossflow_force)
        self._beam = Beam(
            node_number=self._crossflow_force_data.node_number,
            mass_ratio=mass_ratio,
            top_tension=top_tension,
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

    # y
    def plot_crossflow_force(
            self,
            # for force_analysis
            start_time,
            end_time,
            # ,
            medium_step,
            short_step,
            force_min,
            force_max,
            force_delta,
            force_std,
            velocity_min,
            velocity_max,
            curvature_min,
            curvature_max,
            curvature_delta,
            curvature_std,
            frequency_min,
            frequency_max,
            mode_number_min,
            mode_number_max,
            modal_weight_force_delta,
            reduced_velocity_max,
            period,
            make_animation=False,
            contour_curvature=False,
            wavelet_analysis=False):

        crossflow_force_analysis = ForceAnalysis(
            force_data=self._crossflow_force_data,
            beam=self._beam,
            start_time=start_time,
            end_time=end_time,
            frequency_min=frequency_min,
            frequency_max=frequency_max,
            mode_number_min=mode_number_min,
            mode_number_max=mode_number_max,
            modal_analysis=False,
            frequency_domain_analysis=False,
            wavelet_analysis=wavelet_analysis, )

        crossflow_force_visulization = ForceVisualization(
            crossflow_force_analysis)

        start_time = crossflow_force_visulization.force_analysis.start_time
        end_time = crossflow_force_visulization.force_analysis.end_time

        crossflow_force_visulization.plot_force_along_time_function(
            out_filenames=self._append_filetypes(
                'crossflow_force_reduced_velocity_node_40'),
            node_i=40,
            start_time=start_time,
            end_time=end_time,
            time_function=lambda time: reduced_velocity_max*numpy.abs(numpy.cos(2*numpy.pi*time/period)),
            xlabel=r'$U_r$',
            force_label=r'$A_y\cdot D^{-1}$',
            force_min=0, # set as 0
            force_max=force_max,
            #reference_line_factor_tuple=(0.7071, 0.7071),
            grid=True,
            # figsize=(style.SINGLE_COLUMN_WIDTH,
            #          style.SINGLE_COLUMN_SHORT_HEIGHT / 1.5)
            )

        print('Crossflow plotting finished.')
