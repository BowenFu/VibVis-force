'''
Designed to process data.
'''

from __future__ import division, print_function

import logging
import numpy


def nextpow2(number):
    '''
    Return the smallest power more than the input integal number.

    Parameters
    ----------
    number : int
        The integal number for extending to power of 2.

    Returns
    -------
    int
        The next power of two more than the input integal number.

    >>> print(nextpow2(10))
    16

    '''
    index_float = numpy.log2(number)
    index_int = int(numpy.ceil(index_float))
    return 2**index_int


class ForceAnalysis:
    '''
    Do force analysis.
    Provide all kinds of analyses for vibraton visualization.

    Parameters
    ----------
    force_data : ForceData instance
        Dimensionalized force data.
    beam : Beam instance
        Beam property.
    start_time : float
        Start of concerned time range.
    end_time : float
        End of concerned time range.
    frequency_min : float
        Lower limit of concerned frequency range.
    frequency_max : float
        Upper limit of concerned frequency range.
    mode_number_min : int
        Lower limit of concerned mode range.
    mode_number_max : int
        Upper limit of concerned mode range.
    frequency_domain_analysis : bool
        True if frequency domain analysis is required.
    wavelet_analysis : bool
        True if time frequency analysis is required.
    modal_domain_analysis : bool
        True if modal analysis is required.

    Input: force_data, beam
    Output: modal analysis and frequency analysis
    Both mode_number_min and mode_number_max is included in the modal analysis

    '''

    def __init__(self,
                 force_data,
                 beam,
                 fluid_density,
                 fluid_velocity,
                 start_time,
                 end_time,
                 frequency_min,
                 frequency_max,
                 mode_number_min,
                 mode_number_max,
                 frequency_domain_analysis=False,
                 wavelet_analysis=False,
                 modal_analysis=False):
        self._force_data = force_data
        self._beam = beam

        self._fluid_density = fluid_density
        self._fluid_velocity = fluid_velocity

        self._start_time = start_time
        self._end_time = end_time

        self._frequency_min = frequency_min
        self._frequency_max = frequency_max

        self._mode_number_min_minus_one = mode_number_min - 1
        self._mode_number_max = mode_number_max

        self._preprocess_force_data()
        self._do_force_analysis(frequency_domain_analysis,
                                    wavelet_analysis, modal_analysis)

    def _do_force_analysis(self, frequency_domain_analysis,
                               wavelet_analysis, modal_analysis):
        self._frequency_domain_analysis = frequency_domain_analysis
        self._wavelet_analysis = wavelet_analysis
        self._modal_analysis = modal_analysis

        self._do_time_domain_analysis()

        if frequency_domain_analysis:
            self._do_frequency_domain_analysis()
        if wavelet_analysis:
            self._do_wavelet_analysis()
        if modal_analysis:
            self._do_modal_analysis()

    def _preprocess_force_data(self):
        if self.start_time < self._force_data.time[0]:
            self._start_time = self._force_data.time[0]
            logging.warning('''
                    Start time has been set as {start_time}.
                    '''.format(start_time=self.start_time))

        if self.end_time == -1 or self.end_time > self._force_data.time[
                -1]:
            self._end_time = self._force_data.time[-1]
            logging.warning('''
                    End time has been set as {end_time}.
                    '''.format(end_time=self.end_time))

        time_mask = ((self._force_data.time >= self._start_time) &
                     (self._force_data.time <= self._end_time))
        self._time = self._force_data.time[time_mask]
        self._force = self._force_data.force[
            time_mask] / (0.5 * self._fluid_density * self._fluid_velocity**2 * self._beam.diameter)

        # calculate angle and curvature
        angle = numpy.gradient(self.force, self.sampling_period,
                               self._beam.node_spacing /
                               self._beam.diameter)[1]
        # divide this by diameter to make it non-dimenional

        self._curvature = numpy.gradient(angle, self.sampling_period,
                                         self._beam.node_spacing /
                                         self._beam.diameter)[1]
        # divide this by diameter to make it non-dimenional

        # self._start_time = self._time[0]
        # self._end_time = self._time[-1]
        # start = self.time_index(self.start_time)
        # end = self.next_time_index(self.end_time)
        # self._time = self._time[start:end]

        self._fft_length = nextpow2(self._time.size)
        # for frequency domain analysis
        self._fft_frequency = numpy.fft.fftfreq(
            self._fft_length,
            self.sampling_period)[:int(self._fft_length // 2)]

        # desirable frequency interval
        self._frequency_mask = numpy.multiply(
            self.fft_frequency - self.frequency_min,
            self.fft_frequency - self.frequency_max) <= 0

        # we unify modal frequency limit to frequency limit,
        # so crop fft frequency directly
        self._fft_frequency = self.fft_frequency[self._frequency_mask]

    def _do_time_domain_analysis(self):
        '''
        Calculate during the whole period.

        '''

        self._force_mean = numpy.mean(self.force, axis=0)
        self._force_min = numpy.min(self.force, axis=0)
        self._force_max = numpy.max(self.force, axis=0)
        self._force_deviation = self.force - self.force_mean
        self._force_std = numpy.std(self.force, axis=0)

        self._curvature_mean = numpy.mean(self.curvature, axis=0)
        self._curvature_min = numpy.min(self.curvature, axis=0)
        self._curvature_max = numpy.max(self.curvature, axis=0)
        self._curvature_deviation = self.curvature - self.curvature_mean
        self._curvature_std = numpy.std(self.curvature, axis=0)

        # calculate velocity and accelaration
        self._velocity = numpy.gradient(
            self.force,
            self.sampling_period,
            self._beam.node_spacing,
            edge_order=2)[0]
        self._accelaration = numpy.gradient(
            self._velocity,
            self.sampling_period,
            self._beam.node_spacing,
            edge_order=2)[0]

    def _do_frequency_domain_analysis(self):
        '''
        Do frequency domain analysis.
        Get oscillatory frequencies for every node.
        Get power spectra for every node.

        '''
        self._calculate_oscillatory_frequencies()

    def _do_wavelet_analysis(self):
        '''
        Do time-frequency analysis.

        '''

        self._calculate_wavelet()

    def _do_modal_analysis(self):
        '''
        Do modal analysis.

        '''
        self._calculate_modal_weight_forces()
        self._calculate_modal_oscillatory_frequencies()

    def _calculate_oscillatory_frequencies(self):
        '''
        calculate power spectra of the force along the span

        '''
        # no longer use the fft_length-piont fft
        # use force rather than deviation
        # fft = numpy.fft.fft( self.force, axis=0)
        fft = numpy.fft.fft(self.force, n=self._fft_length, axis=0)
        # only plotting the first half since it's mirrored
        # change weight to power spectra
        self._power_spectral_density = numpy.square(
            numpy.abs(fft[0:(self._fft_length // 2)])) / (self._time.size // 2)
        self._power_spectral_density = self._power_spectral_density[
            self._frequency_mask]

        power_spectral_density_max_index = numpy.argmax(
            self._power_spectral_density, axis=0)

        self._oscillatory_frequencies = self.fft_frequency[
            power_spectral_density_max_index]

    def _calculate_wavelet(self):
        '''
        Get a list of WaveletAnalyses, one for every node

        '''

        from wavelets import WaveletAnalysis
        self._wavelet = [
            WaveletAnalysis(
                force, time=self._time, dt=self.sampling_period)
            for force in self.force.T
        ]
        frequencies = self.wavelet[0].fourier_frequencies
        f_mask = numpy.multiply(frequencies - self.frequency_min,
                                frequencies - self.frequency_max) <= 0
        self._wavelet_frequency = frequencies[f_mask]
        wavelet_dominant_frequencies = []
        self._wavelet_power = []
        # for i in range(len(self.wavelet)):
        for wavelet_ in self.wavelet:
            power = wavelet_.wavelet_power
            power = power[f_mask]
            wavelet_power_spectral_density_max_index = numpy.argmax(
                power, axis=0)
            self._wavelet_power.append(power)
            wavelet_dominant_frequencies.append(self._wavelet_frequency[
                wavelet_power_spectral_density_max_index])
        # _wavelet_power[node, time, frequency]
        # self._wavelet_power=numpy.stack(wavelet_power)
        # wavelet_dominant_frequencies[node, time]
        self._wavelet_dominant_frequencies = numpy.stack(
            wavelet_dominant_frequencies)

    def _calculate_modal_weight_forces(self):
        # here we use the least-squares solution
        self._modal_weight_force = numpy.linalg.lstsq(
            self._beam.shapes[1:-1, self._mode_number_min_minus_one:
                              self._mode_number_max],
            self.force[:, 1:-1].T)[0].T

    def _calculate_modal_oscillatory_frequencies(self):
        '''
        Calculate modal frequency of the force along the span

        '''
        modal_fft = numpy.fft.fft(self._modal_weight_force,
                                  n=self._fft_length,
                                  axis=0) / self._time.size
        self._modal_power_spectral_density = numpy.square(
            numpy.abs(modal_fft[0:(self._fft_length // 2)])) / (
                self._time.size // 2)

        # self._modal_power_spectral_density[f_mask] = numpy.nan
        self._modal_power_spectral_density =\
            self._modal_power_spectral_density[self._frequency_mask]
        modal_power_spectral_density_max_index = numpy.argmax(
            self._modal_power_spectral_density, axis=0)
        self._modal_oscillatory_frequencies = self._fft_frequency[
            modal_power_spectral_density_max_index]
        self._modal_weight_force_std = numpy.std(
            self._modal_weight_force, axis=0)

    """
    def time_index(self, time):
        '''
        Return the index of specific time in time sequence.

        Parameters
        ----------
        time : float
            The specific time to be fetched in the time sequence.

        Returns
        -------
        int
            The index of specific time in time sequence.

        '''

        if time == -1:
            return self.time_index(self._end_time)
        if time == 0:
            return 0
        elif time < 0:
            raise ValueError('time cannot be negative.')
        elif time < self._start_time or time > self._end_time:
            raise ValueError(
                'time {:.4f} exceed the time range ({:.4f}:{:4f}).'.format(
                    time, self._start_time, self._end_time
                ))
        else:
            # return
            # int(numpy.around((time-self._time[0])/self.sampling_period))
            return int(
                numpy.around(
                    (time - self._start_time) /
                    self.sampling_period
                ))
            # return
            # int((time-self._time[0])/self.sampling_period)

    def next_time_index(self, time):
        '''
        Return the next index of specific time in time sequence for the use of end index.

        Parameters
        ----------
        time : float
            The specific time to be fetched in the time sequence.

        Returns
        -------
        int
            The next index of specific time in time sequence.

            None for end time, time index + 1 for others.

        '''
        if time == -1 or time == self._time[-1]:
            return None
        else:
            return self.time_index(time) + 1

    """

    @property
    def node_number(self):
        '''
        The node number of the beam.
        '''
        return self._beam.node_number

    @property
    def span(self):
        '''
        The span of the beam.
        '''
        return numpy.linspace(0, 1, num=self.node_number)

    @property
    def time(self):
        '''
        The time sequence of the force.

        '''
        return self._time

    @property
    def start_time(self):
        '''
        The start time of analysis.

        '''
        return self._start_time

    @property
    def end_time(self):
        '''
        The end time of analysis.

        '''
        return self._end_time

    @property
    def force(self):
        '''
        Rank-2 nondimensionalized force [time, node number].
        '''

        return self._force

    @property
    def velocity(self):
        '''
        Rank-2 nondimensionalized velocity [time, node number].
        '''

        return self._velocity

    @property
    def accelaration(self):
        '''
        Rank-2 nondimensionalized accelaration [time, node number].
        '''

        return self._accelaration

    @property
    def force_min(self):
        '''
        Rank-1 nondimensionalized force min value for each node.
        '''
        return self._force_min

    @property
    def force_max(self):
        '''
        Rank-1 nondimensionalized force max value for each node.
        '''
        return self._force_max

    @property
    def force_mean(self):
        '''
        Rank-1 nondimensionalized force mean value for each node.
        '''
        return self._force_mean

    @property
    def force_deviation(self):
        '''
        Rank-2 nondimensionalized force deviation value [time, node number].
        '''
        return self._force_deviation

    @property
    def force_std(self):
        '''
        Rank-1 nondimensionalized force standard deviation value for each node.
        '''
        return self._force_std

    @property
    def curvature(self):
        '''
        Rank-2 nondimensionalized curvature [time, node number].
        '''
        return self._curvature

    @property
    def curvature_min(self):
        '''
        Rank-1 nondimensionalized curvature min value for each node.
        '''
        return self._curvature_min

    @property
    def curvature_max(self):
        '''
        Rank-1 nondimensionalized curvature max value for each node.
        '''
        return self._curvature_max

    @property
    def curvature_mean(self):
        '''
        Rank-1 nondimensionalized curvature mean value for each node.
        '''
        return self._curvature_mean

    @property
    def curvature_deviation(self):
        '''
        Rank-2 nondimensionalized curvature deviation [time, node number].
        '''
        return self._curvature_deviation

    @property
    def curvature_std(self):
        '''
        Rank-1 nondimensionalized curvature standard deviation value for each node.
        '''
        return self._curvature_std

    @property
    def frequency_min(self):
        '''
        Lower limit of frequency analysis.
        '''
        return self._frequency_min

    @property
    def frequency_max(self):
        '''
        Upper limit of frequency analysis.
        '''
        return self._frequency_max

    @property
    def fft_frequency(self):
        '''
        Frequency used for frequency analysis.
        '''
        return self._fft_frequency

    @property
    def wavelet(self):
        '''
        Wavelet.
        '''
        return self._wavelet

    @property
    def wavelet_power(self):
        '''
        Wavelet power.
        '''
        return self._wavelet_power

    @property
    def wavelet_frequency(self):
        '''
        Frequency range used for wavelet.
        '''
        return self._wavelet_frequency

    @property
    def wavelet_dominant_frequencies(self):
        '''
        Wavelet dominant frequency varying along the time.
        '''
        return self._wavelet_dominant_frequencies

    @property
    def power_spectral_density(self):
        '''
        Power spectral density for each node [frequency, node number].
        '''
        return self._power_spectral_density

    @property
    def modal_power_spectral_density(self):
        '''
        modal power spectral density for each mode [frequency, mode number]
        '''

        return self._modal_power_spectral_density

    @property
    def mode_number_min(self):
        '''
        Lower limit for modal analysis.
        '''

        return self._mode_number_min_minus_one + 1

    @property
    def mode_number_max(self):
        '''
        Upper limit for modal analysis.
        '''

        return self._mode_number_max

    @property
    def modal_natural_frequencies(self):
        '''
        Modal natural frequency of the beam for each mode.
        '''

        return self._beam.natural_frequencies

    @property
    def oscillatory_frequencies(self):
        '''
        Oscillatory frequency for each node.
        '''
        return self._oscillatory_frequencies

    @property
    def modal_weight_force(self):
        '''
        Modal weight of force force for concerned modes.
        '''
        return self._modal_weight_force

    @property
    def modal_weight_force_std(self):
        '''
        Standard deviations of modal weight of force force for concerned modes.
        '''
        return self._modal_weight_force_std

    @property
    def modal_oscillatory_frequencies(self):
        '''
        Modal oscillatory frequency for each mode.
        '''
        return self._modal_oscillatory_frequencies

    @property
    def time_domain_analysis(self):
        '''
        True if time domain analysis has been done.
        '''
        return self._time_domain_analysis

    @property
    def frequency_domain_analysis(self):
        '''
        True if frequency domain analysis has been done.
        '''
        return self._frequency_domain_analysis

    @property
    def wavelet_analysis(self):
        '''
        True if time frequency analysis has been done.
        '''
        return self._wavelet_analysis

    @property
    def modal_analysis(self):
        '''
        True if modal analysis has been done.
        '''
        return self._modal_analysis

    @property
    def sampling_period(self):
        '''
        Sampling period of force data.
        '''
        return self._force_data.sampling_period
