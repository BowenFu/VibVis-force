'''
Designed to read in file and provide data to ForceAnalysis.
'''

import numpy


class ForceData:
    '''
    Preprocess (reshape and nondimensionize) force data.

    Parameters
    ----------
    time : rank-1 array.
        Time sequence, repeating node number times for each time step.
    force : rank-1 array.
        Force sequence, force for each node for each time step.

    time and force share the same length.

    '''

    DECIMALS = 6
    '''Global precision'''

    def __init__(self, time, force):

        # [time axis, node axis]
        if time.size != force.size:
            raise ValueError(
                'The length of time is not equal to that of force!')

        # check if appropriate to reshape
        self._time, self._node_number = numpy.unique(time, return_counts=True)
        self._node_number = numpy.unique(self._node_number)
        if self._node_number.size != 1:
            print(self._node_number)
            raise ValueError(
                'Force data is not complete! Node number is not unique!')
        self._node_number = self._node_number[0]
        self._sampling_period, count = numpy.unique(
            numpy.around(
                numpy.diff(self._time), decimals=self.DECIMALS),
            return_counts=True)
        if count.size != 1:
            print(self._sampling_period)
            print(count)
            print(numpy.diff(self._time))
            raise ValueError(
                'Force data should be cropped so to make sampling period unique!'
            )
        self._sampling_period = self._sampling_period[0]
        self._force = force.reshape((-1, self._node_number))

    @property
    def time(self):
        '''
        Time array in which duplicate has been eliminated.
        '''
        return self._time

    @property
    def node_number(self):
        '''
        Node number of force.
        '''
        return self._node_number

    @property
    def sampling_period(self):
        '''
        Sampling period of force.
        '''
        return self._sampling_period

    @property
    def force(self):
        '''
        Rank-2 force [time, node number].
        '''
        return self._force
