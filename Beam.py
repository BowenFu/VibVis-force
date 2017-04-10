'''
Designed to calculate force modes and corresponding natural frequencies.
'''
import numpy
import numpy.matlib
import style
import matplotlib.pyplot
WATER_DENSITY = 1000
GRAVITATIONAL_ACCELERATION = 9.8


class Beam:
    '''
    Get natural frequencies and shapes of modes
    '''

    def __init__(self, node_number, mass_ratio, top_tension,
                 length, diameter, bending_stiffness, horizontal, tension_include_boyancy):
        self._node_number = node_number
        self._mass_ratio = mass_ratio
        self._top_tension = top_tension

        self._horizontal = horizontal
        self._tension_include_boyancy = tension_include_boyancy
        self._length = length
        self._diameter = diameter
        self._bending_stiffness = bending_stiffness
        self._calculate_modes()

    def _calculate_modes(self):
        '''
        To obtain the modal natural frequency
        '''
        # element length
        self._node_spacing = self._length / (self._node_number - 1)
        l = self._node_spacing
        heights = numpy.linspace(
            self._length, 0, num=self._node_number, endpoint=True)
        heights = (heights[:-1] + heights[1:]) / 2

        water_density_per_length = numpy.pi * self._diameter**2 / 4 * WATER_DENSITY
        density_minus_water = (self._mass_ratio - 1) * water_density_per_length
        density = self._mass_ratio * water_density_per_length
        n_dof = self._node_number * 2
        # mass matrix
        M = numpy.matlib.zeros((n_dof, n_dof))
        # spngness matrix plus
        K = numpy.matlib.zeros((n_dof, n_dof))

        ke = numpy.matrix([
            [12, 6 * l, -12, 6 * l], [6 * l, 4 * l**2, -6 * l, 2 * l**2],
            [-12, -6 * l, 12, -6 * l], [6 * l, 2 * l**2, -6 * l, 4 * l**2]
        ]) * self._bending_stiffness / l**3

        me = numpy.matrix([
            [156, 22 * l, 54, -13 * l],
            [22 * l, 4 * l**2, 13 * l, -3 * l**2],
            [54, 13 * l, 156, -22 * l],
            [-13 * l, -3 * l**2, -22 * l, 4 * l**2],
        ]) * density * l / 420

        # geometric spngness matrix element
        kge_ = numpy.matrix(
            [[36, 3 * l, -36, 3 * l], [3 * l, 4 * l**2, -3 * l, -l**2],
             [-36, -3 * l, 36, -3 * l], [3 * l, -l**2, -3 * l, 4 * l**2]])

        for i, height in enumerate(heights):
            T = self._top_tension
            if self._horizontal is False:
                if self._tension_include_boyancy:
                    T = self._top_tension - density_minus_water * \
                        GRAVITATIONAL_ACCELERATION * height
                else:
                    T = self._top_tension - density * \
                        GRAVITATIONAL_ACCELERATION * height
            kge = kge_ * T / l / 30
            S = numpy.matlib.zeros((4, 2 * self._node_number))
            for j in range(4):
                S[j, 2 * i + j] = 1
            K += S.T * (ke + kge) * S
            M += S.T * me * S

        self.K = K[1:, 1:]
        self.K = numpy.delete(self.K, numpy.s_[-2], 0)
        self.K = numpy.delete(self.K, numpy.s_[-2], 1)

        self.M = M[1:, 1:]
        self.M = numpy.delete(self.M, numpy.s_[-2], 0)
        self.M = numpy.delete(self.M, numpy.s_[-2], 1)

        invM = numpy.linalg.inv(self.M)
        invMK = invM * self.K
        w, shapes = numpy.linalg.eig(invMK)
        zero_r = numpy.zeros((1, shapes.shape[1]))
        shapes = numpy.r_[zero_r, shapes[:-1, :], zero_r, shapes[-1, :]]
        order = numpy.argsort(w)
        w = numpy.sqrt(w[order].real)
        shapes = shapes[:, order]

        # python slice from 0
        # so the first one must be included
        shapes = shapes[2:-3:2, :]
        shapes = shapes[:, :shapes.shape[0]]
        shapes_max = numpy.amax(numpy.abs(shapes), axis=0)
        zero_r = numpy.zeros((1, shapes.shape[1]))
        shapes = numpy.r_[zero_r, shapes, zero_r]
        self._shapes = shapes / shapes_max

        self._natural_frequencies = w / numpy.pi / 2

    @property
    def shapes(self):
        '''
        Modal shapes of forces [mode_number, node_number].
        '''
        return self._shapes

    @property
    def natural_frequencies(self):
        '''
        Natural frequencies of forces (rank-1 array).
        '''
        return self._natural_frequencies

    @property
    def diameter(self):
        '''
        Diameter of the beam.
        '''
        return self._diameter

    @property
    def node_number(self):
        '''
        Node number of the beam.
        '''
        return self._node_number

    @property
    def node_spacing(self):
        '''
        Node spacing of the beam.
        '''
        return self._node_spacing

    def plot_modal_shapes(self, out_filenames, max_order):
        '''
        Plot modal shapes up to max_order.

        Parameters
        ----------
        out_filenames : string list.
            A list of filenames for figure saving.
        max_order : int.
            The max order up to which the figure is plotted.
        '''
        print('plot_modal_shapes')
        matplotlib.pyplot.clf()
        matplotlib.pyplot.gcf().set_size_inches(
            style.SINGLE_COLUMN_WIDTH, style.SINGLE_COLUMN_SHORT_HEIGHT)
        matplotlib.pyplot.xlabel(r'$z\cdot L^{-1}$')
        matplotlib.pyplot.ylabel(r'$\phi^m$')
        matplotlib.pyplot.grid()
        z = numpy.linspace(0, 1, num=self._node_number)
        matplotlib.pyplot.plot(z, self.shapes[:, :max_order])
        matplotlib.pyplot.tight_layout()
        for out_filename in out_filenames:
            matplotlib.pyplot.savefig(out_filename)
