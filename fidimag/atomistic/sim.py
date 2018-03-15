from .llg import LLG
from .sllg import SLLG
from .llg_stt import LLG_STT
from .llg_stt_cpp import LLG_STT_CPP
from .steepest_descent import SteepestDescent
from .minimiser import Minimiser

import fidimag.common.skyrmion_number
import fidimag.common.helper as helper
import fidimag.extensions.clib as clib
import numpy as np

from fidimag.common.sim_base import SimBase

KNOWN_DRIVERS = {'llg': LLG,
                 'sllg': SLLG,
                 'llg_stt': LLG_STT,
                 'llg_stt_cpp': LLG_STT_CPP,
                 'steepest_descent': SteepestDescent,
                 'minimiser': Minimiser
                 }


class Sim(SimBase):
    """

    The Simulation class for an atomistic simulation

    ARGUMENTS:

    mesh    :: For a square mesh, use an instance of CuboidMesh,
               from fidimag.common

               For a hexagonal mesh, use an instance of Hexagonalmesh
               from fidimag.atomistic

    OPTIONAL ARGUMENTS:

    name    :: A string with the simulation name
    driver  :: A string with one of the following drivers to evolve
               the micromagnetic system:

                    llg         - (default) the Landau-Lifshitz-Gilbert
                                  equation
                    llg_stt     - LLG w. spin transfer torque
                    llg_stt_cpp - LLG w. spin transfer torque perpendicular
                                  to the plane (?)

    ** Most of the methods inherit from the Base Simulation class that
    can be found in the common folder

    """

    def __init__(self, mesh, name='unnamed', driver='llg',
                 # Integrator arguments:
                 integrator='sundials',
                 use_jac=False
                 ):

        super(Sim, self).__init__(mesh, name)

        # We should change this:
        self._micromagnetic = False

        # Magnetic moments definitions:
        # self._mu_s = np.zeros(self.n, dtype=np.float)
        self._mu_s = self._magnetisation
        self._mu_s_inv = np.zeros(self.n, dtype=np.float)

        # This is only for old C files using the xperiodic variable
        (self.xperiodic,
         self.yperiodic) = mesh.periodicity[0], mesh.periodicity[1]

        # Here we link one of the drivers to evolve the LLG equation
        if driver not in KNOWN_DRIVERS:
            raise NotImplementedError("""Driver '{}' is not implemented.
                                      Valid choices: one of '{}'.""".format(driver,
                                                                            KNOWN_DRIVERS.keys())
                                      )

        self.driver = KNOWN_DRIVERS[driver](self.mesh,
                                            self.spin,
                                            self._mu_s,
                                            self._mu_s_inv,
                                            self.field,
                                            self._pins,
                                            self.interactions,
                                            self.name,
                                            self.data_saver,
                                            integrator=integrator,
                                            use_jac=use_jac
                                            )

        # Some references to functions in the corresponding driver classes
        # that can be accessed through the Simulation class
        self.relax = self.driver.relax
        self.compute_effective_field = self.driver.compute_effective_field
        self.save_vtk = self.driver.save_vtk
        self.save_m = self.driver.save_m
        self.save_skx = self.driver.save_skx

    # -------------------------------------------------------------------------

    def get_mu_s(self):
        """
        Returns the array with the magnetic moments values per lattice site
        """
        return self._mu_s

    def set_mu_s(self, value):
        """

        Set the magnetic moments values of the system as a uniform or spatially
        dependent scalar field, where values are in J / T. It is recommended to
        set values in Bohr magneton units and to use the
        fidimag.common.constant.mu_B magnitude (see example). Mesh sites with
        no material can be specified with a magnetic moment of zero magnitude.
        Samples with different materials can be specified with different
        magnetic moments in specific regions of the system.

        ARGUMENTS:

        value     :: * For a homogeneous single material sample, specify a
                       float with a magnitude in J / T. It is recommended to
                       use Bohr magneton units taking the constant from the
                       `constant` library in fidimag.common. For example:

                            import fidimag.common.constant as const

                            mu_s = 2 * const.mu_B

                     * In addition, you can specify a function that returns
                       values in J / T, which depends on the spatial
                       coordinates. For example, a 2 nm wide cylinder centered
                       at (x, y) = (1, 1) can be specified with (if you set the
                       unit_length to 1e-9 in the mesh):

                            import fidimag.common.constant as const
                            mu_s = 2 * const.mu_B

                            def mu_s_profile(r):
                                for (r[0] - 1) ** 2 + (r[1] - 1) ** 2 <= 1 ** 2:
                                    return mu_s
                                else:
                                    return 0

                            Sim.set_mu_s(mu_s_profile)

                     * You can also manually specify an array with n values
                     with the magnetisation values, in the same order than the
                     mesh coordinates array.

                     * Alternatively, if you previously saved the magnetic
                     moments array to a numpy file, you can load it using
                     numpy.load(my_array).

        """

        self._mu_s[:] = helper.init_scalar(value, self.mesh)
        nonzero = 0
        for i in range(self.n):
            if self._mu_s[i] > 0.0:
                self._mu_s_inv[i] = 1.0 / self._mu_s[i]
                nonzero += 1

        # We moved this variable to the micro_driver class
        self.n_nonzero = nonzero

        for i in range(len(self._mu_s)):
            if self._mu_s[i] == 0.0:
                self._pins[i] = 1

                # Set the neighbour index to -1 for sites with mu_s = 0
                self.mesh.neighbours[self.mesh.neighbours == i] = -1

        # TODO: Check if this is necessary here, it is only defined
        # for the LLG STT in the drivers
        self.driver.mu_s_const = np.max(self._mu_s)

    mu_s = property(get_mu_s, set_mu_s)

    def skyrmion_number(self, method='FiniteSpinChirality'):
        """
        Calculate the skyrmion number using different methods:

        method      :: 'FiniteSpinChirality', 'BergLuscher'

        which are specifically defined for discrete spin lattices.

        Calling this function will fill the _skx_number array with the skyrmion
        number density per lattice site.

        It is important to mention that this calculation only works for a
        2D layer in the XY plane.

        See the corresponding C code at fidimag/atomistic/lib/util.c for a
        detailed documentation about the methods.

        """
        nx = self.mesh.nx
        ny = self.mesh.ny
        nz = self.mesh.nz

        if method == 'FiniteSpinChirality':
            if self.mesh.mesh_type == 'cuboid':
                number = clib.compute_skyrmion_number(
                    self.spin, self._skx_number, nx, ny, nz,
                    self.mesh.neighbours, self.mesh.n_ngbs)
            else:
                raise ValueError('FiniteSpinChirality method only'
                                 ' defined for cuboid meshes')

        elif method == 'BergLuscher':
            number = clib.compute_skyrmion_number_BergLuscher(
                self.spin, self._skx_number, nx, ny, nz,
                self.mesh.neighbours, self.mesh.n_ngbs)
        else:
            raise ValueError('Specify a valid method')

        return number

    skyrmion_number_slice = fidimag.common.skyrmion_number.skyrmion_number_slice
    skyrmion_number_lee = fidimag.common.skyrmion_number.skyrmion_number_lee
