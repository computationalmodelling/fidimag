import fidimag

from . import llg
from . import llg_stt
from . import llg_stt_cpp
from . import baryakhtar
from fidimag.common import steepest_descent
from fidimag.common import hubert_minimiser

import fidimag.extensions.micro_clib as micro_clib
import fidimag.common.helper as helper

from fidimag.common.sim_base import SimBase

import numpy as np

KNOWN_DRIVERS = {'llg': llg.LLG,
                 'llg_stt': llg_stt.LLG_STT,
                 'llg_stt_cpp': llg_stt_cpp.LLG_STT_CPP,
                 'llbar': baryakhtar.LLBar,
                 'llbar_full': baryakhtar.LLBarFull,
                 'steepest_descent': steepest_descent.SteepestDescent,
                 'hubert_minimiser': hubert_minimiser.HubertMinimiser,
                 }


class Sim(SimBase):
    """

    The Simulation class for a micromagnetic simulation.

    ARGUMENTS:

    mesh    :: A CuboidMesh instance, from fidimag.common

    OPTIONAL ARGUMENTS:

    name    :: A string with the simulation name
    driver  :: A string with one of the following drivers to evolve
               the micromagnetic system:

                llg               - (default) the Landau-Lifshitz-Gilbert
                                    equation

                llg_stt           - LLG w. spin transfer torque

                llg_stt_cpp       - LLG w. spin transfer torque perpendicular
                                    to the plane (?)

                llbar             - Landau-Lifshitz-Baryakhtar equation
                llbar_full

                steepest_descent  - Optimised steepest descent minimisation
                                    [JAP 115, 17D118 (2014)]

    ** Most of the methods inherit from the Base Simulation class that
    can be found in the common folder

    """

    def __init__(self, mesh, name='unnamed', driver='llg',
                 # Integrator arguments:
                 integrator='sundials', use_jac=False
                 ):

        super(Sim, self).__init__(mesh, name)

        self._micromagnetic = True

        # Saturation magnetisation definitions:
        # self._Ms = np.zeros(self.n, dtype=np.float)
        # David: be careful to change these references to the common mag array
        self._Ms = self._magnetisation
        # Remember this is a 3 * n array:
        self._Ms_inv = self._magnetisation_inv

        # Here we link one of the drivers to evolve the LLG equation
        if driver not in KNOWN_DRIVERS:
            raise NotImplementedError("""Driver '{}' is not implemented.
                                      Valid choices: one of '{}'.""".format(driver, KNOWN_DRIVERS.keys()))

        self.driver = KNOWN_DRIVERS[driver](self.mesh,
                                            self.spin,
                                            self._Ms,
                                            self._Ms_inv,
                                            self.field,
                                            self._pins,
                                            self.interactions,
                                            self.name,
                                            self.data_saver,
                                            integrator=integrator,
                                            use_jac=use_jac
                                            )

        # For the SD in the micromagnetic class we need to scale the mxmxH factor
        # by mu0 to leave the field in Tesla units
        if self.driver.__class__.__name__ == 'SteepestDescent':
            self.driver.scale = fidimag.common.constant.mu_0

        # Some references to functions in the corresponding driver classes
        # that can be accessed through the Simulation class
        self.relax = self.driver.relax
        self.compute_effective_field = self.driver.compute_effective_field
        self.save_vtk = self.driver.save_vtk
        self.save_m = self.driver.save_m
        self.save_skx = self.driver.save_skx

    def get_Ms(self):
        """
        Returns the array with the saturation magnetisation values per
        mesh site
        """
        return self._Ms

    def set_Ms(self, value):
        """

        Set the saturation magnetisation of the system as a uniform or
        spatially dependent scalar field, where values are in A /m.  Mesh sites
        with no material can be specified with a magnetisation of zero
        magnitude. Samples with different materials can be specified with
        different magnetisation values in specific regions of the system.

        ARGUMENTS:

        value     :: * For a homogeneous single material sample, specify a
                       float with a magnitude in A /m

                     * Alternatively, you can specify a function that returns
                       values in A /m, which depends on the spatial
                       coordinates. For example, a 2 nm wide cylinder centered
                       at (x, y) = (1, 1) can be specified with (if you set the
                       unit_length to 1e-9 in the mesh):

                            Ms = 1e6  # A / m

                            def Ms_profile(r):
                                for (r[0] - 1) ** 2 + (r[1] - 1) ** 2 <= 1 ** 2:
                                    return Ms
                                else:
                                    return 0

                            Sim.set_Ms(Ms_profile)

                     * You can also manually specify an array with n values
                     with the magnetisation magnitudes, in the same order than
                     the mesh coordinates array.

                     * Alternatively, if you previously saved the magnetisation
                     array to a numpy file, you can load it using
                     numpy.load(my_array).

        """

        self._Ms[:] = helper.init_scalar(value, self.mesh)
        nonzero = 0
        for i in range(self.n):
            if self._Ms[i] > 0.0:
                self._Ms_inv[i] = 1.0 / self._Ms[i]
                nonzero += 1

        # We moved this variable to the micro_driver class
        self.driver.n_nonzero = nonzero

        for i in range(len(self._Ms)):
            if self._Ms[i] == 0.0:
                self._pins[i] = 1

                # Set the neighbour index to -1 for sites with Ms = 0
                self.mesh.neighbours[self.mesh.neighbours == i] = -1

        # TODO: Check if this is necessary here, it is only defined
        # for the LLG STT in the drivers
        self.driver.Ms_const = np.max(self._Ms)

    Ms = property(get_Ms, set_Ms)

    def skyrmion_number(self):
        """
        Returns the skyrmion number from the first layer in the
        XY plane of the mesh
        """
        nx = self.mesh.nx
        ny = self.mesh.ny
        nz = self.mesh.nz
        number = micro_clib.compute_skyrmion_number(
            self.spin, self._skx_number, nx, ny, nz, self.mesh.neighbours)
        return number

    # We need to change these in the future so we can document them:
    skyrmion_number_slice = fidimag.common.skyrmion_number.skyrmion_number_slice
    skyrmion_number_lee = fidimag.common.skyrmion_number.skyrmion_number_lee
