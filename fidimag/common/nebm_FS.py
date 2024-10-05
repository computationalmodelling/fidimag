import numpy as np
import fidimag.extensions.nebm_clib as nebm_clib

from .chain_method_tools import spherical2cartesian, cartesian2spherical, compute_norm
from .chain_method_tools import linear_interpolation_spherical
from .chain_method_tools import interpolation_Rodrigues_rotation
from .chain_method_tools import m_to_zero_nomaterial
from .chain_method_base import ChainMethodBase

from .chain_method_integrators import FSIntegrator

import scipy.integrate as spi

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="fidimag")


# TODO: we can just inherit from the geodesic nebm class!
class NEBM_FS(ChainMethodBase):
    """
    ARGUMENTS -----------------------------------------------------------------

    sim                 :: An instance of a micromagnetic or an atomistic
                           simulation. Every image in the band will be a copy
                           of this simulation.

    initial_images      :: A list containing numpy arrays or space dependent
                           functions to set the magnetisation fields for every
                           image in the band. It is suggested that the first
                           and last elements define stable states of the
                           magnetic system. The arrays and functions are used
                           to load the magnetisation/spin fields through the
                           sim.set_m method from the Simulation object.

    interpolations      :: A list with 1 element less than the initial_images
                           list, where every entry is an integer indicating the
                           number of interpolations between consecutive images.
                           For example, if we defined initial_images as
                           [state_1, state_2, state_3], and we want 10
                           interpolations between state_1 and state_2 and 5
                           interpolations between state_2 and state_3, we set
                           interpolations as [10, 5], making an energy band of
                           17 images. If we do not want any interpolation, we
                           leave this list as None or empty.

    interpolation_method:: In case that a number of interpolations were
                           defined, it is possible to specify how the
                           interpolation is performed using any of these
                           methods:

                                'linear'   : A linear interpolation of the spin
                                             directions using spherical
                                             coordinates

                                'rotation' : Interpolation of the spin
                                             directions using Rodrigue's
                                             rotation formulae

    spring_constant     :: The spring constant magnitude

    name                :: The NEBM simulation name. Folders for VTK and NPY
                           files, and data tables are named according to this
                           string.

    openmp              :: Set this as True to use the parallelised version of
                           CVODE, which is the integrator used to evolve the
                           NEBM minimisation equation.

    ---------------------------------------------------------------------------

    The NEB Method (NEBM) class to find minimum energy paths between two stable
    states in a given magnetic system. This class works both for atomistic and
    micromagnetic simulations. The Geodesic NEBM describes the spins /
    magnetisation vectors in Cartesian coordinates (m_x, m_y, m_z) and
    distances in the energy landscape are measured by a geodesic distance (see
    below). The NEBM is based on the definition of a so called band, which is
    just a sequence of replicas (in terms of geometry and magnetic parameters)
    of the same system (given by our simulation), called images, in different
    magnetic configurations, and where the first and last state are the stable
    states used to find a minimum energy transition between them. Calling the
    images as Y_i, after relaxation an energy band of N+1 images usually looks
    like:


                 Energy                ...
                   ^                                           _
                   |         Y_2   , - ~ ~ ~ - ,   Y_(N-1)    |
                   |           O '               O ,          |  Energy barrier
                   |    Y_1  ,                       ,        |  with respect
                   |        O                         O Y_N   |_ to Y_N
                   |       ,
                   |       O
                   |    Y_0
                   ________________________
                   Distance

    where Y_0 and Y_N are the stable states.

    For more details about the definition of the forces involved in the NEBM,
    see the following papers:

        - Suess et al., Physical Review B 75, 174430 (2007)
        - Henkelman et al., Journal of Chemical Physics 113, 22 (2000)

    """

    def __init__(self, sim,
                 initial_images,
                 interpolations=None,
                 interpolation_method='rotation',
                 spring_constant=1e5,
                 name='unnamed',
                 climbing_image=None,
                 openmp=False,
                 # integrator='sundials'  # or scipy
                 ):

        super(NEBM_FS, self).__init__(sim,
                                      initial_images,
                                      interpolations=interpolations,
                                      spring_constant=spring_constant,
                                      name=name,
                                      climbing_image=climbing_image,
                                      dof=3,
                                      openmp=openmp
                                      )

        # We need the gradient norm to calculate the action
        self.gradientENorm = np.zeros(self.n_images)

        # Initialisation ------------------------------------------------------
        # See the NEBMBase class for details

        self.generate_initial_band(method=interpolation_method)

        self.initialise_energies()

        self.initialise_integrator()

        self.create_tablewriter()

        # ---------------------------------------------------------------------

    def initialise_energies(self):
        # Energy of the images
        self.band = self.band.reshape(self.n_images, -1)
        for i in range(self.n_images):
            self.sim.set_m(self.band[i])
            self.sim.compute_effective_field(t=0)
            self.energies[i] = self.sim.compute_energy()
        self.band = self.band.reshape(-1)

    def initialise_integrator(self):
        self.t = 0
        self.iterations = 0
        self.ode_count = 1

        # Use default integrator parameters. Pass this class object to the integrator
        self.integrator = FSIntegrator(self)
 
    def generate_initial_band(self, method='linear'):
        """
        method      :: linear, rotation

        """

        # Every row will be an image of the band, i.e. the i-th row is
        # the i-t image
        self.band = self.band.reshape(self.n_images, -1)

        # Indexes indicating the image number (position of the images) in the
        # band, for the specified initial images
        i_initial_images = [0]
        for i in range(1, len(self.initial_images)):
            i_initial_images.append(i_initial_images[-1]
                                    + self.interpolations[i - 1]
                                    + 1
                                    )

        for i, m_field in enumerate(self.initial_images[:-1]):

            # Copy the magnetisation field from the i-th and (i + 1)-th to the
            # corresponding rows of the nebm band array To do this, we need to
            # know in which positions are these images in the band, which
            # change according to the number of interpolations. Accordingly,
            # we use the list with the indexes of the initial images
            self.sim.set_m(self.initial_images[i])
            self.band[i_initial_images[i]] = np.copy(self.sim.spin)

            self.sim.set_m(self.initial_images[i + 1])
            self.band[i_initial_images[i + 1]] = np.copy(self.sim.spin)

            # interpolation is an array with *self.interpolations[i]* rows
            # We copy these rows to the corresponding images in the energy
            # band array
            if self.interpolations[i] != 0:
                if method == 'linear':
                    interpolation = linear_interpolation_spherical(
                        cartesian2spherical(self.band[i_initial_images[i]]),
                        cartesian2spherical(self.band[i_initial_images[i + 1]]),
                        self.interpolations[i],
                        self.sim._pins
                        )

                    interpolation = np.apply_along_axis(spherical2cartesian,
                                                        axis=1,
                                                        arr=interpolation)
                elif method == 'rotation':
                    interpolation = interpolation_Rodrigues_rotation(
                        self.band[i_initial_images[i]],
                        self.band[i_initial_images[i + 1]],
                        self.interpolations[i],
                        self.sim._pins
                        )

                interpolation = np.apply_along_axis(lambda m: m_to_zero_nomaterial(m, self.sim),
                                                    axis=1,
                                                    arr=interpolation)

                # We then set the interpolated spins fields at once
                self.band[i_initial_images[i] + 1:
                          i_initial_images[i + 1]] = interpolation

        # expand the energy band array
        self.band = self.band.reshape(-1)

    def compute_effective_field_and_energy(self, y):
        """
        Compute the effective field and the energy of every image in the band,
        using the array `y` as the degrees of freedom of the system (i.e. the
        one that contains all the spin directions of the images in the band).

        The local copy of the `y` array for this NEBM class is the self.band
        array, which we update at the end of every call to the integrator in
        the relaxation function
        """

        self.gradientE.shape = (self.n_images, -1)

        y.shape = (self.n_images, -1)

        # Do not update the extreme images
        for i in range(1, len(y) - 1):

            self.sim.set_m(y[i])
            # elif self.coordinates == 'Cartesian':
            #     self.sim.set_m(self.band[i])

            self.sim.compute_effective_field(t=0)

            self.gradientE[i][:] = -self.sim.field

            self.energies[i] = self.sim.compute_energy()

        y.shape = (-1)
        self.gradientE.shape = (-1)

    def compute_tangents(self, y):
        nebm_clib.compute_tangents(self.tangents, y, self.energies,
                                   self.n_dofs_image, self.n_images
                                   )
        nebm_clib.project_images(self.tangents, y,
                                 self.n_images, self.n_dofs_image
                                 )
        # nebm_clib.normalise_images(self.tangents,
        #                            self.n_images, self.n_dofs_image
        #                            )

    def compute_spring_force(self, y):
        """
        For variable spring constant (which is more effective if we have
        a saddle point), see:
             J. Chem. Phys. 113, 9901 (2000);
        Seems to work when we only have a single saddle point
        (TESTING functionality)
        """
        if self.variable_k:
            E_max = np.max(self.energies)
            E_i = np.maximum(self.energies[1:-1], self.energies[:-2])
            E_ref = max(self.energies[0], self.energies[-1])

            k = np.copy(self.k)
            f = E_i > E_ref
            k[1:-1][f] = k[1:-1][f] - self.dk * ((E_max - E_i[f]) / (E_max - E_ref))
            f = E_i <= E_ref
            k[1:-1][f] = k[1:-1][f] - self.dk

        else:
            k = self.k

        # Compute the distances
        nebm_clib.image_distances_GreatCircle(self.distances,
                                              self.path_distances,
                                              y,
                                              self.n_images,
                                              self.n_dofs_image,
                                              self._material_int,
                                              self.n_dofs_image_material
                                              )

        nebm_clib.compute_spring_force(self.spring_force, y,
                                       self.tangents,
                                       k, self.n_images,
                                       self.n_dofs_image,
                                       self.distances
                                       )

    def compute_action(self):
        """
        """
        # Check that action is computed AFTER calculating and projecting the forces

        # nebm_clib.image_distances_GreatCircle(self.distances,
        #                                       self.path_distances,
        #                                       y,
        #                                       self.n_images,
        #                                       self.n_dofs_image,
        #                                       self._material_int,
        #                                       self.n_dofs_image_material
        #                                       )
        
        # NOTE: Gradient here is projected in the S2^N tangent space
        self.gradientE.shape = (self.n_images, -1)
        # NOTE: HEre we have to divide by the number of spins per image,not n_images:
        Gnorms2 = np.sum(self.gradientE**2, axis=1) / self.n_spins

        # Compute the root mean square per image
        self.gradientENorm[:] = np.sqrt(Gnorms2)
        self.gradientE.shape = (-1)

        # DEBUG:
        # print('gradEnorm', self.gradientENorm)

        # TODO: we can use a better quadrature such as Gaussian
        # notice that the gradient norm here is using the RMS
        action = spi.simpson(self.gradientENorm, x=self.path_distances)
        # print('E', self.energies / (self.mesh.dx * self.mesh.dy * self.mesh.dz * self.mesh.unit_length**3))
        # print('gradE norm', self.gradientENorm)
        # print('Path distance', self.path_distances)
        print('Images', self.band.reshape(-1, 3).reshape(self.n_images, -1))

        # DEBUG:
        # print('action from gradE', action)

        # The spring term in the action is added as |F_k|^2 / (2 * self.k) = self.k * x^2 / 2
        # (CHECK) This assumes the spring force is orthogonal to the force gradient (after projection)
        # Knorms2 = np.sum(self.spring_force.reshape(-1, 3)**2, axis=1)
        # # Compute the root mean square per image
        # springF_norms = np.sqrt(np.mean(Knorms2.reshape(self.n_images, -1), axis=1))

        # Norm of the spring force per image (assuming tangents as unit vectors)
        # These are the norms of the inner images
        dist_plus_norm = self.distances[1:]
        dist_minus_norm = self.distances[:-1]
        # dY_plus_norm = distances[i];
        # dY_minus_norm = distances[i - 1];
        springF2 = 0.5 * self.k[1:-1] * ((dist_plus_norm - dist_minus_norm)**2)
        # CHECK: do we need to scale?
        # action += np.sum(springF2) / (self.n_images - 2)

        # DEBUG:
        # print('action spring term', np.sum(springF2) / (self.n_images - 2))

        return action

    def compute_min_action(self):
        dE = self.energies[-1] - self.energies[0]
        minAction = np.sum(np.abs(self.energies[1:] - self.energies[:-1]))
        return 2 * (dE + minAction) / (self.mesh.dx * self.mesh.dy * self.mesh.dz * self.mesh.unit_length**3) / self.path_distances[-1]

    def nebm_step(self, y, ensure_zero_extrema=False):

        self.compute_effective_field_and_energy(y)
        nebm_clib.project_images(self.gradientE, y, self.n_images, self.n_dofs_image)
        self.compute_tangents(y)
        nebm_clib.normalise_spins(self.tangents, self.n_images, self.n_dofs_image)
        self.compute_spring_force(y)

        nebm_clib.compute_effective_force(self.G,
                                          self.tangents,
                                          self.gradientE,
                                          self.spring_force,
                                          self._climbing_image,
                                          self.n_images,
                                          self.n_dofs_image
                                          )

        # The effective force at the extreme images should already be zero, but
        # we will manually remove any value
        if ensure_zero_extrema:
            self.G[:self.n_dofs_image] = 0
            self.G[-self.n_dofs_image:] = 0

        # Should be the same if we project the gradient before, instead
        # of the total force
        # nebm_clib.project_images(self.G, y,
        #                          self.n_images, self.n_dofs_image
        #                          )

    # -------------------------------------------------------------------------
    # Methods -----------------------------------------------------------------
    # -------------------------------------------------------------------------

    def compute_distances(self):
        """
        Compute the distance between corresponding images of self.band

                A                   B
            [ [image_0]         [ [image_0]
              [image_1]     -     [image_1]
              ...                 ...
            ]                     ]

        """

        nebm_clib.image_distances_GreatCircle(self.distances,
                                              self.path_distances,
                                              self.band,
                                              self.n_images,
                                              self.n_dofs_image,
                                              self._material_int,
                                              self.n_dofs_image_material
                                              )

    # def compute_maximum_dYdt(self, A, B, dt):
    #     """
    #     In case we want to use a Geodesic distance instead of the scaled
    #     norm between corresponding images
    #     """
    #     # # We will not consider the images at the extremes to compute dY
    #     band_no_extremes = slice(self.n_dofs_image, -self.n_dofs_image)
    #     dYdt = self.compute_distances(A[band_no_extremes],
    #                                   B[band_no_extremes])
    #     dYdt /= dt
    #     if np.max(dYdt) > 0:
    #         return np.max(dYdt)
    #     else:
    #         return 0

    def step_RHS(self, band):
        """

        This function is called on every iteration of the integrators in
        chain_method_integrators.py

        """

        self.ode_count += 1

        # Update the effective field, energies, spring forces and tangents
        # using the *y* array
        self.nebm_step(band)

        # Now set the RHS of the equation as the effective force on the energy
        # band, which is stored on the self.G array
        # ydot = self.G[:]

        # The effective force at the extreme images should already be zero, but
        # we will manually remove any value
        self.G[:self.n_dofs_image] = 0
        self.G[-self.n_dofs_image:] = 0

        return 0

    def Sundials_RHS(self, t, y, ydot):
        """

        This function is called on every iteration of the integrator (CVODE
        solver). ydot refers to the Right Hand Side of the equation, since
        we are solving dy/dt = 0

        """

        self.ode_count += 1

        # Update the effective field, energies, spring forces and tangents
        # using the *y* array
        self.nebm_step(y)

        # Now set the RHS of the equation as the effective force on the energy
        # band, which is stored on the self.G array
        ydot[:] = self.G[:]

        # Update the step with the optimisation algorithm, in this
        # case we use: dY /dt = Y x Y x D - correction-factor
        # (check the C code in common/)
        nebm_clib.compute_dYdt(
            y, self.G, ydot, self.sim._pins, self.n_images, self.n_dofs_image)

        # The effective force at the extreme images should already be zero, but
        # we will manually remove any value
        ydot[:self.n_dofs_image] = 0
        ydot[-self.n_dofs_image:] = 0

        return 0

