from __future__ import print_function
from __future__ import division
import numpy as np
import fidimag.extensions.nebm_clib as nebm_clib

from .chain_method_tools import spherical2cartesian, cartesian2spherical, compute_norm
from .chain_method_tools import linear_interpolation_spherical
from .chain_method_tools import interpolation_Rodrigues_rotation
from .chain_method_tools import m_to_zero_nomaterial

from .chain_method_base import ChainMethodBase

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="fidimag")


class NEBM_Geodesic(ChainMethodBase):
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

    The NEBM evolves an energy band [Y_0, Y_1, ... , Y_N]  according to the
    equation
                                     _______
                                    /     2
        dY                         /( dY )           2
        --- =  -Y x Y x G + c *   / ( -- )   * (1 - Y ) Y
        dt                      \/  ( dt )

    for every image Y of the energy band, which is the same equation than the
    Cartesian NEBM code (since we describe the spin field in the same fashion).
    The G vector is the effective force on the image, which is defined in terms
    of tangent vectors (tangent along the band) in the energy landscape, whose
    perpendicular components follow directions of largest energy changes to
    find the minimum energy transition.  Additionally, the effective force
    includes a spring force that keeps images equally spaced along the band to
    avoid clustering around minima or saddle points. The Geodesic NEBM requires
    to project firstly the spring force, and then the total effective force, on
    the tangent space defined by the spin/magnetisation field (see references
    below).

    The spacing between images needs a definition of DISTANCE. In this code we
    use an Geodesic distance, normalised by the number of degrees of freedom,
    which is the sum of all spin components, so if we have P spins in the
    system, the number of dofs is 3 * P, i.e. the 3 directions per spin.  The
    distance is defined as:

                                      ___________________
                                     / P
                                    /  __             2
        distance(Y_i, Y_j) =       /  \   ( L(i,j)_a )
                               \  /   /__
                                \/    a=1

    where
                                 ->       ->        ->       ->
         L(i,j)_a  =  arctan2( | m(i)_a x m(j)_a |, m(i)_a o m(j)_a )

    where m(i)_a is the a-th spin of the image Y_i. This is known as Vincenty's
    formula.

    The second term in the minimisation equation is to correct the
    magnetisation/spin vectors length, since at 0 K this length is fixed. For
    now, we use the *c* factor as 6.

    For more details about the definition of the forces involved in the NEBM,
    see the following papers:

        - Suess et al., Physical Review B 75, 174430 (2007)
        - Henkelman et al., Journal of Chemical Physics 113, 22 (2000)
        - Bessarab et al., Computer Physics Communications 196 (2015) 335-347

    """

    def __init__(self, sim,
                 initial_images,
                 interpolations=None,
                 interpolation_method='rotation',
                 spring_constant=1e5,
                 name='unnamed',
                 climbing_image=None,
                 openmp=False,
                 integrator='sundials'  # or scipy
                 ):

        super(NEBM_Geodesic, self).__init__(sim,
                                            initial_images,
                                            interpolations=interpolations,
                                            spring_constant=spring_constant,
                                            name=name,
                                            climbing_image=climbing_image,
                                            dof=3,
                                            openmp=openmp
                                            )

        # Initialisation ------------------------------------------------------
        # See the NEBMBase class for details

        self.generate_initial_band(method=interpolation_method)

        self.initialise_energies()

        self.initialise_integrator(integrator=integrator)

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
        using the array *y* as the degrees of freedom of the system (i.e. the
        one that contains all the spin directions of the images in the band).

        The local copy of the *y* array for this NEBM class is the self.band
        array, which we update at the end of every call to the integrator in
        the relaxation function

        """

        self.gradientE = self.gradientE.reshape(self.n_images, -1)

        y = y.reshape(self.n_images, -1)

        # Only update the extreme images
        for i in range(1, len(y) - 1):

            self.sim.set_m(y[i])
            # elif self.coordinates == 'Cartesian':
            #     self.sim.set_m(self.band[i])

            self.sim.compute_effective_field(t=0)

            self.gradientE[i][:] = -self.sim.field

            self.energies[i] = self.sim.compute_energy()

        y = y.reshape(-1)
        self.gradientE = self.gradientE.reshape(-1)

    def compute_tangents(self, y):
        nebm_clib.compute_tangents(self.tangents, y, self.energies,
                                   self.n_dofs_image, self.n_images
                                   )
        nebm_clib.project_images(self.tangents, y,
                                 self.n_images, self.n_dofs_image
                                 )
        nebm_clib.normalise_images(self.tangents,
                                   self.n_images, self.n_dofs_image
                                   )

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

        # GreatCircle calculation of the geodesic seems to be more stable
        # in complex landscapes than Vincenty's formula
        # nebm_geodesic.image_distances_Vincenty(self.distances,
        #                                        self.path_distances,
        #                                        y,
        #                                        self.n_images,
        #                                        self.n_dofs_image,
        #                                        self._material_int,
        #                                        self.n_dofs_image_material
        #                                        )

    def nebm_step(self, y):

        self.compute_effective_field_and_energy(y)
        nebm_clib.project_images(self.gradientE, y,
                                 self.n_images, self.n_dofs_image
                                 )
        self.compute_tangents(y)
        self.compute_spring_force(y)

        nebm_clib.compute_effective_force(self.G,
                                          self.tangents,
                                          self.gradientE,
                                          self.spring_force,
                                          self._climbing_image,
                                          self.n_images,
                                          self.n_dofs_image
                                          )

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

    def step_RHS(self, t, y):
        """

        This function is called on every iteration of the integrators in
        chain_method_integrators.py

        """

        self.ode_count += 1

        # Update the effective field, energies, spring forces and tangents
        # using the *y* array
        self.nebm_step(y)

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
