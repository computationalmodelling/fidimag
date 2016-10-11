from __future__ import print_function
from __future__ import division
import numpy as np

import fidimag.extensions.nebm_cartesian_clib as nebm_cartesian
import fidimag.extensions.nebm_clib as nebm_clib

from .nebm_tools import spherical2cartesian, cartesian2spherical, compute_norm
from .nebm_tools import linear_interpolation_spherical
from .nebm_tools import interpolation_Rodrigues_rotation

from .nebm_base import NEBMBase

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="fidimag")


class NEBM_Cartesian(NEBMBase):

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

    k                   :: The spring constant magnitude

    name                :: The NEBM simulation name. Folders for VTK and NPY
                           files, and data tables are named according to this
                           string.

    openmp              :: Set this as True to use the parallelised version of
                           CVODE, which is the integrator used to evolve the
                           NEBM minimisation equation.

    ---------------------------------------------------------------------------

    The NEB Method (NEBM) class to find minimum energy paths between two stable
    states in a given magnetic system. This class works both for atomistic and
    micromagnetic simulations. The NEBM in Cartesian coordinates describes the
    spins / magnetisation vectors in Cartesian coordinates (m_x, m_y, m_z). The
    NEBM is based on the definition of a so called band, which is just a
    sequence of replicas (in terms of geometry and magnetic parameters) of the
    same system (given by our simulation), called images, in different magnetic
    configurations, and where the first and last state are the stable states
    used to find a minimum energy transition between them. Calling the images
    as Y_i, after relaxation an energy band of N+1 images usually looks like:


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

    for every image Y of the energy band. The G vector is the effective force
    on the image, which is defined in terms of tangent vectors (tangent along
    the band) in the energy landscape, whose perpendicular components follow
    directions of largest energy changes to find the minimum energy transition.
    Additionally, the effective force includes a spring force that keeps images
    equally spaced along the band to avoid clustering around minima or saddle
    points. The spacing between images needs a definition of DISTANCE. In this
    code we use an Euclidean distance, normalised by the number of degrees of
    freedom, which is the sum of all spin components, so if we have P spins in
    the system, the number of dofs is 3 * P, i.e. the 3 directions per spin.
    The distance is defined as:

                                             _________________________________
                                            /  P
                                1          /  __   __                        2
        distance(Y_i, Y_j) =   ---        /  \    \   [ m(i,a)_b - m(j,a)_b ]
                               3*P    \  /   /__  /__
                                       \/    a=1  b=x,y,z

    where m(i,a)_b is the b-component of the a-th spin in the image Y_i.

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
                 interpolation_method='linear',
                 spring_constant=1e5,
                 name='unnamed',
                 climbing_image=None,
                 openmp=False
                 ):

        super(NEBM_Cartesian, self).__init__(sim,
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

            self.sim.spin[:] = y[i][:]
            # elif self.coordinates == 'Cartesian':
            #     self.sim.set_m(self.band[i])

            self.sim.compute_effective_field(t=0)

            self.gradientE[i][:] = -self.sim.field[:]

            self.energies[i] = self.sim.compute_energy()

        y = y.reshape(-1)
        self.gradientE = self.gradientE.reshape(-1)

    def compute_tangents(self, y):
        nebm_clib.compute_tangents(self.tangents, y, self.energies,
                                   self.n_dofs_image, self.n_images
                                   )

        nebm_cartesian.normalise_images(self.tangents,
                                        self.n_images, self.n_dofs_image
                                        )

    def compute_spring_force(self, y):
        nebm_cartesian.compute_spring_force(self.spring_force, y,
                                            self.tangents,
                                            self.k, self.n_images,
                                            self.n_dofs_image,
                                            self._material_int,
                                            self.n_dofs_image_material
                                            )

    def nebm_step(self, y):

        # The convergence of the algorithm depends on how we redefine the
        # angles: Redefining the tangents and spring force helps a little
        self.compute_effective_field_and_energy(y)
        self.compute_tangents(y)
        self.compute_spring_force(y)

        nebm_clib.compute_effective_force(self.G,
                                          self.tangents,
                                          self.gradientE,
                                          self.spring_force,
                                          self.climbing_image,
                                          self.n_images,
                                          self.n_dofs_image
                                          )

    # -------------------------------------------------------------------------
    # Methods -----------------------------------------------------------------
    # -------------------------------------------------------------------------

    def compute_distances(self, A, B):
        """
        Compute the distance between corresponding images of the bands A and B

                A                   B
            [ [image_0]         [ [image_0]
              [image_1]     -     [image_1]
              ..self.               ...
            ]                     ]

            We discard sites without material

        """

        A_minus_B = A - B

        A_minus_B.shape = (-1, self.n_dofs_image)

        A_minus_B = np.apply_along_axis(
            lambda y: compute_norm(y[self._material],
                                   scale=True
                                   ),
            axis=1,
            arr=A_minus_B
            )

        return A_minus_B.reshape(-1)

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
        nebm_cartesian.compute_dYdt(
            y, self.G, ydot, self.sim._pins, self.n_images, self.n_dofs_image)

        # The effective force at the extreme images should already be zero, but
        # we will manually remove any value
        ydot[:self.n_dofs_image] = 0
        ydot[-self.n_dofs_image:] = 0

        return 0
