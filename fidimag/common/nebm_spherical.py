from __future__ import print_function
from __future__ import division
import numpy as np

import fidimag.extensions.nebm_spherical_clib as nebm_clib
from .nebm_tools import spherical2cartesian, cartesian2spherical, compute_norm
from .nebm_tools import linear_interpolation_spherical

from .nebm_base import NEBMBase

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="fidimag")


class NEBM_Spherical(NEBMBase):
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
    micromagnetic simulations. The NEBM in Spherical coordinates describes the
    spins / magnetisation vectors using the polar and azimuthal angles in
    spherical coordinates

        (m_x, m_y, m_z) - > (theta, phi)

    with m_x = sin(theta)cos(phi), m_y=... etc.

    The NEBM is based on the definition of a so called band, which is just a
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

                             dY
                             --  =  G
                             dt

    for every image Y of the energy band. It is worth noticing that the
    spin/magnetisation length is implicitly fixed, thus we do not need to
    correct it in the minimisation equation. The G vector is the effective
    force on the image, which is defined in terms of tangent vectors (tangent
    along the band) in the energy landscape, whose perpendicular components
    follow directions of largest energy changes to find the minimum energy
    transition.  Additionally, the effective force includes a spring force that
    keeps images equally spaced along the band to avoid clustering around
    minima or saddle points. The spacing between images needs a definition of
    DISTANCE. In this code we use an Euclidean distance, normalised by the
    number of degrees of freedom, which is the sum of all spin components, so
    if we have P spins in the system, the number of dofs is 2 * P, i.e. the 2
    spherical angles per spin.  The distance is defined as:

                                ______________________________________________
                               /  P
                      1       /  __                    2                      2
   dist(Y_i, Y_j) =  ---     /  \   [ t(i,a) - t(j,a) ]  + [ p(i,a) - p(j,a) ]
                     2*P \  /   /__
                          \/    a=1

    where t(i,a) is the theta (polar) angle of the a-th spin in the image Y_i,
    and p refers to the phi (azimuthal) angle. Notice that in spherical
    coordinates, the azimuthal angle gets undefined as the polar angle is
    closer to the poles (theta=0, PI), hence it is recommended to avoid
    defining systems where the large majority of spins prefer to point in the z
    direction. In addition, the distances should be always in the range [0,PI],
    thus we do some scaling on the phi angles (which are in the [0,2*PI] range)
    when computing a difference between two spins or when computing the
    distance. This way the differences for both angles have the same range.

    For more details about the definition of the forces involved in the NEBM,
    see the following papers:

        - Dittrich et al., JMMM 250 (2002) L12–L19
        - Henkelman et al., Journal of Chemical Physics 113, 22 (2000)
        - Bessarab et al., Computer Physics Communications 196 (2015) 335-347

    """

    def __init__(self, sim,
                 initial_images, interpolations=None,
                 spring_constant=1e5,
                 name='unnamed',
                 openmp=False
                 ):

        super(NEBM_Spherical, self).__init__(sim,
                                             initial_images,
                                             interpolations=interpolations,
                                             spring_constant=spring_constant,
                                             name=name,
                                             dof=2,
                                             openmp=openmp
                                             )

        # Since we use Spherical coordinates for the energy band (dof=2), when
        # saving files we convert the coordinates to Cartesian
        self.files_convert_f = spherical2cartesian

        # Initialisation ------------------------------------------------------
        # See the NEBMBase class for details

        self.generate_initial_band()

        self.initialise_energies()

        self.initialise_integrator()

        self.create_tablewriter()

        # ---------------------------------------------------------------------

    def initialise_energies(self):
        # Energy of the images
        self.band = self.band.reshape(self.n_images, -1)
        for i in range(self.n_images):
            self.sim.set_m(spherical2cartesian(self.band[i]))
            self.sim.compute_effective_field(t=0)
            self.energies[i] = self.sim.compute_energy()
        self.band = self.band.reshape(-1)

    def generate_initial_band(self):
        """

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
            self.band[i_initial_images[i]] = cartesian2spherical(self.sim.spin)

            self.sim.set_m(self.initial_images[i + 1])
            self.band[i_initial_images[i + 1]] = cartesian2spherical(self.sim.spin)

            # interpolation is an array with *self.interpolations[i]* rows
            # We copy these rows to the corresponding images in the energy
            # band array
            if self.interpolations[i] != 0:
                interpolation = linear_interpolation_spherical(
                    self.band[i_initial_images[i]],
                    self.band[i_initial_images[i + 1]],
                    self.interpolations[i],
                    self.sim.pins
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

            self.sim.set_m(spherical2cartesian(y[i]))
            # elif self.coordinates == 'Cartesian':
            #     self.sim.set_m(self.band[i])

            self.sim.compute_effective_field(t=0)

            self.gradientE[i][:] = energygradient2spherical(self.sim.field,
                                                            y[i]
                                                            )
            # elif self.coordinates == 'Cartesian':
            #     self.H_eff[i][:] = self.sim.spin

            self.energies[i] = self.sim.compute_energy()

        y = y.reshape(-1)
        self.gradientE = self.gradientE.reshape(-1)

    def compute_tangents(self, y):
        nebm_clib.compute_tangents(self.tangents, y, self.energies,
                                   self.n_dofs_image, self.n_images
                                   )

    def compute_spring_force(self, y):
        nebm_clib.compute_spring_force(self.spring_force, y, self.tangents,
                                       self.k, self.n_images, self.n_dofs_image,
                                       self._material_int,
                                       self.n_dofs_image_material
                                       )
        nebm_clib.normalise_images(self.tangents,
                                   self.n_images, self.n_dofs_image
                                   )

    def nebm_step(self, y):

        # The convergence of the algorithm depends on how we redefine the
        # angles: Redefining the tangents and spring force helps a little
        correct_angles(y)
        self.compute_effective_field_and_energy(y)
        self.compute_tangents(y)
        self.compute_spring_force(y)

        nebm_clib.compute_effective_force(self.G,
                                          self.tangents,
                                          self.gradientE,
                                          self.spring_force,
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

        """

        A_minus_B = A - B
        redefine_angles(A_minus_B)

        A_minus_B.shape = (-1, self.n_dofs_image)
        A_minus_B = np.apply_along_axis(
            lambda y: compute_norm(y[self._material],
                                   scale=self.n_dofs_image),
            axis=1,
            arr=A_minus_B
            )

        return A_minus_B.reshape(-1)


# -----------------------------------------------------------------------------

def correct_angles(A):
    """
    Correct THETA and PHI angles
    """
    A[::2][A[::2] > np.pi] = 2 * np.pi - A[::2][A[::2] > np.pi]
    A[::2][A[::2] < 0] = np.abs(A[::2][A[::2] < 0])
    # A[::2][A[::2] > np.pi] = np.pi
    # A[::2][A[::2] < 0] = 0

    A[1::2][A[1::2] > np.pi] = A[1::2][A[1::2] > np.pi] - 2 * np.pi
    A[1::2][A[1::2] < -np.pi] = 2 * np.pi + A[1::2][A[1::2] < -np.pi]


def redefine_angles(A):
    """
    This redefines PHI angles to lie in the [-PI, PI] range
    """
    A[1::2][A[1::2] > np.pi] = 2 * np.pi - A[1::2][A[1::2] > np.pi]
    A[1::2][A[1::2] < -np.pi] = 2 * np.pi + A[1::2][A[1::2] < -np.pi]


def energygradient2spherical(Hxyz, spin):
    """

    Transform the gradient of the energy (with respect to the
    magnetisation) in Cartesian coordinates, into the gradient in spherical
    coordinates (r, t, p) using the transformation matrix:

    | sin t cos p  | sin t sin p | cos t  | | dE/dm_x |   | dE/dm_r |
    | cos t cos p  | cos t sin p | -sin t | | dE/dm_y | = | dE/dm_t |
    | -sin p sin t | cos p sin t |   0    | | dE/dm_z |   | dE/dm_p |

    Notice that the gradient is the negative of the effective field, i.e.
    dE/dm_i = - Heff_i , thus we can use the effective field components in
    Cartesian coordinates as input instead of computing the gradient
    numerically

    This formula can be derived from the chain rule applied to the energy
    derivative to obtain the effective field: d E / d M, with M = M_s (sin
    t cos p, sin t sin p, cos t)

    (see Suss et al., Magnetics, IEEE Transactions 36 (5) pp.3282-3284,
    2000)

    The function only returns the (t, p) = (theta, phi) coordinates of
    dE/dm since we asume that the r component is fixed


    INPUTS:

    Hxyz    :: Effective field in Cartesian coordinates
    spin    :: The magnetisation/spin field in Spherical coordinates

    """
    Hxyz = Hxyz.reshape(-1, 3)

    # theta_phi = cartesian2spherical(spin)
    # t, p = theta_phi[::2], theta_phi[1::2]
    t, p = spin[::2], spin[1::2]

    # The gradient is just a vector field, thus we use the spin field size
    # in spherical coordinates
    gradientE = np.zeros_like(spin).reshape(-1, 2)

    # Theta components
    gradientE[:, 0] = (np.cos(t) * np.cos(p) * (-Hxyz[:, 0]) +
                       np.cos(t) * np.sin(p) * (-Hxyz[:, 1]) +
                       np.sin(t) * Hxyz[:, 2]
                       )

    # Phi components
    gradientE[:, 1] = (np.sin(p) * Hxyz[:, 0] +
                       np.cos(p) * (-Hxyz[:, 1])
                       )

    Hxyz = Hxyz.reshape(-1)

    return gradientE.reshape(-1)
