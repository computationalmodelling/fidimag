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

    NEB Method equations:

        G =

    ARGUMENTS

    sim                 :: An instance of a micromagnetic or an atomistic
                           simulation

    initial_images      :: A sequence of arrays or functions to set up the
                           magnetisation field profile

    inteprolations      ::

    """
    def __init__(self, sim,
                 initial_images, interpolations=None,
                 k=1e5,
                 name='unnamed',
                 openmp=False
                 ):

        super(NEBM_Spherical, self).__init__(sim,
                                             initial_images,
                                             interpolations=interpolations,
                                             k=k,
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
    Correct PI and PHI angles
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
