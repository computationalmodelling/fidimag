from __future__ import print_function
from __future__ import division
import numpy as np


class NEBMethod(object):
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
                 coordinates='Spherical'
                 ):

        self.coordinates = coordinates
        # Degrees of Freedom per spin
        if self.coordinates == 'Spherical':
            self.dof = 2
        # elif self.coordinates == 'Cartesian':
        #     self.dof = 3

        self.sim = sim
        self.mesh = self.sim.mesh
        # Number of spins
        self.n_spins = len(self.mesh)

        # Initial states ------------------------------------------------------

        # We assume the extremes are fixed
        self.initial_images = initial_images

        if self.interpolations:
            self.interpolations = interpolations
        else:
            self.interpolations = [0 for i in len(initial_images) - 1]

        # Number of images with/without the extremes
        self.n_band = len(self.images) + len(self.interpolations)
        self.n_inner_band = self.n_band - 2

        # NEBM Arrays ---------------------------------------------------------
        # We will initialise every array using the total number of images,
        # but we must have in mind that the images at the extrema of the band
        # are kept fixed, so we do not compute their gradient, tangents, etc.
        # This might be not memory optimised but the code is understood better
        # when we perform the loops when calculating the effective fields
        # and forces

        # The array containing every degree of freedom
        self.band = np.zeros(self.n_band * (self.dof * self.n_spins))

        # The gradient with respect to the magnetisation (effective field)
        self.H_eff = np.zeros_like(self.band)

        # The effective force
        self.G = np.zeros_like(self.band)

        self.tangents = np.zeros_like(self.band)
        self.energies = np.zeros(self.n_band)
        self.spring_force = np.zeros_like(self.band)

        # Initialisation ------------------------------------------------------

        self.generate_initial_band()

        # Energy of the images at the extremes, which are kept fixed
        self.band = self.band.reshape(-1, 3 * self.n_spins)
        for i in [0, self.n_band - 1]:
            self.sim.set_m(self.band[i])
            self.energy[i] = self.sim.compute_energy()

        self.band = self.band.reshape(-1)

        # ---------------------------------------------------------------------

    def cartesian2spherical(self, y_cartesian):
        """
        y_cartesian     :: [y_x0 y_y0 y_z0 y_x1 y_y1 ...]
        """
        theta_phi = np.zeros((len(y_cartesian.reshape(-1, 3)), 2))

        # r = sqrt (m_x ** 2 + m_y ** 2)
        # r = np.sqrt(y_cartesian[::3] ** 2 + y_cartesian[1::3])

        theta_phi[:, 0] = np.arccos(y_cartesian[2::3])  # theta
        theta_phi[:, 1] = np.arctan2(y_cartesian[1::3],
                                     y_cartesian[::3]
                                     )                  # phi

        return theta_phi.reshape(-1)

    def spherical2cartesian(self, y_spherical):
        y_cartesian = np.zeros((len(y_spherical.reshape(-1, 2)), 3))

        theta, phi = y_spherical[::2], y_spherical[1::2]
        y_cartesian[:, 0] = np.sin(theta) * np.cos(phi)
        y_cartesian[:, 1] = np.sin(theta) * np.sin(phi)
        y_cartesian[:, 2] = np.cos(theta)

        return y_cartesian.reshape(-1)

    def effectivefield2spherical(self):
        pass

    def linear_interpolation(self, y_initial, y_final, n):
        """
        y_initial, y_final      :: In Spherical coordinates
        """
        interpolations = np.tile(y_initial, (n, 1))

        # We will not interpolate pinned spins
        _filter = self.sim.pins == 0

        # y_initial_spherical = self.cartesian2spherical(y_initial)
        # y_final_spherical = self.cartesian2spherical(y_final)

        # dy_spherical = ((y_final_spherical - y_initial_spherical) / (n + 1))
        dy = (y_final - y_initial) / (n + 1)

        for i in range(1, n + 1):
            interpolations[i - 1][_filter] = (y_initial + i * dy)[_filter]

        return interpolations

    def generate_initial_band(self):

        # Every row will be an image of the band, i.e. the i-th row is
        # the i-t image
        self.band = self.band.reshape(-1, self.dof * self.n_spins)

        # Indexes indicating the image number (position of the images) in the
        # band, for the specified initial images
        i_initial_images = [0]
        for i in range(1, len(self.initial_images)):
            i_initial_images.append(i + i_initial_images[i - 1])

        for i, m_field in enumerate(self.initial_images[:-1]):

            # Copy the magnetisation field from the i-th and (i + 1)-th to the
            # corresponding rows of the nebm band array To do this, we need to
            # know in which positions are these images in the band, which
            # change according to the number of interpolations. Accordingly,
            # we use the list with the indexes of the initial images
            self.sim.set_m(self.initial_images[i])
            self.band[i_initial_images[i]] = self.cartesian2spherical(self.sim.spin)

            self.sim.set_m(self.initial_images[i + 1])
            self.band[i_initial_images[i + 1]] = self.cartesian2spherical(self.sim.spin)

            # interpolation is an array with *self.interpolations[i]* rows
            # We copy these rows to the corresponding images in the energy
            # band array
            if self.interpolations[i] != 0:
                interpolation = self.linear_interpolation(
                    self.band[i_initial_images[i]],
                    self.band[i_initial_images[i + 1]],
                    self.interpolations[i]
                    )

                # We then set the interpolated spins fields at once
                self.band[i_initial_images[i] + 1:
                          i_initial_images[i + 1]] = interpolation

        # expand the energy band array
        self.band = self.band.reshape(-1)

    def compute_effective_field_and_energy(self):
        """
        """

        self.H_eff = self.H_eff.reshape(self.n_band, -1)
        self.band = self.band.reshape(self.n_band, -1)

        # Only update the extreme images
        for i in range(1, len(self.band) - 1):

            self.sim.set_m(self.spherical2cartesian(self.band[i]))
            # elif self.coordinates == 'Cartesian':
            #     self.sim.set_m(self.band[i])

            self.sim.compute_effective_field(t=0)

            self.H_eff[i][:] = self.effectivefield2spherical(self.sim.field)
            # elif self.coordinates == 'Cartesian':
            #     self.H_eff[i][:] = self.sim.spin

            self.energy[i] = self.sim.compute_energy()

    def compute_tangents(self):
        # neb_clib.compute_tangents(self.band, self.effectie_field,
        #                           self.n_band, self.dof
        #                           )
        pass

    def compute_spring_force(self):
        pass

    def nebm_step(self):

        self.compute_effective_field_and_energy()
        self.compute_tangents()
        self.compute_spring_force()

        self.G = (self.H_eff
                  - np.dot(self.H_eff, self.tangents) * self.tangents
                  + self.spring_force
                  )
