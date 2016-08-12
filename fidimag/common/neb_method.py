from __future__ import print_function
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
                 initial_images, interpolations=None
                 ):

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

        # The array containing every degree of freedom
        self.band = np.zeros(self.n_band * (3 * self.n_spins))

        # The gradient
        self.G = np.zeros_like(self.images)

        self.tangents = np.zeros(self.n_inner_band * (3 * self.n_spins))

        # ---------------------------------------------------------------------

    def cartesian2spherical(self, y_cartesian):
        pass

    def spherical2cartesian(self, y_spherical):
        pass

    def linear_interpolation(self, y_initial, y_final, n):
        interpolations = np.zeros((self.n, 3 * self.n_spins))
        pass

    def generate_initial_state_by_interpolation(self):

        self.band = self.band.reshape(-1, 3 * self.n_spins)

        i_initial_images = [0]
        for i in range(1, len(self.initial_images)):
            i_initial_images.append(i + i_initial_images[i - 1])

        for i, m_field in enumerate(self.initial_images[:-1]):
            self.sim.set_m(self.initial_images[i])
            self.band[i_initial_images[i]] = self.sim.spin[:]

            self.sim.set_m(self.initial_images[i + 1])
            self.band[i_initial_images[i + 1]] = self.sim.spin[:]

            # interpolation is a (n_interpolations, n_spins) array
            interpolation = self.linear_interpolation(
                self.band[i_initial_images[i]],
                self.band[i_initial_images[i + 1]],
                self.interpolations[i]
                )

            # We then set the interpolated spins fields at once
            self.band[i_initial_images[i] + 1:
                      i_initial_images[i + 1]] = interpolation

