from __future__ import print_function
from __future__ import division
import numpy as np

import fidimag.extensions.cvode as cvode
import fidimag.extensions.neb_method_clib as nebm_clib
from fidimag.common.vtk import VTK

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="fidimag")


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
                 coordinates='Spherical', k=1e5
                 ):

        self.coordinates = coordinates
        # Degrees of Freedom per spin
        if self.coordinates == 'Spherical':
            self.dof = 2
        # elif self.coordinates == 'Cartesian':
        #     self.dof = 3

        self.sim = sim
        self.mesh = self.sim.mesh

        # Number of spins in the system
        self.n_spins = len(self.mesh)

        self.k = k

        # Initial states ------------------------------------------------------

        # We assume the extremes are fixed
        self.initial_images = initial_images

        if self.interpolations:
            self.interpolations = interpolations
        else:
            self.interpolations = [0 for i in len(initial_images) - 1]

        # Number of images with/without the extremes
        self.n_images = len(self.initial_images) + np.sum(self.interpolations)
        self.n_images_inner_band = self.n_images - 2

        # Number of degrees of freedom per image
        self.n_dofs_image = (self.dof * self.n_spins)

        # Total number of degrees of freedom in the NEBM band
        self.n_band = self.n_images * self.n_spins_image

        # NEBM Arrays ---------------------------------------------------------
        # We will initialise every array using the total number of images,
        # but we must have in mind that the images at the extrema of the band
        # are kept fixed, so we do not compute their gradient, tangents, etc.
        # This might be not memory efficient but the code is understood better
        # when we perform the loops when calculating the effective fields
        # and forces

        # The array containing every degree of freedom
        self.band = np.zeros(self.n_band)

        # The gradient with respect to the magnetisation (effective field)
        self.H_eff = np.zeros_like(self.band)

        # The effective force
        self.G = np.zeros_like(self.band)

        self.tangents = np.zeros_like(self.band)
        self.energies = np.zeros(self.n_images)
        self.spring_force = np.zeros_like(self.band)
        self.distances = np.zeros(self.n_images - 1)

        # Initialisation ------------------------------------------------------

        self.generate_initial_band()

        # Energy of the images at the extremes, which are kept fixed
        self.band = self.band.reshape(self.n_images, -1)
        for i in [0, self.n_images - 1]:
            self.sim.set_m(self.band[i])
            self.energy[i] = self.sim.compute_energy()

        self.band = self.band.reshape(-1)

        self.initialise_integrator()

        # ---------------------------------------------------------------------

    def initialise_integrator(self, rtol=1e-6, atol=1e-6):
        self.t = 0
        self.iterations = 0
        self.ode_count = 1

        self.integrator = cvode.CvodeSolver(self.coords, self.sundials_rhs)
        self.integrator.set_options(rtol, atol)

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

    def energygradient2spherical(self, Hxyz, spin):
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
        spin    :: The magnetisation/spin field in Cartesian coordinates

        """
        Hxyz = Hxyz.reshape(-1, 3)

        theta_phi = self.cartesian2spherical(spin)
        t, p = theta_phi[::2], theta_phi[1::2]

        # The gradient is just a vector field, thus we use the spin field size
        # in spherical coordinates
        gradientE = np.zeros_like(theta_phi).reshape(-1, 2)

        # Theta components
        gradientE[:, 0] = (np.cos(t) * np.cos(p) * (-Hxyz[:, 0]) +
                           np.cos(t) * np.sin(p) * (-Hxyz[:, 1]) +
                           np.sin(t) * Hxyz[:, 2]
                           )

        # Phi components
        gradientE[:, 1] = np.sin(t) * (np.sin(p) * Hxyz[:, 0] +
                                       np.cos(p) * (-Hxyz[:, 1])
                                       )

        Hxyz = Hxyz.reshape(-1)

        return gradientE.reshape(-1)

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
        self.band = self.band.reshape(self.n_images, -1)

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

    def compute_effective_field_and_energy(self, y):
        """

        Compute the effective field and the energy of every image in the band,
        using the array *y* as the degrees of freedom of the system (i.e. the
        one that contains all the spin directions of the images in the band).

        The local copy of the *y* array for this NEBM class is the self.band
        array, which we update at the end of every call to the integrator in
        the relaxation function

        """

        self.H_eff = self.H_eff.reshape(self.n_images, -1)

        y = y.reshape(self.n_images, -1)

        # Only update the extreme images
        for i in range(1, len(y) - 1):

            self.sim.set_m(self.spherical2cartesian(y[i]))
            # elif self.coordinates == 'Cartesian':
            #     self.sim.set_m(self.band[i])

            self.sim.compute_effective_field(t=0)

            self.gradE[i][:] = self.energygradient2spherical(self.sim.field)
            # elif self.coordinates == 'Cartesian':
            #     self.H_eff[i][:] = self.sim.spin

            self.energy[i] = self.sim.compute_energy()

        y = y.reshape(-1)

    def compute_tangents(self, y):
        nebm_clib.compute_tangents(self.tangents, y, self.energies,
                                   self.n_dofs_image, self.n_images
                                   )

    def compute_spring_force(self):
        pass

    def nebm_step(self, y):

        self.compute_effective_field_and_energy(y)
        self.compute_tangents()
        self.compute_spring_force()

        self.G = (-self.gradE
                  + np.dot(self.gradE, self.tangents) * self.tangents
                  + self.spring_force
                  )

    # -------------------------------------------------------------------------
    # CVODE solver ------------------------------------------------------------
    # -------------------------------------------------------------------------

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

        # The effective force at the extreme images should already be zero, but
        # we will manually remove any value
        ydot[:self.n_dofs_image] = 0
        ydot[-self.n_dofs_image:] = 0

        return 0

    def compute_norm(self, y):
        """

        Compute the norm of the *y* array, which contains spin directions in
        Spherical coordinates,

        y = [ y_theta0 y_phi0 y_theta1 y_phi1 ... y_thetaN y_phiN]

        If the absolute value of a component is larger than PI, we redefine the
        difference to be smaller than PI

        We scale the norm by the array size

        """

        y[y > np.pi] = 2 * np.pi - y[y > np.pi]
        y[y < -np.pi] = 2 * np.pi + y[y < -np.pi]
        y = np.sqrt(np.sum(y ** 2.)) / len(y)
        return y

    def compute_maximum_dYdt(self, y, dt):
        """

        Compute the maximum difference from the images of the *y* array with
        the images of the *self.band* array, divided by the last time step
        size, i.e.  (t - self.t)

        The differences are not computed for the images at the extremes, since
        these images are fixed and do not change with the iterator

        For instance, in spherical coordinates, if we have a band of (N + 1)
        images, labeled from 0 to N, we start by

        dY = [y1_theta0 y1_phi0 y1_theta1 ... y(N-1)_theta0 y(N-1)_phi0 y(N-1)_theta1 ... ]
            - [G1_theta0 G1_phi0 G1_theta1 ... G(N-1)_theta0 G(N-1)_phi0 G(N-1)_theta1 ... ]

        where y(i)_theta(j) is the theta componenet of the j-th spin of the
        i-th image in the band

        Then we calculate the norm of every difference:

        ||dY|| =  [ || dY1_theta0 dY1_phi0 dY1_theta1         ...       ||
                                   ...
                    || dY(N-1)_theta0 dY(N-1)_phi0 dY(N-1)_theta1  ...  ||
                  ]

        we divide by dt:

        ||dY|| -- > || dY || / dt = dYdt

        and finally take the maximum value of dYdt

        """

        # We will not consider the images at the extremes to compute dY
        # Since we removed the extremes, we only have *n_images_inner_band*
        # images
        band_no_extremes = slice(self.n_dofs_image, -self.n_dofs_image)
        dYdt = (y[band_no_extremes] -
                self.band[band_no_extremes]
                ).reshape(self.n_images_inner_band, -1)

        # Compute the norm of the difference dY for every image in the array
        dYdt = np.apply_along_axis(self.compute_norm, dYdt, axis=1)
        dYdt /= dt

        return np.max(dYdt)

    def relax(self, dt=1e-8, stopping_dYdt=1, max_iterations=1000,
              save_npys_every=100, save_vtks_every=100
              ):

        """

        """

        log.debug("Relaxation parameters: "
                  "stopping_dmdt={} (degrees per nanosecond), "
                  "time_step={} s, max_iterations={}.".format(stopping_dYdt,
                                                              dt,
                                                              max_iterations))

        for i in range(max_iterations):

            # Update the iterations number counter
            self.iterations += 1

            # Get the current size of the time discretisation from the
            # integrator (variable step size)
            cvode_dt = self.integrator.get_current_step()

            # If the step size of the integrator is larger than the specified
            # discretisation, use the current integrator step size to
            # compute the next iteration. Otherwise, just stick to the
            # specified step
            if cvode_dt > dt:
                increment_dt = cvode_dt
            else:
                increment_dt = dt

            self.integrator.run_until(self.t + increment_dt)

            # Compute the maximum change in the integrator step
            max_dYdt = self.compute_maximum_dYdt(self.integrator.y,
                                                 increment_dt)

            # Copy the updated energy band to our local array
            self.band[:] = self.integrator.y[:]

            # Stop criteria:
            if max_dYdt < stopping_dYdt:
                break

            # Update the current step
            self.t = self.t + increment_dt

        log.info("Relaxation finished at time step = {:.4g}, "
                 "t = {:.2g}, call rhs = {:.4g} "
                 "and max_dYdt = {:.3g}".format(self.step,
                                                self.t,
                                                self.ode_count,
                                                max_dYdt)
                 )
