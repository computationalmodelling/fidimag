import numpy as np
import fidimag.extensions.nebm_clib as nebm_clib
import fidimag.common.constant as const

from .chain_method_tools import spherical2cartesian, cartesian2spherical, compute_norm
from .chain_method_tools import linear_interpolation_spherical
from .chain_method_tools import interpolation_Rodrigues_rotation
from .chain_method_tools import m_to_zero_nomaterial

from .chain_method_base import ChainMethodBase

import scipy.interpolate as si

import logging
log = logging.getLogger(name="fidimag")


class StringMethod(ChainMethodBase):
    """
    The String Method class, to find minimum energy paths between two stable
    states in a given magnetic system.

    Unlike the NEBM, the images are not kept apart by a spring force: they
    are evolved by the energy gradient alone and then redistributed along
    the string by spline interpolation, in ``run_until``.

    Parameters
    ----------
    sim
        An instance of a micromagnetic or an atomistic simulation.
        Every image in the band will be a copy of this simulation.
    initial_images
        A list containing numpy arrays or space dependent functions
        to set the magnetisation fields for every image in the
        band. It is suggested that the first and last elements
        define stable states of the magnetic system. The arrays and
        functions are used to load the magnetisation/spin fields
        through the sim.set_m method from the Simulation object.
    interpolations
        A list with 1 element less than the initial_images list,
        where every entry is an integer indicating the number of
        interpolations between consecutive images. For example, if
        we defined initial_images as [state_1, state_2, state_3],
        and we want 10 interpolations between state_1 and state_2
        and 5 interpolations between state_2 and state_3, we set
        interpolations as [10, 5], making an energy band of 17
        images. If we do not want any interpolation, we leave this
        list as None or empty.
    interpolation_method
        If a number of interpolations was defined, how those
        interpolations are computed:

        - ``'linear'``: a linear interpolation of the spin directions
          using spherical coordinates.
        - ``'rotation'``: an interpolation of the spin directions using
          Rodrigues' rotation formula.
    name
        The simulation name. Folders for VTK and NPY files, and
        data tables are named according to this string.
    openmp
        Set this as True to use the parallelised version of CVODE,
        which is the integrator used to evolve the minimisation
        equation.
    integrator
        The integrator evolving the string, passed on to
        ``initialise_integrator``. Defaults to ``'verlet'``, since the
        variable step integrator from Sundials does not work well with the
        String Method (see ``Sundials_RHS``).
    """

    def __init__(self, sim,
                 initial_images,
                 interpolations=None,
                 interpolation_method='rotation',
                 name='unnamed',
                 openmp=False,
                 integrator='verlet'  # or scipy
                 ):

        super().__init__(sim,
                                           initial_images,
                                           interpolations=interpolations,
                                           name=name,
                                           dof=3,
                                           openmp=openmp
                                           )

        # Initialisation ------------------------------------------------------
        # See the ChainMethodBase class for details

        self.generate_initial_band(method=interpolation_method)

        self.initialise_energies()

        self.initialise_integrator(integrator=integrator)

        self.create_tablewriter()

        # ---------------------------------------------------------------------

    def initialise_energies(self):
        """
        Populate the ``self.energies`` array with the energy of every image
        of the band, for the 0th step of the algorithm.
        """
        # Energy of the images
        self.band = self.band.reshape(self.n_images, -1)
        for i in range(self.n_images):
            self.sim.set_m(self.band[i])
            self.sim.compute_effective_field(t=0)
            self.energies[i] = self.sim.compute_energy()
        self.band = self.band.reshape(-1)

    def generate_initial_band(self, method='linear'):
        """
        Generate the initial string from the initial images, interpolating
        between them as requested by the ``interpolations`` list.

        Parameters
        ----------
        method
            How the interpolations between consecutive initial images are
            computed. Either ``'linear'``, for a linear interpolation of the
            spin directions in spherical coordinates, or ``'rotation'``, for
            an interpolation using Rodrigues' rotation formula.
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
        Compute the effective field and the energy of every image in the
        band.

        Parameters
        ----------
        y
            The degrees of freedom of the system, i.e. the array containing
            all the spin directions of the images in the band. The local
            copy of this array for this String class is ``self.band``, which
            is updated at the end of every call to the integrator in the
            relaxation function.
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

    def string_method_step(self, y):
        """
        Update the effective field and the energies from the band *y*, and
        set the effective force ``self.G`` from the energy gradient alone,
        i.e. a steepest descent step.

        Parameters
        ----------
        y
            The degrees of freedom of the system, i.e. the array containing
            all the spin directions of the images in the band.
        """
        # Use the projection mehtod from the NEBM C libs
        self.compute_effective_field_and_energy(y)
        nebm_clib.project_images(self.gradientE, y,
                                 self.n_images, self.n_dofs_image
                                 )

        # The string method only requires to update the images using the
        # gradient (steepest descent)

        # Is it necessary to rescale the effective field?
        # scale = np.tile(np.repeat(const.mu_0 * self.sim.Ms, 3), self.n_images)
        self.G[:] = -self.gradientE[:]
        # print(self.G)

    # -------------------------------------------------------------------------
    # Methods -----------------------------------------------------------------
    # -------------------------------------------------------------------------

    def compute_distances(self):
        """
        Compute the distances between consecutive images of the string, and
        store them in ``self.distances`` and ``self.path_distances``.

        Notes
        -----
        The path distances are rescaled to lie between 0 and 1. The
        distances are computed with the Geodesic library, whereas the
        original String Method uses a Euclidean norm.
        """

        nebm_clib.image_distances_GreatCircle(self.distances,
                                              self.path_distances,
                                              self.band,
                                              self.n_images,
                                              self.n_dofs_image,
                                              self._material_int,
                                              self.n_dofs_image_material
                                              )

        self.path_distances[1:] = np.cumsum(self.distances)
        self.distances /= self.path_distances[-1]
        self.path_distances[1:] /= self.path_distances[-1]

    def step_RHS(self, t, y):
        """
        Right hand side of the String Method equation, as called on every
        iteration of the step integrators in
        ``chain_method_integrators.py``.

        The effective force is left in ``self.G``, with the entries of the
        images at the extremes of the string set to zero, since those images
        are kept fixed.

        Parameters
        ----------
        t
            Current time of the integrator. Unused, since the equation does
            not depend explicitly on time.
        y
            The band at which the right hand side is evaluated.

        Returns
        -------
        int
            Always 0, as expected by the integrators.
        """

        self.ode_count += 1

        # Update the effective field, energies, spring forces and tangents
        # using the *y* array
        self.string_method_step(y)

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
        Right hand side of the String Method equation, as called on every
        iteration of the CVODE integrator.

        Parameters
        ----------
        t
            Current time of the integrator. Unused, since the equation does
            not depend explicitly on time.
        y
            The band at which the right hand side is evaluated.
        ydot
            Output array where the right hand side is stored, since we are
            solving ``dy/dt = 0``. Its entries for the images at the
            extremes of the string are set to zero, since those images are
            kept fixed.

        Returns
        -------
        int
            Always 0, as expected by the integrator.

        Warnings
        --------
        The variable step integrator from Sundials does not work well with
        the String Method: it makes the algorithm overshoot the solutions
        for large time steps, driving the images towards the images at the
        extremes. This could possibly be fixed by redefining the positions
        of the images after every integrator step, instead of after a
        certain number of steps, but that requires tuning the Sundials
        Python wrapper, and the stopping criteria of the algorithm would
        have to be checked as well.
        """

        self.ode_count += 1

        # Update the effective field, energies, spring forces and tangents
        # using the *y* array
        self.string_method_step(y)

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

    def run_until(self, t):
        """
        Evolve the string until the time *t*, then redistribute the images
        along it.

        Parameters
        ----------
        t
            Time up to which the integrator is run. If it is not larger than
            the current ``self.t``, nothing is done.

        Returns
        -------
        float or None
            The maximum rate of change of the string, as computed by
            ``compute_maximum_dYdt``, or ``None`` if the integrator was not
            run.

        Notes
        -----
        After the number of integration steps given by
        ``self.integrator.run_until(t)``, the positions of the images are
        redefined using splines, which are computed with Scipy: not the most
        efficient method for now, but a very accurate one. The image
        positions are normalised in the distance calculation method.
        """

        if (t) <= self.t:
            return

        self.integrator.run_until(t)

        band = self.integrator.y
        # Compute distances and SCALE
        self.compute_distances()

        new_dist = np.linspace(self.path_distances[0],
                               self.path_distances[-1],
                               self.path_distances.shape[0]
                               )
        # print(self.path_distances)
        # print(new_dist)

        # Restructure the string by interpolating every spin component
        # print(self.integrator.y[self.n_dofs_image:self.n_dofs_image + 10])
        band = self.integrator.y.reshape(self.n_images, self.n_dofs_image)
        for i in range(self.n_dofs_image):

            cs = si.CubicSpline(self.path_distances, band[:, i])
            band[:, i] = cs(new_dist)

        # print('LAST', self.integrator.y[self.n_dofs_image:self.n_dofs_image + 10])

        # Copy the updated energy band to our local array
        self.band[:] = self.integrator.y[:]

        # Update fields with the new relocation of images
        # self.compute_effective_field_and_energy(self.band)
        self.string_method_step(self.band)

        # Compute the maximum change in the integrator step
        max_dYdt = self.compute_maximum_dYdt(self.integrator.y, self.last_Y,
                                             t - self.t)

        self.last_Y[:] = self.band[:]

        # Update the current step
        self.t = t

        return max_dYdt

    def compute_tangents(self, y):
        """
        Compute the tangents to the string, using the Geodesic NEBM
        definition, so that the interpolation methods for the energy band
        can be used.

        Parameters
        ----------
        y
            The degrees of freedom of the system, i.e. the array containing
            all the spin directions of the images in the band.

        Notes
        -----
        An alternative approximation of the band could be made with simple
        spline interpolations, although the GNEBM interpolation carries more
        information about the gradients.
        """
        nebm_clib.compute_tangents(self.tangents, y, self.energies,
                                   self.n_dofs_image, self.n_images
                                   )
        nebm_clib.project_images(self.tangents, y,
                                 self.n_images, self.n_dofs_image
                                 )
        nebm_clib.normalise_images(self.tangents,
                                   self.n_images, self.n_dofs_image
                                   )
