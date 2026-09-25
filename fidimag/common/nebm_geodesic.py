import numpy as np
from scipy.interpolate import CubicHermiteSpline
import fidimag.extensions.nebm_clib as nebm_clib

from .chain_method_tools import spherical2cartesian, cartesian2spherical, compute_norm
from .chain_method_tools import linear_interpolation_spherical
from .chain_method_tools import interpolation_Rodrigues_rotation
from .chain_method_tools import m_to_zero_nomaterial

from .chain_method_base import ChainMethodBase

import logging
log = logging.getLogger(name="fidimag")


class NEBM_Geodesic(ChainMethodBase):
    r"""
    The NEB Method (NEBM) class to find minimum energy paths between two
    stable states in a given magnetic system.

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
    spring_constant
        The spring constant magnitude
    name
        The NEBM simulation name. Folders for VTK and NPY files,
        and data tables are named according to this string.
    climbing_image
        Any iterable with the indexes of the climbing images, which are
        driven towards a saddle point. A negative index sets a falling
        image instead.
    openmp
        Set this as True to use the parallelised version of CVODE,
        which is the integrator used to evolve the NEBM
        minimisation equation.
    integrator
        The integrator evolving the band, passed on to
        ``initialise_integrator``.

    Notes
    -----
    This class works both for atomistic and micromagnetic simulations. The
    Geodesic NEBM describes the spins / magnetisation vectors in Cartesian
    coordinates (m_x, m_y, m_z) and distances in the energy landscape are
    measured by a geodesic distance (see below). The NEBM is based on the definition of a so called band, which is
    just a sequence of replicas (in terms of geometry and magnetic parameters)
    of the same system (given by our simulation), called images, in different
    magnetic configurations, and where the first and last state are the stable
    states used to find a minimum energy transition between them. Calling the
    images as Y_i, after relaxation an energy band of N+1 images usually looks
    like::

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
    equation::

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
    distance is defined as::

                                      ___________________
                                     / P
                                    /  __             2
        distance(Y_i, Y_j) =       /  \   ( L(i,j)_a )
                               \  /   /__
                                \/    a=1

    where::

                                 ->       ->        ->       ->
         L(i,j)_a  =  arctan2( | m(i)_a x m(j)_a |, m(i)_a o m(j)_a )

    where m(i)_a is the a-th spin of the image Y_i. This is known as Vincenty's
    formula.

    The second term in the minimisation equation is to correct the
    magnetisation/spin vectors length, since at 0 K this length is fixed. For
    now, we use the *c* factor as 6.

    References
    ----------
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
                 integrator='cvode_bdf'  # or scipy
                 ):

        super().__init__(sim,
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
        """
        Populate the ``self.energies`` array with the energy of every image
        of the band, for the 0th step of the algorithm.
        """
        # Energy of the images
        self.band.shape = (self.n_images, -1)
        for i in range(self.n_images):
            self.sim.set_m(self.band[i])
            self.sim.compute_effective_field(t=0)
            self.energies[i] = self.sim.compute_energy()
        self.band.shape = (-1)

    def generate_initial_band(self, method='linear'):
        """
        Generate the initial band from the initial images, interpolating
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
        self.band.shape = (self.n_images, -1)

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
        self.band.shape = (-1)

    def compute_effective_field_and_energy(self, y):
        """
        Compute the effective field and the energy of every image in the
        band.

        Parameters
        ----------
        y
            The degrees of freedom of the system, i.e. the array containing
            all the spin directions of the images in the band. The local
            copy of this array for this NEBM class is ``self.band``, which is
            updated at the end of every call to the integrator in the
            relaxation function.
        """

        self.gradientE.shape = (self.n_images, -1)
        y.shape = (self.n_images, -1)

        # Only update the extreme images
        for i in range(1, len(y) - 1):

            # self.sim.set_m(y[i])  # -> memory leak if sim has a driver
            self.sim.spin[:] = y[i]

            self.sim.compute_effective_field(t=0)

            self.gradientE[i][:] = -self.sim.field

            self.energies[i] = self.sim.compute_energy()

        y.shape = (-1)
        self.gradientE.shape = (-1)

    def compute_tangents(self, y):
        """
        Compute the tangents to the band, project them on the tangent space
        of the spin field and normalise them, storing the result in
        ``self.tangents``.

        Parameters
        ----------
        y
            The degrees of freedom of the system, i.e. the array containing
            all the spin directions of the images in the band.
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

    def compute_spring_force(self, y):
        """
        Compute the spring force keeping the images equally spaced along the
        band.

        Parameters
        ----------
        y
            The degrees of freedom of the system, i.e. the array containing
            all the spin directions of the images in the band.

        Notes
        -----
        The variable spring constant, enabled with ``self.variable_k``, is
        more effective when the band has a saddle point, and seems to work
        when there is a single one. Experimental functionality.

        References
        ----------
        J. Chem. Phys. 113, 9901 (2000)
        """
        if self.variable_k:
            E_max = np.max(self.energies)
            # E_i is the energy used to decide the i-th image's spring
            # constant: the largest of its two *neighbours'* energies
            # (images i-1 and i+1), following the reference above. Do NOT
            # compare against the image's own energy (self.energies[1:-1])
            # -- that was a past bug here: it made the scheme only look at
            # the left neighbour and ignore the right one, breaking the
            # left/right symmetry the method relies on.
            E_i = np.maximum(self.energies[:-2], self.energies[2:])
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

        # Optionally combine the geodesic distance (self.path_distances)
        # with the energy (self.spring_weighting == 'energy', default) or
        # its curvature (self.spring_weighting == 'curvature') to decide
        # the image spacing used in the spring force -- see
        # compute_energy_weighted_spring_lengths and
        # compute_curvature_weighted_spring_lengths for what each of these
        # actually does and does not achieve.
        if self.spring_force_ratio > 0:
            if self.spring_weighting == 'curvature':
                spring_distances = self.compute_curvature_weighted_spring_lengths()
            else:
                spring_distances = self.compute_energy_weighted_spring_lengths()
        else:
            spring_distances = self.distances

        nebm_clib.compute_spring_force(self.spring_force, y,
                                       self.tangents,
                                       k, self.n_images,
                                       self.n_dofs_image,
                                       spring_distances
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

    def compute_energy_weighted_spring_lengths(self, n_interpolations=20):
        """
        Combine the geodesic path distance (``self.path_distances``) and the
        energy (``self.energies``) into a single, energy-weighted arc length
        per band segment, instead of spacing the images purely by geodesic
        distance.

        Parameters
        ----------
        n_interpolations
            Number of sub-segments used to discretise every image-to-image
            segment when integrating the weighted arc length.

        Returns
        -------
        numpy.ndarray
            An array of length ``n_images - 1``, where entry ``i`` is the
            energy-weighted length of the segment between image ``i`` and
            image ``i + 1``. This is the same layout as ``self.distances``,
            so it can be used as a drop-in replacement in
            ``nebm_clib.compute_spring_force``.

        Notes
        -----
        A cubic Hermite spline is built through the points
        ``(path_distances_i, energies_i, dE/d(path_distance)_i)`` and used
        to measure a distance combining both axes, so images bunch up more
        strongly wherever the energy changes fast (e.g. on the flanks
        approaching a saddle point) rather than being purely equidistant
        along the path. ``self.spring_force_ratio``, in ``(0, 1]``, sets the
        weight given to the energy axis, and the rest,
        ``1 - spring_force_ratio``, is given to the path distance axis.

        This weighting by ``dE/d(path_distance)`` is not just a heuristic:
        for a converged MEP the perpendicular component of the energy
        gradient vanishes by construction, so ``dE/d(path_distance)`` there
        equals the full gradient norm ``||grad E||``. Fabian & Shcherbakov
        derive that ``||grad E||`` is exactly the natural weight (a
        Fermat's-principle-like "refractive index") for a
        thermodynamic-action-minimising arc-length parametrisation of the
        path, i.e. this scheme is the discretised, image-spacing analogue of
        their continuum result. It is therefore *not* expected to
        concentrate images at a critical point itself, since ``||grad E||``
        (and ``dE/d(path_distance)``) vanishes there too, at minima and
        maxima alike. See ``compute_curvature_weighted_spring_lengths`` for
        a spacing scheme that does target critical points specifically.

        References
        ----------
        - Bessarab, Uzdin, Jonsson, Computer Physics Communications 196
          (2015) 335-347, for the path-length weighting (also cited in this
          class' docstring).
        - Fabian & Shcherbakov, Geophysical Journal International 215(1),
          314-324 (2018).
        """
        path_distances, energies, spline = self._energy_hermite_spline()

        ratio_E = min(1.0, self.spring_force_ratio)
        ratio_x = 1.0 - ratio_E

        x_fine = self._fine_path_grid(path_distances, n_interpolations)
        e_fine = spline(x_fine)

        range_x = path_distances[-1] - path_distances[0]
        range_E = np.max(e_fine) - np.min(e_fine)
        if range_E < 1e-30:
            range_E = 1.0

        dx = ratio_x * np.diff(x_fine) / range_x
        dE = ratio_E * np.diff(e_fine) / range_E
        sub_lengths = np.sqrt(dx ** 2 + dE ** 2) * range_x

        # Sum the sub-segments back into one length per image-to-image
        # segment (same layout as self.distances)
        return sub_lengths.reshape(self.n_images - 1, n_interpolations).sum(axis=1)

    def compute_curvature_weighted_spring_lengths(self, n_interpolations=20):
        """
        Weight the spring-force segment spacing by the local energy
        curvature ``|d^2E/d(path_distance)^2|`` instead of by the energy
        value, so that images concentrate around the critical points of the
        path.

        Parameters
        ----------
        n_interpolations
            Number of sub-segments used to discretise every image-to-image
            segment when integrating the weighted length.

        Returns
        -------
        numpy.ndarray
            An array of length ``n_images - 1``, where entry ``i`` is the
            curvature-weighted length of the segment between image ``i`` and
            image ``i + 1``. This is the same layout as ``self.distances``,
            so it can be used as a drop-in replacement in
            ``nebm_clib.compute_spring_force``.

        See Also
        --------
        compute_energy_weighted_spring_lengths

        Notes
        -----
        ``dE/d(path_distance)`` vanishes at every critical point of the path
        -- minima *and* maxima alike -- so the energy-value weighting of
        ``compute_energy_weighted_spring_lengths`` refines the monotonic
        flanks between critical points, not the critical points themselves.
        Curvature does the opposite: it is largest exactly at a sharp,
        well-resolved extremum (e.g. a saddle point pinned by a climbing
        image) and near zero on the gently-curved stretches in between, so
        weighting by it concentrates images around critical points instead.

        This uses the adaptive-mesh/equidistribution approach: a monitor
        (weight) function
        ``w(x) = (1 - ratio_C) + ratio_C * |curvature(x)| / max|curvature|``,
        in ``[1 - ratio_C, 1]``, is integrated along the path, with
        ``self.spring_force_ratio`` setting ``ratio_C``, the same convention
        as ``compute_energy_weighted_spring_lengths``. Where ``w(x)`` is
        large, a given physical path-distance step accumulates more of this
        weighted length, so a *smaller* step is needed to fill out an equal
        share of it, concentrating images there.

        Unlike ``compute_energy_weighted_spring_lengths``, this is an ad hoc
        heuristic, not something derived from a reference. Weighting the
        discrete band spacing by curvature is not, by itself, how the GNEB
        literature addresses poor resolution around sharply peaked barriers:
        Schrautzer et al. note this exact limitation of GNEB and instead
        resolve it with a *separate* minimum-mode-following saddle-point
        solver (inverting the gradient along the lowest Hessian eigenmode,
        using only its lowest two eigenpairs), run after GNEB/climbing image
        gets close, rather than by reshaping the band's spacing. That is the
        properly-motivated way to get real resolution at a saddle point;
        this method is a cheaper, purely local stand-in that concentrates
        *existing* images there instead of adding genuinely new information
        about the landscape.

        References
        ----------
        Schrautzer, Sallermann, Bessarab, Jonsson, arXiv:2403.11799 (2024)
        """
        path_distances, energies, spline = self._energy_hermite_spline()
        # Exact (not finite-difference) 2nd derivative of the piecewise
        # cubic spline. CubicHermiteSpline is only C1-continuous (it
        # matches value + 1st derivative at each node), so this curvature
        # can jump discontinuously right at the image nodes; evaluating it
        # at the fine sub-segment midpoints below (x_mid) rather than at
        # the nodes themselves avoids having to pick a side at that jump.
        curvature_spline = spline.derivative(nu=2)

        ratio_C = min(1.0, self.spring_force_ratio)

        x_fine = self._fine_path_grid(path_distances, n_interpolations)
        x_mid = 0.5 * (x_fine[:-1] + x_fine[1:])
        curvature_mid = np.abs(curvature_spline(x_mid))

        curvature_scale = np.max(curvature_mid)
        if curvature_scale < 1e-30:
            curvature_scale = 1.0

        weight = (1.0 - ratio_C) + ratio_C * (curvature_mid / curvature_scale)
        sub_lengths = weight * np.diff(x_fine)

        # Sum the sub-segments back into one length per image-to-image
        # segment (same layout as self.distances)
        return sub_lengths.reshape(self.n_images - 1, n_interpolations).sum(axis=1)

    def _energy_hermite_spline(self):
        """
        Build a cubic Hermite spline through
        ``(path_distance_i, energy_i, dE/d(path_distance)_i)`` for every
        image, including the ones at the extremes of the band.

        Shared by ``compute_energy_weighted_spring_lengths`` and
        ``compute_curvature_weighted_spring_lengths``.

        Returns
        -------
        path_distances : numpy.ndarray
            The distance of every image measured from the 0th image.
        energies : numpy.ndarray
            The energy of every image.
        spline : scipy.interpolate.CubicHermiteSpline
            The spline of the energy as a function of the path distance.

        Notes
        -----
        For the inner images ``self.gradientE`` is the raw effective field,
        not yet a true energy derivative, so it must be converted with the
        same ``self.scale`` factor used in ``compute_polynomial_factors``
        (``mu_0 * Ms * dV`` per degree of freedom for micromagnetics,
        ``mu_s`` per degree of freedom for atomistic simulations) before
        taking the dot product with the tangent. The images at the extremes
        are fixed throughout the relaxation and never get an effective
        field or a tangent computed for them, so a cheap finite-difference
        estimate of the local slope is enough to seed the boundary
        derivative of the spline, with no extra field evaluation required at
        the fixed endpoints.
        """
        path_distances = self.path_distances
        energies = self.energies

        dE_dx = np.zeros(self.n_images)
        gradientE = self.gradientE.reshape(self.n_images, -1)
        tangents = self.tangents.reshape(self.n_images, -1)
        for i in range(1, self.n_images - 1):
            dE_dx[i] = np.dot(self.scale * gradientE[i], tangents[i])
        dE_dx[0] = ((energies[1] - energies[0])
                    / (path_distances[1] - path_distances[0]))
        dE_dx[-1] = ((energies[-1] - energies[-2])
                     / (path_distances[-1] - path_distances[-2]))

        spline = CubicHermiteSpline(path_distances, energies, dE_dx)
        return path_distances, energies, spline

    def _fine_path_grid(self, path_distances, n_interpolations):
        """
        Sub-divide every segment ``[x_i, x_{i+1}]`` of the path into
        ``n_interpolations`` pieces, with the segment boundaries falling
        exactly on the path distance nodes.

        Parameters
        ----------
        path_distances
            The distance of every image measured from the 0th image.
        n_interpolations
            Number of sub-segments per image-to-image segment.

        Returns
        -------
        numpy.ndarray
            The sub-divided path distances, of length
            ``(n_images - 1) * n_interpolations + 1``.
        """
        return np.concatenate(
            [np.linspace(path_distances[i], path_distances[i + 1],
                         n_interpolations, endpoint=False)
             for i in range(self.n_images - 1)] + [[path_distances[-1]]]
            )

    def nebm_step(self, y):
        """
        Update the effective field, the energies, the tangents, the
        distances and the spring force from the band *y*, and combine them
        into the effective force stored in ``self.G``.

        Parameters
        ----------
        y
            The degrees of freedom of the system, i.e. the array containing
            all the spin directions of the images in the band.
        """
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
        Compute the geodesic distance between consecutive images of
        ``self.band``, and store them in ``self.distances`` and
        ``self.path_distances``::

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
        Right hand side of the NEBM equation, as called on every iteration
        of the step integrators in ``chain_method_integrators.py``.

        The effective force is left in ``self.G``, with the entries of the
        images at the extremes of the band set to zero, since those images
        are kept fixed.

        Parameters
        ----------
        t
            Current time of the integrator. Unused, since the NEBM equation
            does not depend explicitly on time.
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
        Right hand side of the NEBM equation, as called on every iteration
        of the CVODE integrator.

        Parameters
        ----------
        t
            Current time of the integrator. Unused, since the NEBM equation
            does not depend explicitly on time.
        y
            The band at which the right hand side is evaluated.
        ydot
            Output array where the right hand side is stored, since we are
            solving ``dy/dt = 0``. Its entries for the images at the
            extremes of the band are set to zero, since those images are
            kept fixed.

        Returns
        -------
        int
            Always 0, as expected by the integrator.
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
