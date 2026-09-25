import numpy as np
import os
import time

import fidimag.extensions.cvode as cvode
from .chain_method_integrators import VerletIntegrator, StepIntegrator
from fidimag.common.vtk import VTK
from .driver_base import DriverBase
from .chain_method_tools import compute_norm
# from .chain_method_tools import linear_interpolation_spherical
from .fileio import DataSaver

import fidimag.common.constant as const
import scipy.interpolate as si

import logging
log = logging.getLogger(name="fidimag")


class ChainMethodBase:
    """
    Base class for chain methods, such as the NEBM or String Method codes.

    This class only sets up the arrays describing an energy band and the
    machinery shared by every chain method: the integrators, the data
    writers, the VTK/npy savers and the interpolation of the energy band.
    The physics -- effective field, tangents, distances and the right hand
    side of the evolution equation -- is left to the subclasses, which must
    define the abstract methods listed in the Notes section.

    Parameters
    ----------
    sim
        An instance of a micromagnetic or an atomistic simulation. Its
        driver is removed, since the band is evolved by this class instead.
    initial_images
        A sequence of arrays or functions setting the magnetisation field
        profile of the images that define the initial band. The images at
        the extremes of the band are kept fixed.
    interpolations
        A list of integers with the number of images to interpolate between
        every pair of consecutive ``initial_images``, hence of length
        ``len(initial_images) - 1``. If not specified, no images are
        interpolated.
    spring_constant
        Magnitude of the spring constant, which keeps the images from
        drifting towards the minima of the energy band.
    name
        Name of the chain method simulation, used as the prefix of every
        output file.
    climbing_image
        Any iterable with the indexes of the climbing images, i.e. images
        whose spring force is removed and whose energy gradient is inverted
        along the tangent, so they are driven towards a saddle point. A
        negative index sets a falling image instead, whose total force is
        the energy gradient only.
    dof
        Degrees of freedom per spin. Spherical coordinates have ``dof=2``
        and Cartesian coordinates have ``dof=3``.
    openmp
        Set to ``True`` to evolve the band with the OpenMP version of the
        CVODE integrator.

    Attributes
    ----------
    dof : int
        Degrees of freedom of the coordinates used in the band (e.g.
        spherical coordinates have ``dof=2``).
    sim : object
        Fidimag atomistic or micromagnetic simulation object.
    mesh : object
        Fidimag simulation mesh object.
    name : str
        Name of the chain method simulation.
    n_spins : int
        Number of spins per image.
    k : numpy.ndarray
        Spring constant of every image.
    variable_k : bool
        Set to ``True`` to update the ``k`` values according to the
        energies, which also requires setting ``dk``. Experimental, and
        ``False`` by default.
    VTK : object
        Fidimag VTK object used to save VTK files.
    files_convert_f : callable or None
        Function converting the coordinates of the band into Cartesian
        coordinates. ``None`` means the band is already Cartesian.
    initial_images : list
        The initial images of the band, as Numpy arrays or as space
        dependent functions of the magnetisation/spin field.
    interpolations : list of int
        Number of images interpolated between every pair of consecutive
        initial images.
    n_images : int
        Number of images in the band, computed from the number of initial
        images and the interpolations between them.
    n_images_inner_band : int
        Number of images excluding the two at the extremes of the band.
    n_dofs_image : int
        Number of degrees of freedom per image, i.e. the number of spins
        times ``dof``.
    n_band : int
        Total number of degrees of freedom in the whole band.
    band : numpy.ndarray
        Every degree of freedom (spin direction) of the band, ordered per
        image and, within an image, in the XYZ format.
    gradientE : numpy.ndarray
        Components of the energy gradient.
    G : numpy.ndarray
        Effective force, as defined by the chain method in use.
    tangents : numpy.ndarray
        Tangents to the band at every image.
    energies : numpy.ndarray
        Energy of every image in the band, of length ``n_images``.
    spring_force : numpy.ndarray
        Components of the spring force.
    distances : numpy.ndarray
        Distances between adjacent images, ``[0-1, 1-2, 2-3, ...]``.
    path_distances : numpy.ndarray
        Distance of every image measured from the 0th image, so its first
        element is always zero.
    last_Y : numpy.ndarray
        The band as computed in the previous step of the integrator, with
        the same layout as ``band``.
    _material : numpy.ndarray of bool
        Array of size ``dof * n_spins``, i.e. the number of degrees of
        freedom of a single image, which is ``True`` where ``Ms`` or
        ``mu_s`` are larger than zero. It filters out the spins that should
        not be counted in, for example, a scaled norm.
    n_dofs_image_material : int
        Number of degrees of freedom of a single image where ``mu_s`` or
        ``Ms`` are larger than zero.
    scale : numpy.ndarray
        Factor rescaling the energy gradient into the right units:
        ``mu_0 * Ms * cell_volume`` for micromagnetics and ``mu_s`` for
        atomistic simulations.
    G_log : list
        The maximum force norms logged during ``relax``, also written to
        ``<name>_G_log.txt``.

    Notes
    -----
    Subclasses must define the following abstract methods:

    ``compute_distances``
        Compute the distances between corresponding images of two bands.
        The inputs are two arrays with at least one full image and the
        output is a 1D array with the distances, so an input with ``x``
        images must return an array with ``x`` entries.
    ``compute_effective_field_and_energy``
        Compute the effective field and the energies of the images of the
        band, according to the number of degrees of freedom (i.e. the total
        number of spin components). The effective field is stored in
        ``gradientE`` and the energies in ``energies``.
    ``initialise_energies``
        Populate the ``energies`` array with the energies of every image at
        the 0th step of the algorithm.
    ``generate_initial_band``
        Use the initial states and the interpolations to generate the
        initial band.
    ``Sundials_RHS``
        Right hand side of the chain method equation, as called by Sundials.
    """
    def __init__(self, sim,
                 initial_images, interpolations=None,
                 spring_constant=1e5,
                 name='unnamed',
                 climbing_image=None,
                 dof=2,
                 openmp=False
                 ):

        self.openmp = openmp

        # Whether relax() has been called on this object before. Used to
        # only save the pre-loop initial-state entries once per object
        # lifetime, since initialise_integrator() (e.g. when switching
        # integrator mid-relaxation) resets self.iterations to 0 without
        # the band actually returning to its original initial state.
        self._relax_called = False

        # Degrees of Freedom per spin
        self.dof = dof

        self.sim = sim
        self.mesh = self.sim.mesh
        self.name = name

        # Number of spins in the system
        self.n_spins = len(self.mesh.coordinates)

        # We will use this filter to know which sites of the system has
        # material, i.e. M_s or mu_s > 0 and norm(m) = 1
        if self.sim._micromagnetic:
            self._material = np.repeat(self.sim.Ms, self.dof) > 1e-10
        else:
            # We will assume, for now, that the magnetic moment in atomistic
            # simulations is in units of mu_B
            self._material = np.repeat(self.sim.mu_s / const.mu_B,
                                       self.dof) > 1e-10

        # self._material = self._material
        # For C, we use 1 and 0s
        self._material_int = np.copy(self._material).astype(np.int32)
        self.n_dofs_image_material = np.sum(self._material)

        # VTK saver for the magnetisation/spin field --------------------------
        self.VTK = VTK(self.mesh, directory='vtks', filename='image')

        # Functions to convert the energy band coordinates to Cartesian
        # coordinates when saving VTK and NPY files We assume Cartesian
        # coordinates by default, i.e. we do not transform anything
        self.files_convert_f = None

        # Initial states ------------------------------------------------------

        # We assume the extremes are fixed
        self.initial_images = initial_images

        if interpolations:
            self.interpolations = interpolations
        else:
            self.interpolations = [0 for i in range(len(initial_images) - 1)]

        # Number of images with/without the extremes
        self.n_images = len(self.initial_images) + np.sum(self.interpolations)
        self.n_images_inner_band = self.n_images - 2

        # Number of degrees of freedom per image
        self.n_dofs_image = (self.dof * self.n_spins)

        # Total number of degrees of freedom in the string/band
        self.n_band = self.n_images * self.n_dofs_image

        # Spring constant -----------------------------------------------------

        # Spring constant (we could use an array in the future)
        self.k = spring_constant * np.ones(self.n_images)

        # Set to True to update spring constant values relative to the energies
        # (TESTING)
        self.variable_k = False
        self.dk = 1

        # Weight (in (0, 1]) given to the energy axis when combining the
        # path distance and the energy into a single spring-force spacing
        # metric (see e.g. NEBM_Geodesic.compute_energy_weighted_spring_lengths).
        # 0 (default) disables it and uses the plain path distance, as before.
        self.spring_force_ratio = 0

        # Which quantity spring_force_ratio weights the spacing by:
        # 'energy' (default) uses NEBM_Geodesic.compute_energy_weighted_spring_lengths
        # (dE/d(path_distance), refines the flanks between critical points),
        # 'curvature' uses NEBM_Geodesic.compute_curvature_weighted_spring_lengths
        # (d^2E/d(path_distance)^2, refines around critical points themselves).
        # Only used when spring_force_ratio > 0.
        self.spring_weighting = 'energy'

        # Divides the max|G|/max|gradE|/max|F_k| values in the relax()
        # debug log, purely for display. Those are reported in the raw units
        # of the effective field (A/m for micromagnetics, Tesla for
        # atomistic simulations), whose magnitude is already readable and is
        # what stopping_max_force is compared against, so
        # this is left at 1 by default. Set it only if some other unit is
        # more convenient for a particular system, keeping in mind that it
        # then no longer matches the stopping_max_force threshold.
        self.log_energy_scale = 1.0

        # Climbing Image ------------------------------------------------------

        # Set a list with the images where 1 is for climbing image and 0 for
        # normal
        self._climbing_image = np.zeros(self.n_images, dtype=np.int32)
        if climbing_image is not None:
            self.climbing_image = climbing_image

        # Chain Method Arrays -------------------------------------------------
        # We will initialise every array using the total number of images,
        # but we must have in mind that the images at the extrema of the band
        # are kept fixed, so we do not compute their gradient, tangents, etc.
        # This might be not memory efficient but the code is understood better
        # when we perform the loops when calculating the effective fields
        # and forces

        # The array containing every degree of freedom
        self.band = np.zeros(self.n_band)

        # The gradient with respect to the magnetisation (effective field)
        self.gradientE = np.zeros_like(self.band)

        # The effective force
        self.G = np.zeros_like(self.band)

        self.tangents = np.zeros_like(self.band)
        self.energies = np.zeros(self.n_images)
        self.spring_force = np.zeros_like(self.band)
        self.distances = np.zeros(self.n_images - 1)
        # Total distance starting from image_0
        # (first element shoud always be zero)
        self.path_distances = np.zeros(self.n_images)

        self.last_Y = np.zeros_like(self.band)

        # ---------------------------------------------------------------------

        # If the integrator uses an LLG-like equation to relax the energy band
        # we need to set this variable
        # This variable only affects the StepIntegrators, NOT Sundials
        self._llg_evolve = False

        # ---------------------------------------------------------------------
        # Factors for interpolating the energy band
        # For now we only have a 3rd order polynomial interp, thus we set
        # 4 factors
        self.interp_factors = np.zeros((4, self.n_images))

        # Somehow we need to rescale the gradient by the right units. In the
        # case of micromag, we use mu0 * Ms, and for the atomistic case we
        # simply use mu_s. This must be related to the way we derive the
        # effective field to calculate the negative energy gradient, which is
        # the functional derivative of the energy
        if self.sim._micromagnetic:
            self.scale = np.repeat(self.mesh.dx * self.mesh.dy * self.mesh.dz *
                                   (self.mesh.unit_length ** 3.) *
                                   const.mu_0 * self.sim.Ms, 3)
        else:
            self.scale = np.repeat(self.sim.mu_s, 3)

        # ---------------------------------------------------------------------

        self.G_log = []

        # Make sure the sim object does not have a driver: (see sim_base.py -> set_m)
        sim.driver = None

    # TODO: Move this property to the NEBM classes because they are only
    # relevant to the NEBM and not the string method
    @property
    def climbing_image(self):
        """
        The climbing images of the band, as an array of length
        ``n_images`` whose entries are ``1`` for a climbing image, ``-1``
        for a falling image and ``0`` otherwise.

        Set it with an image index, or with any iterable of image indexes,
        of the images that will climb towards a saddle point. A negative
        index sets a falling image instead, whose total force is the energy
        gradient only. Indexes of the images at the extremes of the band are
        not allowed, since those images are kept fixed. Deleting it resets
        every image back to a normal one.

        Raises
        ------
        Exception
            If an index does not belong to an image of the inner band.
        """
        return self._climbing_image

    @climbing_image.setter
    def climbing_image(self, climbing_image_list):
        self._climbing_image[:] = 0
        images = range(self.n_images)[1:-1]
        for ci in np.array([climbing_image_list]).flatten():
            if abs(ci) in images:
                if ci > 0:
                    self._climbing_image[ci] = 1
                # Falling images are specified with negative index
                elif ci < 0:
                    self._climbing_image[abs(ci)] = -1
            else:
                raise Exception('Cannot set image={} as climbing image'.format(ci))

    @climbing_image.deleter
    def climbing_image(self):
        self._climbing_image[:] = 0

    def initialise_energies(self):
        pass

    def save_VTKs(self, coordinates_function=None):
        """
        Save a VTK file per image of the band.

        Files are saved in a different folder per simulation name and step,
        as ``vtks/simname_simstep/image_00000x.vti`` (or ``.vtp`` for a
        hexagonal mesh).

        Parameters
        ----------
        coordinates_function
            A function transforming the coordinates of the band into
            Cartesian coordinates. For example, for a band in spherical
            coordinates this is the ``spherical2cartesian`` function from
            ``chain_method_tools``. If ``None``, the band is saved as it is.
        """
        # Create the directory
        directory = "vtks/{}_{:05d}".format(self.name, self.iterations)
        self.VTK.directory = directory

        self.band.shape = (self.n_images, -1)

        # We use Ms from the simulation assuming that all the images are the
        # same
        for i in range(self.n_images):
            self.VTK.reset_data()
            # We will try to save for the micromagnetic simulation (Ms) or an
            # atomistic simulation (mu_s) TODO: maybe this can be done with an:
            # isinstance
            if self.sim._micromagnetic:
                self.VTK.save_scalar(self.sim.Ms, name='M_s')
            else:
                self.VTK.save_scalar(self.sim.mu_s, name='mu_s')

            if coordinates_function:
                self.VTK.save_vector(
                    coordinates_function(self.band[i]).reshape(-1, 3),
                    name='spins'
                    )
            else:
                self.VTK.save_vector(
                    self.band[i].reshape(-1, 3),
                    name='spins'
                    )

            self.VTK.write_file(step=i)

        self.band.shape = (-1, )

    def save_npys(self, coordinates_function=None):
        """
        Save a npy file per image of the band.

        Files are saved in a different folder per simulation name and step,
        as ``npys/simname_simstep/image_x.npy``.

        Parameters
        ----------
        coordinates_function
            A function transforming the coordinates of the band into
            Cartesian coordinates. If ``None``, the band is saved as it is.
        """
        # Create directory as simname_simstep
        directory = 'npys/%s_%d' % (self.name, self.iterations)

        if not os.path.exists(directory):
            os.makedirs(directory)

        # Save the images with the format: 'image_{}.npy'
        # where {} is the image number, starting from 0
        self.band.shape = (self.n_images, -1)
        for i in range(self.n_images):
            name = os.path.join(directory, 'image_{:06}.npy'.format(i))
            if coordinates_function:
                np.save(name, coordinates_function(self.band[i, :]))
            else:
                np.save(name, self.band[i])
        self.band.shape = (-1)

    def initialise_integrator(self, integrator='cvode_bdf', rtol=1e-6, atol=1e-6,
                              linear_solver='spgmr', maxl=30, maxrs=10):
        """
        Start the integrator that evolves the band, and set it in
        ``self.integrator``.

        Parameters
        ----------
        integrator
            One of ``'cvode_bdf'`` (the CVODE solver from Sundials, and the
            default), ``'rk4'``, ``'euler'`` or ``'verlet'``.
        rtol
            Relative tolerance of the CVODE integrator.
        atol
            Absolute tolerance of the CVODE integrator.
        linear_solver
            Only used by the CVODE integrator. ``'spgmr'`` (default) is
            restarted GMRES(``maxl``) with at most ``maxrs`` restarts, and
            ``'diag'`` is CVODE's diagonal approximate Jacobian solver.
        maxl
            Only used by the CVODE integrator with the ``'spgmr'`` linear
            solver. Dimension of the Krylov basis.
        maxrs
            Only used by the CVODE integrator with the ``'spgmr'`` linear
            solver. Maximum number of GMRES restarts.

        Raises
        ------
        Exception
            If the specified integrator is not one of the valid options.

        Notes
        -----
        GMRES keeps a Krylov basis of ``maxl + 1`` copies of the whole band,
        so ``maxl`` directly sets the memory footprint of the integrator:
        ``(maxl + 1) * n_images * n_dofs_image * 8`` bytes. SUNDIALS reserves
        that basis up front, but its pages are only faulted in as GMRES
        actually uses the vectors, so an oversized ``maxl`` does not show up
        as one large allocation: it shows up as resident memory that keeps
        creeping upwards during a relaxation, which is easily mistaken for a
        memory leak. The restarts keep the total iteration reach at
        ``maxl * (1 + maxrs)`` while bounding the basis.
        """
        self.t = 0
        self.iterations = 0
        self.ode_count = 1

        integrator = DriverBase._canonical_integrator(integrator)

        if integrator == 'cvode_bdf':
            if not self.openmp:
                self.integrator = cvode.CvodeSolver(self.band, self.Sundials_RHS,
                                                    linear_solver=linear_solver,
                                                    maxl=maxl, maxrs=maxrs)
                self.integrator.set_options(rtol, atol)
            else:
                self.integrator = cvode.CvodeSolver_OpenMP(self.band, self.Sundials_RHS,
                                                           linear_solver=linear_solver,
                                                           maxl=maxl, maxrs=maxrs)
                self.integrator.set_options(rtol, atol)
        # elif integrator == 'scipy':
        #     self.integrator = ScipyIntegrator(self.band, self.step_RHS)
        #     self.integrator.set_options()
        elif integrator == 'rk4' or integrator == 'euler':
            self.integrator = StepIntegrator(self.band, self.step_RHS,
                                             step=integrator,
                                             stepsize=1e-3)
            self.integrator.set_options()
            self._llg_evolve = True
        elif integrator == 'verlet':
            self.integrator = VerletIntegrator(self.band,    # y
                                               self.G,       # forces
                                               self.step_RHS,
                                               self.n_images,
                                               self.n_dofs_image,
                                               mass=1,
                                               stepsize=1e-4)
            self.integrator.set_options()
            # In Verlet algorithm we only use the total force G and not YxYxG:
            self._llg_evolve = False
        else:
            raise Exception(
                'No valid integrator specified. Available: "cvode_bdf", '
                '"euler", "rk4", "verlet"')

    def create_tablewriter(self):
        entities_energy = {
            'step': {'unit': '<1>',
                     'get': lambda sim: sim.iterations,
                     'header': 'iterations'},
            'energy': {'unit': '<J>',
                       'get': lambda sim: sim.energies,
                       'header': ['image_%d' % i
                                  for i in range(self.n_images)]}
        }

        self.tablewriter = DataSaver(
            self, '%s_energy.ndt' % (self.name),  entities=entities_energy)

        entities_dm = {
            'step': {'unit': '<1>',
                     'get': lambda sim: sim.iterations,
                     'header': 'iterations'},
            'dYs': {'unit': '<1>',
                    'get': lambda sim: sim.distances,
                    'header': ['image_%d_%d' % (i, i + 1)
                               for i in range(self.n_images - 1)]}
        }

        self.tablewriter_dm = DataSaver(
            self, '%s_dYs.ndt' % (self.name), entities=entities_dm)

        # ---------------------------------------------------------------------

    def generate_initial_band(self):
        pass

    def compute_effective_field_and_energy(self, y):
        pass

    # NEBM only:
    # def compute_tangents(self, y):
    #     pass

    # def compute_spring_force(self, y):
    #     pass

    def compute_distances(self):
        pass

    # -------------------------------------------------------------------------
    # CVODE solver ------------------------------------------------------------
    # -------------------------------------------------------------------------

    def compute_norms(self, A, B):
        """
        Compute the norms of the difference between corresponding images of
        the bands *A* and *B*.

        Every norm is scaled by the number of degrees of freedom, and is
        computed using only the mesh/lattice sites with material, i.e. those
        where ``mu_s`` or ``Ms`` are larger than zero.

        Parameters
        ----------
        A
            The degrees of freedom of a band, with at least one full image.
        B
            The degrees of freedom of a band, with at least one full image.

        Returns
        -------
        numpy.ndarray
            A 1D array with one norm per image.
        """

        A_minus_B = A - B

        A_minus_B.shape = (-1, self.n_dofs_image)
        A_minus_B = np.apply_along_axis(
            lambda y: compute_norm(y[self._material], scale=True),
            axis=1,
            arr=A_minus_B
            )

        return A_minus_B.reshape(-1)

    def step_RHS(self, t, y):
        """
        Right hand side of the ODE solved by the chain method, as called by
        the step integrators (Euler, Runge-Kutta and Verlet).

        Parameters
        ----------
        t
            Current time of the integrator.
        y
            The band at which the right hand side is evaluated.
        """
        pass

    def Sundials_RHS(self, t, y, ydot):
        """
        Right hand side of the ODE solved by the chain method, as called on
        every iteration of the CVODE integrator.

        Parameters
        ----------
        t
            Current time of the integrator.
        y
            The band at which the right hand side is evaluated.
        ydot
            Output array where the right hand side is stored, since we are
            solving ``dy/dt = 0``.
        """

        pass

    def compute_maximum_dYdt(self, A, B, dt):
        """
        Compute the maximum difference between the images of the *A* array
        and the images of the *B* array, divided by *dt*.

        Parameters
        ----------
        A
            The degrees of freedom of a band, typically the band from the
            last step of the integrator.
        B
            The degrees of freedom of a band, typically the band from the
            previous step of the integrator.
        dt
            Time step separating the two bands.

        Returns
        -------
        float
            The largest rate of change of the band, or zero if no image
            changed.

        Notes
        -----
        The differences are not computed for the images at the extremes,
        since these images are fixed and do not change with the integrator.

        For instance, in spherical coordinates, if we have a band of
        ``N + 1`` images labelled from 0 to ``N``, we start by::

            dY = [A1_theta0 A1_phi0 A1_theta1 ... A(N-1)_theta0 ... ]
                 - [B1_theta0 B1_phi0 B1_theta1 ... B(N-1)_theta0 ... ]

        where ``A(i)_theta(j)`` is the theta component of the j-th spin of
        the i-th image in the band. Then we calculate the norm of every
        difference, using only the mesh/lattice sites with material, i.e.
        those where ``mu_s`` or ``Ms`` are larger than zero::

            ||dY|| = [ || dY1_theta0 dY1_phi0 dY1_theta1 ... ||,
                                        ...
                       || dY(N-1)_theta0 dY(N-1)_phi0 ...  || ]

        Finally we divide by *dt*, so ``||dY|| -> ||dY|| / dt = dYdt``, and
        take the maximum value of ``dYdt``.
        """

        # We will not consider the images at the extremes to compute dY.
        # Since we removed the extremes, we only have *n_images_inner_band*
        # images
        band_no_extremes = slice(self.n_dofs_image, -self.n_dofs_image)
        dYdt = self.compute_norms(
            A[band_no_extremes],
            B[band_no_extremes]).reshape(self.n_images_inner_band, -1)

        dYdt /= dt

        if np.max(dYdt) > 0:
            return np.max(dYdt)
        else:
            return 0

    def run_until(self, t):

        if (t) <= self.t:
            return

        self.integrator.run_until(t)

        # Copy the updated energy band to our local array
        self.band[:] = self.integrator.y[:]

        # CVODE's CV_NORMAL mode returns the solution at *t* by
        # interpolating its internal (adaptive) steps, which can land past
        # *t*. self.energies/self.gradientE (and, for NEBM subclasses that
        # define it, self.tangents/self.spring_force/self.G) were last set
        # by Sundials_RHS at that internal, generally different, y -- so
        # they are stale here and must be refreshed against the actual
        # self.band being returned, otherwise everything computed from them
        # afterwards (logging, climbing-image selection, stopping_max_force)
        # is off by a fraction of a step. This gap grows as the integrator
        # takes larger internal steps (e.g. later in the relaxation), which
        # is why it is barely noticeable early on but visibly compounds
        # after several relax() calls.
        nebm_step = getattr(self, 'nebm_step', None)
        if nebm_step is not None:
            nebm_step(self.band)
        else:
            self.compute_effective_field_and_energy(self.band)

        # Compute the maximum change in the integrator step
        max_dYdt = self.compute_maximum_dYdt(self.integrator.y, self.last_Y,
                                             t - self.t)

        self.last_Y[:] = self.band[:]

        # Update the current step
        self.t = t

        return max_dYdt

    def relax(self, dt=1e-8, stopping_dYdt=1, max_iterations=1000,
              save_npys_every=100, save_vtks_every=100,
              save_initial_state=True, stopping_max_force=None
              ):

        """
        Relax the energy band for a given number of steps of the integrator
        set in ``initialise_integrator``.

        Parameters
        ----------
        dt
            Time step of the relaxation. Its meaning depends on the
            integrator in use:

            - CVODE: the initial step size, which is then updated by CVODE
              itself, so the number of evaluations is decided by CVODE.
            - Step integrators (Euler, Runge-Kutta, Verlet): the relaxation
              step size. These integrators evolve the band using an internal
              ``stepsize``, so the number of evaluations is
              ``dt / stepsize``. The internal step can be updated with
              ``self.integrator.stepsize = 1e-4``.
        stopping_dYdt
            Stop the relaxation once the largest rate of change of the band
            (see ``compute_maximum_dYdt``) drops below this value.
        max_iterations
            Maximum number of steps of the integrator.
        save_npys_every
            Save a npy file per image every this number of steps.
        save_vtks_every
            Save a VTK file per image every this number of steps.
        save_initial_state
            Save the VTK/npy files and the table entries of the initial band
            before the relaxation starts. They are only saved on the first
            ``relax`` call of an object, since any later call continues from
            the band that the previous one already saved.
        stopping_max_force
            Optional. If set, also stop once the largest force norm on the
            band (``max|G|``, over the inner images) drops below this value.
            This is a step-size-independent convergence check, unlike
            ``stopping_dYdt``, which only measures how much the Cartesian
            coordinates changed in the last step and can be small simply
            because the integrator took a tiny internal step, even far from
            a true force equilibrium. Off (``None``) by default, so the
            previous ``stopping_dYdt``-only behaviour is unchanged unless
            this is explicitly requested.

        Notes
        -----
        The ``stopping_max_force`` comparison uses the largest ``|G|`` over
        the sites that carry material, in the raw units of the effective
        field: A/m for micromagnetics and Tesla for atomistic simulations.
        For a micromagnetic system this is the same kind of quantity as
        OOMMF's ``|m x H x m|``, so a value around ``1e-2`` means here what
        it means there; for an atomistic system the threshold is a field in
        Tesla and has to be chosen on that scale. In both cases it is
        exactly the ``max|G|`` written to the debug log, so a threshold can
        be read straight off a previous run. The same values are collected
        in ``self.G_log`` and written to ``<name>_G_log.txt`` at the end of
        ``relax``.

        Cells with no material are excluded on purpose: they keep feeling
        the stray field of the magnetised region while nothing relaxes them,
        so including them would hold ``max|G|`` at a fixed floor and this
        criterion would never be met.
        """

        # Units of the max|G| / max|gradE| / max|F_k| values reported below,
        # and of the stopping_max_force threshold they are compared against:
        # the effective field is in A/m for micromagnetics and in Tesla
        # (energy per magnetic moment) for atomistic simulations.
        force_unit = "A/m" if self.sim._micromagnetic else "T"
        if self.log_energy_scale != 1.0:
            force_unit += f"/{self.log_energy_scale:g}"

        log.debug("Relaxation parameters: " +
                  f"stopping_dYdt={stopping_dYdt}, " +
                  f"time_step={dt} s, " +
                  f"max_iterations={max_iterations}, " +
                  f"stopping_max_force={stopping_max_force}, " +
                  f"forces in [{force_unit}] over material sites")

        # Do not re-save the initial VTK/npy/table log entries on a
        # restart (i.e. any relax() call after the first one on this
        # object), since they were already saved as the last entries of
        # the previous relax() call. Note self.iterations is not a
        # reliable "first call" check on its own, since
        # initialise_integrator() may reset it to 0 between relax() calls
        # (e.g. when switching integrator or enabling climbing images)
        # without the band returning to its original initial state.
        if save_initial_state and not self._relax_called:
            self.save_VTKs(coordinates_function=self.files_convert_f)
            self.save_npys(coordinates_function=self.files_convert_f)

        # Save the initial state i=0 in the data table
        # Update self.distances and self.path_distances:
        self.compute_distances()

        if not self._relax_called:
            self.tablewriter.save()
            self.tablewriter_dm.save()

        self._relax_called = True

        INNER_DOFS = slice(self.n_dofs_image, -self.n_dofs_image)
        # Per-spin mask selecting the sites that carry material, repeated for
        # every inner image, to pick out the meaningful entries of the force
        # norms reported below. self._material holds one entry per degree of
        # freedom and the dof components of a spin share the same site.
        INNER_MATERIAL = np.tile(self._material.reshape(-1, self.dof)[:, 0],
                                 self.n_images - 2)

        for i in range(max_iterations):

            # Update the iterations number counter
            self.iterations += 1

            # Get the current size of the time discretisation from the
            # integrator (variable step size)
            try:
                cvode_dt = self.integrator.get_current_step()
            except AttributeError:
                cvode_dt = dt

            # If the step size of the integrator is larger than the specified
            # discretisation, use the current integrator step size to
            # compute the next iteration. Otherwise, just stick to the
            # specified step
            if cvode_dt > dt:
                increment_dt = cvode_dt
            else:
                increment_dt = dt

            # Integrator steps: ***********************************************
            # This can be redefined according to the chosen chain method
            max_dYdt = self.run_until(self.t + increment_dt)
            # *****************************************************************

            # Update the current step
            self.t = self.t + increment_dt

            # Save data -------------------------------------------------------

            if self.iterations % save_vtks_every == 0:
                self.save_VTKs(coordinates_function=self.files_convert_f)
            if self.iterations % save_npys_every == 0:
                self.save_npys(coordinates_function=self.files_convert_f)

            self.compute_distances()
            self.tablewriter.save()
            self.tablewriter_dm.save()

            # Print information about the simulation and the forces.
            # The last two terms are the largest gradient and spring
            # force norms from the spins (not counting the extrema)
            #
            # Only the sites that carry material are taken into account. In a
            # patterned sample the surrounding cells have Ms (or mu_s) equal
            # to zero, but they still see the stray field of the magnetised
            # region: the exchange and anisotropy fields vanish there because
            # they are scaled by the magnetisation, while the demagnetising
            # field does not. Nothing relaxes at those cells, so their |G|
            # stays at whatever the stray field imposes for the whole
            # relaxation, which would put a floor under max|G| that the
            # stopping_max_force criterion could never get below.
            #
            # The norms are left in the raw units of the effective field:
            # A/m for micromagnetics, and Tesla for atomistic simulations,
            # where the field is an energy per magnetic moment, which is why
            # self.scale is mu_s there rather than mu_0 * Ms * dV. For a
            # micromagnetic system this is the same kind of quantity as
            # OOMMF's |m x H x m|, so the values printed here mean the same as
            # the mxHxm figures OOMMF reports and can be used as they are to
            # choose stopping_max_force; for an atomistic system the numbers
            # are fields in Tesla and live on a different scale. self.G,
            # self.gradientE and self.spring_force are all in these units by
            # construction, since G is assembled from the other two, so the
            # three printed numbers can always be compared with one another.
            G_norms = np.linalg.norm(
                self.G[INNER_DOFS].reshape(-1, self.dof), axis=1)[INNER_MATERIAL]
            gradE_norms = np.linalg.norm(
                self.gradientE[INNER_DOFS].reshape(-1, self.dof), axis=1)[INNER_MATERIAL]
            Fk_norms = np.linalg.norm(
                self.spring_force[INNER_DOFS].reshape(-1, self.dof), axis=1)[INNER_MATERIAL]

            # For DEBUGGING purposes: -----------------------------------------
            # mean_G_norms_per_image = np.mean(G_norms.reshape(self.n_images - 2, -1), axis=1)
            # print(mean_G_norms_per_image)
            # gradE_dot_t = np.einsum('ij,ij->i',
            #                         self.gradientE.reshape(self.n_images, -1),
            #                         self.tangents.reshape(self.n_images, -1))
            # gradE_perp = (self.gradientE.reshape(self.n_images, -1)
            #               - np.einsum('i,ij->ij', gradE_dot_t,
            #                           self.tangents.reshape(self.n_images, -1))
            #               )
            # print(np.max(gradE_perp))
            # print(np.max(gradE_norms.reshape(self.n_images - 2, -1), axis=1))
            # print(np.max(Fk_norms.reshape(self.n_images - 2, -1), axis=1))
            # -----------------------------------------------------------------

            log.debug(time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime()) +
                      f"step: {self.iterations:>6d}, " +
                      f"step_size: {increment_dt:>8.3g}, " +
                      f"max dYdt: {max_dYdt:>8.3g} " +
                      "max|G|: {:>8.3g} ".format(np.max(G_norms) / self.log_energy_scale) +
                      "max|gradE|: {:>8.3g} ".format(np.max(gradE_norms) / self.log_energy_scale) +
                      "and max|F_k|: {:>8.3g}".format(np.max(Fk_norms) / self.log_energy_scale)
                      )

            self.G_log.append(np.max(G_norms))

            # -----------------------------------------------------------------

            # Stop criteria:
            if max_dYdt < stopping_dYdt:
                break
            if (stopping_max_force is not None
                    and np.max(G_norms) < stopping_max_force):
                break

        np.savetxt(self.name + '_G_log.txt', np.array(self.G_log))

        log.info("Relaxation finished at time step = {:.4g}, "
                 "t = {:.2g}, call rhs = {:.4g} "
                 "and max_dYdt = {:.3g}".format(self.iterations,
                                                self.t,
                                                self.ode_count,
                                                max_dYdt)
                 )

        self.save_VTKs(coordinates_function=self.files_convert_f)
        self.save_npys(coordinates_function=self.files_convert_f)

    # -------------------------------------------------------------------------
    # Interpolations

    def compute_polynomial_factors(self, compute_fields=True):

        """
        Compute a smooth approximation of the energy band, using a third
        order polynomial, and store its prefactors in the
        ``self.interp_factors`` array.

        The approximation uses the tangents and the derivatives of every
        image of the band as the information to estimate the curvatures.

        Parameters
        ----------
        compute_fields
            Update the effective field, the tangents and the distances
            before computing the prefactors. Necessary if the band has not
            been relaxed before calling this method.

        References
        ----------
        Bessarab et al., Computer Physics Communications 196 (2015) 335-347
        """

        # To be sure, update the effective field and tangents when calling
        # this function (this is necesary if we call the function without
        # relaxing the band before)
        if compute_fields:
            self.compute_effective_field_and_energy(self.band)
            self.compute_tangents(self.band)
            self.compute_distances()

        self.gradientE.shape = (self.n_images, -1)
        self.tangents.shape = (self.n_images, -1)

        deltas = np.zeros(self.n_images)
        for i in range(self.n_images):
            deltas[i] = np.dot(self.scale * self.gradientE[i], self.tangents[i])

        self.gradientE.shape = (-1)
        self.tangents.shape = (-1)

        ds = self.path_distances
        E = self.energies

        # The coefficients for the polynomial approximation
        self.interp_factors[2][:] = deltas
        self.interp_factors[3][:] = E

        # Populate the a and b coefficients for every image
        for i in range(self.n_images - 1):
            # a factor
            self.interp_factors[0][i] = (deltas[i + 1] + deltas[i]) / (ds[i + 1] - ds[i]) ** 2.
            self.interp_factors[0][i] -= 2 * (E[i + 1] - E[i]) / (ds[i + 1] - ds[i]) ** 3.

            # b factor
            self.interp_factors[1][i] = -(deltas[i + 1] + 2 * deltas[i]) / (ds[i + 1] - ds[i])
            self.interp_factors[1][i] += 3 * (E[i + 1] - E[i]) / (ds[i + 1] - ds[i]) ** 2.

    def compute_polynomial_approximation_energy(self, n_points):

        """
        Compute a smooth approximation of the energy band, using the third
        order polynomial whose prefactors are computed by the
        ``compute_polynomial_factors`` method.

        Parameters
        ----------
        n_points
            Number of points of the interpolation.

        Returns
        -------
        x : numpy.ndarray
            The distance of every data point measured from the 0th image.
        E_interp : numpy.ndarray
            An ``n_points`` long array with the interpolated energy band.
        """

        ds = self.path_distances
        # The arrays with the data points and the interpolated energy values
        x = np.linspace(0, ds[-1], n_points)
        E_interp = np.array([self._compute_polynomial_approximation_energy(i) for i in x])

        return x, E_interp

    def _compute_polynomial_approximation_energy(self, x):

        """
        Return the polynomial interpolation of the energy at the point *x*.

        Parameters
        ----------
        x
            Distance from the 0th image, within the path distances of the
            band.

        Returns
        -------
        float
            The interpolated energy.

        Raises
        ------
        Exception
            If *x* lies outside the path distances of the band.
        """

        ds = self.path_distances
        if x < 0.0 or x > ds[-1]:
            raise Exception('x lies outside the valid interpolation range')
        # Find index of the ds array for the value that is closest to x
        ds_idx = np.abs(x - ds).argmin()
        # If x is smaller than the given ds, use the previous ds value so
        # that we use ds(i) when x lies in the interval ds(i) < x < ds(i+1)
        if x < ds[ds_idx]:
            ds_idx -= 1

        E_interp = (self.interp_factors[0][ds_idx] * ((x - ds[ds_idx]) ** 3.) +
                    self.interp_factors[1][ds_idx] * ((x - ds[ds_idx]) ** 2.) +
                    self.interp_factors[2][ds_idx] * ((x - ds[ds_idx])) +
                    self.interp_factors[3][ds_idx]
                    )

        return E_interp

    # -------------------------------------------------------------------------

    def compute_Bernstein_polynomials(self, compute_fields=True):

        """
        Compute the Bernstein polynomials approximating the energy curve,
        and store them in the ``self.Bernstein_polynomials`` list.

        These polynomials are used by the
        ``compute_Bernstein_approximation_energy`` method.

        Parameters
        ----------
        compute_fields
            Update the effective field, the tangents and the distances
            before computing the polynomials. Necessary if the band has not
            been relaxed before calling this method.
        """

        # To be sure, update the effective field and tangents when calling
        # this function (this is necesary if we call the function without
        # relaxing the band before)
        if compute_fields:
            self.compute_effective_field_and_energy(self.band)
            self.compute_tangents(self.band)
            self.compute_distances()

        derivatives = np.zeros(self.n_images)
        for i in range(self.n_images):
            derivatives[i] = np.dot(
                self.scale * (self.gradientE).reshape(self.n_images, -1)[i],
                self.tangents.reshape(self.n_images, -1)[i])

        E = self.energies

        # The coefficients for the polynomial approximation
        # self.interp_factors[0][:] = E
        # self.interp_factors[1][:] = deltas

        # Store the polynomial functions
        self.Bernstein_polynomials = []
        for i, ds in enumerate(self.distances):
            self.Bernstein_polynomials.append(
                si.BPoly.from_derivatives(
                    [self.path_distances[i], self.path_distances[i + 1]],
                    [[E[i], derivatives[i]],
                     [E[i + 1], derivatives[i + 1]]]
                )
            )

    def _compute_Bernstein_approximation_energy(self, x):

        """
        Return the Bernstein interpolation of the energy at the point *x*.

        Parameters
        ----------
        x
            Distance from the 0th image, within the path distances of the
            band.

        Returns
        -------
        float
            The interpolated energy.

        Raises
        ------
        Exception
            If *x* lies outside the path distances of the band.
        """

        ds = self.path_distances
        if x < 0.0 or x > ds[-1]:
            raise Exception('x lies outside the valid interpolation range')
        elif x == 0.0:
            return self.distances[0]
        elif x == ds[-1]:
            return self.distances[-1]
        # Find index of the ds array for the value that is closest to x
        ds_idx = np.abs(x - ds).argmin()
        # If x is smaller than the given ds, use the previous ds value so
        # that we use ds(i) when x lies in the interval ds(i) < x < ds(i+1)
        if x < ds[ds_idx]:
            ds_idx -= 1

        return self.Bernstein_polynomials[ds_idx](x)

    def compute_Bernstein_approximation_energy(self, n_points):

        """
        Compute a smooth approximation of the energy band, using the
        Bernstein polynomials computed by the
        ``compute_Bernstein_polynomials`` method.

        Parameters
        ----------
        n_points
            Number of points of the interpolation.

        Returns
        -------
        x : numpy.ndarray
            The distance of every data point measured from the 0th image.
        E_interp : numpy.ndarray
            An ``n_points`` long array with the interpolated energy band.
        """

        ds = self.path_distances
        # The arrays with the data points and the interpolated energy values
        x = np.linspace(0, ds[-1], n_points)
        E_interp = np.zeros_like(x)
        E_interp[0] = self.energies[0]
        E_interp[-1] = self.energies[-1]
        E_interp[1:-1] = np.array(
            [self._compute_Bernstein_approximation_energy(i) for i in x[1:-1]]
            )

        return x, E_interp

        return E_interp
