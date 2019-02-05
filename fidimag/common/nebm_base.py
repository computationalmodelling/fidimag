from __future__ import print_function
from __future__ import division
import numpy as np
import os
import time

import fidimag.extensions.cvode as cvode
from .integrators import ScipyIntegrator, StepIntegrator
from .nebm_integrators import VerletIntegrator
from fidimag.common.vtk import VTK
from .nebm_tools import compute_norm
# from .nebm_tools import linear_interpolation_spherical
from .fileio import DataSaver

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="fidimag")

import fidimag.common.constant as const


class NEBMBase(object):
    """

    Base Class for the NEBM codes.

    Abstract Methods: ---------------------------------------------------------

        compute_distances          :: Function to compute the distances between
                                      corresponding images of two bands. So the
                                      inputs are two arrays with at least one
                                      full image. The output is a 1D array with
                                      the distances. So if the inputs have x
                                      images we must return an array with x
                                      entries.

        compute_effective_field_and_energy   :: Calculate effective field and
                                                the energies of the images for
                                                the energy band, according to
                                                the number of degrees of
                                                freedom (total number of spin
                                                components). Effective field is
                                                stored in the self.gradient and
                                                the enrgies in self.energies

        initialise_energies        :: Populate the self.energies array with the
                                      energies of every image for the 0th step
                                      of the algorithm

        Sundials_RHS               :: Right hand side of the NEB equation for
                                      Sundials

    Methods -------------------------------------------------------------------

        compute_spring_force       :: Compute the spring force for every degree
                                      of freedom of the band (we usually do
                                      this using the C library)
        compute_tangents           :: Compute the tangents of the system (use
                                      the function defined  in
                                      neb_method/nebm_clib.c). This tangents
                                      need to be normalised / projected
                                      according to the variation of the NEB
                                      method.

        compute_maximum_dYdt       :: Compute the maximum distance between the
                                      images (not counting the extremes) of the
                                      last computed band in the evolver
                                      (self.integrator.y) with the band of the
                                      previous step (stored in self.last_Y)
                                      divided by the last time step of the
                                      evolver.

        create_tablewriter         :: Start the data frameworks to output the
                                      energies ( *_energy.ndt) and distances
                                      (*_dYs.ndt) for every step of the
                                      integrator
        generate_initial_band      :: Use the initial states and inteprolations
                                      lists to generate the initial band. The
                                      interpolations are linearly made in
                                      spherical coordinates, using the
                                      nebm_tools.linear_interpolation_spherical
                                      function. We may change this in the
                                      future to generate an inteprolation using
                                      a geodesic path

        initialise_integrator      :: Start the CVODE integrator and define it
                                      in self.integrator
        relax                      :: Relax the band for a specific number of
                                      steps. Arguments of this function
                                      includes saving VTK and NPY files every
                                      certain number of steps and some criteria
                                      for stopping the integrator
        run_until                  :: This function is called from the relax
                                      function (maybe is not very useful, we
                                      may check this in the future)
        save_VTKs                  :: Function to save a VTK of every image for
                                      the current step of the integrator, which
                                      is defined in the self.iterations
                                      variable
        save_npys                  :: Same for NPY files

    ARGUMENTS -----------------------------------------------------------------

    sim                 :: An instance of a micromagnetic or an atomistic
                           simulation

    initial_images      :: A sequence of arrays or functions to set up the
                           magnetisation field profile

    interpolations      ::

    dof                 :: Degrees of freedom per spin. Spherical coordinates
                           have dof=2 and Cartesian have dof=3

    VARIABLES -----------------------------------------------------------------

        self.dof              :: Degree of freedom for the coordinates used in
                                 the band (e.g. Spherical has self.dof=2)
        self.sim              :: Fidimag atomistic or micromagnetic simulation
                                 object
        self.mesh             :: Fidimag simulation mesh object
        self.name             :: Name of the NEBM simulation
        self.n_spins          :: Number of spins per image
        self.k                :: Spring constant
        self.variable_k(TEST) :: Set True to update k values acc to energy.
                                 Need to specify self.dk var.
                                 Default value: False
        self.VTK              :: Fidimag VTK object to save VTK files
        self.files_convert_f  :: Function to convert the coordinates from the
                                 band to Cartesian coordinates
        self.initial_images   :: List with the initial images for the band
                                 (Numpy arrays or space functions with the
                                 magnetisation/spin field)
        self.interpolations   :: Optional List with integers indicating
                                 interpolations between images
        self.n_images         :: Number of images in the band calculated from
                                 the number of initial images and
                                 interpolations between them
        self.n_images_inner_band :: Number of images without considering the
                                    images at the extremes
        self.n_dofs_image     :: Number of degrees of freedom per image, which
                                 is the total number of spins multiplied by the
                                 self.dof
        self.n_band           :: Total number of degrees of freedom in the
                                 whole band
        self.band             :: The array containing all the degrees of
                                 freedom (spin directions). It is ordered in
                                 the XYZ format
        self.gradientE        :: Array with the components of the energy
                                 gradient
        self.G                :: Array with the effective force from the NEB
                                 method definition
        self.tangents         :: Array with the NEBM tangents
        self.energies         :: Array with the energies of every image in the
                                 band (length = self.n_images)
        self.spring_force     :: Array with the spring force components
        self.distances        :: Array with the distances between adjacent
                                 images: [0-1 1-2 2-3 .. etc ]
        self.last_Y           :: Array with the last computed band (like
                                 self.band) from the integrator
        self._material        :: Array of Booleans of size (self.dof*n_spins),
                                 i.e. the number of dofs in an image.
                                 It is True where Ms or mu_s are larger than
                                 zero. This is necessary to filter spins that
                                 not need to be counted in, for example, a
                                 scaled norm.
        self.n_dofs_image_material :: Number of dofs where mu_s or Ms > 0
                                      (in a single image)
        climbing_image        :: Any iterable with the indexes of the climbing
                                 images. Climbing images are stored in the
                                 self._climbing_image array, which the
                                 self.climbing_image decorator populates

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
        self.VTK = VTK(self.mesh,
                       directory='vtks'.format(self.name),
                       filename='image'
                       )

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

        # Total number of degrees of freedom in the NEBM band
        self.n_band = self.n_images * self.n_dofs_image

        # Spring constant -----------------------------------------------------

        # Spring constant (we could use an array in the future)
        self.k = spring_constant * np.ones(self.n_images)

        # Set to True to update spring constant values relative to the energies
        # (TESTING)
        self.variable_k = False
        self.dk = 1

        # Climbing Image ------------------------------------------------------

        # Set a list with the images where 1 is for climbing image and 0 for
        # normal
        self._climbing_image = np.zeros(self.n_images, dtype=np.int32)
        if climbing_image is not None:
            self.climbing_image = climbing_image

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

    @property
    def climbing_image(self):
        return self._climbing_image

    @climbing_image.setter
    def climbing_image(self, climbing_image_list):
        self._climbing_image[:] = 0
        images = range(self.n_images)[1:-1]
        for ci in np.array([climbing_image_list]).flatten():
            if ci in images:
                self._climbing_image[ci] = 1
            else:
                raise Exception('Cannot set image={} as climbing image'.format(ci))

    @climbing_image.deleter
    def climbing_image(self):
        self._climbing_image[:] = 0

    def initialise_energies(self):
        pass

    def save_VTKs(self, coordinates_function=None):
        """

        Save VTK files in different folders, according to the simulation name
        and step. Files are saved as vtks/simname_simstep_vtk/image_00000x.vtk

        coordinates_function    :: A function to transform the coordinates of
                                   the band to Cartesian coordinates. For
                                   example, in spherical coordinates we need
                                   the spherical2cartesian function from
                                   nebm_tools

        """
        # Create the directory
        directory = 'vtks/%s_%d' % (self.name, self.iterations)
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
        Save npy files in different folders according to
        the simulation name and step
        Files are saved as: npys/simname_simstep/image_x.npy
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

    def initialise_integrator(self,
                              integrator='sundials',
                              rtol=1e-6, atol=1e-6):
        self.t = 0
        self.iterations = 0
        self.ode_count = 1

        if integrator == 'sundials':
            if not self.openmp:
                self.integrator = cvode.CvodeSolver(self.band,
                                                    self.Sundials_RHS)
                self.integrator.set_options(rtol, atol)
            else:
                self.integrator = cvode.CvodeSolver_OpenMP(self.band,
                                                           self.Sundials_RHS)
                self.integrator.set_options(rtol, atol)
        # elif integrator == 'scipy':
        #     self.integrator = ScipyIntegrator(self.band, self.step_RHS)
        #     self.integrator.set_options()
        elif integrator == 'rk4' or integrator == 'euler':
            self.integrator = StepIntegrator(self.band, self.step_RHS,
                                             step=self.integrator,
                                             stepsize=1e-3)
            self.integrator.set_options()
            self._llg_evolve = True
        elif integrator == 'verlet':
            self.integrator = VerletIntegrator(self.band,
                                               self.step_RHS,
                                               self.n_images,
                                               mass=1,
                                               stepsize=1e-3)
            self.integrator.set_options()
            # In Verlet algorithm we only use the total force G and not YxYxG:
            self._llg_evolve = False
        else:
            raise Exception('No valid integrator specified. Available: '
                            '"sundials", "scipy"')

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

    def compute_tangents(self, y):
        pass

    def compute_spring_force(self, y):
        pass

    def compute_distances(self):
        pass

    # -------------------------------------------------------------------------
    # CVODE solver ------------------------------------------------------------
    # -------------------------------------------------------------------------

    def compute_norms(self, A, B):
        """
        Compute the norms of the difference between corresponding images of the
        bands A and B. Every norm is scaled by the number of degrees of freedom

        The norm is computed only using mesh/lattice sites with material, i.e.
        mu_s or Ms > 0
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
        The RHS of the ODE to be solved for the NEBM
        This function is specified for the Scipy integrator
        """

        self.ode_count += 1

        # Update the effective field, energies, spring forces and tangents
        # using the *y* array
        self.nebm_step(y)
        # Now set the RHS of the equation as the effective force on the energy
        # band, which is stored on the self.G array
        ydot = self.G[:]
        # The effective force at the extreme images should already be zero, but
        # we will manually remove any value
        ydot[:self.n_dofs_image] = 0
        ydot[-self.n_dofs_image:] = 0

        return ydot

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

    def compute_maximum_dYdt(self, A, B, dt):
        """

        Compute the maximum difference from the images of the *A* array with
        the images of the *B* array, divided by dt

        The differences are not computed for the images at the extremes, since
        these images are fixed and do not change with the iterator

        For instance, in spherical coordinates, if we have a band of (N + 1)
        images, labeled from 0 to N, we start by

        dY = [A1_theta0 A1_phi0 A1_theta1 ... A(N-1)_theta0 A(N-1)_phi0 A(N-1)_theta1 ... ]
            - [B1_theta0 B1_phi0 B1_theta1 ... B(N-1)_theta0 B(N-1)_phi0 B(N-1)_theta1 ... ]

        where A(i)_theta(j) is the theta componenet of the j-th spin of the
        i-th image in the band

        Then we calculate the norm of every difference (the norm is computed
        only using mesh/lattice sites with material, i.e.  mu_s or Ms > 0):

        ||dY|| =  [ || dY1_theta0 dY1_phi0 dY1_theta1         ...       ||
                                   ...
                    || dY(N-1)_theta0 dY(N-1)_phi0 dY(N-1)_theta1  ...  ||
                  ]

        We divide by dt:

        ||dY|| -- > || dY || / dt = dYdt

        and finally take the maximum value of dYdt

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

        # Compute the maximum change in the integrator step
        max_dYdt = self.compute_maximum_dYdt(self.integrator.y, self.last_Y,
                                             t - self.t)

        self.last_Y[:] = self.band[:]

        # Update the current step
        self.t = t

        return max_dYdt

    def relax(self, dt=1e-8, stopping_dYdt=1, max_iterations=1000,
              save_npys_every=100, save_vtks_every=100
              ):

        """
        Relax the energy band according to the specified integrator

            Sundials:        dt is the initial stepsize which is updated
                             by CVODE. Number of calls is given by CVODE
            StepIntegrators: dt is the relaxation stepsize. These integrators
                             evolve using an internal `stepsize`, thus the
                             number of evaluations is dt / stepsize.
                             You can update the integrator evolve step using:
                                self.integrator.stepsize = 1e-4

        """

        log.debug("Relaxation parameters: "
                  "stopping_dYdt={}, "
                  "time_step={} s, "
                  " max_iterations={}.".format(stopping_dYdt,
                                               dt,
                                               max_iterations))

        self.save_VTKs(coordinates_function=self.files_convert_f)
        self.save_npys(coordinates_function=self.files_convert_f)

        # Save the initial state i=0
        self.distances = self.compute_distances(self.band[self.n_dofs_image:],
                                                self.band[:-self.n_dofs_image]
                                                )
        # Compute total distance relative to Y_0
        self.path_distances[1:] = np.cumsum(self.distances)

        self.tablewriter.save()
        self.tablewriter_dm.save()

        INNER_DOFS = slice(self.n_dofs_image, -self.n_dofs_image)

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

            # max_dYdt = self.run_until(self.t + increment_dt)

            if self.t + increment_dt <= self.t:
                break

            self.integrator.run_until(self.t + increment_dt)

            # Copy the updated energy band to our local array
            self.band[:] = self.integrator.y[:]

            # Compute the maximum change in the integrator step
            max_dYdt = self.compute_maximum_dYdt(self.band,
                                                 self.last_Y,
                                                 increment_dt)

            self.last_Y[:] = self.band[:]

            # Update the current step
            self.t = self.t + increment_dt

            # Save data -------------------------------------------------------

            if self.iterations % save_vtks_every == 0:
                self.save_VTKs(coordinates_function=self.files_convert_f)
            if self.iterations % save_npys_every == 0:
                self.save_npys(coordinates_function=self.files_convert_f)

            self.distances = self.compute_distances(
                self.band[self.n_dofs_image:],
                self.band[:-self.n_dofs_image]
                )
            self.path_distances[1:] = np.cumsum(self.distances)
            self.tablewriter.save()
            self.tablewriter_dm.save()

            # Print information about the simulation and the NEBM forces.
            # The last two terms are the largest gradient and spring
            # force norms from the spins (not counting the extrema)
            G_norms = np.linalg.norm(self.G[INNER_DOFS].reshape(-1, 3), axis=1)
            # mean_G_norms_per_image = np.mean(G_norms.reshape(self.n_images - 2, -1), axis=1)
            # print(mean_G_norms_per_image)
            Fk_norms = np.linalg.norm(self.spring_force[INNER_DOFS].reshape(-1, 3), axis=1)
            log.debug(time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime()) +
                      "step: {:.3g}, step_size: {:.3g}, "
                      "max dYdt: {:.3g} "
                      "max|G|: {:.3g} "
                      "and max|F_k|: {:.3g}".format(self.iterations,
                                                    increment_dt,
                                                    max_dYdt,
                                                    np.max(G_norms),
                                                    np.max(Fk_norms)
                                                    )
                      )

            self.G_log.append(np.max(G_norms))

            # -----------------------------------------------------------------

            # Stop criteria:
            if max_dYdt < stopping_dYdt:
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
        Compute a smooth approximation for the band, using a third order
        polynomial approximation. Prefactors are stored in the
        self.interp_factors array

        This approximation uses the tangents and derivatives from each image of
        the band as information to estimate the curvatures. The formula can be
        found in

        - Bessarab et al., Computer Physics Communications 196 (2015) 335-347

        """

        # To be sure, update the effective field and tangents when calling
        # this function (this is necesary if we call the function without
        # relaxing the band before)
        if compute_fields:
            self.compute_effective_field_and_energy(self.band)
            self.compute_tangents(self.band)
            self.distances = self.compute_distances(self.band[self.n_dofs_image:],
                                                    self.band[:-self.n_dofs_image]
                                                    )
            self.path_distances[1:] = np.cumsum(self.distances)

        deltas = np.zeros(self.n_images)
        for i in range(self.n_images):
            deltas[i] = np.dot(self.scale * self.gradientE.reshape(self.n_images, -1)[i],
                               self.tangents.reshape(self.n_images, -1)[i]
                               )

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

        Compute a smooth approximation of the energy band, using a third order
        polynomial approximation. The pre-factors in the polynomial are
        computed from the self.compute_polynomial_factors function

        This function returns a tuple with two elements:

            0. An array with the distance of every data point from the 0th
            image.

            1. A n_points long array, with the cubic interpolated energy band

        ARGUMENTS

        n_points    :: The number of points for the interpolation

        """

        ds = self.path_distances
        # The arrays with the data points and the interpolated energy values
        x = np.linspace(0, ds[-1], n_points)
        E_interp = np.array([self._compute_polynomial_approximation_energy(i)
                             for i in x])

        return x, E_interp

    def _compute_polynomial_approximation_energy(self, x):

        """
        Return interpolated energy value for a point x
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
