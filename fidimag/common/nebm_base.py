from __future__ import print_function
from __future__ import division
import numpy as np
import os

import fidimag.extensions.cvode as cvode
from fidimag.common.vtk import VTK
# from .nebm_tools import compute_norm
# from .nebm_tools import linear_interpolation_spherical
from .fileio import DataSaver

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="fidimag")


class NEBMBase(object):
    """

    NEB Method equations:

        G =

    ARGUMENTS

    sim                 :: An instance of a micromagnetic or an atomistic
                           simulation

    initial_images      :: A sequence of arrays or functions to set up the
                           magnetisation field profile

    inteprolations      ::

    dof                 :: Degrees of freedom per spin. Spherical coordinates
                           have dof=2 and Cartesian have dof=3

    """
    def __init__(self, sim,
                 initial_images, interpolations=None,
                 k=1e5,
                 name='unnamed',
                 dof=2
                 ):

        # Degrees of Freedom per spin
        self.dof = dof

        self.sim = sim
        self.mesh = self.sim.mesh
        self.name = name

        # Number of spins in the system
        self.n_spins = len(self.mesh.coordinates)

        # Spring constant (we could use an array in the future)
        self.k = k

        # VTK saver for the magnetisation/spin field --------------------------
        self.VTK = VTK(self.mesh,
                       directory='vtks'.format(self.name),
                       filename='m'
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
            self.interpolations = [0 for i in len(initial_images) - 1]

        # Number of images with/without the extremes
        self.n_images = len(self.initial_images) + np.sum(self.interpolations)
        self.n_images_inner_band = self.n_images - 2

        # Number of degrees of freedom per image
        self.n_dofs_image = (self.dof * self.n_spins)

        # Total number of degrees of freedom in the NEBM band
        self.n_band = self.n_images * self.n_dofs_image

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

        self.last_Y = np.zeros_like(self.band)

        # ---------------------------------------------------------------------

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
            try:
                self.VTK.save_scalar(self.sim.Ms, name='M_s')
            except:
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
            name = os.path.join(directory, 'image_%d.npy' % i)
            if coordinates_function:
                np.save(name, coordinates_function(self.band[i, :]))
            else:
                np.save(name, self.band[i])
        self.band.shape = (-1)

    def initialise_integrator(self, rtol=1e-6, atol=1e-6):
        self.t = 0
        self.iterations = 0
        self.ode_count = 1

        self.integrator = cvode.CvodeSolver(self.band, self.Sundials_RHS)
        self.integrator.set_options(rtol, atol)

    def create_tablewriter(self):
        entities_energy = {
            'step': {'unit': '<1>',
                     'get': lambda sim: sim.iterations,
                     'header': 'iterations'},
            'energy': {'unit': '<J>',
                       'get': lambda sim: sim.energies,
                       'header': ['image_%d' % i for i in range(self.n_images)]}
        }

        self.tablewriter = DataSaver(
            self, '%s_energy.ndt' % (self.name),  entities=entities_energy)

        entities_dm = {
            'step': {'unit': '<1>',
                     'get': lambda sim: sim.iterations,
                     'header': 'iterations'},
            'dms': {'unit': '<1>',
                    'get': lambda sim: sim.distances,
                    'header': ['image_%d_%d' % (i, i + 1) for i in range(self.n_images -1)]}
        }

        self.tablewriter_dm = DataSaver(
            self, '%s_dms.ndt' % (self.name), entities=entities_dm)

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
        dYdt = self.compute_distances(
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

        """

        log.debug("Relaxation parameters: "
                  "stopping_dmdt={} (degrees per nanosecond), "
                  "time_step={} s, max_iterations={}.".format(stopping_dYdt,
                                                              dt,
                                                              max_iterations))

        self.save_VTKs(coordinates_function=self.files_convert_f)
        self.save_npys(coordinates_function=self.files_convert_f)

        # Save the initial state i=0
        self.distances = self.compute_distances(self.band[self.n_dofs_image:],
                                                self.band[:-self.n_dofs_image]
                                                )
        self.tablewriter.save()
        self.tablewriter_dm.save()

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

            max_dYdt = self.run_until(self.t + increment_dt)

            # Save data -------------------------------------------------------

            if self.iterations % save_vtks_every == 0:
                self.save_VTKs(coordinates_function=self.files_convert_f)
            if self.iterations % save_npys_every == 0:
                self.save_npys(coordinates_function=self.files_convert_f)

            self.distances = self.compute_distances(
                self.band[self.n_dofs_image:],
                self.band[:-self.n_dofs_image]
                )
            self.tablewriter.save()
            self.tablewriter_dm.save()
            log.debug("step: {:.3g}, step_size: {:.3g}"
                      " and max_dYdt: {:.3g}.".format(self.iterations,
                                                      increment_dt,
                                                      max_dYdt)
                      )

            # -----------------------------------------------------------------

            # Stop criteria:
            if max_dYdt < stopping_dYdt:
                break

        log.info("Relaxation finished at time step = {:.4g}, "
                 "t = {:.2g}, call rhs = {:.4g} "
                 "and max_dYdt = {:.3g}".format(self.iterations,
                                                self.t,
                                                self.ode_count,
                                                max_dYdt)
                 )

        self.save_VTKs()

    # -------------------------------------------------------------------------
