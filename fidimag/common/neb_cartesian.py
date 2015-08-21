import os
import fidimag.extensions.cvode as cvode
import numpy as np
from fileio import DataSaver, DataReader
from save_vtk import SaveVTK
import fidimag.extensions.neb_clib as neb_clib

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="fidimag")


def cartesian2spherical(xyz):
    """
    Transform cartesian coordinates (x, y, z)
    in spherical coordinates. The function only returns
    the (theta, phi) pair since
    the magnetisation is fixed at zero Temperature
    (the r-component is constant) and is
    fully characterised by two degrees of freedom.
    (We use this to specifically transform M coordinates)
    Are we asuming? that xyz is:
     [x1, y1, z1, x2, y2, z2, x3 ... ]
    """
    # Transform to a -- x 3 array
    xyz.shape = (-1, 3)
    r_xy = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2)
    theta = np.arctan2(r_xy, xyz[:, 2])
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    xyz.shape = (-1,)

    theta_phi = np.concatenate((theta, phi))
    # Return [theta1, theta2, ... , phi1, phi2, ... ]
    return theta_phi


def spherical2cartesian(theta_phi):
    """
    Returns the (x, y, z) cartesian components
    from spherical coordinates (theta, phi)
    for a r-component equal to 1
    (thus, (x,y,z) is normalised)
    We are assuming that theta_phi has the form

    [theta1, theta2, ... , phi1, phi2, ...]

    as in the previous function
    """
    theta_phi.shape = (2, -1)
    theta = theta_phi[0]
    phi = theta_phi[1]
    # print 't=', theta
    # print 'phi=', phi
    mxyz = np.zeros(3 * len(theta))
    mxyz.shape = (-1, 3)
    mxyz[:, 0] = np.sin(theta) * np.cos(phi)
    mxyz[:, 1] = np.sin(theta) * np.sin(phi)
    mxyz[:, 2] = np.cos(theta)
    mxyz.shape = (-1, )
    # Return [mx1, my1, mz1, mx2, ...]
    theta_phi.shape = (-1, )
    return mxyz


def linear_interpolation_two(m0, m1, n, pin_ids):
    """
    Define a linear interpolation between
    two states of the energy band (m0, m1) to get
    an initial state. The interpolation is
    done in the magnetic moments that constitute the
    magnetic system.

    ** The interpolation has phi undefined when x=y=0 and
    the spin direction is in a pole (theta=0 or pi), thus it
    is suggested to give a very small deviation

    """
    theta_phi0 = cartesian2spherical(m0)
    theta_phi1 = cartesian2spherical(m1)
    m0.shape = (-1, 3)

    dtheta = (theta_phi1 - theta_phi0) / (n + 1)
    coords = []
    for i in range(n):
        theta_phi = theta_phi0 + (i + 1) * dtheta
        new_m = spherical2cartesian(theta_phi)
        new_m.shape = (-1, 3)
        new_m[pin_ids, :] = m0[pin_ids, :]
        new_m.shape = (-1,)
        coords.append(new_m)
    m0.shape = (-1,)
    return coords


def normalise_m(a):
    """
    Normalise the magnetisation array.
    We asume:

    ******* FIX
    a = [mx1, mx2, ..., my1, my2, ..., mz1, mz2, ...]
    to transform this into

    [ [mx1, mx2, ...],
      [my1, my2, ...],
      [mz1, mz2, ...]
    ]

    *******

    normalise the matrix, and return again a  1 x -- array
    """
    # Transform to matrix
    a.shape = (-1, 3)
    # Compute the array 'a' length
    lengths = np.sqrt(a[:, 0] * a[:, 0] + a[:, 1] * a[:, 1] + a[:, 2] * a[:, 2])
    # Normalise all the entries
    # a[:] /= lengths
    a.T[:, 0] /= lengths
    a.T[:, 1] /= lengths
    a.T[:, 2] /= lengths
    # Return to original shape
    a.shape = (-1, )


def compute_dm(m0, m1):

    dm = m0 - m1
    length = len(dm)
    dm = np.sqrt(np.sum(dm ** 2)) / length
    return dm


class NEB_Sundials(object):

    """
    Nudged elastic band method by solving the differential equation using Sundials.
    """

    def __init__(self, sim,
                 initial_images,
                 climbing_image=None,
                 interpolations=None,
                 spring=5e5, name='unnamed'):
        """
          *Arguments*

              sim: the Simulation class

              initial_images: a list contain the initial value, which can have
              any of the forms accepted by the function 'finmag.util.helpers.
              vector_valued_function', for example,

                    initial_images = [(0,0,1), (0,0,-1)]

              or with given defined function

                    def init_m(pos):
                        x=pos[0]
                        if x<10:
                            return (0,1,1)
                        return (-1,0,0)

                  initial_images = [(0,0,1), (0,0,-1), init_m ]

              are accepted forms.


              climbing_image : An integer with the index (from 1 to the total
              number of images minus two; it doesn't have any sense to use the
              extreme images) of the image with the largest energy, which will
              be updated in the NEB algorithm using the Climbing Image NEB
              method (no spring force and "with the component along the elastic
              band inverted" [*]). See: [*] Henkelman et al., The Journal of
              Chemical Physics 113, 9901 (2000)

              interpolations : a list only contain integers and the length of
              this list should equal to the length of the initial_images minus
              1, i.e., len(interpolations) = len(initial_images) - 1 ** THIS IS
              not well defined in CARTESIAN coordinates**

              spring: the spring constant, a float value

              disable_tangent: this is an experimental option, by disabling the
              tangent, we can get a rough feeling about the local energy minima
              quickly.

        """

        self.sim = sim
        self.name = name
        self.spring = spring

        # We set a minus one because the *sundials_rhs* function
        # only uses an array without counting the extreme images,
        # whose length is self.image_num (see below)
        if climbing_image is not None:
            self.climbing_image = climbing_image - 1
        else:
            self.climbing_image = climbing_image

        # Dolfin function of the new _m_field (instead of _m)
        #self.effective_field = sim.llg.effective_field

        if interpolations is None:
            interpolations = [0 for i in range(len(initial_images) - 1)]

        self.initial_images = initial_images
        self.interpolations = interpolations

        if len(interpolations) != len(initial_images) - 1:
            raise RuntimeError("""The length of interpolations should be equal to
                the length of the initial_images array minus 1, i.e.,
                len(interpolations) = len(initial_images) - 1""")

        if len(initial_images) < 2:
            raise RuntimeError("""At least two images must be provided
                               to create the energy band""")

        # the total image number including two ends
        self.total_image_num = len(initial_images) + sum(interpolations)
        self.image_num = self.total_image_num - 2

        self.nxyz = sim.nxyz

        self.coords = np.zeros(3 * self.nxyz * self.total_image_num)
        self.last_m = np.zeros(self.coords.shape)

        self.Heff = np.zeros(self.coords.shape)
        self.Heff.shape = (self.total_image_num, -1)

        self.tangents = np.zeros(3 * self.nxyz * self.image_num)
        self.tangents.shape = (self.image_num, -1)

        self.energy = np.zeros(self.total_image_num)

        self.springs = np.zeros(self.image_num)

        self.pin_ids = np.array(
            [i for i, v in enumerate(self.sim.pins) if v > 0], dtype=np.int32)

        self.t = 0
        self.step = 0
        self.ode_count = 1
        self.integrator = None

        self.initial_image_coordinates()
        self.create_tablewriter()

    def create_tablewriter(self):
        entities_energy = {
            'step': {'unit': '<1>',
                     'get': lambda sim: sim.step,
                     'header': 'steps'},
            'energy': {'unit': '<J>',
                       'get': lambda sim: sim.energy,
                       'header': ['image_%d' % i for i in range(self.image_num + 2)]}
        }

        self.tablewriter = DataSaver(
            self, '%s_energy.ndt' % (self.name),  entities=entities_energy)

        entities_dm = {
            'step': {'unit': '<1>',
                     'get': lambda sim: sim.step,
                     'header': 'steps'},
            'dms': {'unit': '<1>',
                    'get': lambda sim: sim.distances,
                    'header': ['image_%d_%d' % (i, i + 1) for i in range(self.image_num + 1)]}
        }
        self.tablewriter_dm = DataSaver(
            self, '%s_dms.ndt' % (self.name), entities=entities_dm)

    def initial_image_coordinates(self):
        """

        Generate the coordinates linearly according to the number of
        interpolations provided.

        Example: Imagine we have 4 images and we want 3 interpolations
        between every neighbouring pair, i.e  interpolations = [3, 3, 3]

        1. Imagine the initial states with the interpolation numbers
           and choose the first and second state

            0          1           2          3
            X -------- X --------- X -------- X
                  3          3           3

            2. Counter image_id is set to 0

            3. Set the image 0 magnetisation vector as m0 and append the
               values to self.coords[0]. Update the counter: image_id = 1 now

            4. Set the image 1 magnetisation values as m1 and interpolate
               the values between m0 and m1, generating 3 arrays
               with the magnetisation values of every interpolation image.
               For every array, append the values to self.coords[i]
               with i = 1, 2 and 3 ; updating the counter every time, so
               image_id = 4 now

            5. Append the value of m1 (image 1) in self.coords[4]
               Update counter (image_id = 5 now)

            6. Move to the next pair of images, now set the 1-th image
               magnetisation values as m0 and append to self.coords[5]

            7. Interpolate to get self.coords[i], for i = 6, 7, 8
               ...
            8. Repeat as before until move to the pair of images: 2 - 3

            9. Finally append the magnetisation of the last image
               (self.initial_images[-1]). In this case, the 3rd image

        Then, for every magnetisation vector values array (self.coords[i])
        append the value to the simulation and store the energies
        corresponding to every i-th image to the self.energy[i] arrays

        Finally, flatten the self.coords matrix (containing the magnetisation
        values of every image in different rows)


        ** Our generalised coordinates in the NEBM are the magnetisation values

        """

        # Initiate the counter
        image_id = 0
        self.coords.shape = (self.total_image_num, -1)

        # For every interpolation between images (zero if no interpolations
        # were specified)
        for i in range(len(self.interpolations)):
            # Store the number
            n = self.interpolations[i]

            # Save on the first image of a pair (step 1, 6, ...)
            self.sim.set_m(self.initial_images[i])
            m0 = self.sim.spin.copy()

            self.coords[image_id][:] = m0[:]
            image_id = image_id + 1

            # Set the second image in the pair as m1 and interpolate
            # (step 4 and 7), saving in corresponding self.coords entries
            self.sim.set_m(self.initial_images[i + 1])
            m1 = self.sim.spin.copy()
            # Interpolations (arrays with magnetisation values)
            coords = linear_interpolation_two(m0, m1, n, self.pin_ids)

            for coord in coords:
                self.coords[image_id][:] = coord[:]
                image_id = image_id + 1

            # Continue to the next pair of images

        # Append the magnetisation of the last image
        self.sim.set_m(self.initial_images[-1])
        m2 = self.sim.spin
        self.coords[image_id][:] = m2[:]

        # Save the energies
        for i in range(self.total_image_num):
            self.sim.spin[:] = self.coords[i][:]
            self.sim.compute_effective_field(t=0)
            self.energy[i] = self.sim.compute_energy()

        # Flatten the array
        self.coords.shape = (-1,)

    def add_noise(self, T=0.1):
        noise = T * np.random.rand(self.total_image_num, 3, self.nxyz)
        noise[:, :, self.pin_ids] = 0
        noise[0, :, :] = 0
        noise[-1, :, :] = 0
        noise.shape = (-1,)
        self.coords += noise

    def save_vtks(self):
        """
        Save vtk files in different folders, according to the
        simulation name and step.
        Files are saved as simname_simstep/image_00000x.vtu
        """

        # Create the directory
        directory = 'vtks/%s_%d' % (self.name, self.step)

        self.vtk = SaveVTK(self.sim.mesh, directory)

        self.coords.shape = (self.total_image_num, -1)

        for i in range(self.total_image_num):
            self.vtk.save_vtk(self.coords[i], step=i, vtkname='m')

        self.coords.shape = (-1, )

    def save_npys(self):
        """
        Save npy files in different folders according to
        the simulation name and step
        Files are saved as: simname_simstep/image_x.npy
        """
        # Create directory as simname_simstep
        directory = 'npys/%s_%d' % (self.name, self.step)

        if not os.path.exists(directory):
            os.makedirs(directory)

        # Save the images with the format: 'image_{}.npy'
        # where {} is the image number, starting from 0
        self.coords.shape = (self.total_image_num, -1)
        for i in range(self.total_image_num):
            name = os.path.join(directory, 'image_%d.npy' % i)
            np.save(name, self.coords[i, :])

        self.coords.shape = (-1, )

    def create_integrator(self, rtol=1e-6, atol=1e-6, nsteps=10000):

        self.integrator = cvode.CvodeSolver(self.coords, self.sundials_rhs)

        self.integrator.set_options(rtol, atol)

    def compute_effective_field(self, y):

        y.shape = (self.total_image_num, -1)

        for i in range(self.image_num):

            self.sim.spin[:] = y[i + 1][:]
            #
            self.sim.compute_effective_field(t=0)
            # Compute effective field, which is the gradient of
            # the energy in the NEB method (derivative with respect to
            # the generalised coordinates)
            h = self.sim.field
            #
            self.Heff[i + 1, :] = h[:]
            # Compute the total energy
            self.energy[i + 1] = self.sim.compute_energy()

            # Compute the 'distance' or difference between neighbouring states
            # around y[i+1]. This is used to compute the spring force
            #
            dm1 = compute_dm(y[i], y[i + 1])
            dm2 = compute_dm(y[i + 1], y[i + 2])
            self.springs[i] = self.spring * (dm2 - dm1)

        # Use the native NEB (C++ code) to compute the tangents according
        # to the improved NEB method, developed by Henkelman and Jonsson
        # at: Henkelman et al., Journal of Chemical Physics 113, 22 (2000)

        #native_neb.compute_tangents(y, self.energy, self.tangents)
        neb_clib.compute_tangents(
            y, self.energy, self.tangents, self.total_image_num, 3 * self.nxyz)
        # native_neb.compute_springs(y,self.springs,self.spring)
        y.shape = (-1, )

    def sundials_rhs(self, time, y, ydot):
        """

        Right hand side of the optimization scheme used to find the minimum
        energy path. In our case, we use a LLG kind of equation:

            d Y / dt = Y x Y x D

            D = -( nabla E + [nabla E * t] t ) + F_spring

        where Y is an image: Y = (M_0, ... , M_N) and t is the tangent vector
        defined according to the energy of the neighbouring images (see
        Henkelman et al publication)

        If a climbing_image index is specified, the corresponding image
        will be iterated without the spring force and with an inversed
        component along the tangent

        """
        # Update the ODE solver
        self.ode_count += 1

        # Compute the eff field H for every image, H = -nabla E
        # (derived with respect to M)
        self.compute_effective_field(y)

        # Reshape y and ydot in a matrix of total_image_num rows
        y.shape = (self.total_image_num, -1)
        ydot.shape = (self.total_image_num, -1)

        # Compute the total force for every image (not the extremes)
        # Rememeber that self.image_num = self.total_image_num - 2
        # The total force is:
        #   D = - (-nabla E + [nabla E * t] t) + F_spring
        # This value is different is a climbing image is specified:
        #   D_climb = -nabla E + 2 * [nabla E * t] t
        for i in range(self.image_num):
            h = self.Heff[i + 1]
            t = self.tangents[i]
            sf = self.springs[i]

            if not (self.climbing_image and i == self.climbing_image):
                h3 = h - np.dot(h, t) * t + sf * t
            else:
                h3 = h - 2 * np.dot(h, t) * t

            h[:] = h3[:]

            #ydot[i+1, :] = h3[:]

        # Update the step with the optimisation algorithm, in this
        # case we use: dY /dt = Y x Y x D
        # (check the C++ code in finmag/native/src/)
        neb_clib.compute_dm_dt(
            y, self.Heff, ydot, self.sim._pins, self.total_image_num, self.nxyz)

        ydot[0, :] = 0
        ydot[-1, :] = 0

        y.shape = (-1,)
        ydot.shape = (-1,)
        return 0

    def compute_distance(self):

        distance = []

        ys = self.coords
        ys.shape = (self.total_image_num, -1)
        for i in range(self.total_image_num - 1):
            dm = compute_dm(ys[i], ys[i + 1])
            distance.append(dm)

        ys.shape = (-1, )
        self.distances = np.array(distance)

    def run_until(self, t):

        if t <= self.t:
            return

        self.integrator.run_until(t)
        self.coords[:] = self.integrator.y[:]

        m = self.coords
        y = self.last_m

        m.shape = (self.total_image_num, -1)
        y.shape = (self.total_image_num, -1)

        max_dmdt = 0
        for i in range(1, self.image_num + 1):
            dmdt = compute_dm(y[i], m[i]) / (t - self.t)
            if dmdt > max_dmdt:
                max_dmdt = dmdt

        m.shape = (-1,)
        y.shape = (-1,)
        self.last_m[:] = m[:]
        self.t = t

        return max_dmdt

    def relax(self, dt=1e-8, stopping_dmdt=1e4,
              max_steps=1000, save_npy_steps=100,
              save_vtk_steps=100):

        if self.integrator is None:
            self.create_integrator()

        log.debug("Relaxation parameters: "
                  "stopping_dmdt={} (degrees per nanosecond), "
                  "time_step={} s, max_steps={}.".format(stopping_dmdt,
                                                         dt, max_steps))
        # Save the initial state i=0
        self.compute_distance()
        self.tablewriter.save()
        self.tablewriter_dm.save()

        for i in range(max_steps):

            if i % save_vtk_steps == 0:
                self.save_vtks()

            if i % save_npy_steps == 0:
                self.save_npys()

            self.step += 1

            cvode_dt = self.integrator.get_current_step()

            increment_dt = dt

            if cvode_dt > dt:
                increment_dt = cvode_dt

            dmdt = self.run_until(self.t + increment_dt)

            self.compute_distance()
            self.tablewriter.save()
            self.tablewriter_dm.save()
            log.debug("step: {:.3g}, step_size: {:.3g}"
                      " and max_dmdt: {:.3g}.".format(self.step,
                                                      increment_dt,
                                                      dmdt))

            if dmdt < stopping_dmdt:
                break

        log.info("Relaxation finished at time step = {:.4g}, "
                 "t = {:.2g}, call rhs = {:.4g} "
                 "and max_dmdt = {:.3g}".format(self.step,
                                                self.t,
                                                self.ode_count,
                                                dmdt))
        self.save_vtks()
        self.save_npys()
