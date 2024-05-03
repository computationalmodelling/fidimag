import numpy as np
import warnings
from .integrators import BaseIntegrator
from .integrators import euler_step, runge_kutta_step

import fidimag.extensions.nebm_clib as nebm_clib

EPSILON = 1e-16


class StepIntegrator(BaseIntegrator):
    """
    A simple integrator where spins are normalised at every integrator step.
    Integrator options are `euler` and `rk4`
    """
    def __init__(self, spins, rhs_fun, step="euler", stepsize=1e-15):
        super(StepIntegrator, self).__init__(spins, rhs_fun)

        self.set_step(step)
        self.stepsize = stepsize

    def run_until(self, t):
        while abs(self.t - t) > EPSILON:
            self.t, self.y, evals = self._step(self.t, self.y,
                                               self.stepsize, self.rhs)
            normalise_spins(self.y)

            self.rhs_evals_nb += evals
            if self.t > t:
                break
        return 0

    def set_options(self, rtol=1e-8, atol=1e-8):
        warnings.warn("Tolerances not available for StepIntegrator")

    def set_step(self, step):
        step_choices = {'euler': euler_step, 'rk4': runge_kutta_step}
        if step not in step_choices:
            raise NotImplementedError("step must be euler or rk4")
        self._step = step_choices[step]


class VerletIntegrator(BaseIntegrator):
    """A quick Verlet integration in Cartesian coordinates

    See: J. Chem. Theory Comput., 2017, 13 (7), pp 3250–3259
    """
    def __init__(self, band, forces, rhs_fun, n_images, n_dofs_image,
                 mass=0.1, stepsize=1e-15):
        super(VerletIntegrator, self).__init__(band, rhs_fun)

        self.n_images = n_images
        self.n_dofs_image = n_dofs_image
        self.mass = mass
        self.stepsize = stepsize
        self.velocity = np.zeros_like(band).reshape(n_images, -1)
        self.velocity_new = np.zeros_like(band).reshape(n_images, -1)
        self.forces_prev = np.zeros_like(band).reshape(n_images, -1)
        # self.G :
        self.forces = forces

    def run_until(self, t):
        while abs(self.t - t) > EPSILON:
            self.t = self._step(self.t, self.y,
                                self.stepsize, self.rhs)

            self.rhs_evals_nb += 0

            # If we could make the C function to work:
            # nebm_clib.step_Verlet(
            #     self.forces,
            #     self.forces_prev,
            #     self.velocity,
            #     self.velocity_new,
            #     self.y,
            #     self.t,
            #     self.stepsize,
            #     self.mass,
            #     self.n_images,
            #     self.n_dofs_image,
            #     self.rhs
            #     )

            if self.t > t:
                break
        return 0

    def set_options(self, rtol=1e-8, atol=1e-8):
        warnings.warn("Tolerances not available for VerletIntegrator")

    def _step(self, t, y, h, f):
        """
        Quick-min Velocity Verlet step
        """

        f(t, y)
        force_images = self.forces
        force_images.shape = (self.n_images, -1)
        # In this case f represents the force: a = dy/dt = f/m
        # * self.m_inv[:, np.newaxis]
        y.shape = (self.n_images, -1)
        # print(force_images[2])

        # Loop through every image in the band
        for i in range(1, self.n_images - 1):

            force = force_images[i]
            velocity = self.velocity[i]
            velocity_new = self.velocity_new[i]

            # Update coordinates from Newton eq, which uses the "acceleration"
            # At step 0 velocity is zero
            y[i][:] = y[i] + h * (velocity + (h / (2 * self.mass)) * force)

            # Update the velocity from a mean with the prev step
            # (Velocity Verlet)
            velocity[:] = velocity_new[:] + (h / (2 * self.mass)) * (self.forces_prev[i] + force)

            # Project the force of the image into the velocity vector: < v | F >
            velocity_proj = np.einsum('i,i', force, velocity)

            # Set velocity to zero when the proj = 0
            if velocity_proj <= 0:
                velocity_new[:] = 0.0
            else:
                # Norm of the force squared <F | F>
                force_norm_2 = np.einsum('i,i', force, force)
                factor = velocity_proj / force_norm_2
                # Set updated velocity: v = v * (<v|F> / |F|^2)
                velocity_new[:] = factor * force

            # New velocity from Newton equation (old Verlet)
            # velocity[:] = velocity_new + (h / self.mass) * force

        # Store the force for the Velocity Verlet algorithm
        self.forces_prev[:] = force_images

        y.shape = (-1,)
        force_images.shape = (-1,)
        normalise_spins(y)
        tp = t + h
        return tp


# -----------------------------------------------------------------------------
# NEBM integrator for F-S algorithm


class FSIntegrator(BaseIntegrator):
    """A step integrator considering the action of the band
    """
    def __init__(self, band, forces, action, rhs_fun, n_images, n_dofs_image,
                 max_steps=1000,
                 maxCreep=5, eta_scale=1.0, stopping_dE=1e-6, dEta=2,
                 etaMin=0.001,
                 # perturbSeed=42, perturbFactor=0.1,
                 nTrail=10, resetMax=20, mXgradE_tol=0.1
                 ):
        super(FSIntegrator, self).__init__(band, rhs_fun)

        self.i_step = 0
        self.n_images = n_images
        self.n_dofs_image = n_dofs_image
        self.forces_prev = np.zeros_like(band).reshape(n_images, -1)
        # self.G :
        self.forces = forces
        self.max_steps = max_steps

        self.y_last = np.zeros_like(self.y)  # y -> band
        self.step = 0
        self.nTrail = nTrail

    def run_until(self, t):
        pass

    def run_for(self, n_steps):

        nStart = 0
        exitFlag = False
        totalRestart = True
        resetCount = 0
        creepCount = 0
        self.trailE = np.zeros(nTrail)
        trailPool = cycle(range(nTrail))  # cycle through 0,1,...,(nTrail-1),0,1,...
        eta = 1.0

        while not exitFlag:

            if totalRestart:
                if self.step > 0:
                    print('Restarting')
                self.y[:] = self.y_last

                # Compute from self.band. Do not update the step at this stage:
                # This step updates the forces in the G array of the nebm module,
                # using the current band state self.y
                self.rhs(t, self.y)


                # self.step += 1
                self.gradE_last[:] = -self.field  # Scale field??
                self.gradE_last[~_material] = 0.0
                self.gradE[:] = self.gradE_last
                self.totalE_last = self.totalE
                self.trailE[nStart] = self.totalE
                nStart = next(trailPool)
                eta = 1.0
                totalRestart = False

            creepCount = 0

            # Creep stage: minimise with a fixed eta
            while creepCount < maxCreep:
                # Update spin. Avoid pinned or zero-Ms sites
                self.y[:] = self.y_last - eta * eta_scale * self.forces_images






        # while abs(self.i_step - steps) > EPSILON:
        st = 1
        while st < n_steps:
            self.i_step = self._step(self.t, self.y, self.stepsize, self.rhs)

            self.rhs_evals_nb += 0

            if self.i_step > self.max_steps:
                break

            st += 1
        return 0

    def set_options(self):
        pass

    def _step(self, t, y, h, f):
        """
        """

        f(t, y)
        force_images = self.forces
        force_images.shape = (self.n_images, -1)
        # In this case f represents the force: a = dy/dt = f/m
        # * self.m_inv[:, np.newaxis]
        y.shape = (self.n_images, -1)
        # print(force_images[2])

        # Loop through every image in the band
        for i in range(1, self.n_images - 1):

            force = force_images[i]
            velocity = self.velocity[i]
            velocity_new = self.velocity_new[i]

            # Update coordinates from Newton eq, which uses the "acceleration"
            # At step 0 velocity is zero
            y[i][:] = y[i] + h * (velocity + (h / (2 * self.mass)) * force)

            # Update the velocity from a mean with the prev step
            # (Velocity Verlet)
            velocity[:] = velocity_new[:] + (h / (2 * self.mass)) * (self.forces_prev[i] + force)

            # Project the force of the image into the velocity vector: < v | F >
            velocity_proj = np.einsum('i,i', force, velocity)

            # Set velocity to zero when the proj = 0
            if velocity_proj <= 0:
                velocity_new[:] = 0.0
            else:
                # Norm of the force squared <F | F>
                force_norm_2 = np.einsum('i,i', force, force)
                factor = velocity_proj / force_norm_2
                # Set updated velocity: v = v * (<v|F> / |F|^2)
                velocity_new[:] = factor * force

            # New velocity from Newton equation (old Verlet)
            # velocity[:] = velocity_new + (h / self.mass) * force

        # Store the force for the Velocity Verlet algorithm
        self.forces_prev[:] = force_images

        y.shape = (-1,)
        force_images.shape = (-1,)
        normalise_spins(y)
        tp = t + h
        return tp


def normalise_spins(y):
    # Normalise an array of spins y with 3 * N elements
    y.shape = (-1, 3)
    n = np.sqrt(y[:, 0] ** 2 + y[:, 1] ** 2 + y[:, 2] ** 2)
    fltr = n != 0.0
    y[fltr] = y[fltr] / n[:, np.newaxis][fltr]
    y.shape = (-1,)
