import numpy as np
import warnings
from .integrators import BaseIntegrator
from .integrators import euler_step, runge_kutta_step
from itertools import cycle

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
    def __init__(self, band, forces, rhs_fun, action_fun, n_images, n_dofs_image,
                 maxSteps=1000,
                 maxCreep=5, actionTol=1e-10, forcesTol=1e-6,
                 etaScale=1.0, dEta=2, minEta=0.001,
                 # perturbSeed=42, perturbFactor=0.1,
                 nTrail=10, resetMax=20
                 ):
        super(FSIntegrator, self).__init__(band, rhs_fun)

        self.action_fun = action_fun

        # Integration parameters
        # TODO: move to run function?
        self.dEta = dEta
        self.minEta = minEta
        self.resetMax = resetMax
        self.maxCreep = maxCreep
        self.etaScale = etaScale
        self.forcesTol = forcesTol
        self.actionTol = actionTol
        self.maxSteps = maxSteps

        self.i_step = 0
        self.n_images = n_images
        self.n_dofs_image = n_dofs_image
        # self.forces_prev = np.zeros_like(band).reshape(n_images, -1)
        # self.G :
        self.forces = forces
        self.forces_old = np.zeros_like(forces)

        # Rename y to band
        self.band = self.y
        self.band_old = np.zeros_like(self.band)  # y -> band

        self.step = 0
        self.nTrail = nTrail

    def run_until(self, t):
        pass

    def run_for(self, n_steps):

        t = 0.0
        nStart = 0
        exitFlag = False
        totalRestart = True
        resetCount = 0
        creepCount = 0
        self.trailAction = np.zeros(self.nTrail)
        trailPool = cycle(range(self.nTrail))  # cycle through 0,1,...,(nTrail-1),0,1,...
        eta = 1.0

        # In __init__:
        # self.band_last[:] = self.band

        INNER_DOFS = slice(self.n_dofs_image, -self.n_dofs_image)

        while not exitFlag:

            if totalRestart:
                if self.step > 0:
                    print('Restarting')
                self.band[:] = self.band_old

                # Compute from self.band. Do not update the step at this stage:
                # This step updates the forces in the G array of the nebm module,
                # using the current band state self.y
                # TODO: remove time from chain method rhs
                #       make a specific function to update G??
                self.rhs(t, self.y)
                self.action = self.action_fun()

                # self.step += 1
                self.forces_old[:] = self.forces  # Scale field??
                # self.gradE_last[~_material] = 0.0
                # self.gradE[:] = self.gradE_last
                self.action_old = self.action
                self.trailAction[nStart] = self.action
                nStart = next(trailPool)
                eta = 1.0
                totalRestart = False

            creepCount = 0

            # Creep stage: minimise with a fixed eta
            while creepCount < self.maxCreep:
                # Update spin. Avoid pinned or zero-Ms sites
                self.band[:] = self.band_old + eta * self.etaScale * self.forces
                normalise_spins(self.band)

                self.trailAction

                self.rhs(t, self.y)
                self.action = self.action_fun()

                self.trailAction[nStart] = self.action
                nStart = next(trailPool)

                # Getting averages of forces from the INNER images in the band (no extrema)
                # (forces are given by vector G in the chain method code)
                # TODO: we might use all band images, not only inner ones
                G_norms = np.linalg.norm(self.forces[INNER_DOFS].reshape(-1, 3), axis=1)
                # Compute the root mean square per image
                rms_G_norms_per_image = np.sqrt(np.mean(G_norms.reshape(self.n_images - 2, -1), axis=1))
                mean_rms_G_norms_per_image = np.mean(rms_G_norms_per_image)

                # Average step difference between trailing action and new action
                deltaAction = (np.abs(self.trailAction[nStart] - self.action)) / self.nTrail

                # 10 seems like a magic number; we set here a minimum number of evaulations
                if (nStart > self.nTrail * 10) and (deltaAction < self.actionTol):
                    print('Change in action is negligible')
                    exitFlag = True
                    break  # creep loop

                if (n_steps >= self.maxSteps):
                    print('Number of steps reached maximum')
                    exitFlag = True
                    break  # creep loop

                # If action increases compared to last step, decrease eta and start creeping again
                if (self.action_old < self.action):
                    print('Action increased. Start new creep stage from last step')
                    creepCount = 0
                    eta = eta / (self.dEta * self.dEta)

                    # If eta is too small, reset and start again
                    if (eta < 1e-3):
                        print('')
                        resetCount += 1
                        bestAction = self.action_old
                        RefinePath(self.band)  # Resets the path to equidistant structures (smoothing  kinks?)
                        # PathChanged[:] = True

                        if resetCount > self.resetMax:
                            print('Failed to converge!')
                    # Otherwise, just 
                    else:
                        self.band[:] = self.band_old
                        self.forces[:] = self.forces_old
                # If action decreases, move to next creep step
                else:
                    creepCount += 1
                    self.action_old = self.action
                    self.band_old[:] = self.band
                    self.forces_old[:] = self.forces

                    if (mean_rms_G_norms_per_image < self.forcesTol):
                        print('Change of mean of the force RMSquares negligible')
                        exitFlag = True
                        break  # creep loop
            # END creep while loop

            # After creep loop:
            self.eta = self.eta * self.dEta
            resetCount = 0

    def set_options(self):
        pass

    def _step(self, t, y, h, f):
        """
        """
        pass

def normalise_spins(y):
    # Normalise an array of spins y with 3 * N elements
    y.shape = (-1, 3)
    n = np.sqrt(y[:, 0] ** 2 + y[:, 1] ** 2 + y[:, 2] ** 2)
    fltr = n != 0.0
    y[fltr] = y[fltr] / n[:, np.newaxis][fltr]
    y.shape = (-1,)
