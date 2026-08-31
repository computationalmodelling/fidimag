import numpy as np
import warnings
from .integrators import BaseIntegrator
from .integrators import euler_step, runge_kutta_step

import fidimag.extensions.nebm_clib as nebm_clib

EPSILON = 1e-16


class StepIntegrator(BaseIntegrator):
    """
    A simple integrator where spins are normalised at every inetegrator step
    Integrator options are Euler and RK4
    """
    def __init__(self, spins, rhs_fun, step="euler", stepsize=1e-15):
        super().__init__(spins, rhs_fun)

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

    def get_current_step(self):
        return self.stepsize

    def set_step(self, step):
        step_choices = {'euler': euler_step, 'rk4': runge_kutta_step}
        if step not in step_choices:
            raise NotImplementedError("step must be euler or rk4")
        self._step = step_choices[step]


class VerletIntegrator(BaseIntegrator):
    """
    Quick-min velocity-projection (VP) integration in Cartesian coordinates, with a *global* velocity
    projection across the whole band (all inner images together) rather than projecting each image
    independently -- same idea as SPIRIT's VP solver. See:

        Bessarab, Uzdin, Jonsson, "Method for finding mechanism and activation energy of magnetic
        transitions, applied to skyrmion and antivortex annihilation", Comp. Phys. Comm. 196, 335 (2015).

    Also based on the quick-min velocity Verlet scheme in J. Chem. Theory Comput., 2017, 13 (7), 3250-3259.

    A single global ratio <v|F> / <F|F>, summed over every degree of freedom of every (non-extreme)
    image, is used to project/reset the velocity of the *whole* band at once, rather than letting each
    image decide independently whether to keep moving. This couples the images' descent direction
    together (as in a single, high-dimensional quick-min minimisation of the whole band) and is
    noticeably less prone to getting stuck oscillating between images with locally inconsistent
    velocity projections than the previous per-image version.
    """
    def __init__(self, band, forces, rhs_fun, n_images, n_dofs_image,
                 mass=0.1, stepsize=1e-15):
        super().__init__(band, rhs_fun)

        self.n_images = n_images
        self.n_dofs_image = n_dofs_image
        self.mass = mass
        self.stepsize = stepsize
        self.velocity = np.zeros_like(band).reshape(n_images, -1)
        self.forces_prev = np.zeros_like(band).reshape(n_images, -1)
        # self.G :
        self.forces = forces

    def run_until(self, t):
        while abs(self.t - t) > EPSILON:
            self.t = self._step(self.t, self.y,
                                self.stepsize, self.rhs)

            self.rhs_evals_nb += 0

            if self.t > t:
                break
        return 0

    def set_options(self, rtol=1e-8, atol=1e-8):
        warnings.warn("Tolerances not available for VerletIntegrator")

    def get_current_step(self):
        return self.stepsize

    def _step(self, t, y, h, f):
        """
        Quick-min velocity-projection (VP) step, with the projection done
        globally over the whole band -- see the class docstring.
        """

        f(t, y)
        force_images = self.forces
        force_images.shape = (self.n_images, -1)
        y.shape = (self.n_images, -1)

        inner = slice(1, self.n_images - 1)

        # Leapfrog velocity update: average of the previous and current
        # force, for every inner image
        self.velocity[inner] += (h / (2 * self.mass)) * (
            self.forces_prev[inner] + force_images[inner])

        # Global projection of the velocity onto the force, summed over
        # every degree of freedom of every inner image (i.e. treating the
        # whole band as a single high-dimensional vector)
        projection_full = np.einsum('ij,ij', self.velocity[inner],
                                    force_images[inner])
        force_norm2_full = np.einsum('ij,ij', force_images[inner],
                                     force_images[inner])

        if projection_full <= 0 or force_norm2_full == 0:
            self.velocity[inner] = 0.0
        else:
            ratio = projection_full / force_norm2_full
            self.velocity[inner] = ratio * force_images[inner]

        # Update coordinates using the (just projected) velocity, plus the
        # usual explicit half-step force term
        y[inner] = (y[inner] + h * self.velocity[inner]
                   + (h ** 2 / (2 * self.mass)) * force_images[inner])

        # Store the force for the next step's leapfrog velocity update
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
