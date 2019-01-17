import numpy as np
import warnings
from .integrators import BaseIntegrator
from .integrators import euler_step, runge_kutta_step

EPSILON = 1e-16


class StepIntegrator(BaseIntegrator):
    """
    A simple integrator where spins are normalised at every inetegrator step
    Integrator options are Euler and RK4
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
            raise NotImplemented("step must be euler or rk4")
        self._step = step_choices[step]


class VerletIntegrator(BaseIntegrator):
    """
    A quick Verlet integration in Cartesian coordinates
    See: J. Chem. Theory Comput., 2017, 13 (7), pp 3250–3259
    """
    def __init__(self, spins, rhs_fun, m=0.1, stepsize=1e-15):
        super(VerletIntegrator, self).__init__(spins, rhs_fun)

        self.m = m
        self.stepsize = stepsize
        self.velocity = np.zeros_like(spins).reshape(-1, 3)
        self.velocity_new = np.zeros_like(spins).reshape(-1, 3)
        # self.velocity_proj = np.zeros(len(spins) // 3)

    def run_until(self, t):
        while abs(self.t - t) > EPSILON:
            self.t, self.y = self._step(self.t, self.y,
                                        self.stepsize, self.rhs)

            self.rhs_evals_nb += 0
            if self.t > t:
                break
        return 0

    def set_options(self, rtol=1e-8, atol=1e-8):
        warnings.warn("Tolerances not available for VerletIntegrator")

    def _step(self, t, y, h, f):
        """
        Quick-min Verlet step
        """
        # In this case f represents the force: a = dy/dt = f/m
        force = f(t, y).reshape(-1, 3)  # * self.m_inv[:, np.newaxis]
        y.shape = (-1, 3)

        velocity_proj = np.einsum('ij,ij->i',
                                  force, self.velocity)
        self.velocity_new[velocity_proj <= 0] = 0.0

        force_norm_2 = np.einsum('ij,ij->i', force, force)
        fltr = np.logical_and(velocity_proj > 0, force_norm_2 > 0)
        factor = np.zeros_like(velocity_proj)
        factor[fltr] = velocity_proj[fltr] / force_norm_2[fltr]
        self.velocity_new[fltr] = np.einsum('i,ij->ij', factor[fltr], force[fltr])

        self.velocity = self.velocity_new + (h / m) * force

        yp = y + h * (self.velocity + (h / (2 * m)) * force)

        y.shape = (-1,)

        yp.shape = (-1,)
        normalise_spins(yp)

        tp = t + h
        return tp, yp


def normalise_spins(y):
    # Normalise an array of spins y with 3 * N elements
    y.shape = (-1, 3)
    n = np.sqrt(y[:, 0] ** 2 + y[:, 1] ** 2 + y[:, 2] ** 2)
    fltr = n != 0.0
    y[fltr] = y[fltr] / n[:, np.newaxis][fltr]
    y.shape = (-1,)
