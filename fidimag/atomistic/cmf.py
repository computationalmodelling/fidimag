import logging
from dataclasses import dataclass

import numpy as np

import fidimag.common.constant as const
import fidimag.extensions.clib as clib
from fidimag.common.minimiser_base import MinimiserBase

from .anisotropy import Anisotropy, CubicAnisotropy

log = logging.getLogger(__name__)


@dataclass
class CMFResult:
    """
    What a classical mean-field relaxation did, returned by `minimise`
    """

    converged: bool

    reason: str

    n_iterations: int

    # Largest change of a component of a spin in the last iteration
    max_change: float

    # Sum of the energies of the interactions, at the mean spins
    total_energy: float


class ClassicalMeanField(MinimiserBase):
    """
    Classical mean-field theory of a lattice of spins at a temperature T
    (Hovorka and Sluckin, arXiv:2007.12777).

    Each spin is treated as a classical spin of unit length and moment
    mu_s in the field of the average of its neighbours. For a spin in the
    effective field B the thermal average is

        <S> = L(x) B / |B|,    x = mu_s |B| / (k_B T),    L(x) = coth(x) - 1/x

    The effective field is that of the interactions of the simulation, evaluated
    at the average spins, so any interaction that is bilinear in the spins or
    linear in a spin (exchange, DMI, Zeeman, dipolar) is included. The spins of
    the simulation become the average spins, whose length is at most 1 and
    shrinks as the temperature goes up, and `minimise` finds the self-consistent
    solution by iterating

        S <- S + mixing * (target(S) - S)

    or with Anderson mixing, using the last `history` iterates.

    A single-ion anisotropy is not a pair interaction, so the average of a spin
    in it is not along its field and the same equation does not hold; these
    interactions are refused.

    The iteration finds a self-consistent solution, which is not always the one
    of lowest free energy. Which one it reaches depends on the starting spins,
    which is what makes hysteresis visible, and near T_c or in a textured phase
    it can take many thousands of iterations: `minimise` reports whether it
    converged.
    """

    def __init__(self, mesh, spin,
                 magnetisation, magnetisation_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac=False,
                 integrator=None
                 ):

        super().__init__(mesh, spin,
                         magnetisation, magnetisation_inv, field, pins,
                         interactions, name, data_saver)

        self.t = 0
        self.step = 0
        self.T = 0.0
        self.resid = np.zeros_like(self.spin)

    def _check_interactions(self):
        for obj in self.interactions:
            if isinstance(obj, (Anisotropy, CubicAnisotropy)):
                raise NotImplementedError(
                    "The mean-field driver does not support the single-ion "
                    "anisotropy '{}': the thermal average of a spin in it is "
                    "not along its field. Remove it, or use a driver that "
                    "does not average over the temperature".format(obj.name))

    def total_energy(self):
        return sum(obj.compute_energy() for obj in self.interactions)

    def iterate(self, T, mixing=0.5):
        """One linear mixing step at the temperature T, returning the largest
        change of a component of a spin"""

        self.compute_effective_field()
        change = clib.compute_cmf_residual(self.spin, self.field,
                                           self._magnetisation, self._pins,
                                           self.resid, T, const.k_B)
        clib.cmf_mix_spins(self.spin, self.resid, mixing)
        return mixing * change

    def minimise(self, T, stopping_dm=1e-8, max_steps=20000, mixing=0.5,
                 history=0, save_m_steps=None, save_vtk_steps=None,
                 log_steps=1000):
        """
        Relax the spins to the self-consistent mean-field solution at the
        temperature `T` (kelvin), starting from the current spins.

        The iteration stops when no component of a spin changes by more than
        `stopping_dm`. `mixing` is the fraction of the new average taken in
        each step (1 is the plain iteration, smaller values are more robust).
        With `history` > 0 Anderson mixing is used, which combines the last
        `history` iterates (up to 32) and often needs far fewer steps, in
        particular close to T_c.

        Returns a `CMFResult`; it does not raise if the iteration did not
        converge in `max_steps`.
        """

        if T < 0:
            raise ValueError("The temperature has to be positive, not {}".format(T))
        if not 0 < mixing <= 1:
            raise ValueError("mixing has to be in (0, 1], not {}".format(mixing))
        if history < 0 or history > 32:
            raise ValueError("history has to be between 0 and 32")
        self._check_interactions()
        self.T = T

        n3 = self.spin.size
        if history:
            X = np.zeros((history, n3))
            R = np.zeros((history, n3))
        nh = 0
        best = np.inf
        change = np.inf
        converged = False

        for it in range(1, max_steps + 1):
            self.compute_effective_field()
            max_resid = clib.compute_cmf_residual(self.spin, self.field,
                                                  self._magnetisation,
                                                  self._pins, self.resid,
                                                  T, const.k_B)
            if max_resid < stopping_dm:
                converged = True
                change = max_resid
                it -= 1
                break

            if history:
                # Start again if the residual rises a lot, which is where
                # the history is misleading
                if max_resid > 10 * best:
                    nh = 0
                    best = np.inf
                best = min(best, max_resid)
                if nh == history:
                    X[:-1], R[:-1] = X[1:].copy(), R[1:].copy()
                    nh -= 1
                X[nh], R[nh] = self.spin, self.resid
                nh += 1
                clib.cmf_anderson_step(self.spin, X, R, nh, mixing,
                                       self._magnetisation, self._pins)
            else:
                clib.cmf_mix_spins(self.spin, self.resid, mixing)
            change = mixing * max_resid
            self.step += 1

            if log_steps and it % log_steps == 0:
                log.info("cmf step {}: max residual {:.3e}".format(it, max_resid))
            if save_m_steps and it % save_m_steps == 0:
                self.save_m()
            if save_vtk_steps and it % save_vtk_steps == 0:
                self.save_vtk()

        if converged:
            reason = "converged: the largest change of a spin component is " \
                     "{:.3e}".format(change)
        else:
            reason = "reached max_steps = {} with a residual of {:.3e}".format(
                max_steps, change / mixing)
            log.warning("Mean-field iteration did not converge: " + reason)

        return CMFResult(converged=converged, reason=reason, n_iterations=it,
                         max_change=change, total_energy=self.total_energy())
