import logging
import numpy as np
# import fidimag.extensions.clib as clib
# import fidimag.common.constant as const
from itertools import cycle

from .minimiser_base import MinimiserBase, MinimiserResult

log = logging.getLogger(name='fidimag')


class HubertMinimiser(MinimiserBase):
    """ Gradient descent energy minimisation algorithm.

    A minimisation, energy gradient descent, algorithm, that dynamically
    updates a scaling factor η for the gradient of the energy::

        m_new = m_old + η * η_s * δE/δm|_old

    where η_s is a fixed scaling factor. The E functional gradient is computed
    from the effective field `δE/δm = - H_eff` (note that we are ignoring
    scaling parameters, such as mu0, Ms or mu_s).

    This algorithm is based on the works [1, 2] and is implemented in [3] as
    *Hubert minimiser*. The modifications include using the E gradient instead
    of the magnetisation torque, and using different criteria for stopping the
    algorithm.

    The energy is stored for `t` number of steps during a creep stage, in the
    `trailE` array `[E0, E1 ... Et]`. If the energy decreases at this stage,
    with respect to the trail step, η is increased to accelerate the descent.
    Otherwise, η is decreased, with a given limit which, if it is reached, the
    minimisation is restarted. The saving of the trail energy is cyclic, i.e.
    if the current step reaches `t`, the next step will save `E` at `0`, and
    the energy difference at the current step is computed as `abs(Et - E0)/t`.
    The energy difference is scaled by the length of the trailing energy array,
    `t`.

    The energy in this minimisation class is scaled by the `self.energyScale`
    parameter, so define the `stopping_dE` argument in `minimise` accordingly.
    A more robust and global parameter to stop the minimisation process is the
    torque with respect to the effective field, which should tend to zero in
    a local energy minimum. In Cartesian coordinates, this is the result of
    minimising the energy functional with the constraint of a fixed
    magnetisation length. This criteria is controlled via the `mXgradE_tol`
    parameter.

    The number of evaluation steps is counted according to how many times the
    effective field is computed in the creep stage. This number is stored in
    `self.step`.

    The creep algorithm described above is the default step control of
    `minimise`. A second one, selected with `stepControl='BB'`, replaces the
    fixed `η * η_s` step by a Barzilai-Borwein (spectral) step length, which
    is estimated from the curvature seen along the previous step rather than
    tuned by hand. Because the BB quotients are invariant under a constant
    rescaling of the gradient, that variant needs no `eta_scale` argument at
    all, which is otherwise the parameter that has to be matched to the units
    of the effective field. See `_minimise_BB` and reference [4].

    Notes
    -----

    [1] Berkov, D. (1998a). Numerical calculation of the energy barrier
    distribution in disordered many-particle systems: The path integral method.
    JMM, 186(1–2), 199–213.

    [2] Berkov, D.V., Ramstöcck, K. and Hubert, A. (1993), Solving
    Micromagnetic Problems. Towards an Optimal Numerical Method. phys. stat.
    sol. (a), 137: 207-225.

    [3] Ó Conbhuí, P., Williams, W., Fabian, K., Ridley, P., Nagy, L., &
    Muxworthy, A. R. (2018). MERRILL: Micromagnetic earth related robust
    interpreted language laboratory. Geochemistry, Geophysics, Geosystems, 19,
    1080– 1106.

    [4] Birgin, E. G., Martínez, J. M., & Raydan, M. (2000). Nonmonotone
    spectral projected gradient methods on convex sets. SIAM Journal on
    Optimization, 10(4), 1196–1211. The step lengths are those of Barzilai,
    J., & Borwein, J. M. (1988), Two-point step size gradient methods, IMA
    Journal of Numerical Analysis, 8(1), 141–148.
    """

    def __init__(self, mesh, spin, magnetisation, magnetisation_inv, field,
                 pins, interactions, name, data_saver,
                 use_jac=False, integrator=None):
        """Hubert minimiser constructor
        """

        # Inherit from the base minimiser class
        super().__init__(mesh, spin,
                                              magnetisation, magnetisation_inv,
                                              field,
                                              pins,
                                              interactions,
                                              name,
                                              data_saver
                                              )
        # TODO: spin_last and gradE_last should only be temporal, not
        # driver variables

        self.t = 0.0
        # Not using DAMPING here:
        # self._alpha_field = self._alpha * np.ones_like(self.spin)
        self.gradE = np.zeros_like(self.field)
        self.gradE_last = np.zeros_like(self.field)
        # Polak-Ribiere conjugate gradient search direction (MERRILL's `S`)
        self.PR_searchDirection = np.zeros_like(self.field)
        # If we use Cartesian coordinates then what is decreasing in gradient search is m X gradE
        # This comes form the fact that we minimize: E(m) - λ * (m^2 - 1)
        # with respect to m, i.e. d(...)/dm = 0, and that leads to m x Heff = 0
        # If we use speherical coordinates, the constraint is implicit and we have: dE/dtheta = 0, dE/dphi = 0
        self.mXgradE = np.zeros(mesh.n)
        self.totalE = 0.0
        self.totalE_last = 0.0
        # Energy of each site, summed over the interactions. `_minimise_BB`
        # takes the change of energy from these rather than from `totalE`
        self.cellE = np.zeros(mesh.n)
        self.energyScale = 1.

    # def run_step(self):
    #     self.spin_last[:] = self.spin[:]
    #     self.update_effective_field()
    #     self._new_spin[self._material] = (self.spin + self.field)[self._material]

    # WARNING: obj.compute_energy() computes the energy by calling
    #          compute_field() first, evaluated at t=0; if the simulation has a
    #          time-dependent-field, the answer might be wrong
    def compute_effective_field(self, t=0):
        """
        Modified version of the `compute_effective_field` function to obtain
        the total energy of the system. The energy is scaled by the parameter
        `energyScale`.
        """

        self.field[:] = 0
        self.totalE = 0
        self.cellE[:] = 0
        for obj in self.interactions:
            obj.compute_energy()  # Using self.spin, and
            self.field += obj.field[:]
            # `energy` is the energy of each site in every interaction class,
            # so the contributions add directly
            self.cellE += obj.energy / self.energyScale
            self.totalE += obj.total_energy / self.energyScale

    def _normalise_spin(self, spin):
        """ Normalize all spins
        WARNING: Pass spin by reference
        """
        spin.shape = (-1, 3)
        spin[self._pins == 0] /= np.linalg.norm(spin[self._pins == 0], axis=1)[:, None]
        # spin[self._pins > 0] = 0.
        spin.shape = (-1)

    def _result(self, reason):
        """
        Describe the minimisation that has just finished

        Parameters
        ----------
        reason : str
            Which criterion ended it.

        Returns
        -------
        MinimiserResult
        """
        converged = reason in ('mXgradE_tol', 'stopping_dE')
        n = self.mXgradE.shape[0]
        return MinimiserResult(converged=converged,
                               reason=reason,
                               n_evaluations=self.step,
                               total_energy=float(self.totalE),
                               max_torque=float(self.mXgradE.max()),
                               mean_torque=float(np.sum(np.abs(self.mXgradE))
                                                 / n))

    def _minimise_hubert(self,
                         max_steps=2000,
                         save_data_steps=10, save_m_steps=None,
                         save_vtk_steps=None,
                         log_steps=1,
                         maxCreep=5, eta_scale=1.0, stopping_dE=1e-6, dEta=2,
                         etaMin=0.001,
                         # perturbSeed=42, perturbFactor=0.1,
                         nTrail=10, resetMax=20, mXgradE_tol=0.1
                         ):
        """Performs the minimisation

        Parameters
        ----------
        max_steps
            Maximum number of evaluation steps, which increase according to
            the number of calls of the `compute_effective_field` function
        save_data_steps, save_m_steps, save_vtk_steps
            Multiple of steps at which data, spin field and VTK file is saved
        log_steps
            Show log info every X steps
        maxCreep
            Maximum number of steps for the creeping stage, which is the
            minimisation with a fixed η value
        eta_scale
            Scaling factor for the gradient at the spin update step
        stopping_dE
            Mean energy difference with the trailing energy. Remember that the
            energy is scaled by the `self.energyScale` parameter
        dEta
            Factor to increase/decrease η according to the energy change
        etaMin
            Minimum value that η can reach, otherwise it is reset at 1.0 and
            the minimisation starts over again
        nTrail
            Number of energy trailing steps
        resetMax
            Maximum number of resets in case eta reaches the minimum value
            (indicating slow convergence)
        mXgradE_tol
            Tolerance for the mean of the squared norm of the m X energy gradient,
            product `||m X gradE||^2`. The average is calculated from all spin sites
            in material sites.
        """

        # rstate = np.random.RandomState(perturbSeed)
        self.spin_last = np.zeros_like(self.spin)
        self.step = 0
        nStart = 0
        exitFlag = False
        reason = 'max_steps'
        totalRestart = True
        resetCount = 0
        creepCount = 0
        self.spin_last[:] = self.spin
        self.trailE = np.zeros(nTrail)
        trailPool = cycle(range(nTrail))  # cycle through 0,1,...,(nTrail-1),0,1,...
        eta = 1.0
        # ||gradE||^2 at the last accepted point, used as the denominator in
        # the Polak-Ribiere beta below (MERRILL's `GO2`)
        GSQUARE = 0.0
        # We might want to change this in the future to save memory:
        # pinsField = np.repeat(self._pins, 3).astype(bool)
        _material = self._material_mask()

        while not exitFlag:

            if totalRestart:
                if self.step > 0:
                    log.debug('Restarting')
                self.spin[:] = self.spin_last
                # Compute from self.spin. Do not update the step at this stage:
                self.compute_effective_field()
                # self.step += 1
                self.gradE_last[:] = -self.field  # Scale field??
                self.gradE_last[~_material] = 0.0
                self.gradE[:] = self.gradE_last
                self.totalE_last = self.totalE
                self.trailE[nStart] = self.totalE
                nStart = next(trailPool)
                eta = 1.0
                # Reset the conjugate direction to plain steepest descent
                # (equivalent to beta = 0 in MERRILL, which happens here
                # naturally since the gradient hasn't changed from the last
                # accepted point)
                self.PR_searchDirection[:] = self.gradE_last
                GSQUARE = np.sum(self.gradE_last[_material] ** 2)
                totalRestart = False

            creepCount = 0

            # Creep stage: minimise with a fixed eta
            while creepCount < maxCreep:
                # Update spin. Avoid pinned or zero-Ms sites
                self.spin[_material] = self.spin_last[_material] - eta * eta_scale * self.PR_searchDirection[_material]
                # self.spin[~pinsField] = self.spin_last[~pinsField] - eta * eta_scale * self.PR_searchDirection[~pinsField]
                # self.spin[:] = self.spin_last - eta * eta_scale * self.PR_searchDirection

                # Normalize spin
                self._normalise_spin(self.spin)
                # print(self.spin.reshape(-1, 3)[self._pins == 0][:10])
                # if creepCount > 4:
                #     exitFlag = True
                #     break

                self.compute_effective_field()  # Compute Heff and E using self.spin
                self.step += 1

                # TODO: No-material sites should have gradE=0.0 but we must
                # be sure not taking them into account, in gradE at least
                self.gradE[:] = -self.field  # Scale field??
                self.gradE[~_material] = 0.0
                # Save the energy and move trail index to next site
                self.trailE[nStart] = self.totalE
                nStart = next(trailPool)
                mXgrad = np.cross(self.spin.reshape(-1, 3), self.gradE.reshape(-1, 3), axis=1)
                # np.einsum('ij,ij->i', mXgrad, mXgrad, out=self.mXgradE2)
                self.mXgradE[:] = np.linalg.norm(mXgrad, axis=1)
                self.mXgradE[~_material[::3]] = 0.0
                # self.gradE2[~_material[::3]] = 0.0
                # Compute E difference of current E (totalE) with the trailing E
                deltaE = abs(self.trailE[nStart] - self.totalE) / nTrail

                # Statistics and saving:
                # `mXgradE.max()` is a pass over the sites, so the message is
                # only built when it will be seen. `log_steps` throttles it
                # further, for a run logged at DEBUG
                if (self.step % log_steps == 0
                        and log.isEnabledFor(logging.DEBUG)):
                    log.debug(f'Step = {self.step:>4} Creep n = {creepCount:>3}  reset = {resetCount:>3}  eta = {eta:>5.4e}  E_new = {self.totalE:.4e}  ΔE = {deltaE:.4e}  max(|mX∇E|) = {self.mXgradE.max():.4e}')
                # Note that step == 0 is never saved
                if self.step % save_data_steps == 0:
                    self.data_saver.save()
                if (save_vtk_steps is not None) and (self.step % save_vtk_steps == 0):
                    self.save_vtk()
                if (save_m_steps is not None) and (self.step % save_m_steps == 0):
                    self.save_m()

                # print(self.trailE)
                # print('Tot E last', self.totalE_last)

                # with np.printoptions(precision=2):
                #     print('Creep: ', self.trailE)
                #     print(f': eta = {eta}  maxgradE = {self.gradE.max()}')

                # Only trust a near-flat trailing energy as convergence if
                # this step was actually accepted (the energy did not
                # increase). trailE is filled on every evaluation, including
                # rejected trial steps, so while eta is being shrunk the band
                # can bounce back close to its nTrail-old energy and give a
                # spuriously small deltaE while the true residual (mX∇E) is
                # still large.
                if self.totalE <= self.totalE_last and deltaE < stopping_dE:
                    log.info(f'Delta E = {deltaE} negligible. Stopping calculation.')
                    reason = 'stopping_dE'
                    exitFlag = True
                    break  # creep loop

                if self.step > max_steps:
                    log.warning(f'N of evaluations = {self.step} reached maximum value. Stopping calculation.')
                    reason = 'max_steps'
                    exitFlag = True
                    break  # creep loop

                if self.totalE > self.totalE_last:  # If E increases, decrease eta so minimise slower
                    # print('Decreasing eta')
                    creepCount = 0
                    eta = eta / (dEta * dEta)

                    if eta < etaMin:
                        log.debug(f'Parameter eta smaller than minimum. Restarting minimisation. resetCount = {resetCount}')
                        resetCount += 1
                        # perturbSpins()  # in case of using Sph coordinates
                        if resetCount > resetMax:
                            reason = 'resets'
                            exitFlag = True
                            log.warning(f'N of resets {resetCount} reached maximum value. Stopping calculation.')
                            break  # creep loop
                        totalRestart = True
                        break  # creep loop
                else:  # if E decreases move to next creep step
                    # print('Total E < total E last')
                    creepCount += 1

                    # Polak-Ribiere conjugate gradient direction update.
                    # Must use the OLD gradE_last (gradient at the point we
                    # are leaving) together with the newly accepted gradE
                    # before gradE_last gets overwritten below.
                    GO2 = GSQUARE
                    GSQUARE = np.sum(self.gradE[_material] ** 2)
                    SPG = np.sum(self.gradE_last[_material] * self.gradE[_material])
                    beta = max(0.0, (GSQUARE - SPG) / GO2) if GO2 != 0.0 else 0.0
                    self.PR_searchDirection[_material] = self.gradE[_material] + beta * self.PR_searchDirection[_material]
                    self.PR_searchDirection[~_material] = 0.0

                    # Update Energy, spin and gradE
                    self.spin_last[:] = self.spin[:]
                    self.gradE_last[:] = self.gradE[:]
                    self.totalE_last = self.totalE

                    avGradE = np.sum(np.abs(self.mXgradE)) / self.mXgradE.shape[0]
                    if avGradE < mXgradE_tol:
                        log.info(f'Average torque length |mX∇E|/N = {avGradE} negligible. Stopping calculation.')
                        reason = 'mXgradE_tol'
                        exitFlag = True

            # Stop while creepCount

            # If creeping went OK (i.e. we did not just break out to restart
            # the minimisation because eta got too small), increase
            # minimisation speed by increasing eta for a next creeping stage,
            # and clear the reset-failure counter. Note: resetCount must NOT
            # be cleared on a restart, otherwise it can never accumulate past
            # 1 and the resetMax limit is never enforced (matches MERRILL,
            # where ResetCnt is only zeroed when the creep loop finishes
            # normally, not on the `goto 10` restart path).
            if not totalRestart:
                eta *= dEta
                resetCount = 0

        return self._result(reason)

    # -------------------------------------------------------------------------
    # Shared helpers

    def _material_mask(self):
        """Boolean mask, shaped like `spin`, of sites that may be updated

        These are the sites which are not pinned and which carry a non-zero
        magnetisation.
        """
        mask = ~(np.repeat(self._pins, 3).astype(bool))
        mask.reshape(-1, 3)[self._magnetisation <= 0.0] = False
        return mask

    def _project_gradient(self, out=None):
        """Tangential (Riemannian) gradient of the energy at `self.spin`

        The Cartesian gradient `δE/δm = -H_eff` has a radial component
        `(m · δE/δm) m` which is the Lagrange multiplier enforcing `|m| = 1`.
        That component is annihilated by the re-normalisation of the spins, so
        it carries no information about the descent direction, but it *does*
        pollute any quantity built from differences of gradients (such as the
        Barzilai-Borwein secant pair). Removing it gives::

            g = δE/δm - (m · δE/δm) m = -m × (m × δE/δm)

        whose norm at every site is `||m × δE/δm||`, i.e. the same residual
        already used as the stopping criterion of this class. The norms are
        stored in `self.mXgradE`.
        """
        g = self.gradE if out is None else out
        g[:] = -self.field
        g3 = g.reshape(-1, 3)
        m3 = self.spin.reshape(-1, 3)
        g3 -= np.einsum('ij,ij->i', m3, g3)[:, None] * m3
        g[~self._material] = 0.0
        self.mXgradE[:] = np.linalg.norm(g3, axis=1)
        return g

    # -------------------------------------------------------------------------

    def _minimise_BB(self,
                     max_steps=2000,
                     save_data_steps=10, save_m_steps=None,
                     save_vtk_steps=None,
                     log_steps=1,
                     stopping_dE=1e-6, dEta=2,
                     nTrail=10, resetMax=20, mXgradE_tol=0.1,
                     maxDeltaM=0.1, gamma=1e-4, BBstep='alternate',
                     maxBacktrack=15
                     ):
        """Spectral (Barzilai-Borwein) gradient descent on the unit sphere

        Same descent problem as `_minimise_hubert`, but the step length is not
        a hand-tuned constant: it is estimated at every iteration from the
        secant pair of the last accepted step,

        ::

            s = m_k - m_{k-1}     y = g_k - g_{k-1}

        with `g` the *tangential* gradient (see `_project_gradient`), through
        the two Barzilai-Borwein quotients::

            η_BB1 = (s · s) / (s · y)      η_BB2 = (s · y) / (y · y)

        Both are inverse Rayleigh quotients of the Hessian along the last
        step, so η carries the units of `[m] / [δE/δm]` and the `eta_scale`
        argument of `_minimise_hubert` -- which exists only to convert the
        field units into a spin displacement -- is no longer needed.

        BB steps are not monotone by construction, so a monotone accept/reject
        test would throw away exactly the behaviour that makes them fast.
        Acceptance therefore uses the non-monotone Grippo-Lampariello-Lucidi
        condition over the trailing energies already kept by this class::

            E(m_new) <= max(trailE) - γ λ ||g||^2

        backtracking with `λ / dEta^2` when it fails. This is the SPG scheme
        of Birgin, Martínez & Raydan, with the sphere as the constraint set
        and the re-normalisation of the spins as the projection onto it.

        Parameters
        ----------
        max_steps
            Maximum number of evaluation steps, counted as calls to
            `compute_effective_field` (backtracking trials included)
        save_data_steps, save_m_steps, save_vtk_steps
            Multiple of steps at which data, spin field and VTK file is saved
        log_steps
            Show log info every X steps
        stopping_dE
            Mean energy difference with the trailing energy. Remember that the
            energy is scaled by the `self.energyScale` parameter
        dEta
            The backtracking factor is `dEta ** 2`, matching the rate at which
            `_minimise_hubert` shrinks η on a rejected step
        nTrail
            Number of energy trailing steps. This is also the width of the
            non-monotone acceptance window: `nTrail = 1` recovers a monotone
            (Armijo) line search
        resetMax
            Maximum number of restarts, where a restart is a step that could
            not be accepted within `maxBacktrack` backtracks
        mXgradE_tol
            Tolerance for the mean of the norm of the `m X gradE` torque,
            averaged over all sites
        maxDeltaM
            Trust region: the trial step is clipped so that no spin moves more
            than this distance (with `|m| = 1`) before re-normalisation. It
            keeps the projection onto the sphere accurate and provides the
            first step length, `maxDeltaM / max||g||`, when no secant pair is
            available yet
        gamma
            Sufficient-decrease parameter of the acceptance test
        BBstep
            `'BB1'`, `'BB2'` or `'alternate'` (default). Alternating the two
            quotients is more robust than either alone, since BB1 overshoots
            and BB2 undershoots on opposite sides of the spectrum
        maxBacktrack
            Maximum number of backtracks before declaring a restart
        """

        self.spin_last = np.zeros_like(self.spin)
        self.step = 0
        exitFlag = False
        reason = 'max_steps'
        resetCount = 0
        self._material = self._material_mask()
        _material = self._material

        # Energy and tangential gradient at the starting point
        self.compute_effective_field()
        self._project_gradient(out=self.gradE)
        self.totalE_last = self.totalE

        # Every energy here is measured from the starting point and summed
        # from the accepted steps. Near a minimum the decrease to be detected
        # falls below the spacing of the doubles around the total energy, so a
        # difference of totals cannot see it; a sum of the per site
        # differences, each small against its own site energy, can. This is
        # what OOMMF's Oxs_CGEvolve does.
        #
        # There is a second route worth recording, since it needs no per site
        # energy at all. Writing the energy of the interactions that are
        # quadratic in the magnetisation as
        # E = -1/2 sum_i w_i m_i . H_i, with w_i = mu_0 V Ms_i in the
        # micromagnetic case and mu_s_i in the atomistic one, and using that
        # the field is linear in m with a symmetric operator, so that
        # m_n . H_r = m_r . H_n, the two mixed terms of
        # (m_n - m_r).(H_n + H_r) cancel and
        #
        #     dE = -1/2 sum_i w_i (m_n - m_r)_i . (H_n + H_r)_i
        #
        # a trapezoid that the reciprocity makes exact rather than first
        # order. It carries the small displacement in every term, so it does
        # not cancel either, and a constant Zeeman field is covered by it
        # because its two fields are equal. It was measured to reach the same
        # convergence as the sum above, but it is only exact while the energy
        # is quadratic, a cubic anisotropy being the case to watch; it needs
        # the field at both ends, so 3n stored rather than n; and it gives the
        # change up to a constant that has to be calibrated. The sum below is
        # simpler on all three counts
        cellE_ref = self.cellE.copy()
        Erel = 0.0

        # Fill the whole trailing window with E0, so that max(trailE) is a
        # valid (and initially tight) reference from the very first step
        self.trailE = np.zeros(nTrail)
        trailPool = cycle(range(nTrail))
        nStart = next(trailPool)

        self.spin_last[:] = self.spin
        self.gradE_last[:] = self.gradE

        maxGrad = self.mXgradE.max()
        if maxGrad == 0.0:
            log.warning('Gradient is zero everywhere. Nothing to minimise.')
            return self._result('zero_gradient')
        # Scale-free first step: move the largest torque by maxDeltaM
        eta0 = maxDeltaM / maxGrad
        eta = eta0
        # `gradE` is built from the effective field, so it is the true energy
        # gradient only up to a constant with the units of the interactions
        # (mu_0 Ms V in the micromagnetic case, mu_s in the atomistic one),
        # further divided by `energyScale`. The BB quotients are invariant
        # under that rescaling -- it cancels between η and the gradient in the
        # spin update, which is precisely why no `eta_scale` is needed here --
        # but the sufficient-decrease term below is not, since it compares a
        # gradient norm against an energy. Rather than hard-coding a unit
        # conversion that differs between the micromagnetic and the atomistic
        # classes, the constant is calibrated from the decrease actually
        # observed along the accepted steps, starting from a plain
        # non-monotone decrease test.
        # Only safeguard against a step length collapsing to nothing; the
        # upper end is handled by the maxDeltaM trust region below
        etaMin = 1e-10 * eta0
        BBcount = 0
        gradScale = 0.0

        while not exitFlag:

            gradNorm2 = np.sum(self.gradE[_material] ** 2)
            if gradNorm2 == 0.0:
                log.warning('Gradient is zero everywhere. Stopping calculation.')
                reason = 'zero_gradient'
                break

            # Trust region: never displace a spin by more than maxDeltaM
            maxGrad = self.mXgradE.max()
            lamb = min(eta, maxDeltaM / maxGrad)
            # Non-monotone reference energy of the trailing window
            Eref = self.trailE.max()

            nBacktrack = 0
            accepted = False
            while not accepted:
                self.spin[_material] = (self.spin_last[_material]
                                        - lamb * self.gradE_last[_material])
                self._normalise_spin(self.spin)

                self.compute_effective_field()
                self.step += 1

                dE = np.sum(self.cellE - cellE_ref)
                Etrial = Erel + dE

                # Grippo-Lampariello-Lucidi non-monotone acceptance
                accepted = (Etrial
                            <= Eref - gamma * gradScale * lamb * gradNorm2)

                if self.step > max_steps:
                    log.warning(f'N of evaluations = {self.step} reached '
                                'maximum value. Stopping calculation.')
                    reason = 'max_steps'
                    exitFlag = True
                    break  # backtracking loop

                if accepted:
                    break  # backtracking loop

                nBacktrack += 1
                lamb = lamb / (dEta * dEta)
                if nBacktrack > maxBacktrack or lamb < etaMin:
                    # No decrease along -g: drop the secant information and
                    # start over from the last accepted point
                    resetCount += 1
                    log.debug('Could not decrease the energy along the '
                              'gradient. Restarting minimisation. '
                              f'resetCount = {resetCount}')
                    if resetCount > resetMax:
                        log.warning(f'N of resets {resetCount} reached '
                                    'maximum value. Stopping calculation.')
                        reason = 'resets'
                        exitFlag = True
                    break  # backtracking loop

            if exitFlag:
                break  # main loop

            if not accepted:
                # Restart: recompute the state at the last accepted point and
                # take a plain steepest descent step next
                self.spin[:] = self.spin_last
                self.compute_effective_field()
                self._project_gradient(out=self.gradE)
                self.gradE_last[:] = self.gradE
                self.totalE_last = self.totalE
                # Back at the last accepted point, which `Erel` already
                # describes; only the site energies have to be re-referenced
                cellE_ref[:] = self.cellE
                eta = eta0
                BBcount = 0
                continue  # main loop

            # ------ Accepted step ------------------------------------------
            self._project_gradient(out=self.gradE)

            # Calibrate the units of the gradient against the energy from the
            # first-order model of the decrease just achieved,
            # E_last - E ~ scale * λ ||g||^2
            gradScale = max(gradScale, -dE / (lamb * gradNorm2))

            Erel = Etrial
            cellE_ref[:] = self.cellE

            self.trailE[nStart] = Erel
            nStart = next(trailPool)
            deltaE = abs(self.trailE[nStart] - Erel) / nTrail

            # Secant pair of the step just taken. Note that `s` is the
            # difference of the *projected* spins, not -lamb * g, which is
            # what the SPG framework requires
            s = self.spin[_material] - self.spin_last[_material]
            y = self.gradE[_material] - self.gradE_last[_material]
            sy = np.dot(s, y)

            if sy > 0.0:
                if BBstep == 'BB1':
                    useBB1 = True
                elif BBstep == 'BB2':
                    useBB1 = False
                else:
                    useBB1 = (BBcount % 2 == 0)
                if useBB1:
                    eta = np.dot(s, s) / sy
                else:
                    yy = np.dot(y, y)
                    eta = sy / yy if yy > 0.0 else eta0
            else:
                # Non-positive curvature along the step: the BB quotients are
                # meaningless there, so fall back to the trust region step
                eta = eta0
            if not np.isfinite(eta) or eta <= 0.0:
                eta = eta0
            eta = max(eta, etaMin)
            BBcount += 1

            self.spin_last[:] = self.spin
            self.gradE_last[:] = self.gradE
            self.totalE_last = self.totalE

            # See the note in `_minimise_hubert`: the message costs a pass
            # over the sites, so it is only built when it will be seen
            if self.step % log_steps == 0 and log.isEnabledFor(logging.DEBUG):
                log.debug(f'Step = {self.step:>4}  '
                          f'backtracks = {nBacktrack:>2}  '
                          f'reset = {resetCount:>3}  eta = {eta:>5.4e}  '
                          f'E_new = {self.totalE:.4e}  ΔE = {deltaE:.4e}  '
                          f'max(|mX∇E|) = {self.mXgradE.max():.4e}')
            if self.step % save_data_steps == 0:
                self.data_saver.save()
            if (save_vtk_steps is not None) and (self.step % save_vtk_steps == 0):
                self.save_vtk()
            if (save_m_steps is not None) and (self.step % save_m_steps == 0):
                self.save_m()

            if deltaE < stopping_dE:
                log.info(f'Delta E = {deltaE} negligible. '
                         'Stopping calculation.')
                reason = 'stopping_dE'
                exitFlag = True

            avGradE = np.sum(np.abs(self.mXgradE)) / self.mXgradE.shape[0]
            if avGradE < mXgradE_tol:
                log.info(f'Average torque length |mX∇E|/N = {avGradE} '
                         'negligible. Stopping calculation.')
                reason = 'mXgradE_tol'
                exitFlag = True

        return self._result(reason)

    def minimise(self, max_steps=2000,
                 save_data_steps=10, save_m_steps=None, save_vtk_steps=None,
                 log_steps=1,
                 stopping_dE=1e-6, dEta=2, nTrail=10, resetMax=20,
                 mXgradE_tol=0.1,
                 stepControl='hubert',
                 # `hubert` step control:
                 maxCreep=5, eta_scale=1.0, etaMin=0.001,
                 # `BB` step control:
                 maxDeltaM=0.1, gamma=1e-4, BBstep='alternate',
                 maxBacktrack=15
                 ):
        """Performs the minimisation

        Two step-length strategies are available through `stepControl`:

        `'hubert'` (default)

            The creep algorithm of [1, 2, 3]: a fixed `η * eta_scale` step
            along a Polak-Ribiere conjugate direction, with η grown or shrunk
            according to the sign of the energy change. See
            `_minimise_hubert` for its parameters. Note that `eta_scale` has
            to be chosen by hand to match the units of the effective field.

        `'BB'`

            Barzilai-Borwein spectral gradient descent with a non-monotone
            line search, which infers the step length from the curvature seen
            along the previous step and so needs no `eta_scale`. See
            `_minimise_BB` for its parameters.

        Parameters shared by both are documented in `_minimise_hubert`.

        Returns
        -------
        MinimiserResult
            Whether a convergence criterion was met, which one, the number of
            effective field evaluations spent, and the energy and torque of
            the final configuration. Progress is reported through the
            `fidimag` logger, at DEBUG for the per step lines and INFO for
            the reason it stopped, so a run that converges is silent unless
            that logger is turned up::

                logging.getLogger('fidimag').setLevel(logging.INFO)

            Failures to converge are logged at WARNING and so are shown.
        """

        if stepControl == 'hubert':
            return self._minimise_hubert(
                max_steps=max_steps,
                save_data_steps=save_data_steps, save_m_steps=save_m_steps,
                save_vtk_steps=save_vtk_steps, log_steps=log_steps,
                maxCreep=maxCreep, eta_scale=eta_scale,
                stopping_dE=stopping_dE, dEta=dEta, etaMin=etaMin,
                nTrail=nTrail, resetMax=resetMax, mXgradE_tol=mXgradE_tol)
        elif stepControl == 'BB':
            return self._minimise_BB(
                max_steps=max_steps,
                save_data_steps=save_data_steps, save_m_steps=save_m_steps,
                save_vtk_steps=save_vtk_steps, log_steps=log_steps,
                stopping_dE=stopping_dE, dEta=dEta,
                nTrail=nTrail, resetMax=resetMax, mXgradE_tol=mXgradE_tol,
                maxDeltaM=maxDeltaM, gamma=gamma, BBstep=BBstep,
                maxBacktrack=maxBacktrack)
        else:
            raise ValueError(f'Unknown stepControl `{stepControl}`. Use '
                             '`hubert` or `BB`.')
