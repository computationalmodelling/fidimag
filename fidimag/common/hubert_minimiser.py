import numpy as np
# import fidimag.extensions.clib as clib
# import fidimag.common.constant as const
from .minimiser_base import MinimiserBase
from itertools import cycle


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
    """

    def __init__(self, mesh, spin, magnetisation, magnetisation_inv, field,
                 pins, interactions, name, data_saver,
                 use_jac=False, integrator=None):
        """Hubert minimiser constructor
        """

        # Inherit from the base minimiser class
        super(HubertMinimiser, self).__init__(mesh, spin,
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
        # If we use Cartesian coordinates then what is decreasing in gradient search is m X gradE
        # This comes form the fact that we minimize: E(m) - λ * (m^2 - 1)
        # with respect to m, i.e. d(...)/dm = 0, and that leads to m x Heff = 0
        # If we use speherical coordinates, the constraint is implicit and we have: dE/dtheta = 0, dE/dphi = 0
        self.mXgradE = np.zeros(mesh.n)
        self.totalE = 0.0
        self.totalE_last = 0.0
        self.energyScale = 1.

    # def run_step(self):
    #     self.spin_last[:] = self.spin[:]
    #     self.update_effective_field()
    #     self._new_spin[self._material] = (self.spin + self.field)[self._material]

    def compute_effective_field(self, t=0):
        """
        Modified version of the `compute_effective_field` function to obtain
        the total energy of the system. The energy is scaled by the parameter
        `energyScale`.
        """

        self.field[:] = 0
        self.totalE = 0
        for obj in self.interactions:
            obj.compute_energy()  # Using self.spin, and
            self.field += obj.field[:]
            self.totalE += obj.total_energy / self.energyScale

    def _normalise_spin(self, spin):
        """ Normalize all spins
        WARNING: Pass spin by reference
        """
        spin.shape = (-1, 3)
        spin[self._pins == 0] /= np.linalg.norm(spin[self._pins == 0], axis=1)[:, None]
        # spin[self._pins > 0] = 0.
        spin.shape = (-1)

    def minimise(self,
                 max_steps=2000,
                 save_data_steps=10, save_m_steps=None, save_vtk_steps=None,
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
        totalRestart = True
        resetCount = 0
        creepCount = 0
        self.spin_last[:] = self.spin
        self.trailE = np.zeros(nTrail)
        trailPool = cycle(range(nTrail))  # cycle through 0,1,...,(nTrail-1),0,1,...
        eta = 1.0
        # We might want to change this in the future to save memory:
        # pinsField = np.repeat(self._pins, 3).astype(bool)
        # Only update site with magnetisation > 0 which are not pinned
        _material = ~(np.repeat(self._pins, 3).astype(bool))
        _material.reshape(-1, 3)[self._magnetisation > 0.0] = True

        while not exitFlag:

            if totalRestart:
                if self.step > 0:
                    print('Restarting')
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
                totalRestart = False

            creepCount = 0

            # Creep stage: minimise with a fixed eta
            while creepCount < maxCreep:
                # Update spin. Avoid pinned or zero-Ms sites
                self.spin[_material] = self.spin_last[_material] - eta * eta_scale * self.gradE_last[_material]
                # self.spin[~pinsField] = self.spin_last[~pinsField] - eta * eta_scale * self.gradE_last[~pinsField]
                # self.spin[:] = self.spin_last - eta * eta_scale * self.gradE_last

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
                if self.step % log_steps == 0:
                    print(f'Step = {self.step:>4} Creep n = {creepCount:>3}  reset = {resetCount:>3}  eta = {eta:>5.4e}  E_new = {self.totalE:.4e}  ΔE = {deltaE:.4e}  max(|mX∇E|) = {self.mXgradE.max():.4e}')
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

                if deltaE < stopping_dE:  # if mean E diff is too small
                    print(f'Delta E = {deltaE} negligible. Stopping calculation.')
                    exitFlag = True
                    break  # creep loop

                if self.step > max_steps:
                    print(f'N of evaluations = {self.step} reached maximum value. Stopping calculation.')
                    exitFlag = True
                    break  # creep loop

                if self.totalE > self.totalE_last:  # If E increases, decrease eta so minimise slower
                    # print('Decreasing eta')
                    self.creepCount = 0
                    eta = eta / (dEta * dEta)

                    if eta < etaMin:
                        print(f'Parameter eta smaller than minimum. Restarting minimisation. resetCount = {resetCount}')
                        resetCount += 1
                        # perturbSpins()  # in case of using Sph coordinates
                        if resetCount > resetMax:
                            exitFlag = True
                            print(f'N of resets {resetCount} reached maximum value. Stopping calculation.')
                            break  # creep loop
                        totalRestart = True
                        break  # creep loop
                else:  # if E decreases move to next creep step
                    # print('Total E < total E last')
                    creepCount += 1
                    # Update Energy, spin and gradE
                    self.spin_last[:] = self.spin[:]
                    self.gradE_last[:] = self.gradE[:]
                    self.totalE_last = self.totalE

                    avGradE = np.sum(np.abs(self.mXgradE)) / self.mXgradE.shape[0]
                    if avGradE < mXgradE_tol:
                        print(f'Average torque length |mX∇E|/N = {avGradE} negligible. Stopping calculation.')
                        exitFlag = True

            # Stop while creepCount

            # If creeping went OK, increase minimisation speed by increasing eta
            # for a next creeping stage
            eta *= dEta
            resetCount = 0
