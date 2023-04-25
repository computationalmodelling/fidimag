import numpy as np
# import fidimag.extensions.clib as clib
# import fidimag.common.constant as const
from .minimiser_base import MinimiserBase
from itertools import cycle


class HubertMinimiser(MinimiserBase):
    """ Gradient descent energy minimisation algorithm.

    A minimisation, energy gradient descent, algorithm, that dynamically
    updates a scaling factor α for the gradient of the energy::

        m_new = m_old + α * α_s * δE/δm|_old

    where α_s is a fixed scaling factor. The E functional gradient is computed
    from the effective field `δE/δm = - H_eff`. 

    This algorithm is based on the works [1, 2] and is implemented in [3] as
    *Hubert minimiser*. The modifications include using the E gradient instead
    of the magnetisation torque, and using different criteria for stopping the
    algorithm.

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

    def __init__(self, mesh, spin,
                 magnetisation, magnetisation_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac=False,
                 integrator=None
                 ):

        # Inherit from the base minimiser class
        super(HubertMinimiser, self).__init__(mesh, spin,
                                              magnetisation, magnetisation_inv,
                                              field,
                                              pins,
                                              interactions,
                                              name,
                                              data_saver
                                              )

        self.t = 1e-4
        # self._alpha_field = self._alpha * np.ones_like(self.spin)
        self._new_spin = np.zeros_like(spin)

        self.totalRestart = False
        self.exit = False
        self.gradE = np.zeros_like(self.field)
        self.gradE_last = np.zeros_like(self.field)
        self.gradE2 = np.zeros(mesh.n)
        self.totalE = 0.0
        self.totalE_last = 0.0
        self.energyScale = 1.

    # def run_step(self):
    #     spin_last[:] = self.spin[:]
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
                 max_steps=1000,
                 # save_data_steps=10, save_m_steps=None, save_vtk_steps=None, log_steps=1000,
                 maxCreep=5, alpha_scale=1.0, stopping_dE=1e-6, dAlpha=2,
                 alphaMin=0.001,
                 # perturbSeed=42, perturbFactor=0.1,
                 nTrail=10, resetMax=20, gradEtol=1e-14
                 ):
        """

        """

        # rstate = np.random.RandomState(perturbSeed)
        spin_last = np.zeros_like(self.spin)
        self.step = 0
        nStart = 0
        exitFlag = False
        totalRestart = True
        resetCount = 0
        creepCount = 0
        spin_last[:] = self.spin
        self.trailE = np.zeros(nTrail)
        trailPool = cycle(range(nTrail))  # cycle through 0,1,...,(nTrail-1),0,1,...
        alpha = 1.0
        # We might want to change this in the future to save memory:
        # pinsField = np.repeat(self._pins, 3).astype(bool)
        # Only update site with magnetisation > 0 which are not pinned
        _material = ~(np.repeat(self._pins, 3).astype(bool))
        _material.reshape(-1, 3)[self._magnetisation > 0.0] = True

        while not exitFlag:

            if totalRestart:
                print('Restarting')
                self.spin[:] = spin_last
                # Compute from self.spin:
                self.compute_effective_field()  # Compute Heff and E using self.spin
                self.gradE_last[:] = -self.field  # Scale field??
                self.gradE[:] = self.gradE_last
                self.totalE_last = self.totalE
                self.trailE[nStart] = self.totalE
                nStart = next(trailPool)
                alpha = 1.0
                totalRestart = False

            creepCount = 0

            while creepCount < maxCreep:
                # Update spin. TODO: avoid pinned sites
                self.spin[_material] = spin_last[_material] - alpha * alpha_scale * self.gradE_last[_material]
                # self.spin[~pinsField] = spin_last[~pinsField] - alpha * alpha_scale * self.gradE_last[~pinsField]
                # self.spin[:] = spin_last - alpha * alpha_scale * self.gradE_last

                # Normalize spin
                self._normalise_spin(self.spin)
                # print(self.spin.reshape(-1, 3)[self._pins == 0][:10])
                # if creepCount > 4:
                #     exitFlag = True
                #     break

                self.compute_effective_field()  # Compute Heff and E using self.spin
                self.gradE[:] = -self.field  # Scale field??
                self.trailE[nStart] = self.totalE
                nStart = next(trailPool)
                # TODO: No-material sites should have gradE=0 but we should make sure not
                # taking them into account:
                gE_reshape = self.gradE.reshape(-1, 3)
                np.einsum('ij,ij->i', gE_reshape, gE_reshape, out=self.gradE2)
                deltaE = abs(self.trailE[nStart] - self.totalE) / nTrail
                
                print(f'Creep n = {creepCount:>3}  reset = {resetCount:>3}  alpha = {alpha:>8}  E_new = {self.totalE:.4e}  ΔE = {deltaE:.4e}  max(∇E)^2 = {self.gradE2.max():.4e}')
                # print(self.trailE)
                # print('Tot E last', self.totalE_last)

                # with np.printoptions(precision=2):
                #     print('Creep: ', self.trailE)
                #     print(f': alpha = {alpha}  maxgradE = {self.gradE.max()}')

                if deltaE < stopping_dE:
                    print(f'Delta E = {deltaE} negligible. Stopping calculation.')
                    exitFlag = True
                    break  # creep loop

                # if nEval > max_steps:
                #     ...

                if self.totalE > self.totalE_last:
                    print('Decreasing alpha')
                    self.creepCount = 0
                    alpha = alpha / (dAlpha * dAlpha)

                    if alpha < alphaMin:
                        print(f'Parameter alpha smaller than minimum. Restarting minimisation. resetCount = {resetCount}')
                        resetCount += 1
                        # perturbSpins()
                        if resetCount > resetMax:
                            exitFlag = True
                            print(f'N of resets {resetCount} reached maximum value. Stopping calculation.')
                            break  # creep loop
                        totalRestart = True
                        break  # creep loop
                else:
                    # print('Total E < total E last')
                    creepCount += 1
                    # Update Energy, spin and gradE
                    spin_last[:] = self.spin[:]
                    self.gradE_last[:] = self.gradE[:]
                    self.totalE_last = self.totalE

                    avGradE = np.sum(self.gradE2) / self.gradE2.shape[0]
                    if avGradE < gradEtol:
                        print(f'Average gradient G^2/N = {avGradE} negligible. Stopping calculation.')
                        exitFlag = True

            # Stop while creepCount

            alpha *= dAlpha
            resetCount = 0


            # if self.step % log_steps == 0:
            #     print("#max_tau={:<8.3g} max_dm={:<10.3g} counter={}".format(
            #         np.max(np.abs(self.tau)),
            #         max_dm, self.step))
