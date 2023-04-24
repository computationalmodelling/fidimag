import numpy as np
import fidimag.extensions.clib as clib
# import fidimag.common.constant as const
from .minimiser_base import MinimiserBase
from itertools import cycle


class HubertMinimiser(MinimiserBase):
    """
    A minimisation algorithm that dynamically updates a scaling factor for
    the

        CHECK:

        h = (m_i + H) / || (m_i + H) ||
        m_i+1 = m_i + alpha * (h_i - m_i)
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
        # self._alpha = 0.1
        # self._alpha_field = self._alpha * np.ones_like(self.spin)
        self.spin_last = np.zeros_like(spin)
        self._new_spin = np.zeros_like(spin)

        self.totalRestart = False
        self.exit = False
        self.gradE = np.zeros_like(self.field)
        self.gradE_last = np.zeros_like(self.field)
        self.gradE2 = np.zeros(mesh.n)
        self.totalE = 0.0
        self.totalE_last = 0.0

        # Only update site with magnetisation > 0 which are not pinned
        self._material = np.logical_and(np.repeat(self._magnetisation, 3) > 0.0,
                                        np.repeat(1 - self._pins, 3).astype(np.bool))

    # def run_step(self):
    #     self.spin_last[:] = self.spin[:]
    #     self.update_effective_field()
    #     self._new_spin[self._material] = (self.spin + self.field)[self._material]

    def compute_effective_field(self, t=0):
        """
        Modified version of the `compute_effective_field` function to obtain
        the total energy of the system
        """

        self.field[:] = 0
        self.totalE = 0
        for obj in self.interactions:
            obj.compute_energy()  # Using self.spin, and
            self.field += obj.field[:]
            self.totalE += obj.total_energy

    def normalise_spin(self, spin):
        """WARNING: Pass spin by reference
        """
        spin.shape = (-1, 3)
        spin[self._pins == 0] /= np.linalg.norm(spin[self._pins == 0], axis=1)[:, None]
        spin.shape = (-1)

    def minimise(self,
                 max_steps=1000,
                 # save_data_steps=10, save_m_steps=None, save_vtk_steps=None, log_steps=1000,
                 maxCreep=5, alpha_scale=0.01, stopping_dE=1e-6, dAlpha=2,
                 alphaMin=0.001,
                 # perturbSeed=42, perturbFactor=0.1,
                 nTrail=10, resetMax=20, gradEtol=1e-14
                 ):
        """

        """

        self.step = 0
        nStart = 0
        exitFlag = False
        totalRestart = True
        # rstate = np.random.RandomState(perturbSeed)
        resetCount = 0
        creepCount = 0
        self.spin_last[:] = self.spin
        self.trailE = np.zeros(nTrail)
        trailPool = cycle(range(nTrail))  # cycle through 0,1,...,(nTrail-1),0,1,...
        alpha = 1.0

        while not exitFlag:

            if totalRestart:
                print('Restarting')
                self.spin[:] = self.spin_last
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
                print('Tot E', self.totalE)
                print('Tot E last', self.totalE_last)
                # Update spin. TODO: avoid pinned sites
                # self.spin[self._material][:] = self.spin_last[self._material][:] - alpha * alpha_scale * self.gradE_last[self._material][:]
                self.spin[:] = self.spin_last[:] - alpha * alpha_scale * self.gradE_last[:]

                # Normalize spin
                self.normalise_spin(self.spin)
                print(self.spin.reshape(-1, 3)[:10])
                if creepCount > 4:
                    exitFlag = True
                    break

                self.compute_effective_field()  # Compute Heff and E using self.spin
                self.gradE[:] = -self.field  # Scale field??
                self.trailE[nStart] = self.totalE
                nStart = next(trailPool)
                gE_reshape = self.gradE.reshape(-1, 3)
                np.einsum('ij,ij->i', gE_reshape, gE_reshape, out=self.gradE2)
                deltaE = abs(self.trailE[nStart] - self.totalE) / nTrail
                print('Tot E', self.totalE)
                print('Tot E last', self.totalE_last)

                with np.printoptions(precision=2):
                    print('Creep: ', self.trailE)
                    print(f': alpha = {alpha}  maxgradE = {self.gradE.max()}')


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
                    print('Total E < total E last')
                    creepCount += 1
                    # Update Energy, spin and gradE
                    self.spin_last[:] = self.spin[:]
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
