import numpy as np
import fidimag.extensions.clib as clib
# import fidimag.common.constant as const
from .minimiser_base import MinimiserBase
from itertools import cycle


class SimpleMinimiser(MinimiserBase):
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
        super(SimpleMinimiser, self).__init__(mesh, spin,
                                              magnetisation, magnetisation_inv,
                                              field,
                                              pins,
                                              interactions,
                                              name,
                                              data_saver
                                              )

        self.t = 1e-4
        self._alpha = 0.1
        self._alpha_field = self._alpha * np.ones_like(self.spin)
        self.spin_last = np.zeros_like(spin)
        self._new_spin = np.zeros_like(spin)

        self.totalRestart = False
        self.exit = False
        self.gradE = np.zeros_like(self.field)
        self.gradE_last = np.zeros_like(self.field)
        self.totalE = 0.0
        self.totalE_last = 0.0

        # Only update site with magnetisation > 0 which are not pinned
        self._material = np.logical_and(np.repeat(self._magnetisation, 3) > 0.0,
                                        np.repeat(1 - self._pins, 3).astype(np.bool))

    # def run_step(self):

    #     self.spin_last[:] = self.spin[:]
    #     self.update_effective_field()

    #     self._new_spin[self._material] = (self.spin + self.field)[self._material]
    #     clib.normalise_spin(self._new_spin, self._pins, self.n)

    #     self._new_spin[self._material] = (self.spin_last +
    #                                       self._alpha_field * (self._new_spin - self.spin_last))[self._material]
    #     self.spin[:] = self._new_spin[:]
    #     clib.normalise_spin(self.spin, self._pins, self.n)

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

    def minimise(self, stopping_dm=1e-2, max_steps=2000,
                 save_data_steps=10, save_m_steps=None, save_vtk_steps=None,
                 log_steps=1000,
                 maxCreep=5, alpha_scale=0.1, stopping_dE=1e-6, dAlpha=2,
                 alphaMin=0.01, perturbSeed=42, perturbFactor=0.1, nTrail=10,
                 resetMax=20
                 ):
        """

        """

        self.step = 0
        nStart = 0
        exitFlag = False
        totalRestart = True
        rstate = np.random.RandomState(perturbSeed)
        resetCount = 0
        creepCount = 0
        self.spin_last[:] = self.spin
        self.trailE = np.zeros(nTrail)
        trailPool = cycle(range(nTrail))  # cycle through 0,1,...,(nTrail-1),0,1,...

        while not exitFlag:

            if totalRestart:
                self.spin[:] = self.spin_last
                # Compute from self.spin:
                self.compute_effective_field()  # Compute Heff and E using self.spin
                self.gradE_last[:] = -self.field  # Scale field??
                self.gradE[:] = self.gradE_last
                self.totalE_last = self.totalE
                self.trailE[nStart] = self.totalE
                nStart = yield(trailPool)
                self._alpha = 1.0
                totalRestart = False

            creepCount = 0

            while self.creepCount < maxCreep:

                self.spin[:] = self.spin_last[:] - self._alpha * alpha_scale * self.gradE_last[:]
                self.compute_effective_field()  # Compute Heff and E using self.spin
                self.gradE[:] = -self.field  # Scale field??
                self.trailE[nStart] = self.totalE
                nStart = yield(trailPool)
                self.gradE2 = np.multiply(self.gradE, self.gradE, out=self.gradE)
                self.gradE2 = np.sum(self.gradE.reshape(-1, 3), axis=1)
                deltaE = abs(self.trailE[nStart] - self.totalE) / nTrail

                if deltaE < stopping_dE:
                    exitFlag = True
                    break  # creep loop

                if self.trailE[nStart] > self.trailE[nStart - 1]:
                    self.creepCount = 0
                    self._alpha = self._alpha / (dAlpha ** 2)

                    if self._alpha < alphaMin:
                        resetCount += 1
                        perturbSpins()
                        if resetCount > resetMax:
                            exitFlag = True
                            print(f'N of resets {resetCount} reached maximum value {resetMax}. Stopping calculation.')
                            break  # creep loop
                        totalRestart = True
                        break  # creep loop
                else:
                    creepCount += 1
                    # Update Energy, spin and gradE
                    self.spin_last[:] = self.spin[:]
                    self.gradE_last[:] = self.gradE[:]
                    self.totalE_last = self.totalE

                        



        self.spin_last[:] = self.spin[:]
        self.compute_effective_field()
        while self.step < max_steps:

            self.run_step()

            max_dm = (self.spin - self.spin_last).reshape(-1, 3) ** 2
            max_dm = np.max(np.sqrt(np.sum(max_dm, axis=1)))

            if self.step % log_steps == 0:
                print("#max_tau={:<8.3g} max_dm={:<10.3g} counter={}".format(
                    np.max(np.abs(self.tau)),
                    max_dm, self.step))

            if max_dm < stopping_dm and self.step > 0:
                print("FINISHED AT: max_tau={:<8.3g} max_dm={:<10.3g} counter={}".format(
                      np.max(np.abs(self.tau)),
                      max_dm, self.step))

                self.compute_effective_field()
                self.data_saver.save()

                break

            if self.step % save_data_steps == 0:
                # update field before saving data
                self.compute_effective_field()
                self.data_saver.save()

            if (save_vtk_steps is not None) and (self.step % save_vtk_steps == 0):
                self.save_vtk()
            if (save_m_steps is not None) and (self.step % save_m_steps == 0):
                self.save_m()

            self.step += 1
