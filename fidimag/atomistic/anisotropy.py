import fidimag.extensions.clib as clib
import numpy as np
from energy import Energy
import fidimag.util.helper as helper


class Anisotropy(Energy):

    """
    Hamiltonian = -Ku_x*S_x^2 - Ku_y*S_y^2 -Ku_z*S_z^2 
    ==>H=2D*S_x
    """

    def __init__(self, Ku, axis=(1, 0, 0), name='anis', direction=None):
        self._Ku = Ku
        self.name = name
        self.axis = axis
        self.jac = True

        if direction is not None:
            self.axis = direction

    def setup(self, mesh, spin, mu_s):
        super(Anisotropy, self).setup(mesh, spin, mu_s)

        self.Ku = np.zeros(3 * self.nxyz, dtype=np.float)
        Ku_scalar = helper.init_scalar(self._Ku, self.mesh)
        self.Ku.shape = (3, -1)
        self.Ku[0, :] = Ku_scalar[:] * self.axis[0]
        self.Ku[1, :] = Ku_scalar[:] * self.axis[1]
        self.Ku[2, :] = Ku_scalar[:] * self.axis[2]
        self.Ku.shape = (-1,)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin
        clib.compute_anisotropy(m,
                                self.field,
                                self.Ku,
                                self.nxyz)

        return self.field * self.mu_s_inv

    def compute_energy(self):
        self.energy = clib.compute_anisotropy_energy(self.spin,
                                                     self.Ku,
                                                     self.nxyz)

        return self.energy
