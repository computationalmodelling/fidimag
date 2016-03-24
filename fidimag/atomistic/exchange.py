import fidimag.extensions.clib as clib
import numpy as np
from .energy import Energy


class UniformExchange(Energy):

    """
    The Hamiltonian is defined as

        Hamiltonian = - J \sum_<i,j> S_i \cdot S_j

    where the brackets represent the nearest neighbours and only evaluate once
    for each pair, which means for two spins case, the total energy is
    -J S_1 S_2. Therefore, the effective field at site i is,

        H_i = J \sum_<i,j> S_j

    notice that there is no factor of 2 associated with J.
    """

    def __init__(self, J=0, name='exch'):
        self.J = J
        self.name = name

        self.Jx = self.J
        self.Jy = self.J
        self.Jz = self.J
        self.jac = True

    def compute_field(self, t=0, spin=None):

        if spin is not None:
            m = spin
        else:
            m = self.spin

        clib.compute_exchange_field(m,
                                    self.field,
                                    self.energy,
                                    self.Jx,
                                    self.Jy,
                                    self.Jz,
                                    self.neighbours,
                                    self.n)

        return self.field * self.mu_s_inv

    def compute_energy_directly(self):

        energy = clib.compute_exchange_energy(self.spin,
                                              self.Jx,
                                              self.Jy,
                                              self.Jz,
                                              self.nx,
                                              self.ny,
                                              self.nz,
                                              self.xperiodic,
                                              self.yperiodic)

        return energy


class Exchange(Energy):

    """
    The Hamiltonian is defined as

        Hamiltonian = -  \sum_<i,j> J_ij S_i \cdot S_j

    where the brackets represent the nearest neighbours and only evaluate once
    for each pair, which means for two spins case, the total energy is
    -J S_1 S_2. Therefore, the effective field at site i is,

        H_i = \sum_<i,j> J_ij S_j

    notice that there is no factor of 2 associated with J_ij.
    """

    def __init__(self, J_fun, name='exch'):
        self.J_fun = J_fun
        self.name = name
        self.jac = False
    
    def setup(self, mesh, spin, mu_s):
        super(Exchange, self).setup(mesh, spin, mu_s)

        self._J = np.zeros(self.neighbours.shape)

        if isinstance(self.J_fun, (int, float)):
            self._J[:,:] = self.J_fun
        elif hasattr(self.J_fun, '__call__'):
            n =  self.mesh.n
            for i in range(n):
                value = self.J_fun(self.coordinates[i])
                if len(value) == 6:
                    self._J[i,:] = value[:]
                else:
                    raise Exception('The given spatial function for J is not acceptable!')
            pass
        
        
    def compute_field(self, t=0, spin=None):

        if spin is not None:
            m = spin
        else:
            m = self.spin

        clib.compute_exchange_field_spatial(m,
                                    self.field,
                                    self.energy,
                                    self._J,
                                    self.neighbours,
                                    self.n)

        return self.field * self.mu_s_inv
