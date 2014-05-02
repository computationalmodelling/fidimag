import clib 
import numpy as np


class DMI(object):
    """
    Hamiltonian = D*[S_i x S_j]
        ==> H = D x S_j
    """
    def __init__(self,D,name='dmi'):
        self.D = D
        self.name=name
        
    def setup(self,mesh,spin,mu_s_inv,pbc=None):
        self.mesh=mesh
        self.spin=spin
        self.mu_s_inv = mu_s_inv
        
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        
        self.nxyz=mesh.nxyz
        self.field=np.zeros(3*self.nxyz)
        self.energy=0
    
        self.xperiodic = 0
        self.yperiodic = 0
        
        if pbc=='1d':
            self.xperiodic = 1
        elif pbc=='2d':
            self.xperiodic = 1
            self.yperiodic = 1

    def compute_field(self):
        
        clib.compute_dmi_field(self.spin,
                                    self.field,
                                    self.D,
                                    self.nx,
                                    self.ny,
                                    self.nz,
                                    self.xperiodic,
                                    self.yperiodic)
        
        return self.field*self.mu_s_inv
    
    def compute_energy(self):
        self.energy=clib.compute_dmi_energy(self.spin,
                                            self.D,
                                            self.nx,
                                            self.ny,
                                            self.nz,
                                            self.xperiodic,
                                            self.yperiodic)
        
        return self.energy
