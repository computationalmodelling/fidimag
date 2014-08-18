import clib 
import numpy as np
from constant import mu_0




class DMI(object):
    """
    Hamiltonian = D*[S_i x S_j]
        ==> H = D x S_j
    """
    def __init__(self,D,name='dmi',Dc=None):
        self.D = D
        self.Dc = Dc
        self.name=name
        
        self.Dx = D
        self.Dy = D
        self.Dz = D
        
    def setup(self,mesh,spin,mu_s_inv,pbc=None):
        self.mesh=mesh
        self.spin=spin
        #self.mu_s_inv = mu_s_inv
        self.mu_s_inv = mu_s_inv.copy()
        
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        
        self.dx=mesh.dx*mesh.unit_length
        self.dy=mesh.dy*mesh.unit_length
        self.dz=mesh.dz*mesh.unit_length
        
        if self.Dc is not None:
            self.Dx = 2*self.Dc*self.dy*self.dz
            self.Dy = 2*self.Dc*self.dx*self.dz
            self.Dz = 2*self.Dc*self.dx*self.dy
            self.mu_s_inv[:] *= (1.0/(mu_0*self.dx*self.dy*self.dz))
        
        self.nxyz=mesh.nxyz
        self.field=np.zeros(3*self.nxyz)
        self.energy=np.zeros(3*self.nxyz)
        self.total_energy=0
    
        self.xperiodic = 0
        self.yperiodic = 0
        
        if pbc=='1d':
            self.xperiodic = 1
        elif pbc=='2d':
            self.xperiodic = 1
            self.yperiodic = 1

    def compute_field(self, t=0):
        
        clib.compute_dmi_field(self.spin,
                                    self.field,
                                    self.energy,
                                    self.Dx,
                                    self.Dy,
                                    self.Dz,
                                    self.nx,
                                    self.ny,
                                    self.nz,
                                    self.xperiodic,
                                    self.yperiodic)
        
        return self.field*self.mu_s_inv
    
    def compute_energy(self):
        
        #since we are not always calling this function, so it's okay to call compute_field again
        self.compute_field()
        
        self.total_energy = np.sum(self.energy)/2.0
        
        return self.total_energy
    
    def compute_energy_direct(self):
        """
        mainly for testing
        """
        energy=clib.compute_dmi_energy(self.spin,
                                            self.D,
                                            self.nx,
                                            self.ny,
                                            self.nz,
                                            self.xperiodic,
                                            self.yperiodic)
        
        
        return energy
