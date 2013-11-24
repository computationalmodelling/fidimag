import clib 
import numpy as np


class DMI(object):
    """
    Hamiltonian = D*[S_i x S_j]
        ==> H = D x S_j
    """
    def __init__(self,D,direction=(1,0,0),name='dmi'):
        self.Dx=D*direction[0]
        self.Dy=D*direction[1]
        self.Dz=D*direction[2]
        self.name=name
        
    def setup(self,mesh,spin,unit_length=1.0,mu_s=1.0):
        self.mesh=mesh
        self.spin=spin
        
        self.Dx_mu_s=self.Dx/mu_s
        self.Dy_mu_s=self.Dy/mu_s
        self.Dz_mu_s=self.Dz/mu_s
        
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        
        self.nxyz=mesh.nxyz
        self.field=np.zeros(3*self.nxyz)
        self.energy=0

    def compute_field(self):
        clib.compute_dmi_field(self.spin,
                               self.field,
                               self.Dx_mu_s,
                               self.Dy_mu_s,
                               self.Dz_mu_s,
                               self.nx,
                               self.ny,
                               self.nz)
        
        return -self.field
    
    def compute_energy(self):
        self.energy=clib.compute_dmi_energy(self.spin,
                                            self.Dx,
                                            self.Dy,
                                            self.Dz,
                                            self.nx,
                                            self.ny,
                                            self.nz)
        
        return self.energy
