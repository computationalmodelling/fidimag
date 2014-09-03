import micro_clib
import numpy as np
from pc.constant import mu_0

class DMI(object):
    """
        compute the DMI field in micromagnetics
    """
    def __init__(self,D,name='dmi', type='bulk'):
        """
        type could be 'interfacial' or 'bulk'
        """
        self.D = D
        self.name=name
        self.type=type
    
        
    def setup(self, mesh, spin, Ms):
        self.mesh=mesh
        self.dx=mesh.dx*mesh.unit_length
        self.dy=mesh.dy*mesh.unit_length
        self.dz=mesh.dz*mesh.unit_length
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.spin=spin
            
        self.field=np.zeros(3*mesh.nxyz, dtype=np.float)
        self.energy=np.zeros(3*mesh.nxyz, dtype=np.float)
        self.total_energy = 0
        self.pbc = mesh.pbc
        self.Ms = Ms
        self.Ms_inv = np.zeros(mesh.nxyz)
            
        for i in range(mesh.nxyz):
            if self.Ms[i] == 0.0:
                self.Ms_inv[i] = 0
            else:
                self.Ms_inv[i] = 1.0/self.Ms[i]
            
        self.xperiodic = mesh.xperiodic
        self.yperiodic = mesh.yperiodic
        
    def compute_field(self, t=0):
        if self.type == 'bulk':
            micro_clib.compute_dmi_field_bulk(self.spin,
                                                self.field,
                                                self.energy,
                                                self.Ms_inv,
                                                self.D,
                                                self.dx,
                                                self.dy,
                                                self.dz,
                                                self.nx,
                                                self.ny,
                                                self.nz,
                                                self.xperiodic,
                                                self.yperiodic)
        elif self.type == 'interfacial':
            micro_clib.compute_dmi_field_interfacial(self.spin,
                                                    self.field,
                                                    self.energy,
                                                    self.Ms_inv,
                                                    self.D,
                                                    self.dx,
                                                    self.dy,
                                                    self.dz,
                                                    self.nx,
                                                    self.ny,
                                                    self.nz,
                                                    self.xperiodic,
                                                    self.yperiodic)
        else:
            raise Exception("Unsppourted dmi type:{}, avaiable type: 'bulk','interfacial'.".format(self.type))
        
        return self.field
    
    def compute_energy(self):
        
        #since we are not always calling this function, so it's okay to call compute_field again
        self.compute_field()
        
        self.total_energy = np.sum(self.energy)
        
        return self.total_energy*self.mesh.cellsize
    

