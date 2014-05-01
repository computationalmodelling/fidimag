import numpy as np

class Zeeman(object):
    def __init__(self,H0,direction=(1,0,0),name='Zeeman'):
        self.Hx=H0*direction[0]
        self.Hy=H0*direction[1]
        self.Hz=H0*direction[2]
        self.name = name
        
        
    def setup(self,mesh,spin,unit_length=1.0,mu_s=1.0):
        self.mesh=mesh
        self.dx=mesh.dx
        self.dy=mesh.dy
        self.dz=mesh.dz
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.spin=spin
        n=self.nx*self.ny*self.nz
        self.field=np.zeros((3,n))
        self.field[:,:]=[[self.Hx],[self.Hy],[self.Hz]]
        self.field.shape=(3*n)
        #self.field[:]*=(4*np.pi*1e-7)
        

    def compute_field(self):            
        return self.field

    def compute_energy(self):
        # will add later
        return 0