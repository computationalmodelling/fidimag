import clib 
import numpy as np


class Demag(object):
    def __init__(self,mu_s=1.0):
        self.mu_s=mu_s
        
    def setup(self, mesh, spin, unit_length=1.0):
        self.mesh = mesh
        self.dx = mesh.dx
        self.dy = mesh.dy
        self.dz = mesh.dz
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = spin
        self.n = self.nx * self.ny * self.nz
        self.field = np.zeros(3 * self.n)
        self.mu_s_modifid=self.mu_s/(unit_length**3)
        
        self.demag=clib.FFTDemag(self.mu_s_modifid,
                                 self.dx,self.dy,self.dz,
                              self.nx,self.ny,self.nz)
        
    def compute_field(self):
        self.demag.compute_field(self.spin,self.field)
        return self.field
    
    def compute_exact(self):
        field = np.zeros(3 * self.n)
        self.demag.compute_exact(self.spin,field)
        return field
    

if __name__=='__main__':
    import pccp
    mesh=pccp.FDMesh(nx=4,ny=3,nz=2)
    sim=pccp.Sim(mesh)
    
    demag=Demag()
    sim.add(demag)
    
    def init_m(pos):
        x,y,z=pos
        if x<-2:
            return (1,0,0)
        elif x>2:
            return (-1,0,0)
        else:
            return (0,1,0)
    
    sim.set_m(init_m)
    fft=demag.compute_field()
    exact=demag.compute_exact()
    print fft
    print np.max(fft-exact)
