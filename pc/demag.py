import clib 
import numpy as np


class Demag(object):
    def __init__(self, name='demag'):
        self.name = name
        
    def setup(self, mesh, spin, mu_s=1,unit_length=1.0):
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
        #note that the mu_s we used is from the class init 
        self.scale = 1e-7 / unit_length**3
        self.demag=clib.FFTDemag(mu_s,
                                 self.dx,self.dy,self.dz,
                                 self.nx,self.ny,self.nz)
        
    def compute_field(self):
        self.demag.compute_field(self.spin,self.field)
        self.field[:]*=self.scale
        return self.field
    
    def compute_exact(self):
        field = np.zeros(3 * self.n)
        self.demag.compute_exact(self.spin,field)
        field[:]*=self.scale
        return field

    def compute_energy(self):
        
        energy=self.demag.compute_energy(self.spin,self.field)
        
        return energy
    

if __name__=='__main__':
    from pccp.pc.fd_mesh import FDMesh
    from pccp.pc.sim import Sim
    
    mesh = FDMesh(nx=6,ny=1,nz=1)
    sim = Sim(mesh)
    
    demag=Demag(mu_s=1e3)
    sim.add(demag)
    
    def init_m(pos):
        x,y,z=pos
        if x<=2:
            return (1,0,0)
        elif x>=4:
            return (0,0,1)
        else:
            return (0,1,0)
    
    sim.set_m(init_m)
    fft=demag.compute_field()
    exact=demag.compute_exact()
    print fft
    print exact
    print np.max(fft-exact)
