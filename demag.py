import clib 
import numpy as np


class Demag(object):
        
    def setup(self, mesh, spin):
        self.mesh = mesh
        self.dx = mesh.dx
        self.dy = mesh.dy
        self.dz = mesh.dz
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = spin
        n = self.nx * self.ny * self.nz
        self.field = np.zeros(3 * n)
        
        self.demag=clib.FFTDemag(self.dx,self.dy,self.dz,
                              self.nx,self.ny,self.nz)
        
    def compute_field(self):
        self.demag.compute_field(self.spin,self.field)
        return self.field
    
    def compute_exact(self):
        self.demag.compute_exact(self.spin,self.field)
        return self.field
    

if __name__=='__main__':
    import pccp
    mesh=pccp.FDMesh(nx=6,ny=1,nz=1)
    sim=pccp.Sim(mesh)
    
    demag=Demag()
    sim.add(demag)
    
    sim.set_m((1,0,0))
    print demag.compute_field()
    print demag.compute_exact()
    demag.demag.print_field()
