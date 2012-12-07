import clib 
import numpy as np


class Demag(object):
    def __init__(self):
        pass
        
        
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
        
        lenx = 2 * self.nx - 1
        leny = n * self.ny - 1
        lenz = 2 * self.nz - 1
        length = lenx * leny * lenz
        
        self.Nxx = np.zeros(length)
        self.Nyy = np.zeros(length)
        self.Nzz = np.zeros(length)
        self.Nxy = np.zeros(length)
        self.Nxz = np.zeros(length)
        self.Nyz = np.zeros(length)
        
    def compute_tensor(self):
        clib.compute_all_tensors(self.Nxx,
                                 self.Nyy,
                                 self.Nzz,
                                 self.Nxy,
                                 self.Nxz,
                                 self.Nyz,
                                 self.dx,self.dy,self.dz,
                                 self.nx,self.ny,self.nz)

    def compute_field(self):
        clib.compute_uniform_exchange(self.spin,
                                      self.field,
                                      self.J,
                                      self.dx,
                                      self.dy,
                                      self.dz,
                                      self.nx,
                                      self.ny,
                                      self.nz)
                                      
        return self.field
