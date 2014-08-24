import numpy as np

import util.helper as helper

class Zeeman(object):
    """
    The time independent external field, can vary with space
    """
    def __init__(self,H0,name='Zeeman'):
        self.H0 = H0
        self.name = name
        
    def setup(self,mesh,spin, Me, pbc=None):
        self.mesh = mesh
        self.spin = spin
        self.nxyz = mesh.nxyz
                
        self.field = np.zeros(3*self.nxyz)
        self.field[:] = helper.init_vector(self.H0, self.mesh)
        #print self.field

    def compute_field(self, t=0):     
        return self.field
    
    #Todo: update it later
    def average_field(self):
        hx = self.field[0]
        hy = self.field[self.nxyz]
        hz = self.field[2*self.nxyz]
        return np.array([hx,hy,hz])

    def compute_energy(self):
        
        return 0.0
    

class TimeZeeman(object):
    """
    The time dependent external field, also can vary with space
    """
    def __init__(self, H0,time_fun, name='TimeZeeman'):
        self.H0 = H0
        self.time_fun = time_fun
        self.name = name
        
    def setup(self,mesh, spin, Ms, pbc=None):
        self.mesh=mesh
        self.spin=spin
        self.nxyz = mesh.nxyz
        
        self.mu_s = np.zeros(3*mesh.nxyz)
        #FIX Ms
        
        self.field=np.zeros(3*self.nxyz)
        self.H_init=np.zeros(3*self.nxyz)
        self.H_init[:]=helper.init_vector(self.H0, self.mesh)
    
    def compute_field(self, t=0):
        self.field[:] = self.H_init[:]*self.time_fun(t)
        return self.field
    
    #Todo: update it later
    def average_field(self):
        hx = self.field[0]
        hy = self.field[self.nxyz]
        hz = self.field[2*self.nxyz]
        return np.array([hx,hy,hz])

    def compute_energy(self):
        sf = self.field*self.spin*self.mu_s
        energy = -np.sum(sf)
        return energy