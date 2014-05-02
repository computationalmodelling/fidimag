import numpy as np

import pccp.util.helper as helper

class Zeeman(object):
    """
    The time independent external field, can vary with space
    """
    def __init__(self,H0,name='Zeeman'):
        self.H0 = H0
        self.name = name
        
    def setup(self,mesh,spin,mu_s_inv, pbc=None):
        self.mesh=mesh
        self.spin=spin
        self.nxyz = mesh.nxyz
        
        self.mu_s = np.zeros(mesh.nxyz)
        
        for i in range(self.nxyz):
            if mu_s_inv[i] == 0.0:
                self.mu_s[i] = 0
            else: 
                self.mu_s[i] = 1.0/mu_s_inv[i]
        
        self.field=np.zeros(3*self.nxyz)
        self.field[:]=helper.init_vector(self.H0, self.mesh)
        

    def compute_field(self, t=0):     
        return self.field

    def compute_energy(self):
        sf = self.field*self.spin
        sf.shape=(3,-1)
        
        energy_d = -np.sum(sf, axis=0)*self.mu_s
        
        sf.shape=(-1,)
        
        return np.sum(energy_d)
    

class TimeZeeman(object):
    """
    The time dependent external field, also can vary with space
    """
    def __init__(self,H0,time_fun, name='TimeZeeman'):
        self.H0 = H0
        self.time_fun = time_fun
        self.name = name
        
    def setup(self,mesh,spin,mu_s_inv, pbc=None):
        self.mesh=mesh
        self.spin=spin
        self.nxyz = mesh.nxyz
        
        self.mu_s = np.zeros(mesh.nxyz)
        
        for i in range(self.nxyz):
            if mu_s_inv[i] == 0.0:
                self.mu_s[i] = 0
            else: 
                self.mu_s[i] = 1.0/mu_s_inv[i]
        
        self.field=np.zeros(3*self.nxyz)
        self.H_init=np.zeros(3*self.nxyz)
        self.H_init[:]=helper.init_vector(self.H0, self.mesh)
    
    def compute_field(self, t=0):
        self.field[:] = self.H_init[:]*self.time_fun(t)
        return self.field

    def compute_energy(self):
        sf = self.field*self.spin
        sf.shape=(3,-1)
        
        energy_d = -np.sum(sf, axis=0)*self.mu_s
        
        sf.shape=(-1,)
        
        return np.sum(energy_d)