import clib
import numpy as np
from scipy.integrate import ode
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy


class Sim(object):
    def __init__(self,mesh):
        self.mesh=mesh
        self.nxyz=mesh.nxyz
        self.spin=np.zeros((3,self.nxyz))
        self.field=np.zeros((3,self.nxyz))
        self.dm_dt=np.zeros((3,self.nxyz))
        self.interactions=[]
        self.vode=ode(self.ode_rhs)
        
    def set_options(self,rtol=1e-8,atol=1e-12):
        self.vode.set_integrator('vode',
                                 rtol=rtol,
                                 atol=atol)
    
    def set_m(self,m0=(1,0,0)):
        if isinstance(m0,list) or isinstance(m0,tuple):
            self.spin[0,:]=m0[0]
            self.spin[1,:]=m0[1]
            self.spin[2,:]=m0[2]
        elif hasattr(m0, '__call__'):
            m0(self.spin)
            
    def add(self,interaction):
        interaction.setup(self.mesh,self.spin)
        self.interactions.append(interaction)

    def run_until(self,t):

        pass

    def ode_rhs(self,t,y):
        self.field=0
        for obj in self.interactions:
            self.field+=obj.compute_field()
        
        clib.compute_field(self.field,
                           self.spin,
                           self.field,
                           self.gamma,
                           self.alpha,
                           self.mu_s)

                        
if __name__=='__main__':
    
    mesh=FDMesh()
    sim=Sim(mesh)
    exch=UniformExchange(1)
    sim.add(exch)
    

    anis=Anisotropy(1)
    sim.add(anis)
    
    sim.set_m((1,0,0))
    
    print exch.compute_field()
    print anis.compute_field()
    
    
