import clib
import numpy as np
from scipy.integrate import ode
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from show_vector import VisualSpin


class Sim(object):
    def __init__(self,mesh):
        self.mesh=mesh
        self.nxyz=mesh.nxyz
        self.spin=np.ones(3*self.nxyz)
        self.field=np.zeros(3*self.nxyz)
        self.dm_dt=np.zeros(3*self.nxyz)
        self.interactions=[]
        self.vode=ode(self.ode_rhs)
        
        self.set_options()
        
    def set_options(self,rtol=1e-8,atol=1e-12):
        self.vode.set_integrator('vode',
                                 rtol=rtol,
                                 atol=atol,
                                 nsteps=100000)
        self.alpha=0.1
        self.gamma=1
        self.mu_s=1
    
    def set_m(self,m0=(1,0,0)):
        if isinstance(m0,list) or isinstance(m0,tuple):
            tmp=np.zeros((self.nxyz,3))
            tmp[:,:]=m0
            tmp=np.reshape(tmp, 3*self.nxyz, order='F')
            self.spin[:]=tmp[:]
        elif hasattr(m0, '__call__'):
            m0(self.spin)
            
        self.vode.set_initial_value(self.spin, 0)
            
    def add(self,interaction):
        interaction.setup(self.mesh,self.spin)
        self.interactions.append(interaction)

    def run_until(self,t):
        ode=self.vode
        while ode.successful() and ode.t<t:
            ode.integrate(t)
            print ode.t,ode.y

    def ode_rhs(self,t,y):
        self.field=0
        self.spin[:]=y[:]
        for obj in self.interactions:
            self.field+=obj.compute_field()
        
        clib.compute_llg_rhs(self.dm_dt,
                           self.spin,
                           self.field,
                           self.gamma,
                           self.alpha,
                           self.mu_s,
                           self.nxyz)
        
        return self.dm_dt

                        
if __name__=='__main__':
    
    mesh=FDMesh()
    sim=Sim(mesh)
    exch=UniformExchange(1)
    sim.add(exch)

    anis=Anisotropy(1)
    sim.add(anis)
    
    zeeman=Zeeman(10,(0,0,1))
    sim.add(zeeman)
    
    sim.set_m((1,0,0))
    vs=VisualSpin(sim)
    vs.init()
    
    print exch.compute_field()
    print anis.compute_field()
    
    for t in range(10):
        sim.run_until(t)
        vs.update()
    
    
