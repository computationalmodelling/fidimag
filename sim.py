import clib
import numpy as np
from scipy.integrate import ode
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from show_vector import VisualSpin
import time


class Sim(object):
    def __init__(self,mesh):
        self.t=0
        self.c=0.1
        self.mesh=mesh
        self.nxyz=mesh.nxyz
        self.spin=np.ones(3*self.nxyz)
        self.field=np.zeros(3*self.nxyz)
        self.dm_dt=np.zeros(3*self.nxyz)
        self.interactions=[]
        self.vode=ode(self.ode_rhs)
        self.pin_fun=None
        self.ode_times=0
        self.set_options()
        
    def set_options(self,rtol=1e-8,atol=1e-12):
        self.vode.set_integrator('vode',
                                 rtol=rtol,
                                 atol=atol,
                                 nsteps=100000)
        self.alpha=0.01
        self.gamma=1
        self.mu_s=1
        
    
    def set_m(self,m0=(1,0,0)):
        if isinstance(m0,list) or isinstance(m0,tuple):
            tmp=np.zeros((self.nxyz,3))
            tmp[:,:]=m0
            tmp=np.reshape(tmp, 3*self.nxyz, order='F')
            self.spin[:]=tmp[:]
        elif hasattr(m0, '__call__'):
            tmp=np.zeros((self.nxyz,3))
            for i in range(self.mesh.nxyz):
                tmp[i]=m0(self.mesh.pos[i])
            tmp=np.reshape(tmp, 3*self.nxyz, order='F')
            self.spin[:]=tmp[:]
        elif isinstance(m0,np.ndarray):
            if m0.shape==self.spin.shape:
                self.spin[:]=m0[:]
            
        self.vode.set_initial_value(self.spin, self.t)
            
    def add(self,interaction):
        interaction.setup(self.mesh,self.spin)
        self.interactions.append(interaction)

    def run_until(self,t):
        ode=self.vode
        while ode.successful() and ode.t<t:
            ode.integrate(t)
            #print ode.t,ode.y

    def ode_rhs(self,t,y):
        self.ode_times+=1
        self.t=t
        self.field=0
        self.spin[:]=y[:]
        
        if self.pin_fun:
            self.pin_fun(self.t,self.mesh,self.spin)
        
        for obj in self.interactions:
            self.field+=obj.compute_field()
        
        clib.compute_llg_rhs(self.dm_dt,
                           self.spin,
                           self.field,
                           self.gamma,
                           self.alpha,
                           self.mu_s,
                           self.nxyz,
                           self.c)
        
        return self.dm_dt

                        
if __name__=='__main__':
    
    mesh=FDMesh(nx=1)
    sim=Sim(mesh)
    sim.alpha=0.00
    exch=UniformExchange(1)
    sim.add(exch)

    anis=Anisotropy(1)
    #sim.add(anis)
    
    zeeman=Zeeman(10,(0,0,1))
    sim.add(zeeman)
    
    sim.set_m((1,0,0))
    #vs=VisualSpin(sim)
    #vs.init()
    
    print exch.compute_field()
    #print anis.compute_field()
    
    ts=np.linspace(0, 50, 5001)
    run_times=[]
    mzs=[]
    sim.c=10
    for t in ts:
        sim.run_until(t)
        spin=sim.spin
        dm=np.linalg.norm(spin)-1.0
        print sim.c,'times',sim.ode_times,dm,spin[2]
        mzs.append(spin[2])
        run_times.append(dm)
        #vs.update()
        #time.sleep(0.01)
    import pylab
    pylab.plot(ts,run_times)
    pylab.show()
    pylab.plot(ts,mzs)
    pylab.show()
    
    
