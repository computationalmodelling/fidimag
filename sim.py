import clib
import time
import numpy as np
from scipy.integrate import ode
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from materials import Nickel
#from show_vector import VisualSpin



class Sim(object):
    def __init__(self,mesh,T=0):
        self.t=0
        self.T=T
        self.mesh=mesh
        self.nxyz=mesh.nxyz
        self.unit_length=mesh.unit_length
        self.spin=np.ones(3*self.nxyz)
        self.field=np.zeros(3*self.nxyz)
        self.dm_dt=np.zeros(3*self.nxyz)
        self.interactions=[]
        
        self.pin_fun=None
        self.ode_times=0
        self.set_options()
        
    def set_options(self,rtol=1e-8,atol=1e-20,mat=Nickel(),dt=1e-14):
        self.mat=mat
        self.mu_s=1 #since we already consider mu_s in fields
        self.c=1e11 
        
        if self.T>0:
            self.vode=clib.RK2S(mat.mu_s,
                                dt,
                                self.nxyz,
                                mat.gamma,
                                mat.alpha,
                                self.T,
                                self.c*10,
                                self.spin,
                                self.field,
                                self.stochastic_update_field)
        else:
            self.vode=ode(self.ode_rhs)
            self.vode.set_integrator('vode',
                                 rtol=rtol,
                                 atol=atol,
                                 nsteps=100000)
        
        
        
        
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
        #self.integrator.init(self.sundials_rhs, 0, self.spin)
            
    def add(self,interaction):
        interaction.setup(self.mesh,self.spin,unit_length=self.unit_length)
        self.interactions.append(interaction)

    def run_until(self,t):
        ode=self.vode
        while ode.successful() and ode.t<t:
            #print t,real_t
            ode.integrate(t)
            self.spin[:]=ode.y[:]
            
    def ode_rhs(self,t,y):
        self.ode_times+=1
        self.t=t
        self.field[:]=0
        self.spin[:]=y[:]
        
        if self.pin_fun:
            self.pin_fun(self.t,self.mesh,self.spin)
        
        for obj in self.interactions:
            self.field+=obj.compute_field()
        
        clib.compute_llg_rhs(self.dm_dt,
                           self.spin,
                           self.field,
                           self.mat.gamma,
                           self.mat.alpha,
                           self.mu_s,
                           self.nxyz,
                           self.c)
        
        return self.dm_dt
    
    def stochastic_update_field(self,y):
        self.field[:]=0
        self.spin[:]=y[:]
        
        if self.pin_fun:
            self.pin_fun(self.t,self.mesh,self.spin)
        
        for obj in self.interactions:
            self.field+=obj.compute_field()
            
        #print "from python",self.spin,self.field
        
                        
if __name__=='__main__':
    
    mesh=FDMesh(nx=1)
    sim=Sim(mesh,T=0.1)
    
    ni=Nickel()
    ni.alpha=0.1
    sim.set_options(mat=ni)
    
    
    exch=UniformExchange(ni.J,mu_s=ni.mu_s)
    sim.add(exch)

    anis=Anisotropy(ni.D,mu_s=ni.mu_s)
    sim.add(anis)
    
    zeeman=Zeeman(1e5,(0,0,1))
    sim.add(zeeman)
    
    sim.set_m((1,0,0))
    print sim.vode.successful()
    print sim.vode.successful()
    #vs=VisualSpin(sim)
    #vs.init()
    
    
    #print exch.compute_field()
    #print anis.compute_field()
    
    ts=np.linspace(0, 1e-10, 100)
    run_times=[]
    mzs=[]
    for t in ts:
        sim.run_until(t)
        spin=sim.spin
        print 'from sim',sim.field,sim.spin
        dm=np.linalg.norm(spin)-1.0
        #print sim.c,'times',sim.ode_times,dm,spin[2]
        mzs.append(spin[2])
        run_times.append(dm)
        #vs.update()
        #time.sleep(0.01)
    
    print mzs
    import pylab
    pylab.plot(ts,run_times)
    pylab.show()
    pylab.plot(ts,mzs)
    pylab.show()
    
    
    
    
