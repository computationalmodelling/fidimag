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
from finmag.native import sundials



class Sim(object):
    def __init__(self,mesh):
        self.t=0
        self.c=0.1
        self.mesh=mesh
        self.nxyz=mesh.nxyz
        self.unit_length=mesh.unit_length
        self.spin=np.ones(3*self.nxyz)
        self.field=np.zeros(3*self.nxyz)
        self.dm_dt=np.zeros(3*self.nxyz)
        self.interactions=[]
        self.vode=ode(self.ode_rhs)
        self.pin_fun=None
        self.ode_times=0
        self.set_up_solver()
        self.set_options()
        
    def set_options(self,rtol=1e-10,atol=1e-20,mat=Nickel()):
        self.vode.set_integrator('vode',
                                 rtol=rtol,
                                 atol=atol,
                                 first_step=1e-15,
                                 nsteps=100000)
        
        self.mat=mat
        self.mu_s=1 #since we already consider mu_s in fields
        
    
    def set_up_solver(self, reltol=1e-10, abstol=1e-20, nsteps=10000):
        integrator = sundials.cvode(sundials.CV_BDF, sundials.CV_NEWTON)
        integrator.init(self.sundials_rhs, 0, self.spin)
        integrator.set_linear_solver_sp_gmr(sundials.PREC_NONE)
        integrator.set_scalar_tolerances(reltol, abstol)
        integrator.set_max_num_steps(nsteps)
        
        self.integrator = integrator
        
    
    def sundials_rhs(self, t, y, ydot):

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
                           self.mat.gamma,
                           self.mat.alpha,
                           self.mu_s,
                           self.nxyz,
                           self.c)
        
        ydot[:] = self.dm_dt[:]
        
        return 0
        
    
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
        interaction.setup(self.mesh,self.spin,unit_length=self.unit_length)
        self.interactions.append(interaction)

    def run_until(self,t):
        ode=self.vode
        while ode.successful() and ode.t<t:
            ode.integrate(t)
            #print ode.t,ode.y
            
    def run_until2(self,t):
        if t <= self.t:
            return
        
        self.integrator.advance_time(t, self.spin)

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
                           self.mat.gamma,
                           self.mat.alpha,
                           self.mu_s,
                           self.nxyz,
                           self.c)
        
        return self.dm_dt

                        
if __name__=='__main__':
    
    mesh=FDMesh(nx=1)
    sim=Sim(mesh)
    
    ni=Nickel()
    sim.set_options(mat=ni)
    ni.alpha=0.01
    
    exch=UniformExchange(ni.J,mu_s=ni.mu_s)
    #sim.add(exch)

    anis=Anisotropy(ni.D,mu_s=ni.mu_s)
    #sim.add(anis)
    
    zeeman=Zeeman(1e5,(0,0,1))
    sim.add(zeeman)
    
    sim.set_m((1,0,0))
    #vs=VisualSpin(sim)
    #vs.init()
    
    #print exch.compute_field()
    #print anis.compute_field()
    
    ts=np.linspace(0, 1e-10, 100)
    run_times=[]
    mzs=[]
    sim.c=100
    for t in ts:
        sim.run_until(t)
        spin=sim.spin
        #dm=np.linalg.norm(spin)-1.0
        #print sim.c,'times',sim.ode_times,dm,spin[2]
        mzs.append(spin[2])
        #run_times.append(dm)
        #vs.update()
        #time.sleep(0.01)
    print mzs
    import pylab
    #pylab.plot(ts,run_times)
    #pylab.show()
    pylab.plot(ts,mzs)
    pylab.show()
    
    
