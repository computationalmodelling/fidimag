import clib
import time
import numpy as np
from scipy.integrate import ode
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from demag import Demag
from materials import Nickel
from materials import UnitMaterial
#from show_vector import VisualSpin



class Sim(object):
    def __init__(self,mesh,T=0,mat=UnitMaterial()):
        self.t=0
        self.T=T
        self.mesh=mesh
        self.nxyz=mesh.nxyz
        self.unit_length=mesh.unit_length
        self.spin=np.ones(3*self.nxyz)
        self.field=np.zeros(3*self.nxyz)
        self.dm_dt=np.zeros(3*self.nxyz)
        self.interactions=[]
        self.mat=mat
        self.pin_fun=None
        self.ode_times=0
        self.set_options()
        
    def set_options(self,rtol=1e-8,atol=1e-20,dt=1e-15):
        
        self.mu_s=1 #since we already consider mu_s in fields
        self.c=1e11 
        
        if self.T>0:
            self.vode=clib.RK2S(self.mat.mu_s,
                                dt,
                                self.nxyz,
                                self.mat.gamma,
                                self.mat.alpha,
                                self.T,
                                self.c,
                                self.spin,
                                self.field,
                                self.stochastic_update_field)
        else:
            self.vode=ode(self.ode_rhs)
            self.vode.set_integrator('vode',
                                 rtol=rtol,
                                 atol=atol,
                                 nsteps=100000)
        self.gamma=self.mat.gamma
        self.alpha=self.mat.alpha
        
        
        
        
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
            print 'before',tmp
            tmp=np.reshape(tmp, 3*self.nxyz, order='F')
            self.spin[:]=tmp[:]
            print 'after',tmp
            
            for i in range(self.nxyz):
                j=i+self.nxyz
                k=i+self.nxyz
                tmp=m0(self.mesh.pos[i])
                self.spin[i]=tmp[0]
                self.spin[j]=tmp[1]
                self.spin[k]=tmp[2]
                
            print 'after2',self.spin
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
        #self.spin[:]=y[:]
        
        if self.pin_fun:
            self.pin_fun(self.t,self.mesh,self.spin)
        
        for obj in self.interactions:
            self.field+=obj.compute_field()
            
    def compute_average(self):
        self.spin.shape=(3,-1)
        average=np.sum(self.spin,axis=1)/self.nxyz
        self.spin.shape=(3*self.nxyz)
        return average
            
        #print "from python",self.spin,self.field
    
    #only used for tests
    def spin_at(self,i,j,k):
        nxyz=self.mesh.nxyz
        nyz=self.mesh.nyz
        nz=self.mesh.nz
        i1=i*nyz+j*nz+k
        i2=i1+nxyz
        i3=i2+nxyz
        return (self.spin[i1],
                self.spin[i2],
                self.spin[i3])
        
                        
if __name__=='__main__':
    
    T=1000
    ni=Nickel()
    ni.alpha=0.1
    
    mesh=FDMesh(nx=1,ny=1,nz=1)
    mesh.set_material(ni)
    
    sim=Sim(mesh,T=T,mat=ni) 
    
    exch=UniformExchange(ni.J,mu_s=ni.mu_s)
    sim.add(exch)

    anis=Anisotropy(ni.D,mu_s=ni.mu_s)
    sim.add(anis)
    
    zeeman=Zeeman(1e5,(0,0,1))
    sim.add(zeeman)
    
    #demag=Demag(mu_s=ni.mu_s)
    #sim.add(demag)
    
    sim.set_m((1,0,0))
    print sim.vode.successful()
    print sim.vode.successful()
    #vs=VisualSpin(sim)
    #vs.init()
    
    
    #print exch.compute_field()
    #print anis.compute_field()
    
    ts=np.linspace(0, 1e-10, 100)
    mxs=[]
    mys=[]
    mzs=[]
    dms=[]
    for t in ts:
        sim.run_until(t)
        spin=sim.spin
        #print 'from sim',sim.field,sim.spin
        dm=np.linalg.norm(spin)-1.0
        dms.append(dm)
        #print sim.c,'times',sim.ode_times,dm,spin[2]
        
        av=sim.compute_average()
        mxs.append(av[0])
        mys.append(av[1])
        mzs.append(av[2])
        print av
        
  
        #vs.update()
        #time.sleep(0.01)
    
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    fig=plt.figure()
    plt.plot(ts,mxs,'^-',label='mx')
    plt.plot(ts,mys,'.-',label='my')
    plt.plot(ts,mzs,'o-',label='mz')
    plt.legend()
    fig.savefig("mxyz_T%g.png"%T)
    fig=plt.figure()
    plt.plot(ts,dms,'^-')
    plt.savefig('dm_%g.png'%T)
    
    
    
    
