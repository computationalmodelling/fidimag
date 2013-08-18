from __future__ import division
import clib
import time
import numpy as np
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from demag import Demag
from save_vtk import SaveVTK
from materials import Nickel
from materials import UnitMaterial


#from show_vector import VisualSpin


class Sim(object):
    def __init__(self,mesh,T=0,mat=UnitMaterial(),name='unnamed'):
        self.t=0
        self.mesh = mesh
        self.nxyz = mesh.nxyz
        self.unit_length=mesh.unit_length
        self._T = np.zeros(self.nxyz)
        self._alpha = np.zeros(self.nxyz)
        self.spin = np.ones(3*self.nxyz)
        self.field = np.zeros(3*self.nxyz)
        self.dm_dt = np.zeros(3*self.nxyz)
        self.interactions=[]
        self.mat=mat
        self.pin_fun=None
        self.ode_times=0

        self._T[:] = T
        self._alpha[:]=self.mat.alpha
        self.set_options()

        self.vtk=SaveVTK(self.mesh,self.spin,name=name)


    def set_options(self,rtol=1e-7,atol=1e-7,dt=1e-15):

        if self.T.any()>0:
            self.vode=clib.RK2S(self.mat.mu_s,dt,
                        self.nxyz,
                        self.mat.gamma,
                        self.alpha,
                        self.spin,
                        self.field,
                        self.T,
                        self.update_effective_field)
        else:
            self.vode=clib.CvodeLLG(self.mat.mu_s,self.nxyz,
                                    self.mat.gamma,
                                    self.alpha,
                                    self.spin,
                                    self.field,
                                    self.update_effective_field)


    def set_m(self,m0=(1,0,0),normalise=True):

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

        if normalise:
            self.normalise()

        self.vode.set_initial_value(self.spin, self.t)

    def get_T(self):
        return self._T

    def set_T(self,T0):
        if  hasattr(T0, '__call__'):
            T = np.array([T0(p) for p in self.mesh.pos])
            self.T[:]=T[:]
        else:
            self.T[:] = T0

        self.set_options()
    T = property(get_T, set_T)

    def get_alpha(self):
        return self._alpha

    def set_alpha(self,alpha0):
        if  hasattr(alpha0, '__call__'):
            alpha = np.array([alpha0(p) for p in self.mesh.pos])
            self._alpha[:] = alpha[:]
        else:
            self._alpha[:] = alpha0

        self.set_options()
    alpha = property(get_alpha, set_alpha)

    def add(self,interaction):
        interaction.setup(self.mesh,self.spin,
                          unit_length=self.unit_length,
                          mu_s=self.mat.mu_s)
        self.interactions.append(interaction)

    def run_until(self,t):
        
        if abs(t)<1e-15:
            return
        
        ode=self.vode
        
        ode.run_until(t)
        
        self.spin[:]=ode.y[:]

        self.t = t


    def update_effective_field(self,y):
        
        self.spin[:]=y[:]
        self.field[:]=0

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

    def normalise(self):
        a=self.spin
        a.shape=(3,-1)
        a/=np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
        a.shape=(3*self.nxyz)
        #print a

    #only used for tests
    def spin_at(self,i,j,k):
        nxyz=self.mesh.nxyz
        nxy=self.mesh.nxy
        nx=self.mesh.nx

        i1=k*nxy+j*nx+i

        i2=i1+nxyz
        i3=i2+nxyz

        #print self.spin.shape,nxy,nx,i1,i2,i3
        return (self.spin[i1],
                self.spin[i2],
                self.spin[i3])


    def save_vtk(self):
        self.vtk.save_vtk(self.spin)



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
