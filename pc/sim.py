from __future__ import division
import clib
import cvode
import time
import numpy as np
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from demag import Demag
from fileio import DataSaver, DataReader
from save_vtk import SaveVTK
from materials import Nickel
from materials import UnitMaterial
from constant import Constant

import pccp.util.helper as helper

const = Constant()
#from show_vector import VisualSpin


class Sim(object):
    
    def __init__(self,mesh,T=0,name='unnamed',driver='llg',pbc=None):
        """Simulation object.

        *Arguments*

          mesh : a dolfin mesh

          name : the Simulation name (used for writing data files, for examples)

          pbc : Periodic boundary type: None, '1d' or '2d'

          driver : 'llg', 'sllg' or 'llg_stt'

        """
        
        self.t=0
        self.mesh = mesh
        self.nxyz = mesh.nxyz
        self.unit_length=mesh.unit_length
        self._T = np.zeros(self.nxyz,dtype=np.float)
        self._alpha = np.zeros(self.nxyz,dtype=np.float)
        self._mu_s = np.zeros(3*self.nxyz,dtype=np.float)
        self.mu_s_inv = np.zeros(3*self.nxyz,dtype=np.float)
        self.spin = np.ones(3*self.nxyz,dtype=np.float)
        self.field = np.zeros(3*self.nxyz,dtype=np.float)
        self.dm_dt = np.zeros(3*self.nxyz,dtype=np.float)
        self.interactions=[]
        self.pin_fun=None
        self.driver = driver
        self.pbc = pbc

        self.saver = DataSaver(self, name+'.txt')
        self.saver.entities['E_total'] = {
            'unit': '<J>',
            'get': lambda sim : sim.compute_energy(),
            'header': 'E_total'}
        self.saver.update_entity_order()

        self._T[:] = T
        self.set_options()

        self.vtk=SaveVTK(self.mesh,self.spin,name=name)


    def set_options(self,rtol=1e-8,atol=1e-14,dt=1e-15):
        
        self._alpha[:] = 0.1
        self.gamma = 1
        self._mu_s[:] = 1
        self.mu_s_inv[:] = 1

        if self.driver == 'sllg':
            self.vode=clib.RK2S(self.mu_s,dt,
                        self.nxyz,
                        self.gamma,
                        self.alpha,
                        self.spin,
                        self.field,
                        self.T,
                        self.update_effective_field)
            
        elif self.driver == 'llg':
            self.vode=cvode.CvodeSolver(self.spin,
                                    rtol,atol,
                                    self.sundials_rhs)

        elif self.driver == 'llg_stt':
            self.vode=cvode.CvodeSolver(self.spin,
                                        rtol,atol,
                                        self.sundials_llg_stt_rhs)
            
            self.field_stt = np.zeros(3*self.nxyz)

            self.jx = 0
            self.jy = 0
            self.jz = 0
            self.p = 0.5
            self.beta = 0
            #a^3/dx ==> unit_length^2
            cell_size=self.mesh.cell_size
            self.u0 = const.g_e*const.mu_B*cell_size/(2*const.e*self.mu_s)*self.unit_length**2
            
        else:
            raise Exception("Unsppourted driver:{},avaiable drivers: sllg, llg, llg_s, llg_stt.".format(self.driver))
                    


    def set_m(self,m0=(1,0,0),normalise=True):
        
        self.spin[:]=helper.init_vector(m0,self.mesh, normalise)

        self.vode.set_initial_value(self.spin, self.t)

    def get_T(self):
        return self._T

    def set_T(self,T0):
        if  hasattr(T0, '__call__'):
            T = np.array([T0(p) for p in self.mesh.pos])
            self._T[:] = T[:]
        else:
            self._T[:] = T0

    T = property(get_T, set_T)
    
    def get_chi(self):
        return self._chi
    
    def set_chi(self,chi0):
        if  hasattr(chi0, '__call__'):
            chi = np.array([chi0(p) for p in self.mesh.pos])
            self._chi[:] = chi[:]
        else:
            self._chi[:] = chi0
        
    chi = property(get_chi, set_chi)

    def get_alpha(self):
        return self._alpha

    def set_alpha(self,value):
        self._alpha[:] = helper.init_scalar(value, self.mesh)
        
    alpha = property(get_alpha, set_alpha)
    
    def get_mu_s(self):
        return self._mu_s

    def set_mu_s(self, value):
        self._mu_s[:] = helper.init_scalar(value, self.mesh)
        for i in range(self.nxyz):
            if self._mu_s[i] == 0.0:
                self.mu_s_inv[i] = 0
            else: 
                self.mu_s_inv[i]=1.0/self._mu_s[i]
        
    mu_s = property(get_mu_s, set_mu_s)

    def add(self,interaction):
        interaction.setup(self.mesh,self.spin,
                          mu_s_inv=self.mu_s_inv,
                          pbc=self.pbc)
        
        #TODO: FIX
        for i in self.interactions:
            if i.name == interaction.name:
                interaction.name=i.name+'_2'
        
        self.interactions.append(interaction)
        
        energy_name = 'E_{}'.format(interaction.name)
        self.saver.entities[energy_name] = {
            'unit': '<J>',
            'get': lambda sim:sim.get_interaction(interaction.name).compute_energy(),
            'header': energy_name}
        
        self.saver.update_entity_order()
        
    def get_interaction(self, name):
        for interaction in self.interactions:
            if interaction.name == name:
                return interaction
        else: 
            raise ValueError("Failed to find the interaction with name '{}', "
                             "available interactions: {}.".format(
                        name, [x.name for x in self.interactions]))
        

    def run_until(self,t):
        
        if abs(t)<1e-15:
            return
        
        ode = self.vode
        
        ode.run_until(t)
        
        self.spin[:]=ode.y[:]

        self.t = t
        
        if self.pin_fun:
            self.pin_fun(self.t,self.mesh,self.spin)
            
        length = abs(self.spin_length()-1)
        if np.max(length)>1e-6:
            print length
            raise Exception("the error of spin length is large than 1e-7!!!! please check your code!!!")
            
        self.saver.save()

    def update_effective_field(self, y):
        
        self.spin[:]=y[:]
        
        self.field[:]=0
        
        if self.pin_fun:
            self.pin_fun(self.t,self.mesh,self.spin)
        
        for obj in self.interactions:
            self.field += obj.compute_field()

    def compute_effective_field(self, t):
        
        self.field[:]=0

        if self.pin_fun:
            self.pin_fun(t,self.mesh,self.spin)

        for obj in self.interactions:
            self.field += obj.compute_field(t)
    
    def sundials_rhs(self, t, y, ydot):
        
        self.t = t
        
        #already synchronized when call this funciton
        #self.spin[:]=y[:]
                
        self.compute_effective_field(t)
        
        clib.compute_llg_rhs(ydot,
                             self.spin,
                             self.field,
                             self.alpha,
                             self.gamma,
                             self.nxyz)
        
                
        #ydot[:] = self.dm_dt[:]
                
        return 0

    def sundials_llg_stt_rhs(self, t, y, ydot):
        
        self.t = t
        
        #already synchronized when call this funciton
        #self.spin[:]=y[:]
        
        self.compute_effective_field()
        
        clib.compute_stt_field(self.spin,
                               self.field_stt,
                               self.jx,
                               self.jy,
                               self.jz,
                               self.mesh.dx,
                               self.mesh.dy,
                               self.mesh.dz,
                               self.mesh.nx,
                               self.mesh.ny,
                               self.mesh.nz)        

        clib.compute_llg_stt_rhs(ydot,
                               self.spin,
                               self.field,
                               self.field_stt,
                               self.alpha,
                               self.beta,
                               self.u0*self.p,
                               self.gamma,
                               self.nxyz)
        
        #print 'field',self.field
        

    def compute_average(self):
        self.spin.shape=(3,-1)
        average=np.sum(self.spin,axis=1)/self.nxyz
        self.spin.shape=(3*self.nxyz)
        return average

    def compute_energy(self):
        
        energy=0
        
        for obj in self.interactions:
            energy+=obj.compute_energy()
            
        return energy

    def normalise(self):
        a=self.spin
        a.shape=(3,-1)
        a/=np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
        a.shape=(3*self.nxyz)
        #print a

    def spin_at(self,i,j,k):
        nxyz=self.mesh.nxyz

        i1=self.mesh.index(i,j,k)
        i2=i1+nxyz
        i3=i2+nxyz

        #print self.spin.shape,nxy,nx,i1,i2,i3
        return np.array([self.spin[i1],
                self.spin[i2],
                self.spin[i3]])
    
    def add_monitor_at(self, i,j,k, name='p'):
        """
        Save site spin with index (i,j,k) to txt file.
        """
    
        self.saver.entities[name] = {
            'unit': '<>',
            'get': lambda sim:sim.spin_at(i,j,k),
            'header': (name+'_x',name+'_y',name+'_z')}
        
        self.saver.update_entity_order()


    def save_vtk(self):
        self.vtk.save_vtk(self.spin)

    def stat(self):
        return self.vode.stat()

    def spin_length(self):
        self.spin.shape=(3,-1)
        length=np.sqrt(np.sum(self.spin**2,axis=0))
        self.spin.shape=(-1,)
        return length


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
