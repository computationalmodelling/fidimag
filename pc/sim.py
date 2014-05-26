from __future__ import division
import os
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
from constant import Constant

import pccp.util.helper as helper

const = Constant()

class Sim(object):
    
    def __init__(self,mesh,T=0,name='unnamed',driver='llg',pbc=None):
        """Simulation object.

        *Arguments*

          mesh : a dolfin mesh

          name : the Simulation name (used for writing data files, for examples)

          pbc : Periodic boundary type: None, '1d' or '2d'

          driver : 'llg', 'sllg' or 'llg_stt'

        """
        
        self.t = 0
        self.name = name
        self.mesh = mesh
        self.nxyz = mesh.nxyz
        self.nxyz_nonzero = mesh.nxyz
        self.unit_length=mesh.unit_length
        self._T = np.zeros(self.nxyz,dtype=np.float)
        self._alpha = np.zeros(self.nxyz,dtype=np.float)
        self._mu_s = np.zeros(self.nxyz,dtype=np.float)
        self.mu_s_inv = np.zeros(3*self.nxyz,dtype=np.float)
        self.spin = np.ones(3*self.nxyz,dtype=np.float)
        self.spin_last = np.ones(3*self.nxyz,dtype=np.float)
        self._pins = np.zeros(self.nxyz,dtype=np.int32)
        self.field = np.zeros(3*self.nxyz,dtype=np.float)
        self.dm_dt = np.zeros(3*self.nxyz,dtype=np.float)
        self._skx_number = np.zeros(self.nxyz,dtype=np.float)
        self.interactions=[]
        self.pin_fun = None
        self.driver = driver
        self.pbc = pbc
        
        self.step = 0
        
        self.saver = DataSaver(self, name+'.txt')
        
        self.saver.entities['E_total'] = {
            'unit': '<J>',
            'get': lambda sim : sim.compute_energy(),
            'header': 'E_total'}
        
        self.saver.entities['m_error'] = {
            'unit': '<>',
            'get': lambda sim : sim.compute_spin_error(),
            'header': 'm_error'}
        
        self.saver.entities['skx_num'] = {
            'unit': '<>',
            'get': lambda sim : sim.skyrmion_number(),
            'header': 'skx_num'}
        
        
        self.saver.update_entity_order()

        self._T[:] = T
        self.set_options()
        

        self.vtk=SaveVTK(self.mesh,self.spin,name=name)


    def set_options(self,rtol=1e-8,atol=1e-12, dt=1e-15, theta=1.0, gamma=const.gamma, k_B=const.k_B):
        """
        theta = 1.0 Heun method
        """
        
        self._alpha[:] = 0.1
        self._mu_s[:] = const.mu_s_1
        self.mu_s_inv[:] = 1
        self.gamma = gamma
        self.k_B = k_B
        self.do_procession = True

        if self.driver == 'sllg':
            self.vode=clib.RK2S(dt,
                        self.nxyz,
                        self.gamma,
                        self.k_B,
                        theta,
                        self.mu_s_inv,
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
            raise Exception("Unsppourted driver:{0},avaiable drivers: sllg, llg, llg_s, llg_stt.".format(self.driver))
                    


    def set_m(self,m0=(1,0,0),normalise=True):
        
        self.spin[:]=helper.init_vector(m0,self.mesh, normalise)
        
        #TODO: carefully checking and requires to call set_mu first
        for i in range(len(self.spin)):
            if self.mu_s_inv[i]==0:
                self.spin[i] = 0

        self.vode.set_initial_value(self.spin, self.t)

    def get_T(self):
        return self._T

    def set_T(self,T0):
        self._T[:]=helper.init_scalar(T0, self.mesh)

    T = property(get_T, set_T)
    
    def get_pins(self):
        return self._pins

    def set_pins(self,pin):
        self._pins[:]=helper.init_scalar(pin, self.mesh)
        
        for i in range(len(self._mu_s)):
            if self._mu_s[i] == 0.0:
                self._pins[i] = 1

    pins = property(get_pins, set_pins)
    

    def get_alpha(self):
        return self._alpha

    def set_alpha(self,value):
        self._alpha[:] = helper.init_scalar(value, self.mesh)
        
    alpha = property(get_alpha, set_alpha)
    
    def get_mu_s(self):
        return self._mu_s

    def set_mu_s(self, value):
        self._mu_s[:] = helper.init_scalar(value, self.mesh)
        self.mu_s_inv.shape=(3,-1)
        nonzero = 0 
        for i in range(self.nxyz):
            if self._mu_s[i] == 0.0:
                self.mu_s_inv[:,i] = 0
            else: 
                self.mu_s_inv[:,i]=1.0/self._mu_s[i]
                nonzero += 1
        
        self.nxyz_nonzero = nonzero
        self.mu_s_inv.shape=(-1,)
        
        for i in range(len(self._mu_s)):
            if self._mu_s[i] == 0.0:
                self._pins[i] = 1
        
    mu_s = property(get_mu_s, set_mu_s)

    def add(self,interaction, save_field=False):
        interaction.setup(self.mesh,self.spin,
                          mu_s_inv=self.mu_s_inv,
                          pbc=self.pbc)
        
        #TODO: FIX
        for i in self.interactions:
            if i.name == interaction.name:
                interaction.name=i.name+'_2'
        
        self.interactions.append(interaction)
        
        energy_name = 'E_{0}'.format(interaction.name)
        self.saver.entities[energy_name] = {
            'unit': '<J>',
            'get': lambda sim:sim.get_interaction(interaction.name).compute_energy(),
            'header': energy_name}
        
        if save_field:
            fn = '{0}'.format(interaction.name)
            self.saver.entities[fn] = {
            'unit': '<>',
            'get': lambda sim:sim.get_interaction(interaction.name).average_field(),
            'header': ('%s_x'%fn, '%s_y'%fn, '%s_z'%fn)}
        
        self.saver.update_entity_order()
        
    def get_interaction(self, name):
        for interaction in self.interactions:
            if interaction.name == name:
                return interaction
        else: 
            raise ValueError("Failed to find the interaction with name '{0}', "
                             "available interactions: {1}.".format(
                        name, [x.name for x in self.interactions]))
        

    def run_until(self,t):
        
        if t <= self.t:
            if t == self.t and self.t==0.0:
                self.compute_effective_field(t)
                self.saver.save()
            return
        
        ode = self.vode
        
        self.spin_last[:] = self.spin[:]
        
        flag = ode.run_until(t)
        
        if flag<0:
            raise Exception("Run cython run_until failed!!!")
        
        self.spin[:]=ode.y[:]

        self.t = t
        self.step += 1
        
        #update field before saving data
        self.compute_effective_field(t)
        self.saver.save()

    def update_effective_field(self, y, t):
        
        self.spin[:]=y[:]
        
        self.field[:]=0
        
        if self.pin_fun:
            self.pin_fun(self.t,self.mesh,self.spin)
        
        for obj in self.interactions:
            self.field += obj.compute_field(t)

    def compute_effective_field(self, t):
        
        self.field[:]=0

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
                             self._pins,
                             self.gamma,
                             self.nxyz,
                             self.do_procession)
        
                
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
        average=np.sum(self.spin,axis=1)/self.nxyz_nonzero
        self.spin.shape=(3*self.nxyz)
        return average

    def compute_energy(self):
        
        energy=0
        
        for obj in self.interactions:
            energy+=obj.compute_energy()
            
        return energy

    def skyrmion_number(self):
        nx = self.mesh.nx
        ny = self.mesh.ny
        nz = self.mesh.nz
        number = clib.compute_skymrion_number(self.spin,self._skx_number,nx,ny,nz)
        return number

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
        self.vtk.save_vtk(self.spin, step=self.step)
    
    def save_m(self):
        if not os.path.exists('%s_npys'%self.name):
            os.makedirs('%s_npys'%self.name)
        name = '%s_npys/m_%g.npy'%(self.name,self.step)
        np.save(name,self.spin)

    def stat(self):
        return self.vode.stat()

    def spin_length(self):
        self.spin.shape=(3,-1)
        length=np.sqrt(np.sum(self.spin**2,axis=0))
        self.spin.shape=(-1,)
        return length
    
    def compute_spin_error(self):
        length = self.spin_length()-1.0
        length[self._pins>0] = 0 
        return np.max(abs(length))
    
    def compute_dmdt(self, dt):
        m0 = self.spin_last
        m1 = self.spin
        dm = (m1 - m0).reshape((3, -1))
        max_dm = np.max(np.sqrt(np.sum(dm**2, axis=0))) 
        max_dmdt = max_dm / dt
        return max_dmdt
    
    def relax(self, dt=1.0, stopping_dmdt=0.01, max_steps=1000, save_m_steps=100, save_vtk_steps=100):
                 
        for i in range(0,max_steps+1):
            
            cvode_dt = self.vode.get_current_step()
            
            increment_dt = dt
            
            if cvode_dt > dt:
                increment_dt = cvode_dt

            self.run_until(self.t+increment_dt)
                
            if i%save_vtk_steps==0:
                self.save_vtk()
                
            if i%save_m_steps==0:
                self.save_m()
            
            dmdt = self.compute_dmdt(increment_dt)
            
            print 'step=%d, time=%g, max_dmdt=%g ode_step=%g'%(self.step, self.t, dmdt, cvode_dt)
            
            if dmdt<stopping_dmdt:
                break
        
        
        self.save_vtk()
        self.save_m()
    


if __name__=='__main__':
    pass