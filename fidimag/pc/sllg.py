from __future__ import division

import fidimag.extensions.clib as clib
import numpy as np

from llg import LLG
import fidimag.util.helper as helper

from fidimag.pc.constant import Constant
const = Constant()

class SLLG(LLG):
    
    def __init__(self, mesh, name='unnamed'):
        """Simulation object.

        *Arguments*

          name : the Simulation name (used for writing data files, for examples)

        """
        super(SLLG, self).__init__(mesh, name=name)

        self._T = np.zeros(self.nxyz,dtype=np.float)
        
        self.set_options()
            

    def get_T(self):
        return self._T

    def set_T(self,T0):
        self._T[:] = helper.init_scalar(T0, self.mesh)
        
    T = property(get_T, set_T)
            
    def set_options(self, dt=1e-15, theta=1.0, gamma=const.gamma, k_B=const.k_B):

        self.gamma = gamma
        self.k_B = k_B
        
        self.vode=clib.RK2S(dt,
                        self.nxyz,
                        self.gamma,
                        self.k_B,
                        theta,
                        self._mu_s_inv,
                        self.alpha,
                        self.spin,
                        self.field,
                        self._T,
                        self._pins,
                        self.update_effective_field)
        
    
    def update_effective_field(self, y, t):
        
        self.spin[:]=y[:]
        
        self.field[:]=0
                
        for obj in self.interactions:
            self.field += obj.compute_field(t)
