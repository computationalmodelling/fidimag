import cvode
import baryakhtar_clib as clib 
from micro import LLG
from relax import Relaxation
from relax import Laplace

class LLBarFull(LLG):
    
    def __init__(self, mesh, chi=1e-3, name='unnamed'):

        self.chi = chi
        super(LLBarFull, self).__init__(mesh, name=name)
        self.lap = Laplace(mesh)
        self.add(Relaxation(chi))
    
    def set_options(self, rtol=1e-8, atol=1e-12, gamma=2.21e5, Ms=8.0e5, alpha=0.1, beta=0.0):
        self._alpha[:] = alpha
        self._Ms[:] = Ms
        self.gamma = gamma
        self.beta =  beta
        self.do_procession = True
        
        self.vode = cvode.CvodeSolver(self.spin,
                                    rtol,atol,
                                    self.sundials_rhs)
    
    
    def sundials_rhs(self, t, y, ydot):
        
        self.t = t
        
        #already synchronized when call this funciton
        #self.spin[:]=y[:]
        
        self.compute_effective_field(t)
        delta_h = self.lap.compute_laplace_field(self.field)
        
        clib.compute_llg_rhs_baryakhtar(ydot,
                                        self.spin,
                                        self.field,
                                        delta_h,
                                        self.alpha,
                                        self.beta,
                                        self._pins,
                                        self.gamma,
                                        self.nxyz,
                                        self.do_procession)
                             
                             
        #ydot[:] = self.dm_dt[:]
                             
        return 0


class LLBar(LLG):
    
    def __init__(self, mesh, name='unnamed'):
        
        super(LLBarFull, self).__init__(mesh, name=name)
        self.lap = Laplace(mesh)
    
        self.field_perp = np.zeros(3*self.nxyz,dtype=np.float)
    
    def set_options(self, rtol=1e-8, atol=1e-12, gamma=2.21e5, Ms=8.0e5, alpha=0.1, beta=0.0):
        self.default_c = 1e11
        self._alpha[:] = alpha
        self._Ms[:] = Ms
        self.gamma = gamma
        self.beta =  beta
        self.do_procession = True
        
        self.vode = cvode.CvodeSolver(self.spin,
                                      rtol,atol,
                                      self.sundials_rhs)
    
    
    def sundials_rhs(self, t, y, ydot):
        
        self.t = t
        
        #already synchronized when call this funciton
        #self.spin[:]=y[:]
        
        self.compute_effective_field(t)
        clib.compute_perp_field(self.m, self.field, self.field_perp, self.nxyz)
        delta_h = self.lap.compute_laplace_field(self.field_perp)
        
        clib.compute_llg_rhs_baryakhtar_reduced(ydot,
                                        self.spin,
                                        self.field,
                                        delta_h,
                                        self.alpha,
                                        self.beta,
                                        self._pins,
                                        self.gamma,
                                        self.nxyz,
                                        self.do_procession,
                                        self.default_c)
                                        
                                        
                                        #ydot[:] = self.dm_dt[:]
                                        
        return 0

if __name__=='__main__':
    pass