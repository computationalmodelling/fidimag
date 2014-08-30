from macro import LLG

class LLBar(LLG):
    
    def set_options(self,rtol=1e-8,atol=1e-12, gamma=2.21e5, Ms=8.0e5, alpha=0.1):
        self.default_c = 1e11
        self._alpha[:] = alpha
        self._Ms[:] = Ms
        self.gamma = gamma
        self.do_procession = True
        
        self.vode=cvode.CvodeSolver(self.spin,
                                    rtol,atol,
                                    self.sundials_rhs)
    
    
    
    def set_m(self,m0=(1,0,0),normalise=False):
        
        self.spin[:]=helper.init_vector(m0,self.mesh, normalise)
        
        #TODO: carefully checking and requires to call set_mu first
        self.spin.shape=(3,-1)
        for i in range(len(self.spin)):
            if self.Ms[i]==0:
                self.spin[:,i] = 0
        self.spin.shape=(-1,)
        
        self.vode.set_initial_value(self.spin, self.t)
    
    
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
                             self.do_procession,
                             self.default_c)
                             
                             
                             #ydot[:] = self.dm_dt[:]
                             
        return 0


class LLBarFull(LLG):
    
    def __init__(self, mesh, chi=1, name='unnamed'):
        """Simulation object.

        *Arguments*

          mesh : a dolfin mesh

          name : the Simulation name (used for writing data files, for examples)

        """
        
        super()


    def set_options(self,rtol=1e-8,atol=1e-12, gamma=2.21e5, Ms=8.0e5, alpha=0.1):
        self.default_c = 1e11
        self._alpha[:] = alpha
        self._Ms[:] = Ms
        self.gamma = gamma
        self.do_procession = True
        
        self.vode=cvode.CvodeSolver(self.spin,
                                    rtol,atol,
                                    self.sundials_rhs)



    def set_m(self,m0=(1,0,0),normalise=True):
        
        self.spin[:]=helper.init_vector(m0,self.mesh, normalise)
        
        #TODO: carefully checking and requires to call set_mu first
        self.spin.shape=(3,-1)
        for i in range(len(self.spin)):
            if self.Ms[i]==0:
                self.spin[:,i] = 0
        self.spin.shape=(-1,)

        self.vode.set_initial_value(self.spin, self.t)


    def compute_effective_field(self, t):
        
        #self.spin[:] = y[:]
        
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
                             self.do_procession,
                             self.default_c)
        
                
        #ydot[:] = self.dm_dt[:]
                
        return 0
    
    def compute_average(self):
        self.spin.shape=(3,-1)
        average=np.sum(self.spin,axis=1)/self.nxyz_nonzero
        self.spin.shape=(3*self.nxyz)
        return average





if __name__=='__main__':
    pass