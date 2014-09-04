from llg import LLG

KNOWN_DRIVERS = {'llg': LLG}

def Sim(*args, **kwargs):
    
    driver = 'llg'
    
    if kwargs.has_key('driver'):
        driver = kwargs['driver']
        kwargs.pop('driver')
    
    if driver not in KNOWN_DRIVERS:
        raise NotImplementedError("""Driver '{}' is not implemented.
                                  Valid choices: one of '{}'.""".format(driver, KNOWN_DRIVERS.keys()))
    
    return KNOWN_DRIVERS[driver](*args, **kwargs)





"""
class Sim2(object):
    
    

    def set_options(self,rtol=1e-8,atol=1e-12, dt=1e-15, theta=1.0, gamma=const.gamma, k_B=const.k_B):

        self.default_c = -1.0 
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
                        self._pins,
                        self.update_effective_field)
            

        elif self.driver == 'llg_stt':
            self.vode=cvode.CvodeSolver(self.spin,
                                        rtol,atol,
                                        self.sundials_llg_stt_rhs)
            
            self.field_stt = np.zeros(3*self.nxyz)

            self.jx = np.zeros(self.nxyz,dtype=np.float)
            self.jy = np.zeros(self.nxyz,dtype=np.float)
            
            self.jx0 = np.zeros(self.nxyz,dtype=np.float)
            self.jy0 = np.zeros(self.nxyz,dtype=np.float)
            
            self.p = 0.5
            self.beta = 0
            #a^3/dx ==> unit_length^2
            cell_size=self.mesh.dx*self.mesh.dy*self.mesh.dz
            #FIXME: change the u0 to spatial 
            self.u0 = const.g_e*const.mu_B*cell_size/(2*const.c_e*self.mu_s[0])*self.unit_length**2
            
        else:
            raise Exception("Unsppourted driver:{0},avaiable drivers: sllg, llg, llg_s, llg_stt.".format(self.driver))
                    


    

    def get_T(self):
        return self._T

    def set_T(self,T0):
        self._T[:]=helper.init_scalar(T0, self.mesh)

    T = property(get_T, set_T)
    
        


    def sundials_llg_stt_rhs(self, t, y, ydot):
        
        self.t = t
        
        #already synchronized when call this funciton
        #self.spin[:]=y[:]
        
        self.compute_effective_field(t)
        
        if self.update_j_fun is not None:
            self.jx[:] = self.jx0[:]*self.update_j_fun(t)
            self.jy[:] = self.jy0[:]*self.update_j_fun(t)
            self.jz[:] = self.jz0[:]*self.update_j_fun(t)
        
        clib.compute_stt_field(self.spin,
                               self.field_stt,
                               self.jx,
                               self.jy,
                               self.mesh.dx,
                               self.mesh.dy,
                               self.mesh.nx,
                               self.mesh.ny,
                               self.mesh.nz,
                               self.xperiodic,
                               self.yperiodic)        

        clib.compute_llg_stt_rhs(ydot,
                               self.spin,
                               self.field,
                               self.field_stt,
                               self.alpha,
                               self.beta,
                               self.u0*self.p,
                               self.gamma,
                               self.nxyz)
        


"""
if __name__=='__main__':
    pass