import numpy as np


class UnitMaterial(object):
    def __init__(self):
        self.a=1
        self.b=1
        self.c=1
        self.J=1
        self.D=0.5
        self.mu_s=1
        self.gamma=1
        self.alpha=0.01
        self.unit_length=1


        
class Nickel(object):
    def __init__(self):
        self.mu_0=4*np.pi*1e-7
        self.mu_B=9.27400949e-24
        self.k_B=1.3806505e-23
        self.mu_s=0.617*self.mu_B
        self.S=1.0/2
        self.unit_length=1e-10
        self.a=3.524e-10
        self.b=self.a
        self.c=self.a
        self.K=4.8e3
        self.D=self.K*self.a**3
        self.Tc=630
        self.J=self.k_B*self.Tc/3.0
        self.gamma=2.210173e5/self.mu_0
        
if __name__=='__main__':
    ni=Nickel()
    print ni.K
    print 'gamma',ni.gamma
    print 'mu_s',ni.mu_s
    print 'J',ni.J
    print 'D',ni.D
    print 'D/J:',ni.D/ni.J