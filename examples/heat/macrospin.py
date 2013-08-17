import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from pc import Anisotropy
from pc import FDMesh
from pc import Sim
from pc import Zeeman
import numpy as np


class Material():
    def __init__(self):
        self.mu_0=4*np.pi*1e-7
        self.mu_s=1.12e-17
        self.K=3.36e-18
        self.unit_length=1e-10
        self.gamma=2.210173e5/self.mu_0
        self.alpha=0.5

def single_spin(alpha=0.01):
    
    mat = Material()
    
    mesh=FDMesh(nx=1,ny=1,nz=1)
    
    mat.alpha=alpha
    sim = Sim(mesh,T=18000,mat=mat)
    sim.set_m((1, 1, 1))
    
    #sim.add(Zeeman(1,(0, 0, 1)))
    
    anis=Anisotropy(mat.K, direction=(0,0,1))
    sim.add(anis)
    
    dt = 0.5e-12; 
    ts = np.linspace(0, 1000 * dt, 1001)
    
    sx = []
    sy = []
    for t in ts:
        sim.run_until(t)
        sx.append(sim.spin[0])
        sy.append(sim.spin[1])
        print t
    
    
    
    plt.plot(sx, sy) 
    plt.xlabel("$S_x$")
    plt.ylabel("$S_y$")
    plt.grid()
    plt.axis((-0.9 , 0.9 , -0.9 , 0.9))
    plt.axes().set_aspect('equal')
    
    plt.savefig("macrospin.pdf")
        
   
if __name__=='__main__':
    single_spin()
