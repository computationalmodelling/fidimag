import numpy as np
from pccp.pc import Sim
from pccp.pc import FDMesh
from pccp.pc import DMI
from pccp.pc import UniformExchange
from pccp.pc import Zeeman


def init_m(pos):
    x,y,z = pos
    
    x0,y0,r  = 166,96,25
    
    x1 = x%x0
    y1 = y%y0
    
    m1 = (0.05,0.01,-1)
    m2 = (0,0,1)
    
    if (x1-r)**2 + (y1-r)**2 < r**2:
        return m1
    elif (x1-x0/2.-r)**2 + (y1-y0/2.-r)**2 < r**2:
        return m1
    else:
        return m2

def random_m(pos):
    return np.random.random(3)-0.5


def relax_system(mesh):
    
    sim=Sim(mesh,name='relax',pbc='2d')
    #sim.set_options(rtol=1e-10,atol=1e-14)
    sim.alpha = 1.0
    sim.gamma = 1.0
    sim.mu_s = 1.0
    
    sim.set_m(init_m)
    #sim.set_m(random_m)
    #sim.set_m(np.load('m_10000.npy'))

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)
    
    D = 0.09
    dmi = DMI(D)
    sim.add(dmi)
    
    zeeman = Zeeman([0,0,3.75e-3])
    sim.add(zeeman)
    
    sim.relax(dt=2.0, stopping_dmdt=1e-6, max_steps=1000, save_m_steps=100, save_vtk_steps=50)
    
    np.save('m0.npy',sim.spin)

                        
if __name__=='__main__':
    
    np.random.seed(11)
    
    #mesh = FDMesh(nx=288,ny=288,nz=1)
    mesh = FDMesh(nx=166,ny=96*2,nz=1)
    
    relax_system(mesh)
    
    print 'relax system done'
    #spin_wave(mesh,m0)
    
