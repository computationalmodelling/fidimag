import matplotlib as mpl
mpl.use("Agg")

import numpy as np
from fidimag.atomistic import Sim, FDMesh, DMI, UniformExchange, Zeeman, TimeZeeman
from fidimag.common import DataReader

from fidimag.atomistic.eigen import EigenProblem

import matplotlib.pyplot as plt

def init_m(pos):
    x,y,z = pos
    
    Rx, Ry, r = 87, 50, 12
    
    x = x%Rx + 20 
    y = y%Ry 
    
    m1 = (0.05,0.01,-1)
    m2 = (0,0,1)
    
    if x**2 + y**2 < r**2:
        return m1
    elif x**2 + (y-Ry)**2 < r**2:
        return m1
    elif (x-Rx)**2 + y**2 < r**2:
        return m1
    elif (x-Rx)**2 + (y-Ry)**2 < r**2:
        return m1
    elif (x-Rx/2.0)**2 + (y-Ry/2.0)**2 < r**2:
        return m1
    else:
        return m2


def relax_system(mesh):
    
    sim=Sim(mesh,name='relax')
    sim.set_options(rtol=1e-12,atol=1e-14)
    sim.do_procession = False
    sim.alpha = 0.5
    sim.gamma = 1.0
    sim.mu_s = 1.0
    
    sim.set_m(init_m)

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)
    
    D = 0.18
    dmi = DMI(D)
    sim.add(dmi)
    
    zeeman = Zeeman([0,0e-3,2e-2],name='H')
    sim.add(zeeman)
    
    sim.relax(dt=2.0, stopping_dmdt=1e-8, max_steps=10000, save_m_steps=None, save_vtk_steps=100)
    
    np.save('m0.npy',sim.spin)


def eigen(mesh):
    m0 = np.load('m0.npy')
    ep = EigenProblem(mesh, m0, H=[0,0,2e-2], J=1.0, D=0.18)
    ep.solve_sparse()
    
              
if __name__=='__main__':
    
    #mesh = FDMesh(nx=174,ny=150, nz=1, pbc='2d')
    mesh = FDMesh(nx=87,ny=100, nz=1, pbc='2d')
    
    relax_system(mesh)
    
    eigen(mesh)
