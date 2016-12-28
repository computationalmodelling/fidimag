import numpy as np
from fidimag.atomistic import Sim, DMI, UniformExchange, Anisotropy
from fidimag.common import CuboidMesh
from fidimag.atomistic import NEB_Sundials

def mu_s(pos):
    x,y,z = pos
    x0,y0,r  = 60,60,60.5

    if (x-x0)**2 + (y-y0)**2 <= r**2:
        return 1
    else:
        return 0

def init_m(pos):
    x,y,z = pos

    x0,y0,r  = 60,60,25

    m1 = (0.05,0.01,-1)
    m2 = (0,0,1)

    if (x-x0)**2 + (y-y0)**2 < r**2:
        return m1
    else:
        return m2

def init_m_down(pos):
    x,y,z = pos

    x0,y0,r  = 60,60,25

    m1 = (0.05,0.01,1)
    m2 = (0,0,-1)

    if (x-x0)**2 + (y-y0)**2 < r**2:
        return m1
    else:
        return m2

def create_sim():

    mesh = CuboidMesh(nx=121,ny=121,nz=1)
    sim=Sim(mesh,name='relax')

    sim.driver.alpha = 1.0
    sim.driver.gamma = 0.5
    sim.mu_s = mu_s

    sim.set_m(init_m)

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.08
    dmi = DMI(D)
    sim.add(dmi)

    K = 4e-3
    anis=Anisotropy(K, direction=(0,0,1),name='Ku')
    sim.add(anis)

    return sim


def relax_system():

    sim = create_sim()

    sim.do_precession = False

    sim.relax(dt=2.0, stopping_dmdt=1e-6, max_steps=5000, save_m_steps=100, save_vtk_steps=100)

    np.save('m_up.npy',sim.spin)

    sim.set_m(init_m_down)

    sim.relax(dt=2.0, stopping_dmdt=1e-6, max_steps=5000, save_m_steps=100, save_vtk_steps=100)

    np.save('m_down.npy',sim.spin)

def run_neb():

    sim = create_sim()

    init_images = [np.load('m_up.npy'), np.load('m_down.npy')]

    neb = NEB_Sundials(sim, init_images, interpolations=[10], name='neb', spring=1.0)

    neb.relax(dt=1, max_steps=5000, save_vtk_steps=500, save_npy_steps=500, stopping_dmdt=1e-5)



if __name__=='__main__':

    #relax_system()
    run_neb()
