import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.atomistic import Sim, DMI, UniformExchange, Zeeman, TimeZeeman
from fidimag.common import CuboidMesh


global_nx = 174
global_ny = 150
global_mesh = CuboidMesh(nx=global_nx,ny=global_ny, pbc='2d')

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

def sinc_fun(t):

    w = 0.1

    return np.sinc(w*t)

def relax_system(mesh, Hy=0):

    sim=Sim(mesh,name='relax')
    sim.set_options(rtol=1e-10,atol=1e-12)
    sim.alpha = 0.5
    sim.gamma = 1.0
    sim.mu_s = 1.0

    sim.do_precession = False

    sim.set_m(init_m)
    #sim.set_m(random_m)
    #sim.set_m(np.load('m_10000.npy'))

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.18
    dmi = DMI(D)
    sim.add(dmi)

    zeeman = Zeeman([0,Hy,2e-2],name='H')
    sim.add(zeeman)

    sim.relax(dt=2.0, stopping_dmdt=1e-8, max_steps=10000, save_m_steps=100, save_vtk_steps=50)

    np.save('m0.npy',sim.spin)

def excite_system(mesh, Hy=0):

    sim=Sim(mesh,name='dyn')

    sim.set_options(rtol=1e-10,atol=1e-12)
    sim.alpha = 0.04
    sim.gamma = 1.0
    sim.mu_s = 1.0

    sim.set_m(np.load('m0.npy'))

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.18
    dmi = DMI(D)
    sim.add(dmi)

    zeeman = Zeeman([0,Hy,2e-2],name='H')
    sim.add(zeeman)

    hx = TimeZeeman([0,0,1e-5], sinc_fun, name='h')
    sim.add(hx, save_field=True)


    dt = 5
    steps = 2001
    for i in range(steps):

        sim.run_until(i*dt)
        #sim.save_m()

        #print 'sim t=%g'%(i*dt)

if __name__=='__main__':

	relax_system(global_mesh)
	excite_system(global_mesh)
