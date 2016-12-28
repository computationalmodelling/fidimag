import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, DMI, Zeeman

mu0 = 4 * np.pi * 1e-7

Rx = 52
Ry = 90

def init_m(pos):
    x,y,z = pos
    
    r = 12
    
    x = x%Rx 
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

def init_m2(pos):

    x, y = pos[0] - 0, pos[1] - 0

    if x**2 + y**2 > 20**2:
        return (0, 0, 1)
    else:
        return (0, 0, -1)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')

    sim.driver.set_tols(rtol=1e-6, atol=1e-6)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_precession = False

    sim.set_m(init_m)

    exch = UniformExchange(A=1.3e-11)
    sim.add(exch)

    dmi = DMI(D=-4e-3)
    sim.add(dmi)

    zeeman = Zeeman((0, 0, 4e5))
    sim.add(zeeman, save_field=True)

    sim.relax(dt=1e-13, stopping_dmdt=1e-2,
              save_m_steps=None, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)

def plot_mz(nx=Rx*2,ny=Ry):
    m = np.load('m0.npy')
    m.shape = (-1,3)
    mz = m[:,2]
    mz.shape =(ny, nx)

    plt.imshow(mz)
    plt.savefig('skx_mz.png')


def do_fft(nx=Rx*2,ny=Ry, dx=1, dy=1):
    m = np.load('m0.npy')
    m.shape = (-1,3)
    mz = m[:,2]
    mz.shape =(ny, nx)
    mz = np.transpose(mz)
    
    mzr = mz - np.average(mz)
    Mz = np.fft.fft2(mzr)
    Mz = np.fft.fftshift(Mz)

    qx = np.fft.fftfreq(nx,dx)*2*np.pi
    qy = np.fft.fftfreq(ny,dy)*2*np.pi

    qx = np.fft.fftshift(qx)
    qy = np.fft.fftshift(qy)

    
    plt.figure(figsize=(5,5))
    #plt.imshow(np.abs(Mz)**2, aspect=1, interpolation='none', origin='lower', extent=(qx.min(), qx.max(), qy.min(), qy.max()))
    plt.imshow(np.abs(Mz)**2, aspect=1, origin='lower', extent=(qx.min(), qx.max(), qy.min(), qy.max()))
    #plt.axis('equal')
    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.5, 0.5)
    plt.xlabel('qx (1/nm)')
    plt.ylabel('qy (1/nm)')
    plt.tight_layout()
    plt.savefig('skx.png')




if __name__ == '__main__':

    mesh = CuboidMesh(nx=Rx*2, ny=Ry, nz=1, dx=1.0, dy=1.0, dz=1.0, unit_length=1e-9, periodicity = (True, True, False))

    relax_system(mesh)

    plot_mz()

    do_fft()
