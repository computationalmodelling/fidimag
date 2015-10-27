import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from micro import Sim, UniformExchange, Demag, UniaxialAnisotropy
from common import CuboidMesh
from pc import NEB_Sundials
from util.helper import plot_energy_2d, plot_energy_3d
#from finmag.sim.neb import plot_energy_2d, plot_energy_3d


def init_dw(pos):
    x = pos[0]

    if x < 200:
        return (-1, 0, 0)
    elif x > 300:
        return (1, 0, 0)
    else:
        return (0, 1, 0)


def create_simulation(mesh):

    sim = Sim(mesh)
    sim.Ms = 8.6e5

    sim.set_m((1, 0, 0))
    sim.add(UniformExchange(A=1.3e-11))
    # sim.add(Demag())
    #sim.add(UniaxialAnisotropy(Kx, (1, 0, 0), name='Kx'))
    anis = UniaxialAnisotropy(1e5, axis=(1, 0, 0))
    sim.add(anis)

    return sim


def relax_system(sim):

    # init_images=[np.load('m_init.npy')]
    """
    n=20
    for i in range(1,n):
        theta = i*np.pi/n
        mx = -np.cos(theta)
        my = np.sin(theta)
        init_images.append((mx,my,0))
    """

    # init_images.append(np.load('m_final.npy'))

    init_images = [(-1, 0, 0), init_dw, (1, 0, 0)]

    neb = NEB_Sundials(
        sim, init_images, interpolations=[6, 6], name='neb', spring=1e8)

    # neb.add_noise(0.1)

    neb.relax(dt=1e-8, max_steps=5000, save_vtk_steps=500,
              save_npy_steps=500, stopping_dmdt=100)

if __name__ == "__main__":

    mesh = CuboidMesh(nx=160, ny=1, nz=1, dx=4.0, dy=4.0, dz=4.0, unit_length=1e-9)

    sim = create_simulation(mesh)
    # relax_two_state(mesh)
    relax_system(sim)
    plot_energy_2d('neb')
    plot_energy_3d('neb')
