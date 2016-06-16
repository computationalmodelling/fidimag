import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.atomistic.mc import MonteCarlo
from fidimag.common import DataReader, CuboidMesh


def init_m(pos):
    x,y,z = pos

    return (0,0,1)


def run(mesh):
    
    mc = MonteCarlo(mesh, name='test1')
    mc.set_m(init_m)
    mc.set_options(H=[0,0,0.5], J=50.0, D=0.27*50, T=5.0)
    mc.run(steps=200000, save_m_steps=None, save_vtk_steps=1000, save_data_steps=10)

if __name__=='__main__':

    #mesh = CuboidMesh(nx=174,ny=150, nz=1, pbc='2d')
    mesh = CuboidMesh(nx=100, ny=100, nz=1, periodicity = (True, True, False))

    run(mesh)