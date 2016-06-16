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
    mc.set_options(H=[0,0,2e-2], J=50.0, D=0.18*50, T=10)
    mc.run(steps=100)


if __name__=='__main__':

    #mesh = CuboidMesh(nx=174,ny=150, nz=1, pbc='2d')
    mesh = CuboidMesh(nx=50, ny=60, nz=1, periodicity = (True, True, False))

    run(mesh)