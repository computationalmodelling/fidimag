import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.atomistic import MonteCarlo
from fidimag.common import DataReader, CuboidMesh
import fidimag.common.constant as const

def random_m(pos):
    return np.random.random(3) - 0.5

def run(mesh):
    mc = MonteCarlo(mesh, name='test1')
    mc.set_m(random_m)
    J = 50*const.k_B
    mc.set_options(H=[0,0,1.0], J=J, D=0.27*J, T=4.0)
    mc.run(steps=100000, save_m_steps=None, save_vtk_steps=50000, save_data_steps=50000)

if __name__=='__main__':
    mesh = CuboidMesh(nx=58, ny=33, nz=1, periodicity = (True, True, False))
    run(mesh)
