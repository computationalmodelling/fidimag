import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import fidimag.extensions.clib as clib
import numpy as np
from fidimag.atomistic import MonteCarlo, HexagonalMesh
import fidimag.common.constant as const

def test_random_sphere(do_plot=False):
    """
    test whether spins are uniformly distributed on the surface of a sphere.
    """
    mt19937 = clib.rng_mt19937()
    mt19937.set_seed(123)
    spin = np.zeros(3 * 10000, dtype=np.float)
    mt19937.fill_vector_uniform_sphere(spin, 10000)
    spin.shape = (-1,3)
    n = spin.shape[0]
    x = np.average([spin[i,0] for i in range(n) if spin[i,0]>0])
    y = np.average([spin[i,1] for i in range(n) if spin[i,1]>0])
    z = np.average([spin[i,2] for i in range(n) if spin[i,2]>0])
    print(x,y,z)
    print(abs(x-y), abs(x-z), abs(y-z))
    assert(abs(x-y)<4e-3)
    assert(abs(x-z)<4e-3)
    assert(abs(y-z)<4e-3)

    if do_plot:
        fig = plt.figure()
        ax3D = fig.add_subplot(111, projection='3d')
        ax3D.scatter(spin[:1000,0], spin[:1000,1], spin[:1000,2], s=30, marker='.')
        plt.savefig('test_random_sphere.png')

def random_m(pos):
    return np.random.random(3) - 0.5

def test_mc_run(Hz=6.0, T=5.0):
    #This fast test just shows whether mc can be run or not.
    np.random.seed(100)
    #mesh = CuboidMesh(nx=28*2, ny=16*3, nz=1, periodicity = (True, True, False))
    mesh = HexagonalMesh(1, 28, 16*2, periodicity=(True, True))
    mc = MonteCarlo(mesh, name='test_mc')
    mc.set_m(random_m)
    J = 50*const.k_B
    mc.set_options(H=[0,0,Hz], J=J, D=0.5*J, T=T)
    mc.run(steps=5000, save_m_steps=None, save_vtk_steps=None, save_data_steps=1000)
    skx_number = mc.skyrmion_number()
    assert(skx_number<-2.5)
    assert(skx_number>-3.5)
    #np.save('m.npy', mc.spin)
    #plot_m(mesh, 'm.npy', comp='z')

if __name__ == '__main__':
    test_random_sphere(do_plot=True)
    test_mc_run()
