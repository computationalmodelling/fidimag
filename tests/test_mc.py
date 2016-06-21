import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import fidimag.extensions.clib as clib
import numpy as np
from fidimag.atomistic import MonteCarlo


def my_average(spin):
    spin.shape = (-1,3)
    n = spin.shape[0]
    x1 = np.average([spin[i,0] for i in range(n) if spin[i,0]>0])
    y1 = np.average([spin[i,1] for i in range(n) if spin[i,1]>0])
    z1 = np.average([spin[i,2] for i in range(n) if spin[i,2]>0])
    #print(x1, y1, z1)
    return x1, y1, z1

def test_random_sphere(do_plot=False):
    mt19937 = clib.rng_mt19937()
    mt19937.set_seed(123)
    spin = np.zeros(3 * 10000, dtype=np.float)
    mt19937.fill_vector_uniform_sphere(spin, 10000)
    spin.shape = (-1,3)
    x,y,z = my_average(spin)
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
    

if __name__ == '__main__':
    test_random_sphere(do_plot=True)