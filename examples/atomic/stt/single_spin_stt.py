import numpy as np
from fidimag.atomistic import Sim, DMI, UniformExchange, Zeeman, Anisotropy
from fidimag.common import CuboidMesh, DataReader
import matplotlib.pyplot as plt

def dynamic(mesh):

    sim = Sim(mesh, name='dyn_spin', driver='slonczewski')
    # sim.set_options(rtol=1e-10,atol=1e-14)
    sim.driver.gamma = 1.0
    sim.mu_s = 1.0

    sim.set_m((0.8,0,-1))

    Kx = Anisotropy(Ku=-0.05, axis=(0, 0, 1), name='Kz')
    sim.add(Kx)

    sim.p = (0,0,1)

    sim.u0 = 0.005
    sim.driver.alpha = 0.1

    ts = np.linspace(0, 1200, 401)
    for t in ts:
        sim.run_until(t)
        #sim.save_vtk()
        print t

def plot_all():

    font = {'family': 'serif',
            'weight': 'normal',
            'size': 12,
            }

    plt.rc('font', **font)

    data = DataReader('dyn_spin.txt')
    ts = data['time']

    fig = plt.figure(figsize=(5, 4))

    mx = data['m_x']
    my = data['m_y']
    mz = data['m_z']
    plt.plot(ts, mx, '-', label='mx', color='b')
    plt.plot(ts, my, '-', label='my', color='g')
    plt.plot(ts[::6], mz[::6],'.-', label='mz',  color='r')

    plt.legend(bbox_to_anchor=[0.8, 0.8], shadow=True, frameon=True)
    #plt.xlim([0, 1.01])

    plt.legend()

    plt.xlabel('Time')
    plt.ylabel('m')
    plt.tight_layout()
    fig.savefig('m_ts.pdf')

if __name__ == '__main__':
    mesh = CuboidMesh(nx=1, ny=1, nz=1)
    dynamic(mesh)
    plot_all()
