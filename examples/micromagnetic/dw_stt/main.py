import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from micro import Sim
from common import CuboidMesh

# The energies, we can use DMI in a future simulation
from micro import UniformExchange
from micro import UniaxialAnisotropy
# from micro import DMI

mu0 = 4 * np.pi * 1e-7

# Initial State, a rough DW ina 1D chain


def init_m(pos):

    x = pos[0]

    if x < 400:
        return (1, 0, 0)
    elif 400 <= x < 500:
        return (0, 1, 1)
    else:
        return (-1, 0, 0)


def relax_system(mesh):

    # Only relaxation
    sim = Sim(mesh, name='relax')

    # Simulation parameters
    sim.driver.set_tols(rtol=1e-8, atol=1e-10)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_precession = False

    # The initial state passed as a function
    sim.set_m(init_m)
    # sim.set_m(np.load('m0.npy'))

    # Energies
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    anis = UniaxialAnisotropy(5e4)
    sim.add(anis)

    # dmi = DMI(D=8e-4)
    # sim.add(dmi)

    # mT = 795.7747154594767
    # ONE_DEGREE_PER_NS = 17453292.52

    # Start relaxation and save the state in m0.npy
    sim.relax(dt=1e-14, stopping_dmdt=0.00001, max_steps=5000,
              save_m_steps=None, save_vtk_steps=None)

    np.save('m0.npy', sim.spin)


def deal_plot(_list, output_file):
    """
    The list contain any number of list with two entries:
    the _file path and the component of the magnetisation
    in the formats:

            mx, my, mz

    e.g.
            _list = [ ['m0.npy', 'mz'], ... ]

    output_file     :: output file name

    """

    plt.figure()

    comp = {'mx': 0, 'my': 1, 'mz': 2}

    for element in _list:
        # element[0]  --> _file
        # element[1]  --> m component: mx, my or mz

        data = np.load(element[0])

        data.shape = (3, -1)

        mx = data[comp[element[1]]]

        plt.plot(mx, '--',
                 label=element[1] + ' / %s' % element[0],
                 dashes=(2, 2))

    plt.legend()
    # plt.xlim([0, 0.012])
    # plt.ylim([-5, 100])
    plt.xlabel(r'x  [ nm ]')
    plt.savefig(output_file)


# THIS NEEDS REVISION
def deal_plot_dynamics(spin, m_component='mz',
                       output_file):
    """
    The dynamics of the m_component of the i-th spin
    of the chain, from the magnetisation snapshot
    m_{i} files, located at dyn_npys

    m_component     :: 'mx', 'my' or 'mz'

    """

    comp = {'mx': 0, 'my': 1, 'mz': 2}

    plt.figure()

    _all = []
    for i in range(500):
        data = np.load('dyn_npys/m_%d.npy' % i)
        data.shape = (3, -1)
        _all.append(data[comp[m_component]][spin])

    plt.plot(_all, '--', label=m_component)

    plt.xlabel(r't')
    # plt.ylabel('Susceptibility')
    plt.savefig(output_file)


# This function excites the system with a
# current in the x-direction using Spin Transfer Torque
# formalism
def excite_system(mesh):

    # Specify the stt dynamics in the simulation
    sim = Sim(mesh, name='dyn', driver='llg_stt')

    sim.driver.set_tols(rtol=1e-12, atol=1e-14)
    sim.alpha = 0.05
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5

    # sim.set_m(init_m)
    sim.set_m(np.load('m0.npy'))

    # Energies
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    anis = UniaxialAnisotropy(5e4)
    sim.add(anis)

    # dmi = DMI(D=8e-4)
    # sim.add(dmi)

    # Set the current in the x direction, in A / m
    # beta is the parameter in the STT torque
    sim.jx = -1e12
    sim.beta = 1

    # The simulation will run for 5 ns and save
    # 500 snapshots of the system in the process
    ts = np.linspace(0, 5e-9, 501)

    for t in ts:
        print 'time', t
        sim.run_until(t)
        sim.save_vtk()
        sim.save_m()

if __name__ == '__main__':

    # We will crate a mesh with 1000 elements of elements
    # in the x direction, and 1 along y and z
    # (so we have a 1D system)
    mesh = CuboidMesh(nx=1000, ny=1, nz=1,
                  dx=2, dy=2, dz=2.0,
                  unit_length=1e-9)

    # Relax the initial state
    relax_system(mesh)
    deal_plot([['m0.npy', 'mx']], 'initial_state.pdf')

    # Now we excite the system with the current
    excite_system(mesh)

    # Plot the dynamics of a single spin
    # deal_plot_dynamics()

    # We can plot the m_z component for a number snapshots
    # to observe the DW motion
    # We will plot the 200th and 499th file (the last one
    # is the system at ~5 ns)
    deal_plot([['m0.npy', 'mx'],
               ['dyn_npys/m_200.npy', 'mx'],
               ['dyn_npys/m_499.npy', 'mx']
               ],
              'dw_motion_stt.pdf')
