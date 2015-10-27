import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fidimag.atomistic import Sim, DMI, UniformExchange, Anisotropy
from fidimag.common import Constant, CuboidMesh


const = Constant()

def m_init_dw(pos):

    x = pos[0]

    if x < 140:
        return (0, 0, 1)
    elif x > 160:
        return (0, 0, -1)
    else:
        return (0, 1, 1)


def analytical(xs, A=1.3e-11, D=4e-4, K=8e4):

    delta = np.sqrt(A / (K - D * D / (4 * A))) 

    phi = D / (2 * A) * xs 

    mx = - np.tanh(xs / delta)
    my = 1.0 / np.cosh(xs / delta) * np.cos(phi)
    mz = 1.0 / np.cosh(xs / delta) * np.sin(phi)
    return mx, my, mz

def relax_system(mesh):

    sim = Sim(mesh, name='relax')
    sim.set_default_options(gamma=const.gamma)
    sim.alpha = 0.5
    sim.mu_s = const.mu_s_1
    sim.do_procession = False

    sim.set_m(m_init_dw)

    J = 50.0 * const.k_B
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.1 * J
    dmi = DMI(D, dmi_type = 'interfacial')
    sim.add(dmi)

    K = 0.02 * J
    anis = Anisotropy(K, axis=[0,0,1])
    sim.add(anis)

    ONE_DEGREE_PER_NS = 17453292.52

    sim.relax(dt=1e-13, stopping_dmdt=0.01 * ONE_DEGREE_PER_NS,
              max_steps=1000, save_m_steps=100, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)

    xs = np.array([p[0] for p in mesh.pos]) - 150

    mx, my, mz = analytical(xs, A=J/2.0, D=-D, K=K)
    mxyz = sim.spin.copy()
    mxyz.shape = (3, -1)

    save_plot(xs, mxyz, mx, my, mz)


def save_plot(xs, mxyz, mx, my, mz):
    fig = plt.figure()
    mxyz.shape = (3, -1)
    plt.plot(xs, mxyz[0], '.', label='mx')
    plt.plot(xs, mxyz[1], '.', label='my')
    plt.plot(xs, mxyz[2], '.', label='mz')

    plt.plot(xs, mx, '-')
    plt.plot(xs, my, '-')
    plt.plot(xs, mz, '-')

    plt.ylabel('mxyz')
    plt.xlabel('x/a')
    plt.legend(loc=1)
    plt.grid()
    #plt.xlim([-150, 150])
    plt.ylim([-1.2, 1.2])
    fig.savefig('dw_dmi.pdf')


if __name__ == '__main__':
    mesh = CuboidMesh(nx=300, ny=40, nz=1)
    relax_system(mesh)
    #excite_system(mesh)
