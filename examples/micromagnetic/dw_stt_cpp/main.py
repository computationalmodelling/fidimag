import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np
from fidimag.micro import Sim, UniformExchange, UniaxialAnisotropy
from fidimag.common import CuboidMesh
import fidimag.common.constant as const

from scipy.optimize import curve_fit


mu0 = 4 * np.pi * 1e-7
global_Ms = 8.6e5
global_A = 1.3e-11
global_K = 0.25*mu0*global_Ms**2

global_P = 0.5
global_d = 2e-9 #2nm
global_const = const.h_bar*const.gamma*global_P/(2*const.c_e*global_d*global_Ms)

def init_m(pos):

    x = pos[0]

    if x < 300:
        return (1, 0, 0)
    elif 300 <= x < 400:
        return (0, 1, 1)
    else:
        return (-1, 0, 0)


def relax_system(mesh):

    # Only relaxation
    sim = Sim(mesh, name='relax')

    # Simulation parameters
    sim.driver.set_tols(rtol=1e-8, atol=1e-10)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 2.211e5
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

    # Start relaxation and save the state in m0.npy
    sim.relax(dt=1e-14, stopping_dmdt=0.00001, max_steps=5000,
              save_m_steps=None, save_vtk_steps=None)

    np.save('m0.npy', sim.spin)



def model(xs, x0):
    delta = 16.1245154966
    return -np.tanh((xs-x0)/delta)

def extract_dw(m):

    m.shape=(-1,3)
    mx = m[:,0]
    my = m[:,1]
    mz = m[:,2]

    xs = np.linspace(0.5, 999.5, 1000)

    id = np.argmin(abs(mx))

    popt, pcov = curve_fit(model, xs, mx, p0=[id*1.0])
    print(popt[0], id)

    theta = np.arctan2(mz[id], my[id])
    return popt[0], theta

def sinc_fun(t):
    return np.sinc(1e9*t)  # sinc(x) = sin(\pi x)/(\pi x)

# This function excites the system with a
# current in the x-direction using Spin Transfer Torque
# formalism
def excite_system(mesh, beta=0.0):

    # Specify the stt dynamics in the simulation
    sim = Sim(mesh, name='dyn_%g'%beta, driver='llg_stt_cpp')

    sim.driver.set_tols(rtol=1e-12, atol=1e-12)
    sim.driver.alpha = 0.1
    sim.driver.gamma = 2.211e5
    sim.Ms = 8.6e5

    # sim.set_m(init_m)
    sim.set_m(np.load('m0.npy'))

    # Energies
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    anis = UniaxialAnisotropy(5e4)
    sim.add(anis)

    # beta is the parameter in the STT torque
    sim.a_J = global_const*1e11
    sim.p = (1,0,0)
    sim.beta = beta

    # The simulation will run for 5 ns and save
    # 500 snapshots of the system in the process
    ts = np.linspace(0, 0.5e-9, 21)
    
    xs=[]
    thetas=[]

    for t in ts:
        print('time', t)
        sim.run_until(t)
        spin = sim.spin.copy()
        x, theta = extract_dw(spin)
        xs.append(x)
        thetas.append(theta)
        sim.save_vtk()

    np.savetxt('dw_%g.txt'%beta,np.transpose(np.array([ts, xs,thetas])))

def analytical(ts):
    delta = 16.1245154966
    alpha = 0.1
    beta = 0.16
    aj = global_const*1e11
    theta = (beta-alpha)/(1+alpha**2)*ts*aj
    z = (1+alpha*beta)/(1+alpha**2)*delta*aj*ts
    return z, theta


def plot():
    data=np.loadtxt('dw_0.16.txt')
    ts = data[:,0]
    xs = data[:,1]-data[0,1]
    theta = data[:,2]-data[0,2]

    z_an, theta_an = analytical(ts)

    fig = plt.figure(figsize=(5, 2.6))
    
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 4])
    ax0 = plt.subplot(gs[0])
    ax0.plot(ts*1e9,xs,'.')
    ax0.plot(ts*1e9, z_an, '-')
    plt.ylabel(r'DW shift')
    plt.xlabel(r'Time (ns)')


    ax1 = plt.subplot(gs[1])
    ax1.plot(ts*1e9,theta,'.')
    ax1.plot(ts*1e9, theta_an, '-')
    plt.ylabel(r'$\phi$')
    plt.xlabel(r'Time (ns)')

    plt.tight_layout()
    analytical(ts)

    plt.savefig('xs_theta.pdf')
    

if __name__ == '__main__':

    # We will crate a mesh with 1000 elements of elements
    # in the x direction, and 1 along y and z
    # (so we have a 1D system)
    mesh = CuboidMesh(nx=1000, ny=1, nz=1,
                  dx=1.0, dy=2, dz=2.0,
                  unit_length=1e-9)

    # Relax the initial state
    relax_system(mesh)

    excite_system(mesh, beta=0.16)
    plot()

