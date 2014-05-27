import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from pc import Anisotropy
from pc import FDMesh
from pc import Sim
from pc import Zeeman
import numpy as np

def init_m(pos):
    x,y,z=pos
    if (x,y,z)==(1,2,3):
        return (1,2,3)
    elif z<1:
        return (0,0,-1)
    else:
        return (0,0,1)
    
def init_T(pos):
    return np.sum(pos)

def pin_fun(pos):
    if pos[0]==0:
        return 1
    else:
        return 0

def test_sim_pin():
    mesh=FDMesh(nx=3,ny=2,nz=1)
    sim = Sim(mesh)
    sim.set_m((0,0.8,0.6))
    sim.alpha = 0.1
    sim.gamma = 1.0
    sim.pins = pin_fun
    
    anis=Anisotropy(1.0, direction=(0,0,1), name='Dx')
    sim.add(anis)
    
    sim.run_until(1.0)
    print sim.spin
    assert sim.spin[0] == 0
    assert sim.spin[2] != 0
    
    

def test_sim_init_m():
    mesh=FDMesh(nx=3,ny=4,nz=5)
    sim=Sim(mesh)
    sim.set_m((0,1,0))
    sim.spin.shape=(3,-1)
    spin_y=sim.spin[1]
    assert(spin_y.any()==1)


def test_sim_init_m_fun():
    mesh=FDMesh(nx=3,ny=4,nz=5)
    sim=Sim(mesh)
    sim.set_m(init_m,normalise=False)
    assert(sim.spin_at(1,2,3)[0]==1)
    assert(sim.spin_at(1,2,3)[1]==2)
    assert(sim.spin_at(1,2,3)[2]==3)


def test_sim_T_fun():
    mesh=FDMesh(nx=3,ny=4,nz=5)
    sim=Sim(mesh)
    sim.set_T(init_T)
    assert(sim.T[0]==0)
    assert(sim.T[-1]==9)

def test_m_average():
    mesh=FDMesh(nx=3,ny=4,nz=5)
    sim=Sim(mesh)
    sim.set_m((0,0,1))
    a = sim.compute_average()
    assert a[2] == 1.0
    

def single_spin(alpha, gamma, H0, ts):
    """
    compute single spin under the external field H
    """
    
    precession = gamma/(1+alpha**2)
    beta = precession * H0 * ts
    
    mx = np.cos(beta)/np.cosh(alpha*beta)
    my = np.sin(beta)/np.cosh(alpha*beta)
    mz = np.tanh(alpha*beta)
    
    return mx,my,mz
    
    

def test_sim_single_spin_vode(do_plot=False):

    mesh=FDMesh(nx=1,ny=1,nz=1)

    sim = Sim(mesh,name='spin')
    
    alpha = 0.1
    gamma = 2.21e5
    sim.alpha = alpha
    sim.gamma = gamma
    sim.mu_s = 1.0
    
    sim.set_m((1, 0, 0))
    
    H0 = 1e5
    sim.add(Zeeman((0, 0, H0)))
    
    ts = np.linspace(0, 1e-9, 101)
    
    mx = []
    my = []
    mz = []
    real_ts=[]
    for t in ts:
        sim.run_until(t)
        real_ts.append(sim.t)
        print sim.t,abs(sim.spin_length()[0]-1)
        mx.append(sim.spin[0])
        my.append(sim.spin[1])
        mz.append(sim.spin[2])
    
    mz=np.array(mz)
    #print mz
    a_mx,a_my,a_mz = single_spin(alpha,gamma,H0, ts)

    print sim.stat()
    
    if do_plot:
        ts_ns = np.array(real_ts) * 1e9
        plt.plot(ts_ns, mx, ".", label="mx", color='DarkGreen')
        plt.plot(ts_ns, my, ".", label="my", color='darkslateblue') 
        plt.plot(ts_ns, mz, ".", label="mz", color='m') 
        plt.plot(ts_ns, a_mx, "--", label="analytical", color='b') 
        plt.plot(ts_ns, a_my, "--",  color='b') 
        plt.plot(ts_ns, a_mz, "--",  color='b') 
        plt.xlabel("time (ns)")
        plt.ylabel("m")
        plt.title("integrating a macrospin")
        plt.legend()
        plt.savefig("single_spin.pdf")
        
    print("Max Deviation = {0}".format(
            np.max(np.abs(mz - a_mz))))
   
    assert np.max(np.abs(mz - a_mz)) < 5e-7
    
def test_sim_single_spin_sllg(do_plot=False):

    mesh=FDMesh(nx=1,ny=1,nz=1)

    sim = Sim(mesh,name='spin',driver='sllg')
   
    alpha = 0.1
    gamma = 2.21e5
    
    sim.set_options(dt=1e-15, gamma=gamma)
    
    sim.alpha = alpha
    sim.mu_s = 1.0
    
    sim.set_m((1, 0, 0))
    
    H0 = 1e5
    sim.add(Zeeman((0, 0, H0)))
    
    ts = np.linspace(0, 1e-10, 101)
        
    
    mx = []
    my = []
    mz = []
    real_ts=[]
    for t in ts:
        sim.run_until(t)
        real_ts.append(sim.t)
        print sim.t,abs(sim.spin_length()[0]-1)
        mx.append(sim.spin[0])
        my.append(sim.spin[1])
        mz.append(sim.spin[2])
    
    mz=np.array(mz)
    
    a_mx,a_my,a_mz = single_spin(alpha,gamma,H0, ts)
    
    if do_plot:
        ts_ns = np.array(real_ts) * 1e9
        plt.plot(ts_ns, mx, ".", label="mx", color='DarkGreen')
        plt.plot(ts_ns, my, ".", label="my", color='darkslateblue') 
        plt.plot(ts_ns, mz, ".", label="mz", color='m') 
        plt.plot(ts_ns, a_mx, "--", label="analytical", color='b') 
        plt.plot(ts_ns, a_my, "--",  color='b') 
        plt.plot(ts_ns, a_mz, "--",  color='b') 
        plt.xlabel("time (ns)")
        plt.ylabel("m")
        plt.title("integrating a macrospin")
        plt.legend()
        plt.savefig("single_spin_sllg.pdf")
        
    print("Max Deviation = {0}".format(
            np.max(np.abs(mz - a_mz))))
   
    assert np.max(np.abs(mz - a_mz)) < 5e-6


def disable_test_sim_single_spin_llg_stt(do_plot=False):
    ni = Nickel()
    mesh=FDMesh(nx=1,ny=1,nz=1)
    mesh.set_material(ni)
    
    ni.alpha=0.1
    sim=Sim(mesh,driver='llg_stt')
    sim.set_m((1, 0, 0))
    
    H0 = 1
    sim.add(Zeeman((0, 0, H0)))
    
    dt = 1e-12;
    ts = np.linspace(0, 200 * dt, 101)
    
    precession = ni.gamma/(1+ni.alpha**2)
    
    mz_ref = []
    mxyz = []
    real_ts=[]
    
    
    for t in ts:
        sim.run_until(t)
        real_ts.append(sim.t)
        print sim.t,abs(sim.spin_length()[0]-1), sim.spin
        mz_ref.append(np.tanh(precession * ni.alpha * H0 * sim.t))
        mxyz.append(np.copy(sim.spin))
    
    mxyz=np.array(mxyz)
    
    if do_plot:
        ts_ns = np.array(real_ts) * 1e9
        
        plt.plot(ts_ns, mxyz[:,0], ".-", label="mx")
        plt.plot(ts_ns, mxyz[:,1], ".-", label="my")
        plt.plot(ts_ns, mxyz[:,2], ".-", label="mz")
        plt.plot(ts_ns, mz_ref, "-", label="analytical")
        plt.xlabel("time (ns)")
        plt.ylabel("mz")
        plt.title("integrating a macrospin")
        plt.legend()
        plt.savefig("test_llg_stt.png")
        
        print("Deviation = {0}".format(np.max(np.abs(mxyz[:,2] - mz_ref))))
    
    assert np.max(np.abs(mxyz[:,2] - mz_ref))<1e-9

  

if __name__=='__main__':
    #test_sim_init_m()
    #test_sim_init_m_fun()
    #test_sim_T_fun()
    #test_sim_single_spin_vode(do_plot=True)
    test_sim_single_spin_sllg(do_plot=True)
    #test_sim_single_spin_llg_stt(do_plot=True)
    #test_sim_single_spin(do_plot=True)
