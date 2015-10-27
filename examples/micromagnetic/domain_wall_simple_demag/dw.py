import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fidimag.common import CuboidMesh, Sim, Zeeman, TimeZeeman, UniformExchange, UniaxialAnisotropy, SimpleDemag
from fidimag.common import BatchTasks, DataReader

def init_dw(pos):
    x = pos[0]

    if x < 400:
        return (1,0,0)
    elif x > 600:
        return (-1,0,0) 
    else:
        return (0,1,0)

def relax_system(mesh):
    A = 1.3e-11
    Ms = 8.6e5
    mu0 = 4*np.pi*1e-7

    Nx = 0
    Ny = 0.4
    Nz = 0.6

    K = 0.5*(Ny-Nx)*mu0*Ms**2
    Kp = -0.5*(Nz-Ny)*mu0*Ms**2     

    sim = Sim(mesh, name='relax')
    
    sim.set_tols(rtol=1e-8,atol=1e-8)
    sim.gamma = 2.211e5
    sim.Ms = Ms
    sim.alpha = 0.5
    sim.do_procession = False

    sim.set_m(init_dw)
    
    exch = UniformExchange(A=A)
    sim.add(exch)
    
    kx = UniaxialAnisotropy(K, axis=(1,0,0))
    sim.add(kx)

    kp = UniaxialAnisotropy(Kp, axis=(0,0,1))
    sim.add(kp)

    sim.relax(dt=1e-14, stopping_dmdt=0.01, max_steps=5000, save_m_steps=None, save_vtk_steps=None)

    np.save('m0.npy',sim.spin)

def excite_system_K(mesh, Hx=2000):
    
    A = 1.3e-11
    Ms = 8.6e5
    mu0 = 4*np.pi*1e-7
    G = A/(mu0*Ms**2)

    Nx = 0
    Ny = 0.4
    Nz = 0.6

    K = 0.5*(Ny-Nx)*mu0*Ms**2
    Kp = -0.5*(Nz-Ny)*mu0*Ms**2 

    sim = Sim(mesh, name='dyn_K')
    
    sim.set_tols(rtol=1e-8,atol=1e-8)
    sim.gamma = 2.211e5
    sim.Ms = Ms
    sim.alpha = 0.005

    sim.set_m(np.load('m0.npy'))
    
    exch = UniformExchange(A=A)
    sim.add(exch)

    kx = UniaxialAnisotropy(K, axis=(1,0,0))
    sim.add(kx)

    kp = UniaxialAnisotropy(Kp, axis=(0,0,1))
    sim.add(kp)
    
    hx = Zeeman((Hx,0,0), name='Hx')
    sim.add(hx, save_field=True)
    

    ts = np.linspace(0, 5e-9, 501)
    for t in ts:
        print 'time', t
        sim.run_until(t)

def excite_system_D(mesh, Hx=2000):
    
    A = 1.3e-11
    Ms = 8.6e5
    mu0 = 4*np.pi*1e-7
    G = A/(mu0*Ms**2)

    Nx = 0
    Ny = 0.4
    Nz = 0.6
    
    sim = Sim(mesh, name='dyn_D')
    
    sim.set_tols(rtol=1e-8,atol=1e-8)
    sim.gamma = 2.211e5
    sim.Ms = Ms
    sim.alpha = 0.005

    sim.set_m(np.load('m0.npy'))
    
    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = SimpleDemag(Nx=Nx, Ny=Ny, Nz=Nz)
    sim.add(demag)
    
    hx = Zeeman((Hx,0,0), name='Hx')
    sim.add(hx, save_field=True)
    
    ts = np.linspace(0, 5e-9, 101)
    for t in ts:
        print 'time', t
        sim.run_until(t)
        #sim.save_vtk()

def save_plot():

    data = DataReader('dyn_K.txt')
    ts = data['time']
    xs = data['m_x']

    fig=plt.figure()
    plt.plot(ts,xs)

    data = DataReader('dyn_D.txt')
    ts = data['time']
    xs = data['m_x']

    plt.plot(ts,xs,'+')

    fig.savefig('test.pdf')


if __name__=="__main__":

    mesh = CuboidMesh(nx = 1000, ny=1, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9)
    relax_system(mesh)

    excite_system_D(mesh)
    excite_system_K(mesh)


    save_plot()
    
    
