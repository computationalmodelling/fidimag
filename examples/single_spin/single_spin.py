"""
    Although we assume that the magnetisation length is one, the magnetisation length actually is not exact one 
    due to the reason that the cartesian ode solvers allow magnetisation length to change. In this script 
    we test how the tolerence setting influence the magnetisation length.
"""


import numpy as np
from pccp.pc import Sim
from pccp.pc import FDMesh
from pccp.pc import Zeeman, TimeZeeman
import matplotlib.pyplot as plt

def relax_system(rtol=1e-10,atol=1e-12):
    mesh = FDMesh(nx=1,ny=1,nz=1)
    sim = Sim(mesh)
    sim.set_options(rtol=rtol,atol=atol)
    sim.alpha = 0.0
    sim.gamma = 2.21e5
    sim.mu_s = 1.0
    
    sim.set_m((0,0.6,0.8))
    
    sim.add(Zeeman((0,0,1e5)))
    
    ts=np.linspace(0, 10e-9, 1001)
    mzs=[]
    lengths=[]
    
    for t in ts:
        sim.run_until(t)
        mzs.append(sim.spin[2]-0.8)
        lengths.append(sim.spin_length()[0]-1.0)
    
    return ts,mzs,lengths


def plot_m(ts, mz, lengths):
    fig = plt.figure()
    plt.plot(ts, lengths, '-', label='length')
    plt.plot(ts, mz, '--', label='mz')
    plt.xlabel('time (ns)')
    plt.ylabel('errors')
    plt.legend()
    fig.savefig('error.pdf')


if __name__=='__main__':
    
    ts,mzs,lengths = relax_system()
    plot_m(ts,mzs,lengths)
    
