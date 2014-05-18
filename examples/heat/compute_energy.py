import numpy as np
import  scipy.integrate as integrate
from pccp.pc import Constant

const = Constant()

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


"""
compute the figure 2 in the paper at dx.doi.org/10.1063/1.2169472
""" 

def compute_energy(mx):
    Ms = 1.42e6
    v = 16e-27
    alpha = 0.005
    T = 300
    nu2 = alpha*2*const.k_B*T/(const.mu_0*Ms**2*v)
    mu = 2*alpha/nu2
    
    Dx = 0.0946
    Dy = 0.4132
    
    fun = lambda t: np.exp(-mu/2.0*(Dx*np.cos(t)**2+Dy*np.sin(t)**2))*np.sin(t)
    res = integrate.quad(fun, 0, np.pi)
    Z = res[0]
    
    gl = np.exp(-mu/2.0*(Dx*mx**2+Dy*(1-mx**2)))
    print Z
    
    return gl/Z
    
    

def plot_energy(name='out.pdf'):
    fig=plt.figure()
    mxs = np.linspace(-1,1,21)
    energy=[]
    for x in mxs:
        en = compute_energy(x)
        energy.append(en)
    
    plt.plot(mxs, energy, '.-')
    plt.xlabel('m_x')
    plt.ylabel('Energy')
    plt.ylim([0, 2])
    fig.savefig(name)





if __name__=='__main__':
    print compute_energy(0)
    plot_energy()
