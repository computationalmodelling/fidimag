import numpy as np
import matplotlib as mpl
mpl.use("Agg")

import matplotlib.pyplot as plt
from fidimag.atomistic import Sim, DMI, UniformExchange, Zeeman, TimeZeeman
from fidimag.common import DataReader, BatchTasks, CuboidMesh
    
def deal_plot(Hy=0):

    data = DataReader('dyn.txt')
    ts = data['time']
    N = len(ts)
    #N = 8000
    
    dt = ts[1] - ts[0]
    print 'dt=',dt
    
    freq = np.fft.fftshift(np.fft.fftfreq(N, dt))
    
    H = data['h_z'][0:N]
    M = data['m_z'][0:N]
    
    fH = np.fft.fftshift(np.fft.fft(H))
    fM = np.fft.fftshift(np.fft.fft(M))
    
    a = fH.real
    b = fH.imag
    c = fM.real
    d = fM.imag
    
    rx = (a*c+b*d)/(a*a+b*b)
    ix = (b*c-a*d)/(a*a+b*b)
    
    ind = np.argmax(ix)
    print ind, freq[ind]*2*np.pi
    w_ix = np.array([freq*2*np.pi, ix])
    #np.savetxt('w_ix.txt',np.transpose(w_ix))
    
    plt.plot(freq*2*np.pi, ix, '.-')
    #plt.legend()
    plt.xlim([0, 0.06])
    #plt.ylim([-5, 100])
    plt.xlabel(r'$w$')
    plt.ylabel('Susceptibility')
    plt.savefig('res.pdf')

def plot_v(x,y,filename):
    fig=plt.figure()
    plt.plot(x,y,'^-',markersize=3)
    fig.savefig(filename)
                        
if __name__=='__main__':
    
    deal_plot()
    
