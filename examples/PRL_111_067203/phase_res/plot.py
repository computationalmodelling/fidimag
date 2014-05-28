import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np

    
def plot_res(Ts, Hs, res):
    mesh_x, mesh_y = np.meshgrid(Ts, Hs)
    res.shape=(len(Ts), len(Hs))
    
    fig = plt.figure(figsize=(6,4))
    
    im = plt.imshow(np.transpose(res), cmap=cm.RdBu, extent=[Ts[0],Ts[-1],Hs[0],Hs[-1]], origin='lower', aspect='auto')
    plt.colorbar(im)
    
    plt.text(0.25, 0.4,'FM')
    plt.text(0.05, 0.15,'SkX')
    plt.text(0.1, 0.005,'HL')
    
    plt.xlabel(r'$k_B T/J$')
    plt.ylabel('H')
    plt.title('Skyrmion number')
    plt.tight_layout()
    
    fig.savefig('phase.pdf')


if __name__=='__main__':
    #np.random.seed(125)
    
    Ts = np.linspace(0,0.5,11)
    Hs = np.linspace(0,0.5,11)
    
    res = np.load('res.npy')
    
    plot_res(Ts, Hs, res)
