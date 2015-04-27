from __future__ import division
#import matplotlib as mpl
#mpl.use("Agg")
import fidimag.extensions.clib
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import colorConverter
from matplotlib.collections import PolyCollection, LineCollection

def normalise(a):
    """
    normalise the given array a
    """
    a.shape=(3,-1)
    b = np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
    ids = (b == 0)
    b[ids] = 1.0
    a /= b
    a.shape=(-1,)

def init_vector(m0,mesh,norm=False):
    
    nxyz = mesh.nxyz
    
    spin = np.zeros((nxyz,3))
    
    if isinstance(m0,list) or isinstance(m0,tuple):    
        spin[:,:] = m0
        spin = np.reshape(spin, 3*nxyz, order='F')
    elif hasattr(m0, '__call__'):
        for i in range(nxyz):
            v = m0(mesh.pos[i])
            if len(v)!=3:
                raise Exception('The length of the value in init_vector method must be 3.')
            spin[i,:] = v[:]
        spin = np.reshape(spin, 3*nxyz, order='F')
    elif isinstance(m0,np.ndarray):
        spin.shape=(-1,)
        if m0.shape == spin.shape:
            spin[:] = m0[:]
    
    spin.shape=(-1,)
    
    if norm:
        normalise(spin)

    return spin

def init_scalar(value,mesh):
    
    nxyz = mesh.nxyz
    
    mesh_v = np.zeros(nxyz)
    
    if isinstance(value,(int, float)):    
        mesh_v[:] = value
    elif hasattr(value, '__call__'):
        for i in range(nxyz):
            mesh_v[i] = value(mesh.pos[i])
        
    elif isinstance(value,np.ndarray):
        
        if value.shape == mesh_v.shape:
            mesh_v[:] = value[:]

    return mesh_v


def extract_data(mesh, npys, pos, comp='x'):
    """
    extract data of special positions for given npy data
    
    npys:
        the names of npys
    
    pos:
         something like [(1,0,0),...,(2,3,4)]
    """
    ids = []
    for p in pos:
        ids.append(mesh.index(p[0],p[1],p[2]))
        
    ids = np.array(ids)
    
    if comp == 'x':
        cmpi = 0
    elif comp == 'y':
        cmpi = 1
    elif comp == 'z':
        cmpi = 2
    else:
        raise Exception('Seems given component is wrong!!!')
        
    ids += cmpi*mesh.nxyz 
    
    all_data = []

    for ny in npys:
        
        all_data.append(np.load(ny)[ids])
    
    return np.array(all_data)


def plot_m(mesh, npy, comp='x', based=None):
    
    if comp == 'x':
        cmpi = 0
    elif comp == 'y':
        cmpi = 1
    elif comp == 'z':
        cmpi = 2
    else:
        raise Exception('Seems the given component is wrong!!!')
    
    data = np.load(npy)
    
    if based is not None:
        data = data - based
    
    data.shape = (3,-1)
    m = data[cmpi]
    
    nx = mesh.nx
    ny = mesh.ny
    
    m.shape=(nx,ny)
    
    fig = plt.figure()
    #norm=color.Normalize(-1,1)
    plt.imshow(np.transpose(m), aspect = 1, cmap = plt.cm.coolwarm, origin='lower', interpolation='none')
    plt.autoscale(False)
    plt.xticks([])
    plt.yticks([])
    fig.savefig('%s_%s.png'%(npy[:-4],comp))

def plot_energy_2d(name, step=-1):
    """
    Plot the energy path at given step.

    name is the simulation name.
    """
    import matplotlib.pyplot as plt

    data = np.loadtxt('%s_energy.ndt' % name)
    dms = np.loadtxt('%s_dms.ndt' % name)

    if data.ndim == 1:
        data.shape = (1, -1)
        dms.shape = (1, -1)

    if step < 0:
        step = data[-1, 0]
        id = -1
    else:
        steps = abs(data[:, 0] - step)
        id = np.argmin(steps)
        step = data[id, 0]

    fig = plt.figure()
    xs = range(1, len(data[0, :]))

    for i in range(len(xs)):
        xs[i] = sum(dms[id, 1:i + 1])

    plt.plot(xs, data[id, 1:], '.-')

    plt.legend()
    plt.grid()
    plt.ylabel('Energy (J)')
    plt.xlabel('Position in path (a.u.)')

    fig.savefig('energy_%d.pdf' % step)


def plot_energy_3d(name, key_steps=50, filename=None):
    
    data = np.loadtxt('%s_energy.ndt' % name)

    if data.ndim == 1:
        data.shape = (1, -1)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # image index
    xs = range(1, data.shape[1])

    steps = data[:, 0]

    each_n_step = int(len(steps) / key_steps)

    if each_n_step < 1:
        each_n_step = 1

    cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
    colors = [cc('r'), cc('g'), cc('b'), cc('y')]
    facecolors = []
    line_data = []
    energy_min = np.min(data[:, 1:])

    zs = []
    index = 0
    for i in range(0, len(steps), each_n_step):
        line_data.append(list(zip(xs, data[i, 1:] - energy_min)))
        facecolors.append(colors[index % 4])
        zs.append(data[i, 0])
        index += 1

    poly = PolyCollection(line_data, facecolors=facecolors, closed=False)
    poly.set_alpha(0.7)

    ax.add_collection3d(poly, zs=zs, zdir='x')

    ax.set_xlabel('Steps')
    ax.set_ylabel('images')
    ax.set_zlabel('Energy (J)')

    ax.set_ylim3d(0, len(xs) + 1)
    ax.set_xlim3d(0, int(data[-1, 0]) + 1)
    ax.set_zlim3d(0, np.max(data[:, 1:] - energy_min))

    if filename is None:
        filename = '%s_energy_3d.pdf' % name

    fig.savefig(filename)


def compute_RxRy(mesh, spin):
    res = clib.compute_RxRy(spin, mesh.nx, mesh.ny, mesh.nz)
    return res
    

if __name__=='__main__':

    T=1000
    
