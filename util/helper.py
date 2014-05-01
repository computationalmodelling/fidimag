from __future__ import division

import numpy as np

def normalise(a):
    """
    normalise the given array a
    """
    a.shape=(3,-1)
    b = np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
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
            spin[i] = m0(mesh.pos[i])
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


if __name__=='__main__':

    T=1000
    