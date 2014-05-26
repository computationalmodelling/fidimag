from __future__ import division

import numpy as np

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
    

if __name__=='__main__':

    T=1000
    