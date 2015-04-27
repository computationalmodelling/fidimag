import numpy as np
import fidimag.extensions.clib as clib

def test_sst_field_1d():
    #suppose the system is 1d with 4 spins.
    nx,ny,nz=4,1,1
    dx,dy=1,1
    
    nxyz=nx*ny*nz*3
    spin=np.array([1.0*i for i in range(nxyz)])
    field = np.zeros(nxyz)
    jx = np.zeros(nxyz)
    jy = np.zeros(nxyz)
    jx[:] = 1
    jy[:] = 2
    
    
    spin[1]=0.8
    spin[2]=2
    spin[3]=2.1
    
    print spin
    
    clib.compute_stt_field(spin,
                           field,
                           jx,jy,
                           dx,dy,
                           nx,ny,nz,
                           0,0)
    
    
    assert field[0]==0.8
    assert field[1]==1
    assert field[2]==0.65
    assert abs(field[3]-0.1)<1e-16
    
    clib.compute_stt_field(spin,
                           field,
                           jx,jy,
                           dx,dy,
                           nx,ny,nz,
                           1,0)
    #print field
    assert field[0]==-0.65
    assert field[3]==-1


if __name__=='__main__':
    test_sst_field_1d()
