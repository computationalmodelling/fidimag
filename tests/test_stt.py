import numpy as np
import clib

def test_sst_field_1d():
    #suppose the system is 1d with 4 spins.
    nx,ny,nz=4,1,1
    jx,jy,jz=1,2,3
    dx,dy,dz=1,1,1
    
    nxyz=nx*ny*nz*3
    spin=np.array([1.0*i for i in range(nxyz)])
    field = np.zeros(nxyz)
    
    spin[1]=0.8
    spin[2]=2
    spin[3]=2.1
    
    print spin
    
    clib.compute_stt_field(spin,
                           field,
                           jx,jy,jz,
                           dx,dy,dz,
                           nx,ny,nz)
    
    
    assert field[0]==0.4
    assert field[1]==1
    assert field[2]==0.65
    assert abs(field[3]-0.05)<1e-16


if __name__=='__main__':
    test_sst_field_1d()