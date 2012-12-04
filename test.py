import numpy as np
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy

if __name__=='__main__':
    mesh=FDMesh()
    spin=np.ones((mesh.nxyz,3))
    print spin
    exch=UniformExchange(mesh)
    exch.set_up(1,spin)
    print exch.compute_field()

    anis=Anisotropy(mesh)
    anis.set_up(1,spin)
    print anis.compute_field()
    pass
