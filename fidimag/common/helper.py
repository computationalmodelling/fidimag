from fidimag.extensions.c_clib import normalise
from fidimag.extensions.c_clib import init_scalar
from fidimag.extensions.c_clib import init_vector
from fidimag.extensions.c_clib import init_vector_func_fast
import fidimag.extensions.a_clib as clib
import numpy as np


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
        ids.append(mesh.index(p[0], p[1], p[2]))

    ids = np.array(ids)

    if comp == 'x':
        cmpi = 0
    elif comp == 'y':
        cmpi = 1
    elif comp == 'z':
        cmpi = 2
    else:
        raise Exception('Seems given component is wrong!!!')

    ids += cmpi * mesh.n

    all_data = []

    for ny in npys:

        all_data.append(np.load(ny)[ids])

    return np.array(all_data)


def compute_RxRy(mesh, spin, nx_start=0, nx_stop=-1, ny_start=0, ny_stop=-1):
    res = clib.compute_RxRy(spin, mesh.nx, mesh.ny,
                            mesh.nz, nx_start, nx_stop, ny_start, ny_stop)
    return res
