import numpy as np
cimport numpy as np

def normalise(a):
    """
    normalise the given array a
    """
    a.shape = (-1, 3)
    b = np.sqrt(a[:, 0] ** 2 + a[:, 1] ** 2 + a[:, 2] ** 2)
    ids = (b == 0)
    b[ids] = 1.0
    a[:, 0] /= b
    a[:, 1] /= b
    a[:, 2] /= b
    a.shape = (-1,)

def init_scalar(value, mesh, *args):

    n = mesh.n

    mesh_v = np.zeros(n)

    if isinstance(value, (int, float)):
        mesh_v[:] = value
    elif hasattr(value, '__call__'):
        for i in range(n):
            mesh_v[i] = value(mesh.coordinates[i], *args)

    elif isinstance(value, np.ndarray):
        if value.shape == mesh_v.shape:
            mesh_v[:] = value[:]
        else:
            raise ValueError("Array size must match the mesh size")

    return mesh_v

def init_vector(m0, mesh, norm=False, *args):

    n = mesh.n

    spin = np.zeros((n, 3))

    if isinstance(m0, list) or isinstance(m0, tuple):
        spin[:, :] = m0
        spin = np.reshape(spin, 3 * n, order='C')

    elif hasattr(m0, '__call__'):
        v = m0(mesh.coordinates[0], *args)
        if len(v) != 3:
            raise Exception(
                'The length of the value in init_vector method must be 3.')
        for i in range(n):
            spin[i, :] = m0(mesh.coordinates[i], *args)
        spin = np.reshape(spin, 3 * n, order='C')

    elif isinstance(m0, np.ndarray):
        if m0.shape == (3, ):
            spin[:] = m0  # broadcasting
        else:
            spin.shape = (-1)
            spin[:] = m0  # overwriting the whole thing

    spin.shape = (-1,)

    if norm:
        normalise(spin)

    return spin


# def init_vector_func_fast(m0, mesh, norm=False, *args):
#     """
#     A fast version of init_vector for functions that does
#     not perform any checking of input arguments.
#     """
#     v = m0(mesh.coordinates[0], *args)
#     for i in range(n):
#         spin[3*i:] = m0(&mesh, *args)
#     spin = np.reshape(spin, 3 * n, order='C')