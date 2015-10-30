import numpy as np


def vector_field(mesh, v):
    """
    Returns a np.array with values specified by `v`, where `v` should
    be a iterable of length 3, or a function that returns an iterable of
    length 3 when getting the coordinates of a cell of `mesh`.

    """
    return field(mesh, v, dim=3)


def scalar_field(mesh, s):
    """
    Returns a np.array with values specified by `s`, where `s` should be
    a scalar, or a function that returns a scalar when getting the coordinates
    of a cell of `mesh`.

    """
    return field(mesh, s, dim=1)


def field(mesh, val, dim):
    values = np.zeros(mesh.vector_shape()) if dim == 3 else np.zeros(mesh.scalar_shape())

    if hasattr(val, '__call__'):
        f = val
    else:
        f = lambda r: val

    for i, r in enumerate(mesh.coordinates):
        values[i] = f(r)
    return values
