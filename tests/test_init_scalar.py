import fidimag
from fidimag.common import init_scalar
import numpy as np
import unittest

def m_init(pos, t):
    x, y, z = pos
    return np.sin(t)

def test_init_scalar_function():
    import fidimag
    import numpy as np
    # Create simulation
    nx = ny = nz = 2
    mesh = fidimag.common.CuboidMesh(nx=nx, ny=ny, nz=nz,
                                     dx=1, dy=1, dz=1)

    test_vec = np.zeros(nx * ny * nz)
    test_vec[:] = init_scalar(m_init, mesh, 0)
    assert test_vec[0] == 0.0
    test_vec[:] = init_scalar(m_init, mesh, np.pi/2.0)
    assert test_vec[0] == 1.0
