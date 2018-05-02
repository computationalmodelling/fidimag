import fidimag
from fidimag.common import init_vector
import numpy as np
import unittest

def m_init(pos, t):
    x, y, z = pos
    return (np.sin(t), np.cos(t), 0.0)

def test_init_scalar_function():
    import fidimag
    import numpy as np
    # Create simulation
    nx = ny = nz = 2
    mesh = fidimag.common.CuboidMesh(nx=nx, ny=ny, nz=nz,
                                     dx=1, dy=1, dz=1)

    test_vec = np.zeros(3 * nx * ny * nz)
    test_vec[:] = init_vector(m_init, mesh, False, 0)
    assert test_vec[0] < 1e-15
    assert test_vec[1] - 1.0 < 1e-15
    test_vec[:] = init_vector(m_init, mesh, False, np.pi/2.0)
    assert test_vec[0] - 1.0 < 1e-15
    assert test_vec[1] < 1e-15
