import numpy as np
from fidimag.common import CuboidMesh
from fidimag.micro import TimeZeeman

def time_fun(t, frequency):
   return 10*np.cos(frequency*t)

def fixture_setup(nx, ny, nz):
    """
    Fixtures for the tests
    """
    dx = 10.0/nx
    spin = np.zeros(3*nx*ny*nz)
    Ms = np.zeros(nx*ny*nz)
    mesh = CuboidMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=1, dz=1, unit_length=1e-9)
    frequency = 10e9
    return mesh, frequency, spin, Ms

def test_TimeZeeman():
    mesh, frequency, spin, Ms = fixture_setup(10, 10, 10)
    zee = TimeZeeman(np.array([0.0, 0.0, 1.0]), time_fun, extra_args=[frequency])
    zee.setup(mesh, spin, Ms)
    field1 = zee.compute_field(t=0)
    assert field1[0] == 0
    assert field1[1] == 0
    assert field1[2] == 10
    field2 = zee.compute_field(t=1)
    assert field2[0] == 0
    assert field2[1] == 0
    assert field2[2] == 10*np.cos(frequency)
