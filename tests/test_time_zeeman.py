import numpy as np
from fidimag.common import CuboidMesh
from fidimag.micro import TimeZeeman, TimeZeemanSimple, TimeZeemanFast
from fidimag.user.example import TimeZeemanFast_test_time_fun


def time_fun_spatial(pos, t, frequency):
   return (0.0, 0.0, 10*np.cos(frequency*t))

def time_fun(t, frequency):
   return (0.0, 0.0, 10*np.cos(frequency*t))



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
    zee = TimeZeeman(time_fun_spatial, extra_args=[frequency])
    zee.setup(mesh, spin, Ms)
    field1 = zee.compute_field(t=0)
    assert field1[0] == 0
    assert field1[1] == 0
    assert field1[2] == 10

    field2 = zee.compute_field(t=1)
    assert field2[0] == 0
    assert field2[1] == 0
    assert field2[2] == 10*np.cos(frequency)

def test_TimeZeemanSimple():
    mesh, frequency, spin, Ms = fixture_setup(10, 10, 10)
    zee = TimeZeemanSimple(time_fun, extra_args=[frequency])
    zee.setup(mesh, spin, Ms)
    field1 = zee.compute_field(t=0)
    assert field1[0] == 0
    assert field1[1] == 0
    assert field1[2] == 10

    field2 = zee.compute_field(t=1)
    assert field2[0] == 0
    assert field2[1] == 0
    assert field2[2] == 10*np.cos(frequency)

def test_TimeZeemanFast():
    mesh, frequency, spin, Ms = fixture_setup(10, 10, 10)
    zee = TimeZeemanFast(TimeZeemanFast_test_time_fun, extra_args=[frequency])
    zee.setup(mesh, spin, Ms)
    field1 = zee.compute_field(t=0)
    assert field1[0] == 0
    assert field1[1] == 0
    assert field1[2] == 10

    field2 = zee.compute_field(t=1)
    assert field2[0] == 0
    assert field2[1] == 0
    assert field2[2] == 10*np.cos(frequency)

if __name__ == "__main__":
    import time
    # Run the speed tests...
    for i in range(1, 5):
        nx = 10**i
        ny = nz = 10
        print('\nN = {} * {} * {} = '.format(nx, ny, nz,  nx*ny*nz))
        mesh, frequency, spin, Ms = fixture_setup(nx, ny, nz)

        zee = TimeZeeman(time_fun_spatial, extra_args=[frequency])
        zee.setup(mesh, spin, Ms)
        t1 = time.time()
        field1 = zee.compute_field(t=0)
        field2 = zee.compute_field(t=1)
        t2 = time.time()
        print('  Timing results for TimeZeeman = {}'.format(t2 - t1))

        zee = TimeZeemanSimple(time_fun, extra_args=[frequency])
        zee.setup(mesh, spin, Ms)
        t3 = time.time()
        field1 = zee.compute_field(t=0)
        field2 = zee.compute_field(t=1)
        t4 = time.time()
        print('  Timing results for TimeZeemanSimple = {}'.format(t4 - t3))

        zee = TimeZeemanFast(TimeZeemanFast_test_time_fun, extra_args=[frequency])
        zee.setup(mesh, spin, Ms)
        t5 = time.time()
        field1 = zee.compute_field(t=0)
        field2 = zee.compute_field(t=1)
        t6 = time.time()
        print('  Timing results for TimeZeemanFast = {}'.format(t6 - t5))
        print('\n  Normalised: TimeZeeman = 1, TimeZeemanSimple = {}, TimeZeemanFast = {}'.format((t4-t3)/(t2-t1), (t6-t5)/(t2-t1)))
