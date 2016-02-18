import numpy as np
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, Demag, DMI
from fidimag.micro import Zeeman, TimeZeeman
from fidimag.common.fileio import DataReader

mu0 = 4 * np.pi * 1e-7


def test_compute_field():

    mesh = CuboidMesh(nx=1, ny=1, nz=1, dx=2.0, dy=2.0, dz=2.0, 
                      unit_length=1e-9, periodicity=(True, True, False))

    sim = Sim(mesh, name='relax')

    sim.set_tols(rtol=1e-10, atol=1e-14)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_procession = False

    sim.set_m((0,0,1))

    A = 1.3e-11
    exch = UniformExchange(A=A)

    demag = Demag(pbc_2d=True)
    sim.add(demag)
    field=demag.compute_field()
    print(1 + field[2] / 8.6e5) 
    assert abs(1 + field[2] / 8.6e5) < 1e-10

    #np.save('m0.npy', sim.spin)


if __name__ == '__main__':

    test_compute_field()
