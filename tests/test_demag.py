from __future__ import print_function
from fidimag.common import CuboidMesh
from fidimag.atomistic import Sim
from fidimag.atomistic import Demag
import numpy as np


def test_demag_fft_exact():
    mesh = CuboidMesh(nx=5, ny=3, nz=4, unit_length=1e-9)
    sim = Sim(mesh)

    demag = Demag()
    sim.add(demag)

    def init_m(pos):
        x = pos[0]
        if x <= 2:
            return (1, 0, 0)
        elif x >= 4:
            return (0, 0, 1)
        else:
            return (0, 1, 0)

    sim.set_m(init_m)
    fft = demag.compute_field()
    exact = demag.compute_exact()

    np.testing.assert_allclose(fft, exact, rtol=1e-10)
    



def test_demag_fft_exact_oommf():
    mesh = CuboidMesh(nx=5, ny=3, nz=2, unit_length=1e-9)
    sim = Sim(mesh)
    demag = Demag()
    sim.add(demag)

    def init_m(pos):
        x = pos[0]
        if x <= 2:
            return (1, 0, 0)
        elif x >= 4:
            return (0, 0, 1)
        else:
            return (0, 1, 0)
    sim.set_m(init_m)
    fft = demag.compute_field()
    exact = demag.compute_exact()
    np.testing.assert_allclose(fft, exact, rtol=1e-10)



def test_demag_two_spin_xx():
    mesh = CuboidMesh(nx=2, ny=1, nz=1)
    sim = Sim(mesh)
    demag = Demag()
    sim.add(demag)
    sim.set_m((1, 0, 0))
    field = demag.compute_field()
    print(field)
    assert(field[0] == 2e-7)
    assert(field[3] == 2e-7)


if __name__ == '__main__':
    test_demag_fft_exact()
    test_demag_two_spin_xx()
    test_demag_fft_exact_oommf()
