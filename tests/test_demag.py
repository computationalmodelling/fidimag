from fidimag.atomistic import FDMesh
from fidimag.atomistic import Sim
from fidimag.atomistic import Demag
import numpy as np


def test_demag_fft_exact():
    mesh = FDMesh(nx=5, ny=3, nz=4)
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
    # print fft,exact
    print np.max(np.abs(fft - exact))

    assert np.max(np.abs(fft - exact)) < 5e-22


def test_demag_fft_exact_oommf():
    mesh = FDMesh(nx=5, ny=3, nz=2)
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
    # print fft,exact
    print np.max(np.abs(fft - exact))

    assert np.max(np.abs(fft - exact)) < 1e-15


def test_demag_two_spin_xx():
    mesh = FDMesh(nx=2, ny=1, nz=1)
    sim = Sim(mesh)

    demag = Demag()
    sim.add(demag)

    sim.set_m((1, 0, 0))
    field = demag.compute_field()
    print field
    assert(field[0] == 2e-7)
    assert(field[1] == 2e-7)


if __name__ == '__main__':
    # test_demag_fft_exact()
    # test_demag_two_spin_xx()

    # test_oommf_coefficient()
    # test_demag_fft_exact()
    test_demag_fft_exact_oommf()
