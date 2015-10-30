from fidimag.atomistic import DMI
from fidimag.common import CuboidMesh
from fidimag.atomistic import Sim
import numpy as np


def test_dmi_1d():

    mesh = CuboidMesh(nx=2, ny=1, nz=1)

    sim = Sim(mesh)
    sim.set_m((1, 0, 0))

    dmi = DMI(D=1)
    sim.add(dmi)

    field = dmi.compute_field()

    expected = np.array([0, 0, 0, 0, 0, 0])

    assert (field == expected).all()

    energy = dmi.compute_energy()
    assert energy == 0


def init_m(pos):
    if pos[0] < 1:
        return (0, 1, 0)
    else:
        return (0, 0, 1)


def test_dmi_1d_field():

    mesh = CuboidMesh(nx=2, ny=1, nz=1)

    sim = Sim(mesh)
    sim.set_m(init_m)

    dmi = DMI(D=1.23)
    sim.add(dmi)

    field = dmi.compute_field()

    expected = np.array([0, -1, 0, 0, 0, -1]) * 1.23

    assert np.allclose(field, expected)

    energy = dmi.compute_energy()

    assert energy == 1.23

if __name__ == '__main__':
    test_dmi_1d()
    test_dmi_1d_field()
