from __future__ import print_function
import numpy as np
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import Demag

mu0 = 4 * np.pi * 1e-7


def test_compute_field():
    """In an infinite film, we expect the demag tensor to be (0, 0, -1), and thus the
    magnetisation, if aligned in 0, 0, 1 direction, to create a demag field pointing
    with equal strength in the opposite direction.
    """

    mesh = CuboidMesh(nx=1, ny=1, nz=1, dx=2.0, dy=2.0, dz=2.0,
                      unit_length=1e-9, periodicity=(True, True, False))

    sim = Sim(mesh, name='relax')

    sim.set_tols(rtol=1e-10, atol=1e-14)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_precession = False

    sim.set_m((0, 0, 1))

    demag = Demag(pbc_2d=True)
    sim.add(demag)
    field=demag.compute_field()
    print((1 + field[2] / 8.6e5))
    assert abs(1 + field[2] / 8.6e5) < 1e-10


if __name__ == '__main__':
    test_compute_field()
