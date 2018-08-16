from fidimag.common import CuboidMesh
from fidimag.atomistic import DMI as atomDMI
from fidimag.atomistic import Sim as atomSim
from fidimag.micro import DMI as microDMI
from fidimag.micro import Sim as microSim
import numpy as np
import fidimag.common.constant as C


def test_atom_dmi_1d():

    mesh = CuboidMesh(nx=2, ny=1, nz=1)

    sim = atomSim(mesh)
    sim.set_m((1, 0, 0))

    dmi = atomDMI(D=1)
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


def test_atom_dmi_1d_field():

    mesh = CuboidMesh(nx=2, ny=1, nz=1)

    sim = atomSim(mesh)
    sim.set_m(init_m)

    dmi = atomDMI(D=1.23)
    sim.add(dmi)

    field = dmi.compute_field()

    expected = np.array([0, -1, 0, 0, 0, -1]) * 1.23

    assert np.allclose(field, expected)

    energy = dmi.compute_energy()

    assert energy == 1.23


def test_micro_dmi_array_bulk():
    """
    Test the DMI in a micromagnetic simulation when manually specifying
    the DMI vectors. We use the z modulation of bulk DMI

    We specify DMI only along the +- z directions by setting the
    magnitude of the DMI vectors along the other directions as 0

    The system consists of 4 spins in the xz plane aligned as:
          2        3
            X    X                z ^  . y
                                    | /
           -->  <--                 |/___> x
          0        1

    i.e. in the +x, -x, +y, +y directions. Thus the DMI field in this case, for
    the 0th node, is only computed using the 2nd node. Similarly for the
    others, because the DMI is only defined for NNs along the z directions.

    """
    D = 1
    # The DMI vectors
    DMI_vector = [0., 0, 0,  # -x
                  0, 0, 0,   # +x
                  0, 0, 0,   # -y
                  0, 0, 0,   # +y
                  0, 0, -1,  # -z
                  0, 0, 1,   # +z
                  ]

    mesh = CuboidMesh(nx=2, ny=1, nz=2, dx=1, dy=1, dz=1)
    sim = microSim(mesh)
    sim.Ms = 1
    sim.set_m(np.array([1, 0, 0,
                        -1, 0, 0,
                        0, 1, 0,
                        0, 1, 0]))

    # sim.add(microDMI(DMInorms, dmi_type='bulk'))
    # With the custom option, we manually pass a DMI vector array
    # made of 18 * n components
    sim.add(microDMI(D, dmi_type='custom', dmi_vector=DMI_vector))
    sim.compute_effective_field(0)

    # The field is computed as: - (1 / mu_0 Ms) D ( r_ij X M )
    # Thus the field at the 0th site should be: -D ( z X (0, 1, 0) )
    # which is (1, 0, 0), since Ms=1 and D is defined only for the
    # neighbours as +- z
    D_field = sim.field.reshape(-1, 3)
    hx, hy, hz = C.mu_0 * D_field[0]
    assert np.abs(hx - 1) < 1e-10
    assert np.abs(hy) < 1e-10
    assert np.abs(hz) < 1e-10


def test_micro_dmi_array_interfacial():
    """
    Test the DMI in a micromagnetic simulation when manually specifying
    the DMI array. Here we use interfacial DMI (only defined for NNs
    in the xy plane)

    We specify DMI only along the +x / +y directions by setting the
    magnitude of the DMI vector along the other directions as 0

    The system consists of 3 spins in the xz plane aligned as:

                                     z ^  . y
           -->   X   |                 | /
                     v                 |/___> x
          0      1   2

    i.e. in the +x, +y, -z directions. Thus the DMI field in this case, for
    the 1st node, is only computed using the 2nd node.

    """
    # This norm is the same for every lattice site
    D = 1
    # The DMI vectors as interfacial DMI with -x and -y vecs as zero
    DMI_vector = [0., 0, 0,  # -x
                  0, 1, 0,   # +x
                  0, 0, 0,   # -y
                  -1, 0, 0,  # +y
                  0, 0, -1,  # -z
                  0, 0, 1,   # +z
                  ]

    mesh = CuboidMesh(nx=3, ny=1, nz=1, dx=1, dy=1, dz=1)
    sim = microSim(mesh)
    sim.Ms = 1
    sim.set_m(np.array([1, 0, 0,
                        0, 1, 0,
                        0, 0, -1]))

    sim.add(microDMI(D, dmi_type='custom', dmi_vector=DMI_vector))
    sim.compute_effective_field(0)

    # The field is computed as: - (1 / mu_0 Ms) D ( (r_ij X z) X M )
    # Thus the field at the middle site should be: -D ( (-y) X (0, 0, -1) )
    # which is (1, 0, 0), since Ms=1 and D is defined only for the
    # neighbours as +- z
    D_field = sim.field.reshape(-1, 3)
    hx, hy, hz = C.mu_0 * D_field[1]
    assert np.abs(hx - 1) < 1e-10
    assert np.abs(hy) < 1e-10
    assert np.abs(hz) < 1e-10

    # At the first site the field is zero
    hx, hy, hz = C.mu_0 * D_field[0]
    assert np.abs(hx) < 1e-10
    assert np.abs(hy) < 1e-10
    assert np.abs(hz) < 1e-10


if __name__ == '__main__':
    # test_atom_dmi_1d()
    # test_atom_dmi_1d_field()
    test_micro_dmi_array_bulk()
    test_micro_dmi_array_interfacial()
