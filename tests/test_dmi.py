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
    the DMI array. We use bulk DMI

    We specify DMI only along the +- z directions by setting the
    magnitude of the DMI vector along the other directions as 0

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

    mesh = CuboidMesh(nx=3, ny=1, nz=1, dx=1, dy=1, dz=1)
    sim = microSim(mesh)
    sim.Ms = 1
    sim.set_m(np.array([0, 0, 1]))


    sim.add(microDMI(D, dmi_type='bulk'))
    sim.compute_effective_field(0)
    D_field = C.mu_0 * sim.field.reshape(-1, 3)
    print(D_field)

    np.testing.assert_allclose(np.array([[ 0.,  1.,  0.],
 										  [ 0.,  0.,  0.],
 										  [ 0., -1.,  0.]]), D_field, rtol=1e-14)


def test_micro_dmi_array_interfacial():
    """
    Test the DMI in a micromagnetic simulation when manually specifying
    the DMI array. Here we use interfacial DMI (only defined for NNs
    in the xy plane)

    We specify DMI only along the +x direction by setting the
    magnitude of the DMI vector along the other directions as 0

    The system consists of 3 spins in the xz plane aligned as:

                                     z ^  . y
           -->   X   |                 | /
                     v                 |/___> x
          0      1   2

    i.e. in the +x, +y, -z directions. Thus the DMI field in this case, for
    the 1st node, is only computed using the 2nd node.

    """
    D = 1


    mesh = CuboidMesh(nx=3, ny=1, nz=1, dx=1, dy=1, dz=1)
    sim = microSim(mesh)
    sim.Ms = 1
    sim.set_m(np.array([0, 0, 1]))


    sim.add(microDMI(D, dmi_type='interfacial'))
    sim.compute_effective_field(0)

    # The field is computed as: - (1 / mu_0 Ms) D ( (r_ij X z) X M )
    # Thus the field at the middle site should be: -D ( (-y) X (0, 0, -1) )
    # which is (1, 0, 0), since Ms=1 and D is defined only for the
    # neighbours as +- z
    D_field = sim.field.reshape(-1, 3)
    D_field = C.mu_0 * sim.field.reshape(-1, 3)
    print(D_field)
    np.testing.assert_allclose(np.array([[ -1.,  0.,  0.],
 										  [ 0.,  0.,  0.],
 										  [ 1.,  0.,  0.]]), D_field, rtol=1e-14)


if __name__ == '__main__':
    # test_atom_dmi_1d()
    # test_atom_dmi_1d_field()
    # test_micro_dmi_array_bulk()
    test_micro_dmi_array_interfacial()
