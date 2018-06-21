# The idea behind these tests is to prevent regression back to a bug which was
# caused by incorrect referencing. This bug caused

import fidimag

def setup_fixture_micro():
    mesh = fidimag.common.CuboidMesh(nx=10, ny=10, nz=10,
                                     dx=1, dy=1, dz=1,
                                     unit_length=1e-9)
    sim = fidimag.micro.Sim(mesh)
    return mesh, sim


def setup_fixture_atomistic():
    mesh = fidimag.common.CuboidMesh(nx=10, ny=10, nz=10,
                                     dx=1, dy=1, dz=1,
                                     unit_length=1e-9)
    sim = fidimag.atomistic.Sim(mesh)
    return mesh, sim


def test_exchange_Ms_regression_micro():
    mesh, sim = setup_fixture_micro()
    A1 = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_m(Ms1)
    exch = fidimag.micro.UniformExchange(A1)
    sim.add(exch)
    sim.set_m(Ms2)
    assert exch.Ms is sim._magnetisation, "The Ms in the Exchange Micro class is not a reference to Ms in the Sim class"
    assert exch.Ms_inv is sim._magnetisation_inv, "The Ms_inv in the Exchange Micro class is not a reference to Ms_inv in the Sim Class"
