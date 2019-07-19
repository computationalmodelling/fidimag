import fidimag
import numpy as np

def setup_fixture_micro(driver='llg'):
    mesh = fidimag.common.CuboidMesh(nx=3, ny=3, nz=3,
                                     dx=1, dy=1, dz=1,
                                     unit_length=1e-9)
    sim = fidimag.micro.Sim(mesh, driver=driver)
    return mesh, sim


def setup_fixture_atomistic(driver='llg'):
    mesh = fidimag.common.CuboidMesh(nx=3, ny=3, nz=3,
                                     dx=1, dy=1, dz=1,
                                     unit_length=1e-9)
    sim = fidimag.atomistic.Sim(mesh, driver=driver)
    return mesh, sim

def set_Ms(pos):
    if pos[0] < 2:
        return 0
    else:
        return 5e8

def test_llg_Ms_inv():
    mesh, sim = setup_fixture_micro()
    sim.set_Ms(set_Ms)
    assert np.isnan(np.dot(sim._magnetisation_inv, sim._magnetisation_inv)) == False

def test_llg_atomistic_mu_s_inv():
    mesh, sim = setup_fixture_atomistic()
    sim.set_mu_s(set_Ms)
    assert np.isnan(np.dot(sim._magnetisation_inv, sim._magnetisation_inv)) == False
