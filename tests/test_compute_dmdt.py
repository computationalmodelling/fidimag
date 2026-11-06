"""compute_dmdt: the largest change of a spin per unit of time, the quantity that relax() compares with stopping_dmdt."""
import numpy as np
import pytest

import fidimag.atomistic as atomistic
from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange


def micro():
    mesh = CuboidMesh(nx=4, ny=3, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9)
    sim = Sim(mesh, name='dmdt_micro')
    sim.set_Ms(8e5)
    sim.set_m((0, 0, 1))
    sim.add(UniformExchange(1.3e-11))
    return sim


def atomistic_sim():
    mesh = CuboidMesh(nx=4, ny=3, nz=1, dx=0.5, dy=0.5, dz=0.5, unit_length=1e-9)
    sim = atomistic.Sim(mesh, name='dmdt_atomistic')
    sim.set_mu_s(5.56e-23)
    sim.set_m((0, 0, 1))
    sim.add(atomistic.Exchange(1.2e-21))
    return sim


@pytest.fixture(params=[micro, atomistic_sim], ids=['micro', 'atomistic'])
def driver(request):
    return request.param().driver


def move(driver, changes):
    """Leave spin_last as the spins are, then change the spins: {site: (dx, dy, dz)}"""
    driver.spin_last[:] = driver.spin
    spins = driver.spin.reshape(-1, 3)
    for site, d in changes.items():
        spins[site] += d


def test_nothing_moved(driver):
    move(driver, {})
    assert driver.compute_dmdt(1e-12) == 0.0


def test_one_spin_moves_in_all_directions(driver):
    """The change of a spin is the length of its change in x, y and z together: a 3-4-5 triangle"""
    move(driver, {0: (0.3, 0.4, 0.0)})
    assert driver.compute_dmdt(1e-12) == pytest.approx(0.5 / 1e-12, rel=1e-14)


@pytest.mark.parametrize('site', [0, 1, 5, 11])
def test_it_does_not_depend_on_where_the_spin_is(driver, site):
    move(driver, {site: (0.0, 0.3, 0.4)})
    assert driver.compute_dmdt(2.0) == pytest.approx(0.25, rel=1e-14)


def test_the_largest_change_of_all_the_spins_counts(driver):
    move(driver, {0: (0.3, 0.4, 0.0), 2: (0.1, 0.1, 0.1), 7: (0.0, 0.0, 0.45), 11: (-0.2, 0.2, 0.0)})
    assert driver.compute_dmdt(1.0) == pytest.approx(0.5, rel=1e-14)


def test_neighbouring_spins_do_not_add_up(driver):
    """Two neighbouring spins each move by 0.3 in x: the largest change of one spin is 0.3, not the length of the two together"""
    move(driver, {3: (0.3, 0.0, 0.0), 4: (0.3, 0.0, 0.0)})
    assert driver.compute_dmdt(1.0) == pytest.approx(0.3, rel=1e-14)
