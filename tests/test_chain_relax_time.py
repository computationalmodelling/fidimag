"""relax() of the chain methods advances the time by one step of the integrator for each iteration, not two."""
import numpy as np
import pytest

from fidimag.common import CuboidMesh
from fidimag.common.nebm_geodesic import NEBM_Geodesic
from fidimag.common.nebm_spherical import NEBM_Spherical
from fidimag.common.string_method import StringMethod
from fidimag.micro import Sim, UniaxialAnisotropy


def two_particles():
    mesh = CuboidMesh(nx=3, ny=1, nz=1, dx=3, dy=3, dz=3, unit_length=1e-9)
    sim = Sim(mesh, name='relax_time')
    sim.set_Ms(lambda pos: 3.8e5 if (pos[0] > 6 or pos[0] < 3) else 0)
    sim.add(UniaxialAnisotropy(1e5, axis=(1, 0, 0)))
    return sim


def mid(r):
    return (0.5, 0, 0.2) if r[0] > 4 else (-0.5, 0, 0.2)


def make(kind, integrator):
    images = [(-1, 0, 0), mid, (1, 0, 0)]
    if kind == 'spherical':
        return NEBM_Spherical(two_particles(), images, interpolations=[2, 2], spring_constant=1e4, name='relax_time')
    if kind == 'geodesic':
        return NEBM_Geodesic(two_particles(), images, interpolations=[2, 2], spring_constant=1e4, name='relax_time',
                             integrator=integrator)
    return StringMethod(two_particles(), images, interpolations=[2, 2], name='relax_time', integrator=integrator)


CASES = [('geodesic', 'cvode_bdf'), ('geodesic', 'verlet'), ('string', 'verlet'), ('spherical', None)]


@pytest.mark.parametrize('kind,integrator', CASES, ids=['geodesic-cvode', 'geodesic-verlet', 'string-verlet', 'spherical'])
def test_time_advances_once_per_iteration(kind, integrator, tmp_path, monkeypatch):
    monkeypatch.chdir(str(tmp_path))
    chain = make(kind, integrator)
    dts = []
    max_dYdt = chain.compute_maximum_dYdt

    def recording_max_dYdt(a, b, dt):
        dts.append(dt)
        return max_dYdt(a, b, dt)

    chain_run_until = chain.run_until
    targets = []

    def run_until_of_chain(t):
        before_chain = chain.t
        before_integrator = chain.integrator.t
        result = chain_run_until(t)
        targets.append((before_chain, before_integrator, t, chain.t, chain.integrator.t))
        return result

    chain.run_until = run_until_of_chain
    chain.compute_maximum_dYdt = recording_max_dYdt
    chain.relax(dt=1e-3, stopping_dYdt=1e-30, max_iterations=5, save_npys_every=10 ** 9, save_vtks_every=10 ** 9)

    assert len(targets) == 5 and len(dts) == 5
    for (t_before, integrator_before, target, t_after, integrator_after), dt in zip(targets, dts):
        assert target > t_before
        assert t_after == target                                   # the time of the object is where the integrator was sent
        assert integrator_after == pytest.approx(target)           # and where it arrived
        assert integrator_before == pytest.approx(t_before)        # it began where the time of the object said
        # dY/dt is divided by the time the integrator advanced
        assert dt == pytest.approx(integrator_after - integrator_before)
    # no step is skipped: each begins where the last ended
    for previous, current in zip(targets, targets[1:]):
        assert current[0] == previous[3]
    assert chain.t == targets[-1][2]
