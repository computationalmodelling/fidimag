"""The string method spaces the images along the string that the integrator produced."""
import numpy as np
import pytest

from fidimag.common import CuboidMesh
from fidimag.common.string_method import StringMethod
from fidimag.micro import Sim, UniaxialAnisotropy


def two_particles(name):
    mesh = CuboidMesh(nx=3, ny=1, nz=1, dx=3, dy=3, dz=3, unit_length=1e-9)
    sim = Sim(mesh, name=name)
    sim.set_Ms(lambda pos: 3.8e5 if (pos[0] > 6 or pos[0] < 3) else 0)
    sim.add(UniaxialAnisotropy(1e5, axis=(1, 0, 0)))
    return sim


def mid(r):
    return (0.5, 0, 0.2) if r[0] > 4 else (-0.5, 0, 0.2)


@pytest.mark.parametrize('integrator', ['verlet', 'cvode_bdf'])
def test_the_distances_are_those_of_the_integrated_band(integrator, tmp_path, monkeypatch):
    """`run_until` took the distances from `self.band` before it was updated from the integrator. With CVODE `self.band` holds the last state
    that the solver tried, which is close to the one it returned but not the same."""
    monkeypatch.chdir(str(tmp_path))
    string = StringMethod(two_particles('string_distances'), [(-1, 0, 0), mid, (1, 0, 0)], interpolations=[2, 2], name='string_distances',
                          integrator=integrator)
    seen = []
    compute_distances = string.compute_distances

    def recording():
        seen.append((string.band.copy(), string.integrator.y.copy()))
        return compute_distances()

    string.compute_distances = recording
    for i in range(1, 4):
        string.run_until(i * 1e-4)
    assert len(seen) == 3
    for band, integrated in seen:
        assert np.array_equal(band, integrated)
