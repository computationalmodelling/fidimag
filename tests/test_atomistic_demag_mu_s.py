"""An atomistic Demag sees a mu_s that is set after it was added."""
import numpy as np

import fidimag.common.constant as const
from fidimag.atomistic import Sim, Demag
from fidimag.common import CuboidMesh


def mu_s(r):
    return const.mu_B * (2.0 + r[0] + 0.5 * r[1])


def test_the_demag_follows_a_later_mu_s():
    def build(demag_first):
        sim = Sim(CuboidMesh(nx=4, ny=3, nz=1, dx=0.3, dy=0.3, dz=0.3, unit_length=1e-9), name='demag_mu_s')
        demag = Demag()
        if demag_first:
            sim.mu_s = const.mu_B          # (the moment that the demag was set up with)
            sim.add(demag)
            sim.mu_s = mu_s                # (and the moment that it should use)
        else:
            sim.mu_s = mu_s
            sim.add(demag)
        sim.set_m(lambda r: (np.sin(4 * r[0]), np.cos(3 * r[1]), 0.5))
        return np.array(demag.compute_field()), demag.compute_energy()

    field_after, energy_after = build(demag_first=True)
    field_before, energy_before = build(demag_first=False)
    assert np.abs(field_before).max() > 0
    assert np.allclose(field_after, field_before, rtol=1e-12, atol=0)
    assert energy_after == energy_before or abs(energy_after - energy_before) <= 1e-12 * abs(energy_before)
