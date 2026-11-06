"""The interfacial atomistic DMI with a strength that differs from bond to bond."""
import numpy as np

from fidimag.atomistic import Sim, DMI
from fidimag.common import CuboidMesh


def test_the_field_uses_the_strength_of_each_bond():
    """H_i = (1 / mu_s) sum_j D_ij (d_j x S_j) over the bonds of the plane, with the strength D_ij of the bond i-j.

    (The kernel read the strength with the offset of the neighbour's spin, 3 * neighbour + j, and not with that of the bond, so that the field was
    wrong as soon as the strengths were not all the same.)"""
    mesh = CuboidMesh(nx=4, ny=3, nz=1, periodicity=(False, False, False))
    sim = Sim(mesh, name='dmi_interfacial_spatial')
    mu_s = 2.0
    sim.mu_s = mu_s

    def strengths(r):
        # a different strength for each of the six bonds of a site, and from site to site
        return tuple(1e-3 * (1.0 + 0.7 * j + 0.3 * r[0] + 0.1 * r[1]) for j in range(6))

    dmi = DMI(strengths, dmi_type='interfacial')
    sim.add(dmi)
    rng = np.random.default_rng(3)
    m = rng.normal(size=(mesh.n, 3))
    m /= np.linalg.norm(m, axis=1)[:, None]
    sim.set_m(m.reshape(-1))
    field = np.array(dmi.compute_field()).reshape(-1, 3)

    vectors = np.asarray(dmi.DMI_vector).reshape(-1, 3)
    expected = np.zeros((mesh.n, 3))
    for i in range(mesh.n):
        D = strengths(mesh.coordinates[i])
        for j in range(dmi.n_ngbs_dmi):
            neighbour = mesh.neighbours[i, j]
            if neighbour >= 0:
                expected[i] += D[j] * np.cross(vectors[j], m[neighbour])
    expected /= mu_s
    assert np.allclose(field, expected, rtol=1e-12, atol=1e-15)
    assert np.abs(expected).max() > 0
