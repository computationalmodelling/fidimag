"""
Draw the image spacing produced by the weighted NEBM spring forces.

This is the plotting that used to live inside
tests/test_two_particles_neb_method.py. It asserts nothing, it only renders
the comparison, so it does not belong in the test suite: the tests there now
check the weighted segment lengths directly against a quadrature
(test_energy_weighted_spring_lengths, test_curvature_weighted_spring_lengths)
and check the resulting spacing statistically in the two slow clustering
tests.

Run it when the clustering needs to be looked at rather than measured::

    python sandbox/nebm_spring_weighting_plots.py

It writes neb_2particles_spring_force_ratio_comparison.png and
neb_2particles_curvature_weighted_comparison.png into the working directory.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from fidimag.micro import Sim, UniaxialAnisotropy
from fidimag.common import CuboidMesh
from fidimag.common.nebm_geodesic import NEBM_Geodesic

Kx = 1e5
Ms = 3.8e5
mesh = CuboidMesh(nx=3, ny=1, nz=1, dx=3, dy=3, dz=3, unit_length=1e-9)


def two_part(pos):
    return Ms if (pos[0] > 6 or pos[0] < 3) else 0


def mid_m(pos):
    return (0.5, 0, 0.2) if pos[0] > 4 else (-0.5, 0, 0.2)


def make_neb(simname):
    sim = Sim(mesh)
    sim.Ms = two_part
    sim.add(UniaxialAnisotropy(Kx, axis=(1, 0, 0)))
    return NEBM_Geodesic(sim, [(-1, 0, 0), mid_m, (1, 0, 0)],
                         interpolations=[10, 10], spring_constant=1e4,
                         name=simname, integrator='sundials')


def relaxed_pair(weighting, ratio=0.9):
    """
    Returns (uniform band, weighted band). The weighted one follows the
    staged protocol described in the clustering tests: relax plainly, pull
    one image onto the saddle with a climbing image, then switch on the
    weighting together with the fixed-step VP integrator.
    """
    uniform = make_neb('plot_%s_uniform' % weighting)
    uniform.relax(max_iterations=2000, save_vtks_every=5000,
                  save_npys_every=5000, stopping_dYdt=1e-4,
                  stopping_max_force=1e-2, dt=1e-6)

    weighted = make_neb('plot_%s_weighted' % weighting)
    weighted.relax(max_iterations=2000, save_vtks_every=5000,
                   save_npys_every=5000, stopping_dYdt=1e-4,
                   stopping_max_force=1e-2, dt=1e-6)

    peak = int(np.argmax(weighted.energies))
    weighted.climbing_image = [peak]
    weighted.initialise_integrator(integrator='sundials')
    weighted.relax(max_iterations=30, save_vtks_every=50000,
                   save_npys_every=50000, stopping_dYdt=1e-6, dt=1e-6,
                   save_initial_state=False)

    del weighted.climbing_image
    weighted.spring_weighting = weighting
    weighted.spring_force_ratio = ratio
    weighted.initialise_integrator(integrator='verlet')
    weighted.relax(max_iterations=20000, save_vtks_every=50000,
                   save_npys_every=50000, stopping_dYdt=1e-8,
                   stopping_max_force=1e-2, dt=1e-4,
                   save_initial_state=False)
    return uniform, weighted


def critical_image_indices(neb):
    """Images where dE/d(path_distance) ~ 0, i.e. every local extremum of
    the band including the two fixed end states."""
    E = neb.energies
    idxs = [0, len(E) - 1]
    for i in range(1, len(E) - 1):
        if (E[i] - E[i - 1]) * (E[i + 1] - E[i]) <= 0:
            idxs.append(i)
    return np.array(sorted(set(idxs)))


def plot(uniform, weighted, weighted_title, filename):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex='col')
    for col, (neb, title) in enumerate(
            [(uniform, 'spring_force_ratio = 0 (uniform spacing)'),
             (weighted, weighted_title)]):
        crit = critical_image_indices(neb)

        ax_E = axes[0, col]
        ax_E.plot(neb.path_distances, neb.energies / 1.602e-19, 'o-')
        ax_E.plot(neb.path_distances[crit], neb.energies[crit] / 1.602e-19,
                  'r*', markersize=14, label='critical images (dE/dx ~ 0)')
        ax_E.set_title(title)
        ax_E.set_ylabel('Energy (eV)')
        ax_E.legend()

        # Spacing plotted at the midpoint of each segment, so it lines up
        # under the energy panel
        seg_mid = 0.5 * (neb.path_distances[:-1] + neb.path_distances[1:])
        ax_d = axes[1, col]
        ax_d.plot(seg_mid, neb.distances, 'o-', color='tab:green')
        ax_d.plot(neb.path_distances[crit],
                  np.interp(neb.path_distances[crit], seg_mid, neb.distances),
                  'r*', markersize=14)
        ax_d.set_xlabel('path distance')
        ax_d.set_ylabel('image spacing')

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print('wrote', filename)


if __name__ == '__main__':
    u, w = relaxed_pair('energy')
    plot(u, w, 'spring_force_ratio = 0.9 (energy-weighted)',
         'neb_2particles_spring_force_ratio_comparison.png')

    u, w = relaxed_pair('curvature')
    plot(u, w, 'spring_force_ratio = 0.9 (curvature-weighted)',
         'neb_2particles_curvature_weighted_comparison.png')
