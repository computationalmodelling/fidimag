from __future__ import print_function
import pytest

# FIDIMAG:
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, UniaxialAnisotropy
from fidimag.common.nebm_spherical import NEBM_Spherical
from fidimag.common.nebm_geodesic import NEBM_Geodesic
import numpy as np

# Material Parameters
# Parameters
A = 1e-12
Kx = 1e5
# Strong anisotropy
Ms = 3.8e5


"""
We will define two particles using a 4 sites mesh, letting the
sites in the middle as Ms = 0

"""


def two_part(pos):

    x = pos[0]

    if x > 6 or x < 3:
        return Ms
    else:
        return 0

# Finite differences mesh
mesh = CuboidMesh(nx=3,
                  ny=1,
                  nz=1,
                  dx=3, dy=3, dz=3,
                  unit_length=1e-9
                  )


# Simulation Function
def relax_neb(k, maxst, simname, init_im, interp, save_every=10000,
              coordinates='Cartesian'):
    """
    Execute a simulation with the NEB function of the FIDIMAG code, for an
    elongated particle (long cylinder)

    The simulations are made for a specific spring constant 'k' (a float),
    number of images 'init_im', interpolations between images 'interp'
    (an array) and a maximum of 'maxst' steps.
    'simname' is the name of the simulation, to distinguish the
    output files.

    --> vtks and npys are saved in files starting with the 'simname' string

    """

    # Prepare simulation
    # We define the cylinder with the Magnetisation function
    sim = Sim(mesh)
    sim.Ms = two_part

    # sim.add(UniformExchange(A=A))

    # Uniaxial anisotropy along x-axis
    sim.add(UniaxialAnisotropy(Kx, axis=(1, 0, 0)))

    # Define many initial states close to one extreme. We want to check
    # if the images in the last step, are placed mostly in equally positions
    init_images = init_im

    # Number of images between each state specified before (here we need only
    # two, one for the states between the initial and intermediate state
    # and another one for the images between the intermediate and final
    # states). Thus, the number of interpolations must always be
    # equal to 'the number of initial states specified', minus one.
    interpolations = interp

    if coordinates == 'Spherical':
        neb = NEBM_Spherical(sim,
                             init_images,
                             interpolations=interpolations,
                             spring_constant=k,
                             name=simname
                             )
    if coordinates == 'Geodesic':
        neb = NEBM_Geodesic(sim,
                            init_images,
                            interpolations=interpolations,
                            spring_constant=k,
                            name=simname,
                            integrator='sundials'
                            )

    neb.relax(max_iterations=2000,
              save_vtks_every=save_every,
              save_npys_every=save_every,
              stopping_dYdt=1e-4,
              dt=1e-6
              )

    # print(neb.G)
    # print(neb.tangents)
    # print(neb.spring_force)


def mid_m(pos):
    if pos[0] > 4:
        return (0.5, 0, 0.2)
    else:
        return (-0.5, 0, 0.2)


def test_energy_weighted_spring_force_clustering():
    """
    Check that the energy-weighted spring force
    (NEBM_Geodesic.spring_force_ratio > 0,
    NEBM_Geodesic.compute_energy_weighted_spring_lengths) produces a
    non-uniform, adaptively refined image distribution -- denser on the
    flanks approaching the critical (highest-energy) images and sparser
    in the flat basins -- compared to a plain, distance-only spaced band
    (spring_force_ratio = 0, the previous/default behaviour, which stays
    uniform in path distance).

    Note the tightest spacing is NOT expected exactly at a critical point
    itself: dE/d(path_distance) -> 0 there by definition, so the blended
    (path_distance, energy) metric locally reduces to the plain path
    distance right at the peak. The refinement shows up on the flanks,
    where the energy changes fastest -- i.e. gaps should grow with
    distance (in image index) from the nearest critical point, which is
    what's tested below.

    The energy-weighted band is obtained with a staged protocol rather
    than turning spring_force_ratio on from a fresh, uninterpolated band:

        1) relax normally (spring_force_ratio = 0) to convergence, with
           CVODE, which reliably converges for the plain, distance-based
           spring force;
        2) run a handful of extra steps with a climbing image on the
           highest-energy image found in (1), so a single image actually
           sits at the true saddle point (dE/d(path_distance) ~ 0
           *exactly* there), rather than two images straddling an
           unresolved peak;
        3) turn spring_force_ratio on and relax further, switching to the
           fixed-step, quick-min velocity-projection (VP) integrator
           (integrator='verlet') instead of CVODE.

    Step 3 needs the VP integrator rather than CVODE because
    compute_energy_weighted_spring_lengths refits a global spline (with a
    range-based renormalisation) on every RHS evaluation, which is not
    smooth enough for CVODE's adaptive step-size control to handle well:
    it was found to either falsely declare convergence (small step-size
    -> small max_dYdt, even while the actual spring force was still
    orders of magnitude away from zero) or, when forced to keep
    iterating, to diverge outright. The fixed-step VP integrator doesn't
    depend on RHS smoothness for stability and converges to the correct,
    physically consistent (near-zero max spring force) equilibrium
    instead. See also test_energy_barrier_2particles_verlet, which
    validates the VP integrator against the known CVODE barrier.

    Also saves a plot of energy and image spacing vs. path distance for
    both cases, so the clustering can be checked visually.
    """
    init_im = [(-1, 0, 0), mid_m, (1, 0, 0)]
    interp = [10, 10]

    def max_spring_force(neb):
        return np.max(np.linalg.norm(neb.spring_force.reshape(-1, 3), axis=1))

    def make_neb(simname):
        sim = Sim(mesh)
        sim.Ms = two_part
        sim.add(UniaxialAnisotropy(Kx, axis=(1, 0, 0)))
        return NEBM_Geodesic(sim, init_im, interpolations=interp,
                             spring_constant=1e4, name=simname,
                             integrator='sundials')

    # --- Baseline: plain, distance-only spring force (spring_force_ratio
    # = 0), CVODE converges this reliably on its own -----------------------
    neb_uniform = make_neb('neb_2particles_spacing_uniform')
    neb_uniform.relax(max_iterations=2000, save_vtks_every=5000,
                      save_npys_every=5000, stopping_dYdt=1e-4,
                      stopping_max_force=1e-2, dt=1e-6)

    # --- Energy-weighted band: staged protocol -----------------------------
    neb_weighted = make_neb('neb_2particles_spacing_weighted')

    # Stage 1: plain relaxation to convergence
    neb_weighted.relax(max_iterations=2000, save_vtks_every=5000,
                       save_npys_every=5000, stopping_dYdt=1e-4,
                       stopping_max_force=1e-2, dt=1e-6)

    # Stage 2: climbing image on the peak, a handful of extra steps to
    # pull it exactly onto the saddle point
    peak = int(np.argmax(neb_weighted.energies))
    neb_weighted.climbing_image = [peak]
    neb_weighted.initialise_integrator(integrator='sundials')
    neb_weighted.relax(max_iterations=30, save_vtks_every=50000,
                       save_npys_every=50000, stopping_dYdt=1e-6,
                       dt=1e-6, save_initial_state=False)

    # Stage 3: turn on the energy-weighted spring force and switch to the
    # fixed-step VP integrator (see docstring for why)
    del neb_weighted.climbing_image
    neb_weighted.spring_force_ratio = 0.9
    neb_weighted.initialise_integrator(integrator='verlet')
    neb_weighted.relax(max_iterations=20000, save_vtks_every=50000,
                       save_npys_every=50000, stopping_dYdt=1e-8,
                       stopping_max_force=1e-2, dt=1e-4,
                       save_initial_state=False)

    print('max|spring_force| (uniform) :', max_spring_force(neb_uniform))
    print('max|spring_force| (weighted):', max_spring_force(neb_weighted))
    # Both bands must have actually reached spring-force equilibrium
    # before drawing any conclusions from their spacing
    assert max_spring_force(neb_uniform) < 0.02
    assert max_spring_force(neb_weighted) < 0.02

    def spacing_variation(neb):
        """Coefficient of variation of the per-segment spacing: ~0 means
        perfectly uniform spacing, larger means adaptively refined."""
        return np.std(neb.distances) / np.mean(neb.distances)

    def slope_vs_gap_correlation(neb):
        """
        Correlation between each inner image's local energy slope
        |dE/d(path_distance)| (estimated by a central finite difference on
        the *relaxed* band, i.e. using the actual, possibly non-uniform,
        image placement) and the average gap of its two neighbouring
        segments. dE/d(path_distance) -> 0 at every critical point of the
        path (minima AND maxima, e.g. the two fixed end states, the
        interior minimum, and the energy maximum), so a NEGATIVE
        correlation here is the real signature of the energy-weighted
        spring force at work: images bunch up (small gap) wherever the
        energy changes fast, and spread out (large gap) near every flat/
        critical region, not only next to the highest-energy image.
        """
        E = neb.energies
        x = neb.path_distances
        slope = np.abs((E[2:] - E[:-2]) / (x[2:] - x[:-2]))  # inner images
        local_gap = 0.5 * (neb.distances[:-1] + neb.distances[1:])
        return np.corrcoef(slope, local_gap)[0, 1]

    cv_uniform = spacing_variation(neb_uniform)
    cv_weighted = spacing_variation(neb_weighted)
    corr_weighted = slope_vs_gap_correlation(neb_weighted)

    print('spacing coefficient of variation (uniform)  :', cv_uniform)
    print('spacing coefficient of variation (weighted) :', cv_weighted)
    print('|dE/dx| vs. gap correlation (weighted)      :', corr_weighted)

    # Plain path-distance spacing (spring_force_ratio = 0, the previous/
    # default behaviour) should stay essentially uniform
    assert cv_uniform < 0.01
    # The energy-weighted spring force should produce clearly non-uniform,
    # adaptively refined spacing
    assert cv_weighted > 0.3
    # ... and that refinement should track the local energy slope: images
    # bunch up (small gap) where the energy changes fast, and spread out
    # near every flat/critical region of the path
    assert corr_weighted < -0.3

    # --- Plot for visual verification ---------------------------------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def critical_image_indices(neb):
        """All local extrema along the band (both minima and maxima,
        including the two fixed end states), i.e. every image where
        dE/d(path_distance) ~ 0."""
        E = neb.energies
        idxs = [0, len(E) - 1]
        for i in range(1, len(E) - 1):
            if (E[i] - E[i - 1]) * (E[i + 1] - E[i]) <= 0:
                idxs.append(i)
        return np.array(sorted(set(idxs)))

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex='col')
    for col, (neb, title) in enumerate(
            [(neb_uniform, 'spring_force_ratio = 0 (uniform spacing)'),
             (neb_weighted, 'spring_force_ratio = 0.9 (energy-weighted)')]
            ):
        crit_idxs = critical_image_indices(neb)

        ax_E = axes[0, col]
        ax_E.plot(neb.path_distances, neb.energies / 1.602e-19, 'o-')
        ax_E.plot(neb.path_distances[crit_idxs],
                 neb.energies[crit_idxs] / 1.602e-19,
                 'r*', markersize=14,
                 label='critical images (dE/dx ~ 0)')
        ax_E.set_title(title)
        ax_E.set_ylabel('Energy (eV)')
        ax_E.legend()

        # Image spacing (self.distances) plotted at the midpoint path
        # distance of each segment, so it lines up under the energy panel
        seg_mid = 0.5 * (neb.path_distances[:-1] + neb.path_distances[1:])
        ax_d = axes[1, col]
        ax_d.plot(seg_mid, neb.distances, 'o-', color='tab:green')
        ax_d.plot(neb.path_distances[crit_idxs],
                 np.interp(neb.path_distances[crit_idxs], seg_mid,
                           neb.distances),
                 'r*', markersize=14)
        ax_d.set_xlabel('path distance')
        ax_d.set_ylabel('image spacing')

    plt.tight_layout()
    plt.savefig('neb_2particles_spring_force_ratio_comparison.png', dpi=150)
    plt.close(fig)


def test_energy_barrier_2particles():
    # Initial images: we set here a rotation interpolating
    init_im = [(-1, 0, 0), mid_m, (1, 0, 0)]
    interp = [6, 6]

    coord_list = ['Geodesic']
    barriers = []

    # Define different ks for multiple simulations
    # krange = ['1e8']

    for coordinates in coord_list:
        relax_neb(1e4, 2000,
                  'neb_2particles_k1e8_10-10int_{}'.format(coordinates),
                  init_im,
                  interp,
                  save_every=5000,
                  coordinates=coordinates
                  )

        _file = np.loadtxt('neb_2particles_k1e8_10-10int_{}_energy.ndt'.format(coordinates))
        barriers.append((np.max(_file[-1][1:]) - _file[-1][1]) / 1.602e-19)

        print('Energy barrier for {} is:'.format(coordinates), barriers[-1])
        assert np.abs(barriers[-1] - 0.016019) < 1e-5

    print(barriers)


def test_energy_barrier_2particles_verlet():
    """
    Same physical system and barrier as test_energy_barrier_2particles,
    but relaxed with the fixed-step, quick-min velocity-projection (VP)
    integrator (integrator='verlet') instead of CVODE, to validate that
    the VP rewrite in chain_method_integrators.py (see its docstring)
    converges to the same physics.
    """
    init_im = [(-1, 0, 0), mid_m, (1, 0, 0)]
    interp = [6, 6]

    sim = Sim(mesh)
    sim.Ms = two_part
    sim.add(UniaxialAnisotropy(Kx, axis=(1, 0, 0)))

    neb = NEBM_Geodesic(sim, init_im, interpolations=interp,
                        spring_constant=1e4, name='neb_2particles_verlet',
                        integrator='verlet')
    neb.integrator.stepsize = 1e-4
    neb.relax(max_iterations=5000, save_vtks_every=50000,
              save_npys_every=50000, stopping_dYdt=1e-8,
              stopping_max_force=1e-2, dt=1e-4)

    barrier = (np.max(neb.energies) - neb.energies[0]) / 1.602e-19
    print('Energy barrier for Geodesic+verlet is:', barrier)
    assert np.abs(barrier - 0.016019) < 1e-3


if __name__ == '__main__':
    test_energy_barrier_2particles()
