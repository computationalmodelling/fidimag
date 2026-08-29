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
                            integrator='cvode_bdf'
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


def _relaxed_band(simname, interp=(3, 3), max_iterations=50):
    """
    A small band, briefly relaxed, for the spring-length unit tests below.
    The weighting functions are pure functions of (path_distances, energies,
    gradientE, tangents), so they only need a band that is populated and
    self-consistent, not one that is converged.
    """
    sim = Sim(mesh)
    sim.Ms = two_part
    sim.add(UniaxialAnisotropy(Kx, axis=(1, 0, 0)))
    neb = NEBM_Geodesic(sim, [(-1, 0, 0), mid_m, (1, 0, 0)],
                        interpolations=list(interp), spring_constant=1e4,
                        name=simname, integrator='cvode_bdf')
    neb.relax(max_iterations=max_iterations, save_vtks_every=10 ** 9,
              save_npys_every=10 ** 9, stopping_dYdt=1e-6, dt=1e-6)
    neb.nebm_step(neb.band)      # refresh gradientE / tangents
    return neb


def test_energy_weighted_spring_lengths():
    """
    compute_energy_weighted_spring_lengths measures arc length in the
    normalised (path_distance, energy) plane, so the length of segment i is

        int_{x_i}^{x_i+1} sqrt( (r_x/range_x)^2 + (r_E E'(x)/range_E)^2 ) dx
        * range_x

    with r_E = spring_force_ratio and r_x = 1 - r_E. Check the implementation
    against a direct quadrature of that integrand, and check that it reduces
    to the plain geodesic spacing when the energy is given no weight.
    """
    from scipy.integrate import quad

    neb = _relaxed_band('neb_2particles_energy_weight_unit')

    def reference(ratio):
        x, _, spline = neb._energy_hermite_spline()
        dspline = spline.derivative()
        range_x = x[-1] - x[0]
        e_fine = spline(neb._fine_path_grid(x, 20))
        range_E = max(np.max(e_fine) - np.min(e_fine), 1e-30)
        integrand = lambda t: np.sqrt(((1.0 - ratio) / range_x) ** 2 +
                                      (ratio * dspline(t) / range_E) ** 2) * range_x
        return np.array([quad(integrand, x[i], x[i + 1], limit=200)[0]
                         for i in range(neb.n_images - 1)])

    for ratio in [0.1, 0.5, 0.9, 1.0]:
        neb.spring_force_ratio = ratio
        lengths = neb.compute_energy_weighted_spring_lengths()
        # One length per segment, same layout as self.distances, so it can
        # be passed straight to nebm_clib.compute_spring_force
        assert lengths.shape == neb.distances.shape
        assert np.all(lengths > 0)
        # 20 sub-intervals per segment, so the piecewise-linear sum sits a
        # little under the exact integral
        assert np.allclose(lengths, reference(ratio), rtol=1e-3)

    # With no weight on the energy axis the metric collapses to the plain
    # geodesic path distance, i.e. the spring force keeps its old behaviour
    neb.spring_force_ratio = 1e-12
    assert np.allclose(neb.compute_energy_weighted_spring_lengths(),
                       neb.distances, rtol=1e-9)


def test_curvature_weighted_spring_lengths():
    """
    compute_curvature_weighted_spring_lengths stretches each segment by

        w(x) = (1 - r_C) + r_C |E''(x)| / max|E''|

    so that regions of high energy curvature are treated as longer and end up
    more finely sampled. Check it against a quadrature of w, and check that
    it too reduces to the plain geodesic spacing at zero weight.
    """
    from scipy.integrate import quad

    neb = _relaxed_band('neb_2particles_curvature_weight_unit')

    def reference(ratio):
        x, _, spline = neb._energy_hermite_spline()
        curvature = spline.derivative(nu=2)
        # The normalisation is the largest |E''| over the sub-segment
        # midpoints, so the reference has to use that same sampling
        x_fine = neb._fine_path_grid(x, 20)
        x_mid = 0.5 * (x_fine[:-1] + x_fine[1:])
        scale = max(np.max(np.abs(curvature(x_mid))), 1e-30)
        integrand = lambda t: (1.0 - ratio) + ratio * abs(curvature(t)) / scale
        return np.array([quad(integrand, x[i], x[i + 1], limit=200)[0]
                         for i in range(neb.n_images - 1)])

    neb.spring_weighting = 'curvature'
    for ratio in [0.1, 0.5, 0.9, 1.0]:
        neb.spring_force_ratio = ratio
        lengths = neb.compute_curvature_weighted_spring_lengths()
        assert lengths.shape == neb.distances.shape
        assert np.all(lengths > 0)
        # |E''| is only piecewise linear across the C1 spline's nodes, so the
        # midpoint rule is a little coarser here than for the energy weighting
        assert np.allclose(lengths, reference(ratio), rtol=1e-2)

    neb.spring_force_ratio = 1e-12
    assert np.allclose(neb.compute_curvature_weighted_spring_lengths(),
                       neb.distances, rtol=1e-9)


@pytest.mark.slow
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
    where the energy changes fastest. What is asserted here is only that
    the spacing becomes markedly non-uniform; the weighted lengths
    themselves are checked directly against a quadrature in
    test_energy_weighted_spring_lengths.

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

    sandbox/nebm_spring_weighting_plots.py draws the energy and image
    spacing against path distance for both cases, if the clustering needs
    to be inspected visually.
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
                             integrator='cvode_bdf')

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
    neb_weighted.initialise_integrator(integrator='cvode_bdf')
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

    cv_uniform = spacing_variation(neb_uniform)
    cv_weighted = spacing_variation(neb_weighted)

    print('spacing coefficient of variation (uniform)  :', cv_uniform)
    print('spacing coefficient of variation (weighted) :', cv_weighted)

    # Plain path-distance spacing (spring_force_ratio = 0, the previous/
    # default behaviour) should stay essentially uniform
    assert cv_uniform < 0.01
    # The energy-weighted spring force should produce clearly non-uniform,
    # adaptively refined spacing
    assert cv_weighted > 0.3

@pytest.mark.slow
def test_curvature_weighted_spring_force_clustering():
    """
    Check that the curvature-weighted spring force
    (NEBM_Geodesic.spring_weighting = 'curvature',
    NEBM_Geodesic.compute_curvature_weighted_spring_lengths) does the
    opposite of the energy-weighted one
    (test_energy_weighted_spring_force_clustering): it concentrates images
    around every sharp critical point of the path (minima AND maxima
    alike, since |d^2E/d(path_distance)^2| is large at both), and spreads
    them out on the smoothly-curving flanks in between -- see that
    method's docstring for why this is an ad hoc heuristic rather than
    something derived from the GNEB literature, unlike the energy-weighted
    default.

    Uses the same staged protocol as
    test_energy_weighted_spring_force_clustering (plain CVODE relax ->
    climbing image on the saddle -> VP + curvature weighting), for the
    same reason: compute_curvature_weighted_spring_lengths also refits a
    global spline every RHS call, which CVODE's adaptive step-size control
    doesn't handle well.

    sandbox/nebm_spring_weighting_plots.py draws the energy and image
    spacing against path distance for both cases, if the clustering needs
    to be inspected visually.
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
                             integrator='cvode_bdf')

    # --- Baseline: plain, distance-only spring force -----------------------
    neb_uniform = make_neb('neb_2particles_curvature_uniform')
    neb_uniform.relax(max_iterations=2000, save_vtks_every=5000,
                      save_npys_every=5000, stopping_dYdt=1e-4,
                      stopping_max_force=1e-2, dt=1e-6)

    # --- Curvature-weighted band: staged protocol ---------------------------
    neb_curv = make_neb('neb_2particles_curvature_weighted')

    # Stage 1: plain relaxation to convergence
    neb_curv.relax(max_iterations=2000, save_vtks_every=5000,
                   save_npys_every=5000, stopping_dYdt=1e-4,
                   stopping_max_force=1e-2, dt=1e-6)

    # Stage 2: climbing image on the peak
    peak = int(np.argmax(neb_curv.energies))
    neb_curv.climbing_image = [peak]
    neb_curv.initialise_integrator(integrator='cvode_bdf')
    neb_curv.relax(max_iterations=30, save_vtks_every=50000,
                   save_npys_every=50000, stopping_dYdt=1e-6,
                   dt=1e-6, save_initial_state=False)

    # Stage 3: turn on curvature weighting, switch to VP
    del neb_curv.climbing_image
    neb_curv.spring_force_ratio = 0.9
    neb_curv.spring_weighting = 'curvature'
    neb_curv.initialise_integrator(integrator='verlet')
    neb_curv.relax(max_iterations=20000, save_vtks_every=50000,
                   save_npys_every=50000, stopping_dYdt=1e-8,
                   stopping_max_force=1e-2, dt=1e-4,
                   save_initial_state=False)

    print('max|spring_force| (uniform)  :', max_spring_force(neb_uniform))
    print('max|spring_force| (curvature):', max_spring_force(neb_curv))
    assert max_spring_force(neb_uniform) < 0.02
    assert max_spring_force(neb_curv) < 0.02

    def spacing_variation(neb):
        return np.std(neb.distances) / np.mean(neb.distances)

    def curvature_vs_gap_correlation(neb):
        """
        Correlation between each inner image's local |d^2E/d(path_distance)^2|
        (central finite difference of dE/dx, itself a central finite
        difference of E, on the relaxed band) and the average gap of its
        two neighbouring segments. A NEGATIVE correlation is the
        signature of curvature-weighted spacing: images bunch up (small
        gap) at sharp critical points (large curvature), and spread out
        on the smoothly-curving flanks (small curvature) in between.
        """
        E = neb.energies
        x = neb.path_distances
        dE_dx = np.gradient(E, x)
        d2E_dx2 = np.abs(np.gradient(dE_dx, x))[1:-1]
        local_gap = 0.5 * (neb.distances[:-1] + neb.distances[1:])
        return np.corrcoef(d2E_dx2, local_gap)[0, 1]

    cv_uniform = spacing_variation(neb_uniform)
    cv_curv = spacing_variation(neb_curv)
    corr_curv = curvature_vs_gap_correlation(neb_curv)

    print('spacing coefficient of variation (uniform)  :', cv_uniform)
    print('spacing coefficient of variation (curvature):', cv_curv)
    print('|d2E/dx2| vs. gap correlation (curvature)   :', corr_curv)

    assert cv_uniform < 0.01
    assert cv_curv > 0.15
    assert corr_curv < -0.3

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
