"""
Testing of the Fast Multipole Method (FMM) and Barnes-Hut (BH) approaches
to compute the demagnetising (dipolar) field for atomistic simulations.

Both approximate methods are compared against the brute-force DemagFull
calculation, which sums the dipolar contribution of every pair of spins and
is therefore used here as the ground truth (the same role DemagFull plays in
test_demag_libraries.py for the FFT-based Demag).

DemagFMM is also compared directly against the FFT-based Demag class, since
that is the interaction most atomistic simulations actually use day to day --
DemagFull is brute force and only practical at the small system sizes used
here, whereas Demag and DemagFMM are the two options meant for production
use. At theta=0.0 the FMM/BH tree traversal falls back to an exact pairwise
sum (no multipole approximation), so it should agree with the FFT result to
near machine precision, the same way it agrees with DemagFull. A further test
at theta > 0 exercises the actual multipole approximation, at a looser
tolerance appropriate to that order/theta combination.

Both the demag field and the total demag energy are compared. DemagFMM
computes the energy lazily in its compute_energy method (like the FFT Demag
class), using the same convention as DemagFull.
"""
import numpy as np
import pytest

from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic import Demag, DemagFull
import fidimag.common.constant as const

# The FMM/BH extension (fidimag.extensions.fmm) is an optional C++ module.
# Skip the whole module gracefully if it (or DemagFMM) was not built.
demag_mod = pytest.importorskip(
    "fidimag.atomistic.demag", reason="FMM extension not built"
)
if not hasattr(demag_mod, "DemagFMM"):
    pytest.skip("DemagFMM not available in this build", allow_module_level=True)
DemagFMM = demag_mod.DemagFMM


def m_init_2Dvortex(pos, centre):
    """A hedgehog-like vortex profile, reused from test_demag_libraries.py."""
    x, y = pos[0] - centre[0], pos[1] - centre[1]

    if np.sqrt(x ** 2 + y ** 2) <= 2.5:
        r = (x ** 2 + y ** 2) ** 0.5
        phi = np.arctan2(y, x)
        k = np.pi / 2.5
        return (np.sin(k * r) * np.cos(phi),
                np.sin(k * r) * np.sin(phi),
                np.cos(k * r))
    else:
        return (0, 0, -1)


def _relative_l2_error(approx, exact):
    """Relative L2 norm of the field difference."""
    return np.linalg.norm(approx - exact) / np.linalg.norm(exact)


def _demag_field_and_energy(mesh_kwargs, m_init, interaction):
    """
    Build a Sim on a fresh mesh, add `interaction`, and return its demag
    field and total energy (in meV).
    """
    mesh = CuboidMesh(**mesh_kwargs)
    sim = Sim(mesh)
    sim.mu_s = 2 * const.mu_B
    sim.set_m(m_init)
    sim.add(interaction)

    field = sim.get_interaction(interaction.name).compute_field().copy()
    energy = sim.compute_energy() / const.meV
    return field, energy


def _run_comparison(mesh_kwargs, m_init, demag_kwargs, tol, reference=None,
                     energy_atol=1e-9):
    """
    Compare the DemagFMM (fmm or bh) field and energy against a reference
    interaction (DemagFull, the brute-force sum, by default) on identical
    simulations.
    """
    if reference is None:
        reference = DemagFull()
    ref_name = type(reference).__name__
    kind = demag_kwargs.get('type', 'fmm')

    ref_field, ref_energy = _demag_field_and_energy(
        mesh_kwargs, m_init, reference)
    fmm_field, fmm_energy = _demag_field_and_energy(
        mesh_kwargs, m_init, DemagFMM(**demag_kwargs))

    field_err = _relative_l2_error(fmm_field, ref_field)
    assert field_err < tol, (
        "DemagFMM (%s) field differs from %s: relative L2 error "
        "%.3e >= tol %.3e" % (kind, ref_name, field_err, tol)
    )

    # energy_atol covers the near-zero-energy uniform states, where a
    # relative comparison would be meaningless. The default is tuned for the
    # small (N=6) meshes most tests here use; larger meshes accumulate more
    # absolute rounding even when, as for a uniform state's demag energy,
    # the true value is exactly zero by symmetry, so need a looser one.
    assert np.isclose(fmm_energy, ref_energy, rtol=tol, atol=energy_atol), (
        "DemagFMM (%s) energy %.6e differs from %s %.6e"
        % (kind, fmm_energy, ref_name, ref_energy)
    )


def test_cuboid_demag_fmm_2D():
    """
    FMM demag field on a 2D cuboid mesh (vortex state), compared against the
    brute-force DemagFull calculation.
    """
    N = 15
    a = 0.4
    mesh_kwargs = dict(dx=a, dy=a, dz=a, nx=N, ny=N, nz=1, unit_length=1e-9)
    xc = yc = N * a * 0.5

    _run_comparison(
        mesh_kwargs,
        lambda pos: m_init_2Dvortex(pos, (xc, yc)),
        demag_kwargs=dict(order=8, ncrit=48, theta=0.0, type='fmm'),
        tol=1e-6,
    )


def test_cuboid_demag_fmm_3D():
    """
    FMM demag field on a 3D cuboid block (uniform state), compared against the
    brute-force DemagFull calculation. A 3D block exercises the full multipole
    hierarchy.
    """
    N = 6
    a = 0.4
    mesh_kwargs = dict(dx=a, dy=a, dz=a, nx=N, ny=N, nz=N, unit_length=1e-9)

    _run_comparison(
        mesh_kwargs,
        (0, 0, 1),
        demag_kwargs=dict(order=8, ncrit=48, theta=0.0, type='fmm'),
        tol=1e-6,
    )


def test_cuboid_demag_bh_3D():
    """
    Barnes-Hut demag field on a 3D cuboid block (uniform state), compared
    against the brute-force DemagFull calculation.
    """
    N = 6
    a = 0.4
    mesh_kwargs = dict(dx=a, dy=a, dz=a, nx=N, ny=N, nz=N, unit_length=1e-9)

    _run_comparison(
        mesh_kwargs,
        (0, 0, 1),
        demag_kwargs=dict(order=8, ncrit=48, theta=0.0, type='bh'),
        tol=1e-6,
    )


def test_cuboid_demag_fmm_vs_fft_2D():
    """
    FMM demag field (theta=0.0, exact tree evaluation) compared directly
    against the FFT-based Demag class, rather than the brute-force
    DemagFull. Demag is what atomistic simulations actually use day to
    day, so this checks the two production paths agree with each other
    and not just each against the (impractical at scale) ground truth.
    """
    N = 15
    a = 0.4
    mesh_kwargs = dict(dx=a, dy=a, dz=a, nx=N, ny=N, nz=1, unit_length=1e-9)
    xc = yc = N * a * 0.5

    _run_comparison(
        mesh_kwargs,
        lambda pos: m_init_2Dvortex(pos, (xc, yc)),
        demag_kwargs=dict(order=8, ncrit=48, theta=0.0, type='fmm'),
        tol=1e-6,
        reference=Demag(),
    )


def test_cuboid_demag_fmm_vs_fft_3D():
    """
    Same as test_cuboid_demag_fmm_vs_fft_2D for a 3D cuboid block.
    """
    N = 6
    a = 0.4
    mesh_kwargs = dict(dx=a, dy=a, dz=a, nx=N, ny=N, nz=N, unit_length=1e-9)

    _run_comparison(
        mesh_kwargs,
        (0, 0, 1),
        demag_kwargs=dict(order=8, ncrit=48, theta=0.0, type='fmm'),
        tol=1e-6,
        reference=Demag(),
    )


def test_cuboid_demag_fmm_vs_fft_theta_2D():
    """
    Same comparison with theta > 0, which actually exercises the multipole
    acceptance criterion -- every other test in this module uses theta=0.0,
    where the tree traversal falls back to an exact pairwise sum and never
    approximates a cell as a single expansion. The tolerance is looser to
    match the accuracy expected at this order/theta/ncrit combination, not
    because the comparison itself is any less exact.
    """
    N = 15
    a = 0.4
    mesh_kwargs = dict(dx=a, dy=a, dz=a, nx=N, ny=N, nz=1, unit_length=1e-9)
    xc = yc = N * a * 0.5

    _run_comparison(
        mesh_kwargs,
        lambda pos: m_init_2Dvortex(pos, (xc, yc)),
        demag_kwargs=dict(order=8, ncrit=8, theta=0.2, type='fmm'),
        tol=1e-3,
        reference=Demag(),
    )


def test_cuboid_demag_fmm_excludes_zero_mu_s_sites():
    """
    Sites with mu_s == 0 carry no moment, so DemagFMM should exclude them
    from its tree entirely: on the active (mu_s != 0) sites the field
    should still agree with the brute-force DemagFull, and on the inactive
    (mu_s == 0) sites DemagFMM should return exactly zero, since it no
    longer evaluates a field there at all.

    DemagFull is not a useful reference on the inactive sites here: it is
    an evaluation point like any other to DemagFull, which still sums the
    (nonzero) stray field reaching it from the active region, so DemagFull
    and DemagFMM are expected to disagree there. Only the active sites, and
    the exact-zero claim on the inactive ones, are checked.
    """
    N = 6
    a = 0.4
    mesh = CuboidMesh(dx=a, dy=a, dz=a, nx=N, ny=N, nz=N, unit_length=1e-9)

    mu_s = np.full(mesh.n, 2 * const.mu_B)
    coords = mesh.coordinates
    inactive = coords[:, 0] < (coords[:, 0].max() / 2)
    mu_s[inactive] = 0.0

    def field_for(interaction):
        sim = Sim(mesh)
        sim.mu_s = mu_s.copy()
        sim.set_m((0, 0, 1))
        sim.add(interaction)
        return sim.get_interaction(interaction.name).compute_field().copy()

    fmm_field = field_for(DemagFMM(order=8, ncrit=8, theta=0.0))
    full_field = field_for(DemagFull())

    active3 = np.repeat(~inactive, 3)
    inactive3 = np.repeat(inactive, 3)

    tol = 1e-6
    field_err = _relative_l2_error(fmm_field[active3], full_field[active3])
    assert field_err < tol, (
        "DemagFMM field on active (mu_s != 0) sites differs from DemagFull: "
        "relative L2 error %.3e >= tol %.3e" % (field_err, tol)
    )

    assert np.all(fmm_field[inactive3] == 0.0), (
        "DemagFMM should return exactly zero field at mu_s == 0 sites"
    )


def test_cuboid_demag_fmm_vs_fft_theta_3D_uniform():
    """
    Regression test for a real bug: DemagFMM(order=8, ...) used to silently
    return a field missing its entire far-field (M2L) contribution, because
    fidimag/atomistic/fmmlib/operators.cpp was only generated through order
    7 (FMMGEN_MAXORDER=8 is an *exclusive* bound) while an off-by-one in
    fmm.pyx's validation ("if order > MAXORDER") let order == MAXORDER
    through instead of rejecting it. The per-order dispatch switches in
    operators.cpp have no case for an ungenerated order, so every M2L pair
    silently contributed nothing - wrong, but not obviously so, since nearby
    (P2P) contributions are still exact.

    Every other theta > 0 test in this module uses a 2D vortex state, where
    the far field is small enough that a missing contribution stayed under
    that test's looser tolerance by chance. A solid, uniformly magnetised 3D
    block is far more exposed: its far field is the dominant part of the
    total, not a small correction, so this combination - 3D, uniform state,
    theta > 0, order == the previous ungenerated boundary - is the one this
    test exists to cover.
    """
    N = 15
    a = 0.4
    mesh_kwargs = dict(dx=a, dy=a, dz=a, nx=N, ny=N, nz=N, unit_length=1e-9)

    _run_comparison(
        mesh_kwargs,
        (0, 0, 1),
        demag_kwargs=dict(order=8, ncrit=8, theta=0.2, type='fmm'),
        tol=1e-5,
        reference=Demag(),
        # A uniform state's demag energy is exactly zero by symmetry, so
        # Demag's reference value is itself just rounding noise (~1e-16);
        # the default energy_atol (tuned for the smaller N=6 meshes most
        # tests here use) is tighter than DemagFMM's own rounding at N=15.
        energy_atol=1e-6,
    )


def test_demag_fmm_order_validation():
    """
    order must be in [fmm.MINORDER, fmm.MAXORDER): outside that range,
    operators.cpp has no generated code for it, and the per-order dispatch
    functions (S2M/M2M/M2L/L2L/L2P) silently do nothing rather than raising.
    See test_cuboid_demag_fmm_vs_fft_theta_3D_uniform for what that silent
    failure looks like when it is not caught here.
    """
    import fidimag.extensions.fmm as fmm_ext

    for bad_order in (fmm_ext.MINORDER - 1, fmm_ext.MAXORDER, fmm_ext.MAXORDER + 5):
        try:
            DemagFMM(order=bad_order, ncrit=8, theta=0.2)
        except AssertionError:
            pass
        else:
            raise AssertionError(
                "DemagFMM(order=%d) should have raised: valid range is "
                "[%d, %d)" % (bad_order, fmm_ext.MINORDER, fmm_ext.MAXORDER)
            )


def test_demag_fmm_2D_uses_planar_variant():
    """
    A 2D mesh (nz == 1) should route DemagFMM through fmm.FMM2D (fmmgen's
    planar operator variant, generate_code(..., planar=True)) rather than
    the general 3D fmm.FMM used for everything else, with no separate
    argument needed to ask for it.
    """
    import fidimag.extensions.fmm as fmm_ext

    N = 15
    a = 0.4
    mesh_2d = CuboidMesh(dx=a, dy=a, dz=a, nx=N, ny=N, nz=1, unit_length=1e-9)
    mesh_3d = CuboidMesh(dx=a, dy=a, dz=a, nx=6, ny=6, nz=6, unit_length=1e-9)

    def build(mesh):
        sim = Sim(mesh)
        sim.mu_s = 2 * const.mu_B
        sim.set_m((0, 0, 1))
        demag = DemagFMM(order=8, ncrit=8, theta=0.0)
        sim.add(demag)
        return demag

    demag_2d = build(mesh_2d)
    demag_3d = build(mesh_3d)

    assert demag_2d.is_2d is True
    assert isinstance(demag_2d.fmm, fmm_ext.FMM2D)
    assert demag_3d.is_2d is False
    assert isinstance(demag_3d.fmm, fmm_ext.FMM)


def test_demag_fmm_2D_planar_matches_full():
    """
    fmm.FMM2D's field should agree with the brute-force DemagFull on the
    same 2D vortex state used elsewhere in this module, the same way the
    general fmm.FMM path already does (test_cuboid_demag_fmm_vs_fft_2D).
    """
    N = 15
    a = 0.4
    mesh_kwargs = dict(dx=a, dy=a, dz=a, nx=N, ny=N, nz=1, unit_length=1e-9)
    xc = yc = N * a * 0.5

    _run_comparison(
        mesh_kwargs,
        lambda pos: m_init_2Dvortex(pos, (xc, yc)),
        demag_kwargs=dict(order=8, ncrit=8, theta=0.3, type='fmm'),
        tol=1e-6,
        reference=Demag(),
    )


def test_demag_fmm_2D_and_3D_interleaved_do_not_interfere():
    """
    Regression test for a real hazard the planar variant introduces:
    fmm_select()'s target is a single process-global pointer that
    compute_field_fmm/bh/exact read at CALL time, not at tree-build time
    (see the note in fmm.pyx). A 2D DemagFMM (fmm.FMM2D, always the planar
    variant) and a 3D one (fmm.FMM, full or compressed) can both be alive
    at once - unremarkable in a test suite, or any script/notebook that
    runs more than one simulation - and without re-selecting the right
    variant on every call, whichever was selected most recently would
    silently apply to both. FMM/FMM2D.compute_field() both re-select their
    own variant immediately before every call for exactly this reason; this
    test interleaves calls to catch a regression of that.
    """
    a = 0.4
    mesh_2d = CuboidMesh(dx=a, dy=a, dz=a, nx=15, ny=15, nz=1, unit_length=1e-9)
    mesh_3d = CuboidMesh(dx=a, dy=a, dz=a, nx=6, ny=6, nz=6, unit_length=1e-9)

    def build(mesh, interaction):
        sim = Sim(mesh)
        sim.mu_s = 2 * const.mu_B
        sim.set_m((0, 0, 1))
        sim.add(interaction)
        return sim

    ref_2d = build(mesh_2d, DemagFull()).get_interaction(
        "DemagFull").compute_field().copy()
    ref_3d = build(mesh_3d, DemagFull()).get_interaction(
        "DemagFull").compute_field().copy()

    sim_2d = build(mesh_2d, DemagFMM(order=8, ncrit=8, theta=0.3))
    sim_3d = build(mesh_3d, DemagFMM(order=8, ncrit=8, theta=0.3))

    tol = 1e-3
    for _ in range(10):
        field_2d = sim_2d.get_interaction("DemagFMM").compute_field()
        field_3d = sim_3d.get_interaction("DemagFMM").compute_field()

        err_2d = _relative_l2_error(field_2d, ref_2d)
        err_3d = _relative_l2_error(field_3d, ref_3d)
        assert err_2d < tol, (
            "2D DemagFMM field diverged after an interleaved 3D call: "
            "relative L2 error %.3e >= tol %.3e" % (err_2d, tol)
        )
        assert err_3d < tol, (
            "3D DemagFMM field diverged after an interleaved 2D call: "
            "relative L2 error %.3e >= tol %.3e" % (err_3d, tol)
        )


if __name__ == '__main__':
    test_cuboid_demag_fmm_2D()
    test_cuboid_demag_fmm_3D()
    test_cuboid_demag_bh_3D()
    test_cuboid_demag_fmm_vs_fft_2D()
    test_cuboid_demag_fmm_vs_fft_3D()
    test_cuboid_demag_fmm_vs_fft_theta_2D()
    test_cuboid_demag_fmm_excludes_zero_mu_s_sites()
    test_cuboid_demag_fmm_vs_fft_theta_3D_uniform()
    test_demag_fmm_order_validation()
    test_demag_fmm_2D_uses_planar_variant()
    test_demag_fmm_2D_planar_matches_full()
    test_demag_fmm_2D_and_3D_interleaved_do_not_interfere()
