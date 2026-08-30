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


def _run_comparison(mesh_kwargs, m_init, demag_kwargs, tol, reference=None):
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

    # atol covers the near-zero-energy uniform states, where a relative
    # comparison would be meaningless.
    assert np.isclose(fmm_energy, ref_energy, rtol=tol, atol=1e-9), (
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


if __name__ == '__main__':
    test_cuboid_demag_fmm_2D()
    test_cuboid_demag_fmm_3D()
    test_cuboid_demag_bh_3D()
    test_cuboid_demag_fmm_vs_fft_2D()
    test_cuboid_demag_fmm_vs_fft_3D()
    test_cuboid_demag_fmm_vs_fft_theta_2D()
