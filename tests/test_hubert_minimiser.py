import fidimag
import pytest
import fidimag.common.constant as C
import numpy as np


Ms = 0.86e6
A = 13e-12
Ku = 0.4e6


def _setup_1D_DW():
    """1D micromagnetic sample with a domain wall as the initial state"""

    nx, ny, nz = 100, 1, 1
    dx, dy, dz = 1, 1, 1

    mesh = fidimag.common.CuboidMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz,
                                     periodicity=(False, False, False),
                                     unit_length=1e-9)

    sim = fidimag.micro.Sim(mesh, name='1Dmicro', driver='hubert_minimiser')

    # Define the magnetisation
    sim.set_Ms(Ms)

    # Add the magnetic interactions
    sim.add(fidimag.micro.UniformExchange(A))
    sim.add(fidimag.micro.UniaxialAnisotropy(Ku, axis=(0, 0, 1)))

    xs = mesh.coordinates[:, 0]
    centre_x = (xs.max() + xs.min()) * 0.5 + xs.min()

    def m_initial(r):
        x, y, z = r[0], r[1], r[2]
        if x < centre_x:
            return (0, 0.1, -.9)
        else:
            return (0, 0.1, .9)
    sim.set_m(m_initial)

    mesh_vol = mesh.n * mesh.dx * mesh.dy * mesh.dz * 1e-27
    Kd = C.mu_0 * (Ms ** 2) * 0.5 * mesh_vol
    # Scale the energy for the minimiser
    sim.driver.energyScale = Kd

    return sim, A, Ku


def _dw_MAE(sim, A, Ku):
    """Mean absolute error against the analytical wall, within the DW width

    The maximum "deviation" a spin might have from the curve is 2, since the
    spin-z is in the [-1, 1] range.
    """

    mz = sim.spin.reshape(-1, 3)[:, 2]

    def mz_dw_analyt(x, xc):
        """Analytical DW function using DW width param"""
        # DW width from theory
        deltaB = np.sqrt(A / Ku)
        return np.tanh((x - xc) / deltaB)

    # Analytical model in nm units
    x = sim.mesh.coordinates[:, 0]
    y = mz_dw_analyt(x * 1e-9, 50 * 1e-9)

    deltaB = np.sqrt(A / Ku) * 1e9
    xc = 50  # sample centre
    ftr = np.logical_and(x <= xc + deltaB, x >= xc - deltaB)

    return np.mean(np.abs(y[ftr] - mz[ftr]))


def test_hubert_minimiser_1D_DW():

    sim, A, Ku = _setup_1D_DW()

    sim.driver.minimise(stopping_dE=1e-6, maxCreep=6,
                        eta_scale=1e-6, log_steps=10)

    MAE_dw = _dw_MAE(sim, A, Ku)
    print('Domain wall MAE', MAE_dw)
    # This is at least 1% the max deviation
    assert MAE_dw < 0.02


def test_hubert_minimiser_1D_DW_BB():
    """Same 1D domain wall, relaxed with the Barzilai-Borwein step control

    The BB step length is inferred from the curvature along the previous step,
    so no `eta_scale` is given here: the point of the test is that the same
    minimum is reached without that hand-tuned unit conversion, and in fewer
    energy evaluations than the creep algorithm needs.
    """

    sim, A, Ku = _setup_1D_DW()
    sim.driver.minimise(stepControl='BB', stopping_dE=1e-14,
                        mXgradE_tol=1e-3, log_steps=50)
    BB_steps = sim.driver.step
    MAE_BB = _dw_MAE(sim, A, Ku)
    print('Domain wall MAE (BB)', MAE_BB, 'steps', BB_steps)
    assert MAE_BB < 0.02

    # The creep algorithm, with an eta_scale tuned for these field units,
    # solves the same problem
    sim_h, _, _ = _setup_1D_DW()
    sim_h.driver.minimise(stopping_dE=1e-14, mXgradE_tol=1e-3, maxCreep=6,
                          eta_scale=1e-6, log_steps=50)
    assert _dw_MAE(sim_h, A, Ku) < 0.02
    # Same minimum ...
    assert abs(sim.driver.totalE - sim_h.driver.totalE) < 1e-6
    # ... at a fraction of the number of effective field evaluations
    assert BB_steps < 0.5 * sim_h.driver.step


def test_hubert_minimiser_BB_unknown_step_control():

    sim, _, _ = _setup_1D_DW()
    with pytest.raises(ValueError):
        sim.driver.minimise(stepControl='not-a-step-control')


if __name__ == "__main__":
    test_hubert_minimiser_1D_DW()
    test_hubert_minimiser_1D_DW_BB()
