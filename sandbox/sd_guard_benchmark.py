"""Effective field evaluations of the steepest descent with the energy guard

Reproduces the energy guard numbers of
doc/physics_num_methods/energy_minimisation.rst.

    OMP_NUM_THREADS=4 python sandbox/sd_guard_benchmark.py
"""
import sys
import fidimag

sys.path.insert(0, 'tests')
from test_steepest_descent import _setup_1D_DW, _dw_MAE  # noqa: E402

NEVER = 10 ** 9  # no data files


def dw(tmax, guard):
    sim = _setup_1D_DW(tmax=tmax)
    sim.driver.minimise(stopping_dm=1e-9, max_steps=20000, printing=False,
                        save_data_steps=NEVER, energy_guard=guard)
    return f'evals {sim.driver.nEval:5d}  err {_dw_MAE(sim):.4f}'


def sp4(tmax, guard):
    Ms = 8e5
    mesh = fidimag.common.CuboidMesh(nx=100, ny=25, nz=1, dx=5, dy=5, dz=3, unit_length=1e-9)
    sim = fidimag.micro.Sim(mesh, name='bench_sd_sp4', driver='steepest_descent')
    sim.set_Ms(Ms)
    sim.add(fidimag.micro.UniformExchange(1.3e-11))
    sim.add(fidimag.micro.Demag())
    sim.set_m((1, 1, 1))
    sim.driver.energyScale = fidimag.common.constant.mu_0 * Ms ** 2 * 0.5 * 500 * 125 * 3 * 1e-27
    sim.driver.tmax = tmax
    sim.driver.minimise(stopping_dm=0.0, stopping_torque=1e-6, max_steps=20000,
                        printing=False, save_data_steps=NEVER, energy_guard=guard)
    return f'evals {sim.driver.nEval:5d}  torque {sim.driver.max_torque():.1e}'


if __name__ == '__main__':
    for tmax in (0.1, 1, 3, 10):
        print(f'1D wall  tmax {tmax:4g}  guard off: {dw(tmax, False)}   '
              f'guard on: {dw(tmax, True)}', flush=True)
    for tmax in (3, 10):
        print(f'SP4 2500 tmax {tmax:4g}  guard on: {sp4(tmax, True)}', flush=True)
