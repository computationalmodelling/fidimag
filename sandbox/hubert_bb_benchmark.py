"""Effective field evaluations of the BB path of the Hubert minimiser

Reproduces the numbers of doc/physics_num_methods/energy_minimisation.rst:
the creep vs BB table and the standard problem 4 tables against OOMMF.

    OMP_NUM_THREADS=4 python sandbox/hubert_bb_benchmark.py
"""
import sys
import numpy as np
import fidimag
import fidimag.common.constant as C

NEVER = 10 ** 9  # no data files


def record_torque(driver):
    """Log (evaluations, max torque) every time the tangential gradient is taken"""
    history = []
    project = driver._project_gradient

    def wrapped(out=None):
        g = project(out=out)
        history.append((driver.step, driver.mXgradE.max()))
        return g

    driver._project_gradient = wrapped
    return history


def dw_1d():
    sys.path.insert(0, 'tests')
    from test_hubert_minimiser import _setup_1D_DW, _dw_MAE
    sim, A, Ku = _setup_1D_DW()
    r = sim.driver.minimise(stepControl='BB', stopping_dE=1e-14, mXgradE_tol=1e-3,
                            save_data_steps=NEVER)
    return r, f'MAE {_dw_MAE(sim, A, Ku):.2e}'


def skyrmion_demag(tol=1e-1):
    """Nanodisk skyrmion of hubert_minimiser_atomistic.ipynb, plus demag"""
    radius, cell = 50, 2
    A, D, Ku, Ms = 13e-12, 3e-3, 0.4e6, 0.86e6
    n = int(2 * radius / cell)
    mesh = fidimag.common.CuboidMesh(dx=cell, dy=cell, dz=cell, nx=n, ny=n, nz=1,
                                     unit_length=1e-9)
    sim = fidimag.micro.Sim(mesh, name='bench_sk', driver='hubert_minimiser')

    def Ms_fun(pos):
        x, y = np.array(pos)[:2] - radius
        return Ms if (x ** 2 + y ** 2) ** 0.5 < radius else 0

    def m_init(pos):
        x, y = np.array(pos)[:2] - radius
        return (0, 0.1, 1) if (x ** 2 + y ** 2) ** 0.5 < radius / 2 else (0, 0.1, -1)

    sim.set_Ms(Ms_fun)
    sim.set_m(m_init)
    sim.add(fidimag.micro.UniformExchange(A=A))
    sim.add(fidimag.micro.UniaxialAnisotropy(Ku, axis=(0, 0, 1)))
    sim.add(fidimag.micro.DMI(D=D, dmi_type='interfacial'))
    sim.add(fidimag.micro.Demag())
    sim.driver.energyScale = C.mu_0 * Ms ** 2 * 0.5 * mesh.n * cell ** 3 * 1e-27
    r = sim.driver.minimise(stepControl='BB', stopping_dE=1e-20, mXgradE_tol=tol,
                            max_steps=6000, save_data_steps=NEVER)
    return r, ''


def skyrmion_atomistic(tol=0.1):
    """2D atomistic disk of hubert_minimiser_atomistic.ipynb"""
    J, D, Ku, mus, B = 5.88 * C.meV, 1.56 * C.meV, 0.41 * C.meV, 3 * C.mu_B, 2
    a, az = 0.2715, 0.408
    mesh = fidimag.common.CuboidMesh(nx=100, ny=100, nz=1, dx=a, dy=a, dz=az,
                                     unit_length=1e-9)
    xs = mesh.coordinates[:, 0]
    cx = (xs.max() + xs.min()) * 0.5 + xs.min()
    sim = fidimag.atomistic.Sim(mesh, name='bench_atsk', driver='hubert_minimiser')
    sim.set_mu_s(lambda r: mus if (r[0] - cx) ** 2 + (r[1] - cx) ** 2 < (xs.max() - cx) ** 2 else 0)
    sim.add(fidimag.atomistic.Exchange(J))
    sim.add(fidimag.atomistic.Anisotropy(Ku, axis=(0, 0, 1)))
    sim.add(fidimag.atomistic.DMI(D, dmi_type='interfacial'))
    sim.add(fidimag.atomistic.Zeeman((0, 0, B)))
    sim.set_m(lambda r: (0, 0, -1) if (r[0] - cx) ** 2 + (r[1] - cx) ** 2 < 1 else (0, 0, 1))
    sim.driver.energyScale = J
    r = sim.driver.minimise(stepControl='BB', stopping_dE=1e-10, mXgradE_tol=tol,
                            save_data_steps=NEVER)
    return r, ''


def sp4(cell, max_steps):
    """Standard problem 4 s-state, 500 x 125 x 3 nm, relaxed from (1, 1, 1)"""
    Ms = 8e5
    mesh = fidimag.common.CuboidMesh(nx=int(500 / cell), ny=int(125 / cell), nz=1,
                                     dx=cell, dy=cell, dz=3, unit_length=1e-9)
    sim = fidimag.micro.Sim(mesh, name='bench_sp4', driver='hubert_minimiser')
    sim.set_Ms(Ms)
    sim.add(fidimag.micro.UniformExchange(1.3e-11))
    sim.add(fidimag.micro.Demag())
    sim.set_m((1, 1, 1))
    sim.driver.energyScale = C.mu_0 * Ms ** 2 * 0.5 * 500 * 125 * 3 * 1e-27
    history = record_torque(sim.driver)
    r = sim.driver.minimise(stepControl='BB', stopping_dE=-1.0, mXgradE_tol=0.0,
                            max_steps=max_steps, save_data_steps=NEVER)
    firsts = []
    for tol in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        hit = [step for step, t in history if t < tol]
        firsts.append(str(hit[0]) if hit else '-')
    best = min(t for _, t in history)
    return r, 'first below 1e-1..1e-6: ' + ' '.join(f'{f:>5s}' for f in firsts) + \
        f'   lowest torque {best:.1e} A/m'


if __name__ == '__main__':
    runs = {'1D domain wall': dw_1d,
            'skyrmion + demag': skyrmion_demag,
            'skyrmion + demag 1e-3': lambda: skyrmion_demag(1e-3),
            'atomistic skyrmion': skyrmion_atomistic,
            'atomistic skyrmion 1e-4': lambda: skyrmion_atomistic(1e-4),
            'SP4 2500 cells': lambda: sp4(5, 3000),
            'SP4 10000 cells': lambda: sp4(2.5, 3000)}
    only = sys.argv[1:]
    for name, run in runs.items():
        if only and not any(o in name for o in only):
            continue
        r, extra = run()
        print(f'{name:24s} {r.reason:12s} evals {r.n_evaluations:5d}  '
              f'E {r.total_energy:.12g}  {extra}', flush=True)
