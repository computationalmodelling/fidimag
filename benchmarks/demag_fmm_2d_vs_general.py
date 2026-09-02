"""
Benchmark: DemagFMM's automatic 2D (planar) routing vs the general 3D
operators, on the same 2D atomistic mesh.

A 2D mesh (nz=1) has every site at the same z, so fmmgen's planar operator
variant (generate_code(..., planar=True)) can be used exactly, not as an
approximation - it simply omits the always-zero terms the general 3D
operators would still carry (see doc/physics_num_methods/demag_fmm.rst,
"2D systems"). This times DemagFMM (which auto-selects the planar variant
for a 2D mesh) against fidimag.extensions.fmm.FMM called directly on the
same 2D mesh's coordinates, forcing the general 3D path fmmgen.extensions.
fmm.FMM2D would otherwise replace.

Usage:
    .venv/bin/python3 benchmarks/demag_fmm_2d_vs_general.py
"""
import argparse
import csv

import numpy as np

import fidimag
import fidimag.extensions.fmm as fmm_ext
import fidimag.common.constant as const

A_LATTICE = 0.2715  # nm, matches demag_fmm_vs_fft.py
MU_S = 3 * const.mu_B


def time_one(fn, repeats):
    import time
    t = 0.0
    for _ in range(repeats):
        start = time.time()
        fn()
        t += (time.time() - start) / repeats
    return t


def time_demag(L, order, theta, ncrit, repeats=3):
    """Build an LxL 2D atomistic mesh and time one compute_field() call for
    the planar (FMM2D) and general 3D (FMM) operators on the same
    coordinates. Returns (n_spins, t_planar, t_general, field_rel_err)."""
    mesh = fidimag.common.CuboidMesh(
        nx=L, ny=L, nz=1, dx=A_LATTICE, dy=A_LATTICE, dz=A_LATTICE,
        periodicity=(False, False, False), unit_length=1e-9,
    )
    n = mesh.n
    coords3 = (mesh.coordinates * mesh.unit_length).astype(np.float64)
    coords2 = np.ascontiguousarray(coords3[:, :2])
    m = np.zeros(3 * n)
    m[2::3] = 1.0
    mu_s = np.full(n, MU_S)

    planar = fmm_ext.FMM2D(n, ncrit, theta, order, coords2, m, mu_s, 0)
    field_planar = np.zeros(3 * n)
    t_planar = time_one(lambda: planar.compute_field(field_planar), repeats)
    field_planar = field_planar.copy() * -1e-7

    general = fmm_ext.FMM(n, ncrit, theta, order, coords3, m, mu_s, 0, True)
    field_general = np.zeros(3 * n)
    t_general = time_one(lambda: general.compute_field(field_general), repeats)
    field_general = field_general.copy() * -1e-7

    field_rel_err = (np.linalg.norm(field_planar - field_general)
                      / np.linalg.norm(field_general))

    return n, t_planar, t_general, field_rel_err


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l-min", type=int, default=10)
    parser.add_argument("--l-max", type=int, default=150)
    parser.add_argument("--l-step", type=int, default=10)
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--ncrit", type=int, default=128)
    parser.add_argument("--thetas", type=float, nargs="+",
                         default=[0.3, 0.5, 0.7, 0.9])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--out", type=str,
        default="benchmarks/results/demag_fmm_2d_vs_general.csv")
    args = parser.parse_args()

    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["L", "N", "order", "theta", "ncrit",
                          "t_planar_s", "t_general_s", "field_rel_err"])

        for L in range(args.l_min, args.l_max + 1, args.l_step):
            for theta in args.thetas:
                n, t_planar, t_general, err = time_demag(
                    L, args.order, theta, args.ncrit, args.repeats)
                writer.writerow([L, n, args.order, theta, args.ncrit,
                                  f"{t_planar:.6f}", f"{t_general:.6f}",
                                  f"{err:.6e}"])
                f.flush()
                print(f"L={L:3d} N={n:6d} theta={theta}: "
                      f"planar={t_planar:.4f}s general={t_general:.4f}s "
                      f"(speedup={t_general / t_planar:.2f}x) "
                      f"field_rel_err={err:.2e}")


if __name__ == "__main__":
    main()
