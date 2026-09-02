"""
Benchmark: DemagFMM vs FFT Demag on a sparse nanodisk array.

demag_fmm_vs_fft.py compares the two methods on solid cuboids, which is the
worst possible case for DemagFMM: the FFT method is designed for exactly
that geometry. This script instead targets the case the mu_s == 0 tree
exclusion (see doc/physics_num_methods/demag_fmm.rst, "Where FMM/Barnes-Hut
is worth it") is actually for: a regular array of magnetic nanodisks in an
otherwise empty bounding box, of the kind used for patterned media, MRAM
bit arrays, or spin-torque oscillator arrays.

The array size (number of disks, disk radius, disk thickness) is fixed, so
the number of active (mu_s != 0) atoms stays fixed, while the centre-to-
centre spacing between disks grows. That grows the bounding box, and with
it the FFT method's grid, while DemagFMM's tree - built only from the
active sites since the mu_s == 0 exclusion - should not have to pay for
the growing empty space between disks. Compare with demag_fmm_vs_fft.py,
where growing N grows both methods' costs together.

Usage:
    .venv/bin/python3 benchmarks/demag_fmm_vs_fft_sparse.py
"""
import argparse
import csv
import time

import numpy as np

import fidimag
import fidimag.common.constant as const

A_LATTICE = 0.2715  # nm, bcc Fe-like lattice constant, matches
                     # demag_fmm_vs_fft.py
MU_S = 3 * const.mu_B


def nanodisk_mu_s(coordinates, n_side, period, radius):
    """mu_s array (see Sim.set_mu_s) for an n_side x n_side square array of
    disks of the given radius, spaced `period` apart (both in the mesh's dx
    units), spanning the mesh's full z extent. Built with plain numpy on
    `coordinates` rather than passed as a per-site callback to set_mu_s,
    since that callback is a pure-Python loop over every grid site and is
    far too slow once the bounding box reaches into the millions of sites
    (~1 minute per call at N_grid ~ 2*10^5 in testing)."""
    x, y = coordinates[:, 0], coordinates[:, 1]
    ix = np.clip((x // period).astype(int), 0, n_side - 1)
    iy = np.clip((y // period).astype(int), 0, n_side - 1)
    cx, cy = (ix + 0.5) * period, (iy + 0.5) * period
    inside = (x - cx) ** 2 + (y - cy) ** 2 <= radius ** 2
    return np.where(inside, MU_S, 0.0)


def time_demag(n_side, period, radius, nz, order, theta, ncrit, repeats=3):
    """Build the nanodisk array Sim and time one compute_field() call for
    DemagFMM and for the FFT-based Demag, each averaged over `repeats`
    calls. Returns (n_grid, n_active, t_fmm, t_fft, field_rel_err), where
    field_rel_err is the relative L2 norm between the two methods' field,
    restricted to the active (mu_s != 0) sites -- DemagFMM returns exactly
    zero at mu_s == 0 sites (see the "Where FMM/Barnes-Hut is worth it"
    section of doc/physics_num_methods/demag_fmm.rst), while Demag's FFT
    grid still carries the real (nonzero) stray field reaching those sites
    from the disks, so comparing there would measure that expected,
    deliberate difference rather than DemagFMM's approximation error."""
    side_len = n_side * period
    nx = ny = max(int(round(side_len / A_LATTICE)), 1)

    mesh = fidimag.common.CuboidMesh(
        nx=nx, ny=ny, nz=nz, dx=A_LATTICE, dy=A_LATTICE, dz=A_LATTICE,
        periodicity=(False, False, False), unit_length=1e-9,
    )

    mu_s = nanodisk_mu_s(mesh.coordinates, n_side, period, radius)
    active3 = np.repeat(mu_s != 0, 3)

    sim = fidimag.atomistic.Sim(mesh, name="nanodisk_bench", driver="llg")
    sim.set_mu_s(mu_s)
    n_active = sim.n_nonzero

    demag_fmm = fidimag.atomistic.DemagFMM(order, ncrit, theta)
    demag_fft = fidimag.atomistic.Demag()
    sim.add(demag_fmm)
    sim.add(demag_fft)
    sim.set_m((0, 0, 1))

    t_fmm = 0.0
    for _ in range(repeats):
        start = time.time()
        fmm_field = demag_fmm.compute_field()
        t_fmm += (time.time() - start) / repeats
    fmm_field = fmm_field.copy()

    t_fft = 0.0
    for _ in range(repeats):
        start = time.time()
        fft_field = demag_fft.compute_field()
        t_fft += (time.time() - start) / repeats
    fft_field = fft_field.copy()

    field_rel_err = (np.linalg.norm(fmm_field[active3] - fft_field[active3])
                      / np.linalg.norm(fft_field[active3]))

    return mesh.n, n_active, t_fmm, t_fft, field_rel_err


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-side", type=int, default=4,
                         help="disks per side of the square array")
    parser.add_argument("--radius", type=float, default=2.5,
                         help="disk radius, nm")
    parser.add_argument("--nz", type=int, default=3,
                         help="disk thickness, atomic layers")
    parser.add_argument("--periods", type=float, nargs="+",
                         default=[6, 8, 10, 14, 18, 24, 32, 42, 55, 70, 90],
                         help="centre-to-centre disk spacing, nm")
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--ncrit", type=int, default=128)
    parser.add_argument("--thetas", type=float, nargs="+",
                         default=[0.3, 0.5, 0.7, 0.9])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--out", type=str,
        default="benchmarks/results/demag_fmm_vs_fft_sparse.csv")
    args = parser.parse_args()

    if any(p <= 2 * args.radius for p in args.periods):
        raise ValueError("every period must exceed the disk diameter")

    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["period_nm", "N_grid", "N_active", "sparsity",
                          "order", "theta", "ncrit", "t_fmm_s", "t_fft_s",
                          "field_rel_err"])

        for period in args.periods:
            for theta in args.thetas:
                n_grid, n_active, t_fmm, t_fft, field_rel_err = time_demag(
                    args.n_side, period, args.radius, args.nz,
                    args.order, theta, args.ncrit, args.repeats)
                sparsity = 1 - n_active / n_grid
                writer.writerow([period, n_grid, n_active,
                                  f"{sparsity:.4f}", args.order, theta,
                                  args.ncrit, f"{t_fmm:.6f}", f"{t_fft:.6f}",
                                  f"{field_rel_err:.6e}"])
                f.flush()
                print(f"period={period:6.1f}nm N_grid={n_grid:9d} "
                      f"N_active={n_active:7d} sparsity={sparsity:.3f} "
                      f"theta={theta}: fmm={t_fmm:.4f}s fft={t_fft:.4f}s "
                      f"(fmm/fft={t_fmm / t_fft:.2f}x) "
                      f"field_rel_err={field_rel_err:.2e}")


if __name__ == "__main__":
    main()
