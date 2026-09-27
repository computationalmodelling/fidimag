"""
Benchmark: DemagFMM vs the FFT-based Demag, on cuboid atomistic systems.

This times a single compute_field() call for both interactions, on cubic
LxLxL atomistic cuboid meshes of increasing size. It is a speed comparison
only - see tests/test_demag_fmm.py for correctness checks.

Adapted from exploratory scripts originally written against an older
fmmgen/fidimag integration (test_fmm.py / test_fidimag.sh), rerun here
against the harmonic-compressed FMM operators.

Usage:
    .venv/bin/python3 benchmarks/demag_fmm_vs_fft.py
    .venv/bin/python3 benchmarks/demag_fmm_vs_fft.py --l-min 5 --l-max 60 \
        --thetas 0.3 0.5 0.7 0.9 --order 5 --ncrit 128 \
        --out benchmarks/results/demag_fmm_vs_fft.csv
"""
import argparse
import csv
import time

import numpy as np

import fidimag
import fidimag.common.constant as const

A_LATTICE = 0.2715  # nm, bcc Fe-like lattice constant, matches the source script
J_EXCHANGE = 5.88 * const.meV
MU_S = 3 * const.mu_B


def time_demag(L, order, theta, ncrit, repeats=3):
    """Build an LxLxL cuboid atomistic Sim and time one compute_field() call
    for DemagFMM and for the FFT-based Demag, each averaged over `repeats`
    calls. Returns (n_spins, t_fmm, t_fft, field_rel_err), where
    field_rel_err is the relative L2 norm between the two methods' field --
    a free byproduct of timing both, since both are already computed on the
    same spin configuration."""
    mesh = fidimag.common.CuboidMesh(
        nx=L, ny=L, nz=L, dx=A_LATTICE, dy=A_LATTICE, dz=A_LATTICE,
        periodicity=(False, False, False), unit_length=1e-9,
    )

    sim = fidimag.atomistic.Sim(mesh, name="demag_bench", driver="llg")
    sim.set_mu_s(MU_S)
    sim.add(fidimag.atomistic.Exchange(J_EXCHANGE))
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

    field_rel_err = (np.linalg.norm(fmm_field - fft_field)
                      / np.linalg.norm(fft_field))

    return mesh.n, t_fmm, t_fft, field_rel_err


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l-min", type=int, default=5)
    parser.add_argument("--l-max", type=int, default=60)
    parser.add_argument("--l-step", type=int, default=1)
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--ncrit", type=int, default=128)
    parser.add_argument("--thetas", type=float, nargs="+",
                         default=[0.3, 0.5, 0.7, 0.9])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--out", type=str,
                         default="benchmarks/results/demag_fmm_vs_fft.csv")
    args = parser.parse_args()

    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["L", "N", "order", "theta", "ncrit",
                          "t_fmm_s", "t_fft_s", "field_rel_err"])

        for L in range(args.l_min, args.l_max + 1, args.l_step):
            for theta in args.thetas:
                n, t_fmm, t_fft, field_rel_err = time_demag(
                    L, args.order, theta, args.ncrit, args.repeats)
                writer.writerow([L, n, args.order, theta, args.ncrit,
                                  f"{t_fmm:.6f}", f"{t_fft:.6f}",
                                  f"{field_rel_err:.6e}"])
                f.flush()
                print(f"L={L:3d} N={n:7d} theta={theta}: "
                      f"fmm={t_fmm:.4f}s fft={t_fft:.4f}s "
                      f"(fmm/fft={t_fmm / t_fft:.1f}x) "
                      f"field_rel_err={field_rel_err:.2e}")


if __name__ == "__main__":
    main()
