"""
Render the mu_s == 0 / mu_s != 0 pattern of the nanodisk array used by
demag_fmm_vs_fft_sparse.py, at a single representative period, so the
geometry behind that benchmark's numbers is visible rather than just its
timings.

Usage:
    .venv/bin/python3 benchmarks/plot_nanodisk_array.py \
        --n-side 4 --radius 2.5 --nz 3 --period 18 \
        --out benchmarks/results/nanodisk_array_geometry.png
"""
import argparse

import matplotlib.pyplot as plt

import fidimag
from demag_fmm_vs_fft_sparse import A_LATTICE, nanodisk_mu_s


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-side", type=int, default=4)
    parser.add_argument("--radius", type=float, default=2.5)
    parser.add_argument("--nz", type=int, default=3)
    parser.add_argument("--period", type=float, default=18)
    parser.add_argument(
        "--out", type=str,
        default="benchmarks/results/nanodisk_array_geometry.png")
    args = parser.parse_args()

    side_len = args.n_side * args.period
    nx = ny = max(int(round(side_len / A_LATTICE)), 1)

    mesh = fidimag.common.CuboidMesh(
        nx=nx, ny=ny, nz=args.nz, dx=A_LATTICE, dy=A_LATTICE, dz=A_LATTICE,
        periodicity=(False, False, False), unit_length=1e-9,
    )
    mu_s = nanodisk_mu_s(mesh.coordinates, args.n_side, args.period,
                          args.radius)

    # A single z-layer is enough: the pattern is the same in every layer.
    layer = mu_s.reshape(args.nz, ny, nx)[0]

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.imshow(layer > 0, origin="lower",
              extent=[0, side_len, 0, side_len], cmap="Greys")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title(
        f"Nanodisk array: {args.n_side}x{args.n_side} disks, "
        f"radius={args.radius}nm,\n"
        f"period={args.period}nm, nz={args.nz} atomic layers")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")
    print(f"n_active={int((mu_s > 0).sum())} n_grid={mesh.n}")


if __name__ == "__main__":
    main()
