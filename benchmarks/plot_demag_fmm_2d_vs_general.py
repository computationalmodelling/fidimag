"""
Plot demag_fmm_2d_vs_general.py's output: field-evaluation time vs N for
the planar (FMM2D) and general 3D (FMM) operators on the same 2D mesh, one
line per theta, plus the speedup panel.

Usage:
    .venv/bin/python3 benchmarks/plot_demag_fmm_2d_vs_general.py \
        benchmarks/results/demag_fmm_2d_vs_general.csv \
        benchmarks/results/demag_fmm_2d_vs_general.png
"""
import csv
import sys
from collections import defaultdict

import matplotlib.pyplot as plt


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else \
        "benchmarks/results/demag_fmm_2d_vs_general.csv"
    out_path = sys.argv[2] if len(sys.argv) > 2 else \
        "benchmarks/results/demag_fmm_2d_vs_general.png"

    rows = list(csv.DictReader(open(csv_path)))
    by_theta = defaultdict(list)
    for r in rows:
        by_theta[float(r["theta"])].append(
            (int(r["N"]), float(r["t_planar_s"]), float(r["t_general_s"])))

    fig, (ax_t, ax_s) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

    for theta in sorted(by_theta):
        pts = sorted(by_theta[theta])
        n = [p[0] for p in pts]
        ax_t.plot(n, [p[1] for p in pts], marker="o", markersize=3,
                  label=f"planar (theta={theta})")
        ax_s.plot(n, [p[2] / p[1] for p in pts], marker="o", markersize=3,
                  label=f"theta={theta}")

    # One representative general-operator line, since it is close to
    # theta-independent for a given N compared to the planar speedup itself.
    theta0 = sorted(by_theta)[0]
    pts0 = sorted(by_theta[theta0])
    ax_t.plot([p[0] for p in pts0], [p[2] for p in pts0], marker="o",
              markersize=3, color="black", linewidth=2,
              label=f"general (theta={theta0})")

    ax_t.set_ylabel("compute_field() time (s)")
    ax_t.set_yscale("log")
    ax_t.set_title("DemagFMM's planar (2D) vs general (3D) operators,\n"
                    "same 2D mesh")
    ax_t.legend()

    ax_s.set_xlabel("N (spins)")
    ax_s.set_ylabel("speedup (general / planar)")
    ax_s.axhline(1.0, color="grey", linestyle=":", linewidth=1)
    ax_s.legend()

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
