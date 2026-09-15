"""
Plot the output of demag_fmm_vs_fft.py: field-evaluation time vs system size,
one line per theta for DemagFMM plus one for the FFT-based Demag.

Usage:
    .venv/bin/python3 benchmarks/plot_demag_fmm_vs_fft.py \
        benchmarks/results/demag_fmm_vs_fft.csv \
        benchmarks/results/demag_fmm_vs_fft.png
"""
import csv
import sys
from collections import defaultdict

import matplotlib.pyplot as plt


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else \
        "benchmarks/results/demag_fmm_vs_fft.csv"
    out_path = sys.argv[2] if len(sys.argv) > 2 else \
        "benchmarks/results/demag_fmm_vs_fft.png"

    fmm_by_theta = defaultdict(list)
    fft = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            n = int(row["N"])
            theta = float(row["theta"])
            fmm_by_theta[theta].append((n, float(row["t_fmm_s"])))
            fft.append((n, float(row["t_fft_s"])))

    fft = sorted(set(fft))

    fig, ax = plt.subplots(figsize=(7, 5))
    for theta in sorted(fmm_by_theta):
        pts = sorted(fmm_by_theta[theta])
        ax.plot([p[0] for p in pts], [p[1] for p in pts],
                marker="o", markersize=3, label=f"DemagFMM (theta={theta})")
    ax.plot([p[0] for p in fft], [p[1] for p in fft],
            marker="o", markersize=3, label="Demag (FFT)", color="black",
            linewidth=2)

    ax.set_xlabel("N (spins)")
    ax.set_ylabel("compute_field() time (s)")
    ax.set_yscale("log")
    ax.set_title("DemagFMM vs FFT-based Demag, cubic atomistic meshes")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
