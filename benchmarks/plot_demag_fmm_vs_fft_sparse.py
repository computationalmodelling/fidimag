"""
Plot the output of demag_fmm_vs_fft_sparse.py: field-evaluation time (top)
and accuracy (bottom) vs N_grid (the FFT-based Demag's bounding-box grid
size), one line per theta for DemagFMM plus one line for the FFT-based
Demag, on a fixed nanodisk array with growing spacing between disks.

N_grid, not "sparsity", is the plotted x-axis on purpose: Demag's cost is a
direct function of the bounding-box grid it convolves over, full stop, with
no notion of how much of that grid is empty. Growing the disk spacing here
grows N_grid while N_active stays fixed, which *also* grows the derived
"sparsity" (1 - N_active / N_grid) - but sparsity is a side effect of that
one knob, not an independent cause, so it would be misleading to plot
Demag's time against emptiness itself as though it were the driver.

field_rel_err is the relative L2 field error against Demag, restricted to
the active (mu_s != 0) sites (see demag_fmm_vs_fft_sparse.py's time_demag
docstring for why inactive sites are excluded from the comparison).

Usage:
    .venv/bin/python3 benchmarks/plot_demag_fmm_vs_fft_sparse.py \
        benchmarks/results/demag_fmm_vs_fft_sparse.csv \
        benchmarks/results/demag_fmm_vs_fft_sparse.png
"""
import csv
import sys
from collections import defaultdict

import matplotlib.pyplot as plt


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else \
        "benchmarks/results/demag_fmm_vs_fft_sparse.csv"
    out_path = sys.argv[2] if len(sys.argv) > 2 else \
        "benchmarks/results/demag_fmm_vs_fft_sparse.png"

    rows = list(csv.DictReader(open(csv_path)))

    by_theta = defaultdict(list)
    fft_by_n = {}
    for r in rows:
        n_grid = int(r["N_grid"])
        theta = float(r["theta"])
        by_theta[theta].append((n_grid, float(r["t_fmm_s"]),
                                 float(r["field_rel_err"])))
        fft_by_n[n_grid] = float(r["t_fft_s"])

    fft_points = sorted(fft_by_n.items())

    fig, (ax_t, ax_e) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

    for theta in sorted(by_theta):
        pts = sorted(by_theta[theta])
        n_grid = [p[0] for p in pts]
        ax_t.plot(n_grid, [p[1] for p in pts], marker="o", markersize=3,
                  label=f"DemagFMM (theta={theta})")
        ax_e.plot(n_grid, [p[2] for p in pts], marker="o", markersize=3,
                  label=f"theta={theta}")

    ax_t.plot([p[0] for p in fft_points], [p[1] for p in fft_points],
               marker="o", markersize=3, label="Demag (FFT)", color="black",
               linewidth=2)

    ax_t.set_ylabel("compute_field() time (s)")
    ax_t.set_yscale("log")
    ax_t.set_title("DemagFMM vs FFT Demag, fixed nanodisk array,\n"
                    "growing spacing between disks")
    ax_t.legend()

    ax_e.set_xlabel("N_grid (Demag's bounding-box grid size)")
    ax_e.set_ylabel("field_rel_err vs Demag\n(active sites only)")
    ax_e.set_xscale("log")
    ax_e.set_yscale("log")
    ax_e.legend()

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
