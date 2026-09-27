# Benchmarks

Speed comparisons for fidimag's interactions, as opposed to the correctness
tests under `tests/`.

## `demag_fmm_vs_fft.py`

Times a single `compute_field()` call for `DemagFMM` against the FFT-based
`Demag`, on cubic `LxLxL` atomistic cuboid meshes, sweeping `L` and `theta`.
Adapted from exploratory scripts originally written against an older
fmmgen/fidimag integration, rerun here against the harmonic-compressed FMM
operators (see `doc/physics_num_methods/demag_fmm.rst`).

```
.venv/bin/python3 benchmarks/demag_fmm_vs_fft.py \
    --l-min 5 --l-max 60 --thetas 0.3 0.5 0.7 0.9 \
    --order 8 --ncrit 128 \
    --out benchmarks/results/demag_fmm_vs_fft.csv
```

Results are written incrementally as CSV (`L, N, order, theta, ncrit,
t_fmm_s, t_fft_s, field_rel_err`), the last a free byproduct of timing both
methods on the same spin configuration. Plot with:

```
.venv/bin/python3 benchmarks/plot_demag_fmm_vs_fft.py \
    benchmarks/results/demag_fmm_vs_fft.csv \
    benchmarks/results/demag_fmm_vs_fft.png
```

## `demag_fmm_vs_fft_sparse.py`

Same comparison, but on a fixed array of nanodisks with growing spacing
between them, rather than a growing solid cuboid. The active (`mu_s != 0`)
site count stays fixed while the FFT method's bounding-box grid grows, which
is the geometry the `mu_s == 0` tree exclusion in `DemagFMM` is for. See
`doc/physics_num_methods/demag_fmm.rst`, "Where FMM/Barnes-Hut is worth it".

```
.venv/bin/python3 benchmarks/demag_fmm_vs_fft_sparse.py \
    --n-side 4 --radius 2.5 --nz 3 \
    --periods 6 8 10 14 18 24 32 42 55 70 90 \
    --order 8 --thetas 0.3 0.5 0.7 0.9 --ncrit 128 \
    --out benchmarks/results/demag_fmm_vs_fft_sparse.csv
```

Results are written incrementally as CSV (`period_nm, N_grid, N_active,
sparsity, order, theta, ncrit, t_fmm_s, t_fft_s, field_rel_err`), the last
restricted to the active (`mu_s != 0`) sites (see `time_demag`'s docstring
in the script). Plot with:

```
.venv/bin/python3 benchmarks/plot_demag_fmm_vs_fft_sparse.py \
    benchmarks/results/demag_fmm_vs_fft_sparse.csv \
    benchmarks/results/demag_fmm_vs_fft_sparse.png
```

`plot_nanodisk_array.py` renders the `mu_s` pattern itself at a single
period, so the geometry behind the numbers is visible:

```
.venv/bin/python3 benchmarks/plot_nanodisk_array.py \
    --n-side 4 --radius 2.5 --nz 3 --period 18 \
    --out benchmarks/results/nanodisk_array_geometry.png
```

## `demag_fmm_2d_vs_general.py`

Times `fidimag.extensions.fmm.FMM2D` (fmmgen's planar operator variant,
which `DemagFMM` selects automatically for a 2D, `nz=1` mesh) against the
general `FMM` called directly on the same 2D mesh's coordinates, forcing
the path a 3D mesh would otherwise take. See
`doc/physics_num_methods/demag_fmm.rst`, "2D systems".

```
.venv/bin/python3 benchmarks/demag_fmm_2d_vs_general.py \
    --l-min 10 --l-max 150 --l-step 10 \
    --order 8 --thetas 0.3 0.5 0.7 0.9 --ncrit 128 \
    --out benchmarks/results/demag_fmm_2d_vs_general.csv
```

Results are written incrementally as CSV (`L, N, order, theta, ncrit,
t_planar_s, t_general_s, field_rel_err`), the last the relative L2
difference between the two operator sets on identical input (both correct,
so this should stay small rather than trend to zero). Plot with:

```
.venv/bin/python3 benchmarks/plot_demag_fmm_2d_vs_general.py \
    benchmarks/results/demag_fmm_2d_vs_general.csv \
    benchmarks/results/demag_fmm_2d_vs_general.png
```
