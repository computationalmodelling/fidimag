Version 4.0
-----------

### Energy minimisation

* **Hubert minimiser**: new `stepControl='BB'` argument to `minimise`, which
  uses Barzilai-Borwein step lengths with a non-monotone line search instead of
  the creep algorithm. It needs no `eta_scale`, which had to be matched by hand
  to the units of the effective field, and reaches the same minimum in fewer
  evaluations (139 against 649 on a 1D domain wall, 950 against more than 6000
  on a skyrmion with demagnetising field). The default remains
  `stepControl='hubert'`
* **Steepest descent**: fixed the Barzilai-Borwein step size never reaching
  Python from `sd_compute_step`, which took it by value and dropped it, so the
  class had been running with a constant step and was sensitive to
  `initial_t_step`
* **Steepest descent**: a non-positive Barzilai-Borwein quotient no longer
  produces a step up the energy, which could leave the iteration stalled in a
  configuration that was not a minimum
* **Steepest descent**: new optional `energy_guard` argument, which accepts a
  step only if it passes a non-monotone sufficient-decrease test. It is what
  makes an aggressive `tmax` usable
* New `tests/test_steepest_descent.py`, and a new documentation section on
  energy minimisation

### Demagnetising field

* **Faster FFT demag**: the `O(N)` glue loops around the FFTs are now
  OpenMP-parallelised, and the spectral tensor multiply is bounded to the
  `lenz * leny * (lenx/2 + 1)` points an r2c transform actually produces
  instead of iterating over the whole padded grid, which was about twice the
  necessary work. The padding now goes to the next even 7-smooth size above
  `2n` rather than exactly `2n`, so awkward mesh sizes avoid FFTW's slow
  large-prime path. Measured per demag field call, 2D and 4 threads: 2.1x at
  `n = 32`, 2.5x at `n = 64`, 2.2x at `n = 128`, and 4.7x at `n = 127`, which
  was large-prime bound. Field values are bit-for-bit unchanged
* **Lower memory**: the twelve spectral arrays were allocated with
  `total_length` entries, roughly twice what an r2c transform can fill. Sizing
  them to the spectrum saves about 24% per `Demag` object, 44.1 -> 33.5 MB for
  a 200x200x1 mesh and 135.6 -> 101.2 MB for 64x64x16. The allocation sizes are
  also computed in `size_t` now: the product crossed `INT_MAX` and was
  truncated for a padded grid of about 1.3e8 points, which a 256^3 mesh reaches
  exactly
* **Fixed the 2d_pbc demag**, which never reproduced an infinite film: for a
  uniformly magnetised one it gave `Nzz = 0.329` and `Nxy = -0.365` where a
  film must give `Nzz = 1` and `Nxy = 0`. The kernel was laid out with its zero
  offset at index `n - 1` while the convolution expects it at index 0, and the
  signs of the odd components were wrong. An out-of-bounds read in the same
  tensor fill is also fixed
* **FMM/BH demag restored**: the port to the uv build had dropped `DemagFMM`
  and its C++ extension. It is back and wired into CMake, `compute_energy` no
  longer returns 0, and `tests/test_demag_fmm.py` checks both the `fmm` and
  `bh` methods against brute-force `DemagFull`

### VTK output

* The pyvtk dependency is gone, replaced by a pure-Python writer for the modern
  XML formats: `.vti` (ImageData) for `CuboidMesh` and `.vtp` (PolyData) for
  `HexagonalMesh`, with cell data inline-encoded as base64 float32. The
  `SaveVTK` wrapper class was folded into direct `VTK` usage
* `VTK.save_as()` writes a single file at an arbitrary path, and `save_vtk()`
  takes an optional filename
* The `Scalars` and `Vectors` attributes of `CellData` name one array each, and
  the writer was joining every name with spaces, so a file with more than one
  scalar named no valid array and readers fell back on a default. The first of
  each is named now; the others are still written and selectable
* Every file records what wrote it in `FieldData`, naming the Fidimag version
  and describing the mesh, e.g. `fidimag 4.0 | CuboidMesh nx=4 ny=3 nz=2 n=24
  dx=1 dy=2 dz=3 | unit_length=1e-09`. This is what makes a `.vti` still
  readable as evidence once it is one of several hundred in a `vtks/` directory

### Python 2 clean-up

* Removed `from __future__` imports, `u''` string prefixes, `# -*- coding:
  utf-8 -*-` headers, `class X(object)` bases and `super(X, self)` calls
  throughout
* Converted the remaining Python 2 `print` statements in `examples/` and
  `tests/jacobian_computation.py`, which could not be run at all before
* Updated the notebooks under `doc/user_guide/ipynb/` that still declared a
  Python 2 kernel

### Documentation

* Read the Docs now builds on `ubuntu-26.04` with `miniforge3`: the configured
  `ubuntu-20.04` is no longer offered and `mambaforge` is deprecated
* The documentation dependencies are pinned in `doc/environment.yml`, and the
  `docs` extra installs what the docs actually import
* The installation instructions describe the `uv sync` build instead of
  `pip install -e .` and `make`
* The NEBM section describes the current code: the `ChainMethodBase` class that
  the methods derive from, `StringMethod`, the removal of `NEBM_Cartesian`, the
  single `nebm_clib` Cython module that replaced the one-per-coordinate-system
  ones, the choice of integrator, and the renaming of `nebm_tools` to
  `chain_method_tools`. The class diagram is regenerated, and references to
  Abert, and to Fabian and Shcherbakov, are added
* The `steepest_descent_atomistic` and `hubert_minimiser_atomistic` notebooks
  are re-executed, since their stored outputs predate the fixes to both
  minimisers. The steepest descent ones show the difference plainly: the one
  dimensional relaxation goes from 711 steps to 176, and the two dimensional
  one from not converging in 6000 steps to converging in 192

Version 3.0
-----------

### Build System Modernization

**Major change**: The build system has been completely modernized from `setup.py` to a modern PEP 517/518 compliant system using `scikit-build-core`, `CMake`, and `uv`.

#### Key Changes

* **Build Backend**: Migrated from setuptools to scikit-build-core + CMake
* **Package Manager**: Now using `uv` with lock files (`uv.lock`) for reproducible builds
* **Cython Code Generation**: Now generates files out-of-source in `build/` directory instead of in-place
* **Configuration**: Moved from `setup.py` to declarative `pyproject.toml` + `CMakeLists.txt`
* **CI/CD**: Updated GitHub Actions workflows to use `uv` and the new build system
* **Dependencies**: All C/C++ libraries now properly managed by CMake (SUNDIALS, FFTW, BLAS, LAPACK, OpenMP)

#### Installation Changes

**Old way**:
```bash
pip install -e .
```

**New way**:
```bash
uv sync
```

#### For Developers

* New `CMakeLists.txt` for build configuration
* New `cmake/FindSUNDIALS.cmake` for dependency management
* Added `clean_cython_files.sh` for cleaning old build artifacts
* Added `BUILD.md` for detailed build documentation
* Added `MIGRATION.md` for migration guide

#### Breaking Changes

* `setup.py` deprecated (retained for backward compatibility, will be removed in v3.1)
* Cython-generated `.c` files no longer in source tree (now in `build/` only)
* Requires CMake 3.18+ for building from source
* Python 3.9+ required (Python 3.14 recommended)
* **Python 2 no longer supported** - removed `docker/minimal-py2/`

#### Migration

See `MIGRATION.md` for detailed migration instructions. Quick steps:

```bash
# Clean old artifacts
bash clean_cython_files.sh

# Install with new build system
uv sync
```

#### Benefits

* 10-100x faster package installation with `uv`
* Reproducible builds via lock files
* Cleaner source tree (no generated files)
* Better cross-platform support
* Improved dependency management

---

Version 3.0 Alpha
------------------
* Changes to the helper functions init_scalar and init_vector (fidimag/common/helper.py)
  which allow users to pass additional parameters. These are then used within the sim
  classes to allow 
  
  For example:

  ```
  
  mesh = CuboidMesh(nx=10, ny=1, nz=1, unit_length=1e-9)
  sim = Sim(mesh, Ms)

  def init_domain_wall(pos, domain_centre)
      x, y, z = pos

      if x < domain_centre:
          return (0, 0, 1)

      else:
          return (0, 0, -1)
   
  # Place a domain wall at 5nm
  sim.set_m(init_domain_wall, 5)
  # Place a domain wall at 3nm
  sim.set_m(init_domain_wall, 3)

  ```

* Setting currents is now more general, and is standardised across the simulation types.
  This allows us to use more general functions for setting the current.
  Previously, the current function was set as:
  ```
  sim(mesh, driver='llg_stt')
  sim.driver.jx = 1e14 # A / m^2
  sim.driver.update_j_fun = lambda t: np.sin(t)
  ```
  with the actual current used being multiplicative:

  $ jx * sin(t) $

  For the current-perpendicular to the plane STT ('llg_stt_cpp') driver
   we would now change this to 

  ```
  sim.driver(mesh, driver='llg_stt_cpp')
  sim.driver.j_function = lambda t: 1e14 * np.sin(t)
  ```
  and for the standard STT driver:

  ```
  sim.drive(mesh, driver='llg_stt')
  sim.driver.jz_function = lambda t: 1e14 * np.sin(t)
  # Can also set:
  # sim.driver.jx_function = ...
  # sim.driver.jy_function = ...

* Similarly, the TimeZeeman interaction is also no longer multiplicative;
  you can have an arbitrary function of the form:
 
  def time_function(pos, t):
      x, y, z = pos
      # some calculation of Bx(pos, t), By(pos, t), Bz(pos, t)
      return (Bx, By, Bz)
  zee = TimeZeeman(np.array([0, 0, 0]), time_function)
  sim.add(zee)

* You can now remove energy classes from the simulation.

  This can be useful in cases where you have an interaction
  but no longer need to calculate it because the simulation has 
  reached a certain point, e.g. an applied field has been turned off.
 
  In the data table, the entries corresponding to a given interaction 
  will be zero everywhere once the interaction is removed.
 
  
  For example:

  ```
  sim.add(Zeeman((0, 0, B), name='Zeeman'))
  
  sim.run_until(1e-9)
  sim.remove('Zeeman')
  ```
