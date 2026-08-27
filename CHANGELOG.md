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
