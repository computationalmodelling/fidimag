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
* **Hubert minimiser**: the Barzilai-Borwein path no longer stops short at a
  torque of around 1e-6 A/m, with `Could not decrease the energy along the
  gradient`. That was not the step rule but the precision of the acceptance
  test: near a minimum the decrease to be detected falls below the spacing of
  the doubles around the total energy, 1.9e-35 J against a spacing of 1.4e-34
  J on the standard problem 4 film, so every trial step looked like a failure.
  The change of energy is now summed over the sites,
  `dE = sum_i (E_i_new - E_i_ref)`, rather than taken as the difference of two
  totals. Each term is a difference of two cell energies rather than of two
  whole-sample energies, so nothing large is subtracted. This is what OOMMF's
  `Oxs_CGEvolve` has done since 2002, and what the per cell energy clean-up
  below makes possible here. The same problem now reaches 2.4e-10 A/m instead
  of 1.5e-6, in fewer field evaluations than OOMMF needs to reach 1e-6, and
  gives the same answer to the last bit for any `energyScale` from 1 to 1e-31.
  Scaling the energy was never a cure: five values of `energyScale` that all
  put the total within an order of magnitude of one stopped between 4.5e-8 and
  9.5e-6 A/m, with no trend. A second route, a trapezoid in the effective
  field that needs no per cell energy at all, reaches the same convergence and
  is recorded in the code and the documentation rather than used
* **Hubert minimiser**: measured against OOMMF's `Oxs_CGEvolve` on the standard
  problem 4 relaxation, counted in effective field evaluations, the
  Barzilai-Borwein path now needs fewer than the conjugate gradient at every
  tolerance: 413 against 614 to reach 1e-6 A/m on 2500 cells, 1038 against
  1278 on 10000. See the energy minimisation page
* **Hubert minimiser**: `minimise` returns a `MinimiserResult`, saying whether
  a convergence criterion was met, which one (`mXgradE_tol`, `stopping_dE`,
  `zero_gradient`, `max_steps` or `resets`), how many effective field
  evaluations it spent, and the energy and torque it finished at. There was
  previously no way to tell a converged run from one that ran out of steps
  except by reading the printed output
* **The chain method classes no longer configure the root logger.**
  `chain_method_base.py`, `nebm_spherical.py`, `nebm_geodesic.py` and
  `string_method.py` each called `logging.basicConfig(level=logging.DEBUG)` at
  import, so importing any of them turned DEBUG on globally, for every library
  in the process, not only for Fidimag. Configuring logging belongs to the
  application, not to a library. The NEBM relaxation output is at DEBUG and
  its closing summary at INFO, so both are now silent by default; a script
  that wants them asks for them::

      logging.basicConfig(level=logging.DEBUG)

* **Steepest descent and `SimpleMinimiser`**: progress goes to the `fidimag`
  logger too, at the same levels as the Hubert minimiser, so no minimiser
  writes to stdout or stderr any more. The steepest descent had also been
  reporting its convergence unconditionally, ignoring `printing`, and warning
  about a failure to converge through `sys.stderr`. `printing` and `log_every`
  are kept, and now throttle the DEBUG lines
* **`MinimiserBase.relax`** raises `NotImplementedError` instead of printing
  that it is not implemented and returning `None`

* **Hubert minimiser**: progress goes to the `fidimag` logger rather than to
  `print`, as the NEBM classes already did. The per step lines are at DEBUG
  and the reason it stopped at INFO, so a converged run is now silent; raise
  the level to see them again::

      logging.getLogger('fidimag').setLevel(logging.INFO)

  Failures to converge are logged at WARNING and are still shown. `log_steps`
  is kept, and now throttles the DEBUG lines
* **Steepest descent**: fixed the Barzilai-Borwein step size never reaching
  Python from `sd_compute_step`, which took it by value and dropped it, so the
  class had been running with a constant step and was sensitive to
  `initial_t_step`
* **Steepest descent**: a non-positive Barzilai-Borwein quotient no longer
  produces a step up the energy, which could leave the iteration stalled in a
  configuration that was not a minimum
* **Steepest descent**: new `stopping_torque` argument, which stops once no
  site has a torque `||m x (m x H)||` larger than the value given, in the
  units of the effective field. Prefer it to `stopping_dm`, which stops on how
  far a spin moved in a step: that distance is the product of the step length
  and the torque, and Barzilai-Borwein step lengths swing over orders of
  magnitude, so the test can pass on a short step while the residual is
  unchanged. On the standard problem 4 s-state `stopping_dm = 1e-9` ended the
  iteration at a torque of 7e-3 A/m. The new criterion is the residual the
  Hubert minimiser already stops on and the one OOMMF reports as
  `Max mxHxm`, so the three are directly comparable, and it costs nothing:
  `max_torque()` reads a quantity that is computed every step anyway.
  `stopping_dm` still works and is unchanged when `stopping_torque` is not
  given
* **Steepest descent**: the default step ceiling `tmax` is 1 rather than 0.1.
  Relaxing the standard problem 4 film to its s-state takes 2498 evaluations
  of the effective field at 1 against 12318 at 0.1, reaching the same state.
  Higher is faster still, 1478 at 3 and 886 at 10, but the ceiling cannot be
  raised far without `energy_guard`: at 10 the one dimensional domain wall of
  `tests/test_steepest_descent.py` collapses, to a profile error of 0.90
  rather than 0.0018. The default of 1 leaves an order of magnitude of margin
  below that. Runs that relied on the old ceiling will reach the same state in
  about a fifth of the steps
* **Steepest descent**: new optional `energy_guard` argument, which accepts a
  step only if it passes a non-monotone sufficient-decrease test. It is what
  makes an aggressive `tmax` usable
* New `tests/test_steepest_descent.py`, and a new documentation section on
  energy minimisation

### Time integration

* **The integrator names now say what they run.** `sundials` named the suite
  rather than the solver, which distinguished nothing once ARKODE was used as
  well, and did not say which method was being run. The names are
  `cvode_bdf` (the default), `cvode_bdf_diag`, `cvode_adams`,
  `cvode_bdf_openmp`, `cvode_bdf_diag_openmp`, `arkode_dopri5`,
  `arkode_rkf45`, and their `_normalised` and `_openmp` variants. `euler` and `rk4` keep
  their bare names, since they are implemented here rather than taken from a
  library. The pre-4.0 names `sundials`, `sundials_diag`, `sundials_openmp`
  and `sundials_diag_openmp` still work and raise a `DeprecationWarning`
  naming their replacement, so scripts written against earlier versions keep
  running

* **`cvode_adams`**: CVODE with Adams-Moulton and a fixed point iteration,
  the non-stiff arm of the solver Fidimag already wraps, with no linear solve
  at all. `CV_ADAMS` had been declared in the Cython wrapper since it was
  written but never used, `CVodeCreate` being called with `CV_BDF`
  unconditionally
* **`arkode_dopri5` and `arkode_rkf45`**: explicit Runge-Kutta 5(4) pairs, through ARKODE's
  `ERKStep`. `arkode_dopri5` uses the Dormand and Prince tableau, which is the method
  OOMMF calls `rkf54m`, so the two codes can now be compared directly;
  `arkode_rkf45` is the genuine Fehlberg one. (OOMMF's own default, `rkf54`, is RK5(4)7FC,
  another member of the Dormand and Prince family rather than the Fehlberg
  tableau the name suggests.) ARKODE is now linked, and the Butcher table is
  selected by name, so any of the tables it ships can be reached
* On muMAG standard problem 4 at rtol = atol = 1e-10, all four agree to within
  2e-8 and put the crossing of `<m_x> = 0` in the same place. Relative to the
  BDF default, Adams took 0.66 of the wall time and the explicit methods 0.47.
  That ordering is particular to a problem of this size: the stiffness of the
  LLG equation comes from exchange and grows as `1 / dx ** 2`, so the implicit
  methods recover their advantage on a finer mesh
* **`arkode_dopri5_normalised` and `arkode_rkf45_normalised`**: the same explicit methods,
  rescaling the spins to unit length after every accepted step through
  ARKODE's post-step hook, rather than relying only on the `c * (1 - m^2) * m`
  correction term in the LLG right hand side, which makes `|m| = 1` an
  attracting solution instead of imposing it. The two are independent, and
  `default_c` still selects the correction term's behaviour (0 turns it off).
  On standard problem 4 the projection left the step count unchanged and cost
  about 2% in wall time, improving the error in `|m|` by two orders of
  magnitude. CVODE offers no equivalent hook, so this is explicit-only
* **`_openmp` variants of the arkode integrators**, using the OpenMP
  N_Vector. This threads only the integrator's own vector arithmetic, the
  stage combinations and the error norms; the right hand side is already
  parallel inside Fidimag's C code whichever vector is used. It is therefore
  worth less to an explicit method than to CVODE, whose Krylov solver does
  much more vector work: on standard problem 4 with four threads,
  `arkode_dopri5` went from 4.05 s to 3.74 s (1.08x) while `cvode_bdf` went
  from 8.63 s to 6.67 s (1.29x). At one thread they are indistinguishable
  from the serial classes, so the vector carries no overhead of its own. The
  fastest combination measured was `arkode_dopri5_openmp`, 3.74 s against
  8.63 s for the old default, a factor of 2.3
* New `tests/test_std_problem_4.py`, marked slow, which runs muMAG standard
  problem 4 against the OOMMF reference in `examples/micromagnetic/std4`. That
  comparison existed only inside the example, so nothing ran it. The deviation
  of the mean magnetisation over the full nanosecond is max 3.5e-05 and rms
  1.3e-05; a second test checks that the explicit method reproduces the same
  reference
* New `tests/test_integrators.py`, checking that the methods agree on
  a precessing macrospin and that the projection conserves the spin length
  better at the same cost

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

### Fixes

* **`fidimag.common.plot` no longer uses `mpl_toolkits`.** The colorbars were
  built with `ImageGrid`, whose only purpose here was to make them scale with
  the subplots, which `Figure.colorbar` does itself when told which axes to
  take the room from. That dependency is what broke the module against modern
  matplotlib. Also fixed while there: `plot_atom_hex(scale_by_mag=True)`
  raised `AttributeError` for every component, on a duplicated fragment that
  called `cbar.ax.cax.colorbar`; `savefig` passed `dpi=1`, which renders a
  handful of pixels; the three plotting functions returned `None` and now
  return the figure; and `plot` dispatched on `type(sim) ==` rather than
  `isinstance`. The data drawn is unchanged, checked array by array over
  every component and mesh type

* **Minimising an atomistic system with a demagnetising field was wrong.**
  The minimisers build the effective field by asking each interaction for its
  energy and then reading its `field`, and neither atomistic demagnetising
  class recomputed the field when asked for the energy, so the field belonged
  to whatever configuration came before. `Demag` is an `Energy` subclass, so
  it also left `total_energy` at zero and contributed nothing to the energy;
  the failure was silent, and on a 24x24 lattice with exchange and a field it
  relaxed to a mean m_z of 0.88 where the answer is 0.00. `DemagHexagonal` is
  not, so it raised `AttributeError` instead and could not be used with a
  minimiser at all. Results of atomistic minimisation with either class will
  change. Time integration was never affected, since it builds the field from
  `compute_field` and never asks for an energy
* **Per cell energies now mean the same thing in every interaction.**
  `energy[i]` is the energy of cell or site `i` in joules and `total_energy`
  is their sum, in both the micromagnetic and the atomistic classes. It used
  to be an energy density in most of the micromagnetic ones, an already
  weighted energy in `micro/demag.py`, and nothing at all in
  `micro/zeeman.py`, which never wrote to the array; `SimpleDemag` returned
  the array in place of a total and scaled it by a `mesh.cellsize` that does
  not exist, so any call raised. Both `Zeeman` classes and `DemagHexagonal`
  reimplemented the base class setup by hand, which is how they drifted, and
  now inherit it. `Relaxation` reports zero rather than nothing. No energy
  that was previously correct has changed value
* **The OOMMF comparison tests run again.** The file was named
  `tes_oommf.py`, so pytest never collected it and it had drifted: `DMI` had
  been renamed to take `dmi_type`, the harness wrote its MIF scripts and OOMMF
  output *inside the installed package*, and it looked for the saved field
  under a fixed stage index, which OOMMF 2.1 numbers from zero rather than
  one. It is now `test_oommf.py`, writes to a temporary directory
  (`OOMMF_WORK_DIR` to override), and takes whatever stage file OOMMF wrote.
  The two tests that call OOMMF but were not marked now carry the
  `run_oommf` mark, so a run without OOMMF still skips them
* **Fidimag and OOMMF cannot agree beyond 1.3e-10 on any field carrying a
  factor of 1/mu_0.** Fidimag uses the pre-2019 exact `mu_0 = 4 pi x 1e-7` and
  OOMMF 2.x uses the CODATA value `12.5663706127e-7`, which differ by
  1.32033e-10 relatively. The ratio of the two exchange fields is measured at
  1 + 1.3203e-10, constant to 1.7e-13 across the mesh, so this accounts for
  the whole disagreement. The exchange and DMI tolerances are set from it; the
  demagnetising field has no such factor and still matches to 1e-11
* **`pytest.ini` is gone.** It shadowed the configuration in `pyproject.toml`,
  so the `run_oommf` marker was reported as unknown. Removing it exposed that
  `pyproject.toml` collected only `test_*.py`, dropping the five `*_test.py`
  files; both patterns are now listed
* **The analytical Jacobian works** (issue #21, open since 2016). Passing
  `use_jac=True` appeared to hang. It did not: the Jacobian-times-vector
  product was wrong, so GMRES never converged and spent its whole iteration
  budget on every solve, half a million evaluations for one picosecond of
  simulated time. Three faults: the routine was handed `dm/dt` where it wanted
  the effective field at m, and the field at m had already been overwritten by
  the field at m'; `alpha` and `pins` were indexed per component rather than
  per site, reading past the end of both arrays; and the precession term
  dropped the factors of `m.m` and one term entirely. Against a finite
  difference of the right hand side the product now agrees to 2e-8, where it
  was out by a factor of 6e4 and did not even point the same way. The example
  the issue cites runs in 0.7 s, the same as the default path
* **The atomistic Jacobian path never ran at all.** It defined
  `sundials_jtn`, which nothing calls, referenced `compute_llg_jtimes` through
  the atomistic `clib`, which does not export it, and carried a leftover debug
  `print`. It is now `sundials_jtimes`, calls `common_clib`, and is tested
* **The preconditioner is gone.** Its solve was the identity, so it only added
  a copy per GMRES iteration, and its setup function was attached with the
  wrong signature: a method whose C signature takes `self` first, cast to a
  plain `CVLsPrecSetupFn`, so CVODE called it with every argument shifted.
  Removing it does not change any result. A real preconditioner is worth
  having and is noted as future work
* **The atomistic exchange now enters the Jacobian.** It set `jac = False`,
  so the stiffest term in the problem was missing from it, and the Newton
  iteration suffered for it: on a 30x30x4 lattice the integration took 150
  steps without it and 18 with it. Nothing had to be derived, the exchange
  field being linear in the spins, so its contribution is the field evaluated
  at the direction of differentiation. The micromagnetic class already did
  this. The demagnetising field stays out of the Jacobian, which is a
  deliberate trade: it costs a transform per evaluation, 2.84 s against
  1.70 s on a 30x30x10 mesh of 1 nm cells, for the same 361 steps
* **The effective field is reused across a Newton solve.** The
  Jacobian-vector product recomputed the field at m on every call, although
  CVODE asks for many products at the same time and state; it is now computed
  once per state, by `DriverBase.effective_field_at`. This is what makes the
  analytical Jacobian worth turning on: 20 ps on the mesh above take 2.21 s
  with the difference quotient, 2.42 s with the analytical product
  recomputing the field, and 1.70 s with it reused. `use_jac` stays `False`
  by default
* New `tests/test_jacobian.py`, comparing the analytical product against a
  finite difference of the right hand side, which is the one reference that
  cannot be wrong in the same way. It also records two gaps that are design
  questions rather than bugs: which interactions the `jac` flag admits, and
  that the stabilisation term is differentiated only for a positive
  `default_c`, while a negative one selects `c = 6|dm/dt|` and is the
  atomistic default
* **`relax()` says when the system has not relaxed** (issue #118). Reaching
  `max_steps` with `dmdt` still above `stopping_dmdt` used to return quietly,
  and the result is indistinguishable from a relaxed one: the caller gets a
  magnetisation either way. It now raises a `RuntimeWarning` naming both
  numbers, which `warnings.simplefilter('error', RuntimeWarning)` turns into
  an exception. A second warning covers the case where `max_steps` is below
  the current step, so the loop never runs at all. This immediately caught
  `test_skx_num_micromagnetic`, which was capped at 1000 steps but needs about
  1123, so it had been measuring an unrelaxed state; its cap is raised and the
  skyrmion number it asserts is unchanged
* **Hexagonal meshes with a square arrangement can be periodic in x again**
  (issue #129, open since 2019). `init_grid` asks which neighbouring hexagons
  have already been built, so that a shared corner is not created twice, and
  those lookups used the periodic index: for the first column they pointed at
  the last one, which is built later in the same row, so the lookup ran off
  the end of the list and raised `IndexError`. The non-periodic lookup is also
  the correct one geometrically, since the two edges of a periodic mesh are
  far apart in space and share no corner, and it is what the diagonal
  arrangement already used. The connectivity still wraps, so the physics is
  unchanged, and non-periodic meshes are unaffected. The repository's own
  `test_hexagonal_mesh_creation_periodic_x`, marked `xfail` as unsupported,
  now passes and the marker is removed

### Removals

* **`fidimag/atomistic/exchange_new.py` and `fidimag/atomistic/field.py` are
  gone.** The first was a second `Exchange` class that nothing imported and
  that could not have worked: `compute` referred twice to a `mesh` not in its
  scope, and it offered `compute` and `in_jacobian` where the interaction
  interface wants `compute_field` and `jac`. The second was a byte for byte
  copy of `fidimag/common/field.py`, imported only by the first; the original
  stays, the tests use it. With these and `sim2fdfield` removed, `ruff`
  reports no undefined names anywhere in the package
* **`fidimag/common/sim2fdfield.py` is gone.** It converted a simulation into
  a `finitedifferencefield.Field`, importing a package that is not a
  dependency and has never been one, and it could not have run in any case:
  `numpy` was never imported and three lines referred to a `mesh` that does
  not exist in their scope. Nothing referenced it
* **`SimpleMinimiser` has moved to `sandbox/minimiser_cg/`.** It was a
  template rather than a method, a fixed multiple of the torque with no step
  length rule and no acceptance test, and it was unreachable anyway, being
  commented out of the driver tables of both simulation classes. Use
  `steepest_descent` or `hubert_minimiser`

* The Python 2 era install scripts are gone: `bin/ubuntu_install_script.sh`,
  which built through `make` and wrote PYTHONPATH into `/etc/profile.d`,
  together with the `bin/install-ubuntu-packages.sh` and
  `bin/install-python-packages.sh` it called, which installed `python-numpy`,
  `python-pyvtk` and the rest of the Python 2 packages.
  `bin/install-scikit-odes.sh`, written for SUNDIALS 2.6.2 and cloning over the
  `git://` protocol GitHub disabled in 2021, and `bin/fix_load_path_mac.py`,
  referenced by nothing, go with them. `bin/install-fftw.sh` and
  `bin/install-sundials.sh` remain, and are still what the Dockerfiles, CI and
  the install instructions use
* `doc/user_guide/develop.rst`, a superseded copy of the installation
  instructions that still described the pip build and listed pyvtk

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
