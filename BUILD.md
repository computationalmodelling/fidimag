# Building Fidimag

This guide covers building Fidimag from source using the modern build system (scikit-build-core + CMake + uv).

## Quick Start

### Prerequisites

- Python 3.9+ (Python 3.14 recommended)
- CMake 3.18+
- C compiler (GCC 11+ or Clang 13+)
- OpenMP support
- SUNDIALS 7.6.0
- FFTW 3.x
- BLAS/LAPACK

### Install Dependencies (macOS)

```bash
# Install build tools
brew install cmake gcc@15

# Install SUNDIALS and FFTW to ./local/
cd bin
bash install-fftw.sh
bash install-sundials.sh
cd ..
```

### Install Dependencies (Ubuntu/Debian)

```bash
# Install build tools
sudo apt-get update
sudo apt-get install -y cmake gcc g++ libatlas-base-dev

# Install SUNDIALS and FFTW to ./local/
# Note: We build FFTW from source to get the OpenMP version
cd bin
bash install-fftw.sh
bash install-sundials.sh
cd ..
```

### Build and Install

Using uv (recommended):

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Build and install in editable mode
uv sync

# Or with dev dependencies (pytest, pytest-cov, nbval)
uv sync --all-extras
```

Using pip (alternative):

```bash
pip install -e .
```

On macOS with Homebrew GCC:

```bash
CC=gcc-15 CXX=g++-15 CMAKE_GENERATOR="Unix Makefiles" uv sync
```

## Build System Overview

Fidimag uses a modern Python build system:

- **Build backend**: `scikit-build-core` (PEP 517/518 compliant)
- **Build system**: CMake for C/Cython extensions
- **Package manager**: `uv` with lock files for reproducibility

### How the Build Works

1. **CMake Configuration**: CMake finds native dependencies (Python, NumPy, Cython, OpenMP, SUNDIALS, FFTW, BLAS)
2. **Cython Code Generation**: Cython compiles `.pyx` files to C code in the `build/` directory (out-of-source)
3. **C Compilation**: GCC/Clang compiles the generated C code and any additional C source files
4. **Extension Linking**: Links against SUNDIALS, FFTW, BLAS, LAPACK, and OpenMP
5. **Installation**: Installs compiled `.so` files to `fidimag/extensions/`

## Development Workflow

### Clean Build

```bash
# Rebuild
uv sync --reinstall-package fidimag
```

### Incremental Build

After modifying Cython or C files:

```bash
uv sync
```

uv automatically detects changes and rebuilds only what's necessary.

### Testing

Using make (recommended):

```bash
# Quick tests (skip slow and OOMMF tests)
make test

# All tests
make test-all

# Basic tests
make test-basic

# See all available commands
make help
```

Using uv directly:

```bash
# Run all tests
cd tests
uv run pytest -v

# Run specific tests
uv run pytest -v test_atomistic.py

# Skip slow tests
uv run pytest -v -m "not slow and not run_oommf"
```

Or activate the environment manually:

```bash
source .venv/bin/activate
cd tests
pytest -v
```

## Environment Variables

### Compiler Selection

- `CC`: C compiler (e.g., `gcc-15`, `clang`)
- `CXX`: C++ compiler (e.g., `g++-15`, `clang++`)

### Build Configuration

- `CMAKE_GENERATOR`: CMake generator (use `"Unix Makefiles"` on macOS if Ninja is not available)
- `CMAKE_BUILD_TYPE`: Build type (`Release`, `Debug`, `RelWithDebInfo`)
- `CMAKE_PREFIX_PATH`: Additional paths for CMake to search for libraries

### Dependency Paths

If SUNDIALS or FFTW are installed in non-standard locations:

- `SUNDIALS_INC`: Path to SUNDIALS include directory
- `FFTW_INC`: Path to FFTW include directory

Example:

```bash
export SUNDIALS_INC=/opt/sundials/include
export FFTW_INC=/usr/local/include
uv sync
```

## Troubleshooting

### OpenMP not found

**Error**: `Could NOT find OpenMP`

**Solution**: Install a compiler with OpenMP support:

```bash
# macOS
brew install gcc@15
export CC=gcc-15 CXX=g++-15
```

### SUNDIALS not found

**Error**: `Could NOT find SUNDIALS`

**Solution**: Install SUNDIALS to `./local/`:

```bash
cd bin
bash install-sundials.sh
```

Or specify the path:

```bash
export SUNDIALS_INC=/path/to/sundials/include
```

### FFTW not found

**Error**: `Could NOT find FFTW3`

**Solution**: Install FFTW with OpenMP support to `./local/`:

```bash
cd bin
bash install-fftw.sh
```

On macOS, you can also use Homebrew:

```bash
brew install fftw
```

Note: The apt package `libfftw3-dev` does not include OpenMP support, so build from source using the install script.

### Ninja not found (macOS)

**Error**: `ninja: error: Makefile:5: expected '=', got ':'`

**Solution**: Use Make instead of Ninja:

```bash
CMAKE_GENERATOR="Unix Makefiles" uv sync
```

### Old Cython files causing issues

**Symptoms**: Import errors, symbol not found errors

**Solution**: Clean old generated files:

```bash
bash clean_cython_files.sh
uv sync --reinstall-package fidimag
```

## CI/CD

The GitHub Actions workflow uses uv for building and testing:

```yaml
- name: Install uv
  run: curl -LsSf https://astral.sh/uv/install.sh | sh

- name: Build and test
  run: |
    uv sync
    uv run pytest tests/
```

## Docker

Fidimag provides Docker images for easy deployment and testing.

### Building Docker Images

```bash
# Testing image (includes full test suite)
docker build -t fidimag:testing -f docker/testing/Dockerfile .

# Minimal Python 3 image
docker build -t fidimag:minimal -f docker/minimal-py3/Dockerfile .

# Jupyter notebook image
docker build -t fidimag:notebook -f docker/notebook/Dockerfile .
```

### Running Docker Containers

```bash
# Interactive shell
docker run -it fidimag:testing

# Run tests
docker run fidimag:testing uv run pytest tests/ -v

# Jupyter notebook (expose port 8888)
docker run -p 8888:8888 fidimag:notebook
```

### Docker Image Details

All Docker images:
- Use Ubuntu 22.04 (LTS)
- Install uv for package management
- Build SUNDIALS 7.6.0 and FFTW with OpenMP from source
- Use the modern scikit-build-core build system
- Include all runtime dependencies

## Advanced Topics

### Custom CMake Options

Edit `pyproject.toml` under `[tool.scikit-build.cmake.define]`:

```toml
[tool.scikit-build.cmake.define]
CMAKE_BUILD_TYPE = "Debug"
CMAKE_C_FLAGS = "-Wall -Wextra"
```

### Adding New Extensions

1. Create `.pyx` file in appropriate directory
2. Add to `CMakeLists.txt`:

```cmake
add_fidimag_extension("fidimag.extensions.new_module" "path/to/new_module.pyx")
```

3. Rebuild:

```bash
uv sync --reinstall-package fidimag
```

### Lock File Management

The `uv.lock` file pins all dependencies for reproducibility.

Update dependencies:

```bash
# Update all packages
uv lock --upgrade

# Update specific package
uv lock --upgrade-package numpy

# Sync to updated lock file
uv sync
```

Commit `uv.lock` to version control to ensure reproducible builds.

## Migration from setup.py

If you're migrating from an older fidimag version that used `setup.py`:

1. **Remove old build artifacts**:
   ```bash
   bash clean_cython_files.sh
   ```

2. **Uninstall old fidimag**:
   ```bash
   pip uninstall fidimag
   ```

3. **Install with new build system**:
   ```bash
   uv sync
   ```

The new build system generates Cython files **out-of-source** in the `build/` directory, whereas the old system generated them in-place. This keeps the source tree clean.

## References

- [scikit-build-core documentation](https://scikit-build-core.readthedocs.io/)
- [Cython documentation](https://cython.readthedocs.io/)
- [CMake documentation](https://cmake.org/documentation/)
- [uv documentation](https://github.com/astral-sh/uv)
