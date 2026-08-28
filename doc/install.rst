Installation
============


Please use these instructions to build and run Fidimag on Linux, OS X, or Windows WSL.
Fidimag can work in Windows using a Python disctribution such as Anaconda or Mamba, with a C/C++ compiler.

Fidimag is built with `scikit-build-core` and CMake, and the Python environment
is managed with `uv <https://docs.astral.sh/uv/>`. A more detailed description
of the build, including how to pass options to CMake, is in ``BUILD.md``.


Linux
-----

Install `uv`, if it is not already available:

.. code-block:: bash

    curl -LsSf https://astral.sh/uv/install.sh | sh

Now clone the repository and `cd` into it.

.. code-block:: bash

    git clone git@github.com:computationalmodelling/fidimag.git
    cd fidimag

Fidimag needs FFTW and SUNDIALS, and you may also need the development versions
of

* BLAS
* LAPACK

though many Linux distributions come with these. You will also need CMake (>=
3.18) and a C/C++ compiler with OpenMP support. Using the scripts provided in
Fidimag to build FFTW and SUNDIALS into ``local/``:

.. code-block:: bash

    cd bin
    bash install-fftw.sh
    bash install-sundials.sh
    cd ..

The installation script will automatically download and build SUNDIALS v7.6.0.

Python library dependencies are specified in the `pyproject.toml` file, and the
versions that are known to work together are recorded in `uv.lock`. To create
the environment and build the C/C++ modules into it:

.. code-block:: bash

    uv sync

This installs Fidimag in editable mode, so we can make any changes to the
Python code and not install the library again. Modifying the C/Cython modules
does require building again, which is the same command:

.. code-block:: bash

    uv sync --reinstall-package fidimag

The `Makefile` has this as ``make build``. Now we can call

.. code-block:: bash

    uv run python -c "import fidimag"

Any command that needs the Fidimag environment is run in the same way, with
``uv run``, and there is no environment to activate. If you prefer to activate
it, ``uv`` creates it in ``.venv``.

Simulation scripts do not have to live in the Fidimag directory. From anywhere,
``--project`` points ``uv`` at the checkout, so a script kept with the rest of
a project can still be run against it:

.. code-block:: bash

    uv --project /path/to/fidimag run python myscript.py

The files that the simulation writes, the ``<name>.txt`` data table and the
``<name>_npys`` and ``<name>_vtks`` directories, are created in the directory
the command is run from, not in the Fidimag one.

To also install the optional dependencies, which are the ones needed to run the
tests (``dev``) and to build this documentation (``docs``):

.. code-block:: bash

    uv sync --all-extras

Note that ``uv sync --reinstall-package fidimag`` only installs the main
dependencies, so run ``uv sync --all-extras`` again after rebuilding if you
need the extras.

If you want to check everything has worked correctly, try the command ``make
test`` from the fidimag directory - if all tests pass, then you have a working
installation!


OS X
----

The version of clang that ships with OS X does not support OpenMP, so a GCC
installation is needed. We advise that you use the brew package manager:

.. code-block:: bash

    brew install cmake gcc@15

Then specify the compilers when building, since the build has to use the same
ones for the C and C++ modules:

.. code-block:: bash

    CC=gcc-15 CXX=g++-15 CMAKE_GENERATOR="Unix Makefiles" uv sync

You can then follow the same installation instructions as for Linux, but don't
worry about BLAS and LAPACK as they are taken care of for you.


Troubleshooting
---------------

If there is a problem with finding C/C++ sundials and fftw libraries, it is necessary to update the corresponding env variable

.. code-block:: bash

    export LD_LIBRARY_PATH=/path/to/fidimag/local/lib:$LD_LIBRARY_PATH

The build looks in ``local/`` first and only then falls back to a system
installation, so a system FFTW is not picked up by mistake when both are
present.

If a build fails after the C/Cython sources have changed, or after switching
branches, the old artefacts can be removed with

.. code-block:: bash

    make clean


OOMMF
-----

Some additional tests check Fidimag against OOMMF. To run these, you need a working OOMMF installation, and you need need to tell the system where to
find it. You can do this by setting the environment variable to the directory containing oommf.tcl

.. code-block:: bash

    export OOMMF_PATH=/path/to/folder/containing/OOMMF
