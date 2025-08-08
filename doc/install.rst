Installation
============


Please use these instructions to build and run Fidimag on Linux, OS X, or Windows WSL. 
Fidimag can work in Windows using a Python disctribution such as Anaconda or Mamba, with a C/C++ compiler. 

Linux
-----

We recommend using ![mamba](https://mamba.readthedocs.io/en/latest/) (or conda). First create an environment with Python >= 3.10.

.. code-block:: bash

    mamba create -n fidimag python=3.13 -c conda-forge
    mamba activate fidimag


Now clone the repository and `cd` into it.

.. code-block:: bash

    git clone git@github.com:computationalmodelling/fidimag.git
    cd fidimag


Install FFTW and Sundials. You will need a relatively recent installation of  (> version 3) to use the Sundials script. You may also need to install development versions of

* BLAS
* LAPACK

though many Linux distributions come with these. Using the scripts provided in Fidimag:

.. code-block:: bash

    cd bin
    bash install-fftw.sh
    bash install-sundials.sh

Python library dependencies are specified in the `pyproject.toml` file. We can install the `fidimag` library in editable mode, using `pip`:

.. code-block:: bash

    pip install -e .

This will build the C/C++ modules and setup `fidimag` in our Python environment. We can make any changes to the Python code and not install the library again, unless we modified the C/C++ modules, which requires building again. Now we can simply call

.. code-block:: bash

    python -c "import fidimag"


If you want to check everything has worked correctly, try the command 'make test' from the fidimag directory - if all tests pass, then you have a working installation!

Alternatively, for development, we can install the C/C++ modules

.. code-block:: bash

    make


and link the Fidimag directory to the Python path

.. code-block:: bash

    export PYTHONPATH=/path/to/fidimag:$PYTHONPATH


Any changes to the C/C++ modules will require building only the modified codes and not all of the modules using `make`.


OS X
----

OS X has not shipped with GCC since the release of OS X Mavericks. You therefore need to install this, as the version of clang which ships does not support OpenMP. We advise that you use the brew package manager, and install gcc5. We also strongly advise that you install the Anaconda Python distribution - we do not test against the version of Python that comes with OS X.

Once you have done this, you need to specify the compiler you are using

.. code-block:: bash

    export CC=gcc-5

You can then follow the same installation instructions as for 'Other Linux', but don't worry about BLAS and LAPACK as Anaconda takes care of these for you.

Then, follow the instructions in 'All Systems' below.

Troubleshooting
---------------

If there is a problem with finding C/C++ sundials and fftw libraries, it is necessary to update the corresponding env variable

.. code-block:: bash

    export LD_LIBRARY_PATH=/path/to/fidimag/local/lib:$LD_LIBRARY_PATH


OOMMF
-----

Some additional tests check Fidimag against OOMMF. To run these, you need a working OOMMF installation, and you need need to tell the system where to
find it. You can do this by setting the environment variable to the directory containing oommf.tcl

.. code-block:: bash

    export OOMMF_PATH=/path/to/folder/containing/OOMMF
