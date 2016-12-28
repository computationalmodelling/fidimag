Installation Instructions
=========================


Please use these instructions to build and run Fidimag on Linux or
OS X. We do not currently support Windows, as none of the developers
use this as their operating system, though please make a pull request
with instructions or let us know if you are able to get Fidimag
working on a Windows machine.

Ubuntu
------

Users can run a quick convenience script in the folder Fidimag/bin
by running the command::
    sudo bash ubuntu_install_script.sh

Then, follow the instructions in 'All Systems' below.


Other Linux
-----------

Please install FFTW and Sundials using the scripts below, or use your package
manager to do so.   

* install-fftw.sh
* install-sundials.sh

We also need a number of Python packages:
  
  * numpy
  * scipy
  * cython
  * pytest
  * matplotlib
  * ipywidgets
  * pyvtk
  * ipython

These can be installed through the pip package manager with::

    pip install numpy scipy cython pytest matplotlib ipywidgets pyvtk ipython

You will need a relatively recent installation of CMake (> version 3) to use the Sundials script. You may also need to install development versions of

* BLAS
* LAPACK

though many Linux distributions come with these.

Then, follow the instructions in 'All Systems' below.


OS X
----

OS X has not shipped with GCC since the release of OS X Mavericks. You therefore need to install this, as the version of clang which ships does not support OpenMP. We advise that you use the brew package manager, and install gcc5. We also strongly advise that you install the Anaconda Python distribution - we do not test against the version of Python that comes with OS X.

Once you have done this, you need to specify the compiler you are using::

    export CC=gcc-5

You can then follow the same installation instructions as for 'Other Linux', but don't worry about BLAS and LAPACK as Anaconda takes care of these for you.

Then, follow the instructions in 'All Systems' below.


All systems
-----------

Once you've built Fidimag, you need to add the libraries to your LD_LIBRARY_PATH, so that Fidimag can find them. If you installed SUNDIALS and FFTW using our scripts, you can do this with the command::

    export LD_LIBRARY_PATH=/path/to/fidimag/local/lib:$LD_LIBRARY_PATH

You may want to add this and another command to the file ~/.bashrc on Linux or ~/.bash_profile on OS X::

    export PYTHONPATH=/path/to/fidimag:$PYTHONPATH

Adding Fidimag to your PYTHONPATH allows fidimag to be imported in Python from any directory.

If you want to check everything has worked correctly, try the command 'make test' from the fidimag directory - if all tests pass, then you have a working installation!

OOMMF
-----

Some additional tests check Fidimag against OOMMF. To run these, you need a working OOMMF installation, and you need need to tell the system where to
find it. You can do this by setting the environment variable to the directory containing oommf.tcl::

    export OOMMF_PATH=/path/to/folder/containing/OOMMF
