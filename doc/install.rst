How to install
===============

Install prerequisites
---------------------

On Ubuntu systems, we need to run the following commands::

  # required to compile fidimag
  apt-get install cython python-numpy
  # required for tests and running it
  apt-get install python-pytest python-pyvtk ipython python-matplotlib mayavi2

These are available in the script `bin/install-ubuntu-packages.sh` for convenience.


Install external libraries (FFTW and Sundials)
----------------------------------------------

Run the install.sh by using ::

   bash install.sh

in the fidimag/bin folder.

Install fidimag
---------------------------------------

In the fidimag folder type ::

   make

to build fidimag.

Add the fidimag to python path
---------------------------------------

Add the following to your .bashrc file::

   export PYTHONPATH=/path/to/fidimag:$PYTHONPATH

for example, suppose fidimag is in the directory of ~/work, then::

   export PYTHONPATH=~/work/fidimag:$PYTHONPATH

Add the library path to LD_LIBRARY_PATH
---------------------------------------

By default, the libraries are installed in fidimag/local, so in order
to run fidimag we need to include the libs path in LD_LIBRARY_PATH, so
please add the following to your .bashrc file::

   export LD_LIBRARY_PATH=/path/to/fidimag/local/lib:$LD_LIBRARY_PATH

for instance::

  export LD_LIBRARY_PATH=~/work/fidimag/local/lib:$LD_LIBRARY_PATH


Adding OOMMF path to the system
-------------------------------

For the tests that call OOMMF (in micro/tests), we need to tell the system where to
find it. The way it is currently set up (in util/oommf.py), we need to
find our ``oommf.tcl`` file, and add the path to it to ``.bashrc`` in this way::

  export OOMMF_PATH=/opt/oommf/oommf-1.2a5bis


Installing on Iridis
====================

Few additional notes for installing on iridis.

Loading Modules
---------------

Need to load the hg and numpy modules. This can be done by ::

    module load hg numpy

These modules will only be loaded for your current session on iridis. To have the modules automatically loaded at login, you can add them to the module initlist by ::

    module initadd hg numpy


Installing PyTest and PyVTK
---------------------------

py.test and PyVTK need to be installed locally (at /home/$USER/.local/bin/) by ::

    pip install --user pytest pyvtk

The following path needs to be added to the .bashrc file ::

    export PATH=~/.local/bin:$PATH

Cloning Repository
------------------

Can only be done with ssh keys.

Quick test
==========

Use `make test-quick` to run a set of quick tests (that don not need
OOMMF installed). The output should look something like::

  fangohr$ make test-quick
  cd tests && py.test -v -m "not slow and not run_oommf"
  =================== test session starts ====================================
  platform darwin -- Python 2.7.11, pytest-2.8.1, py-1.4.30,
  collected 59 items

  field_test.py::test_initialise_scalar PASSED
  field_test.py::test_initialise_vector PASSED
  test_2dpbc_cube.py::test_compute_field PASSED
  test_anis.py::test_anis PASSED
  test_atomistic_zeeman.py::test_zeeman PASSED
  test_citation.py::test_citation PASSED
  test_demag.py::test_demag_fft_exact PASSED
  test_demag.py::test_demag_fft_exact_oommf PASSED
  test_demag.py::test_demag_two_spin_xx PASSED
  test_demag_libraries.py::test_hexagonal_demags_1Dchain PASSED
  test_demag_libraries.py::test_cuboid_demags_1Dchain PASSED
  test_demag_libraries.py::test_cuboid_demags_2D PASSED
  test_demag_libraries.py::test_hexagonal_demags_2D PASSED
  test_dmi.py::test_dmi_1d PASSED
  test_dmi.py::test_dmi_1d_field PASSED
  test_domain_wall_cobalt.py::test_domain_wall_cobalt_fast PASSED
  test_dw_atomistic.py::test_dw_dmi_atomistic PASSED
  test_dw_dmi.py::test_dw_dmi PASSED
  test_energy.py::test_energy PASSED
  test_energy.txt SKIPPED
  test_exch.py::test_exch_1d PASSED
  test_exch.py::test_exch_1d_spatial PASSED
  test_exch_micro.py::test_init PASSED
  test_exch_micro.py::test_exch_1d PASSED
  test_exch_uniform.py::test_exch_1d PASSED
  test_exch_uniform.py::test_exch_1d_pbc PASSED
  test_exch_uniform.py::test_exch_2d PASSED
  test_exch_uniform.py::test_exch_2d_pbc2d PASSED
  test_exch_uniform.py::test_exch_3d PASSED
  test_exch_uniform.py::test_exch_energy_1d PASSED
  test_imports.py::test_has_pyvtk_installed PASSED
  test_imports.py::test_has_fidimag_installed PASSED
  test_imports.py::test_has_pytest_installed PASSED
  test_llg.py::test_sim_pin PASSED
  test_llg.py::test_sim_init_m PASSED
  test_llg.py::test_sim_init_m_fun PASSED
  test_llg.py::test_m_average PASSED
  test_llg.py::test_sim_single_spin PASSED
  test_llg_atomistic.py::test_sim_pin PASSED
  test_llg_atomistic.py::test_sim_init_m PASSED
  test_llg_atomistic.py::test_sim_init_m_fun PASSED
  test_llg_atomistic.py::test_m_average PASSED
  test_llg_atomistic.py::test_sim_single_spin_vode PASSED
  test_llg_atomistic.py::test_sim_spins PASSED
  test_llg_atomistic.py::test_sim_single_spin_sllg PASSED
  test_mesh.py::test_mesh1 PASSED
  test_micromagnetic_zeeman.py::test_H0_is_indexable_or_callable PASSED
  test_micromagnetic_zeeman.py::test_zeeman PASSED
  test_oommf_without_run.py::test_exch_field_oommf PASSED
  test_oommf_without_run.py::test_with_oommf_spatial_Ms PASSED
  test_oommf_without_run.py::test_dmi_field_oommf PASSED
  test_oommf_without_run.py::test_demag_field_oommf_large PASSED
  test_prb88_184422.py::test_prb88_184422 PASSED
  test_sky_number.py::test_skx_num PASSED
  test_stt.py::test_sst_field_1d PASSED
  test_stt_slonczewski.py::test_dynamic PASSED

  ============ 3 tests deselected by "-m 'not slow and not run_oommf'" ===========
  ============ 55 passed, 1 skipped, 3 deselected in 9.88 seconds ================



How to set up a virtual machine via vagrant
-------------------------------------------

- install vagrant on your host machine
- run::

    vagrant init ubuntu/trusty64

  to set up a basic linux machine.

- run::

    vagrant up

  to start the machine.

- ssh into the machine with X-forwarding::

    vagrant ssh -- -X

Then within the virtual machine::

  aptitude install git
  git clone https://github.com/fangohr/fidimag.git
  cd fidimag/bin
  sudo sh install-ubuntu-packages.sh
  sh install.sh
  cd ..
  make

To run the tests::

  cd /home/vagrant/fidimag/tests
  py.test

Notes:

- some tests will fail as OOMMF is not installed
- it seems that we need an active X server, on OS X, one may need to
  install XQuartz before the tests can pass (even 'import fidimag'
  failed without a working X server).

Install on OS X
===============

The inbuilt OS X gcc compiler (actually clang) doesn't have OpenMP support. A workaround is to

- install gcc5 (via homebrew, for example: ``brew install gcc --without-multilib``)
- set CC environment variable to point to that compiler: ``export CC=gcc-5``


Once this is done, run ``bin/install.sh`` which will compile fftw3 and
sundials (in a local subdirectory) using this compiler.

Also install pytest (``conda install pytest`` if using conda) and
``pyvtk`` via pip (``pip install pyvtk``).

Then run ``make``.

Set the Pythonpath so that the fidimag source is in the path.
