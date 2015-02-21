

How to install
===============

Install external libraries.
---------------------------------------
Run the install.sh by using ::

   bash install.sh

in the fidimag/libs folder.

Install fidimag
---------------------------------------
In the fidimag folder type ::

   make

to install the fidimag.

Add the fidimag to python path
---------------------------------------
Add the following to your .bashrc file::

   export PYTHONPATH="path_of_fidimag_folder:$PYTHONPATH"
   export PYTHONPATH="path_of_fidimag_parent_folder:$PYTHONPATH"

for example, suppose fidimag is in the directory of ~/work, then::

   export PYTHONPATH="~/work/fidimag:$PYTHONPATH"
   export PYTHONPATH="~/work:$PYTHONPATH"

Add the library path to LD_LIBRARY_PATH
---------------------------------------
In default, the libraries are installed in fidimag/libs, so in order
to run fidimag we need to include the libs path in LD_LIBRARY_PATH, so
please add the following to your .bashrc file::

   export LD_LIBRARY_PATH="path_of_fidimag_folder/libs:$LD_LIBRARY_PATH"

for instance::

  export LD_LIBRARY_PATH="~/work/fidimag/libs:$LD_LIBRARY_PATH"


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
===============
Go to the fidimag/tests folder and type "py.test", a similar result is expected ::

   ============================= test session starts ==============================
   platform linux2 -- Python 2.7.5 -- pytest-2.4.2
   collected 20 items

   test_anis.py .
   test_demag.py ..
   test_dmi.py ..
   test_energy.py .
   test_exch.py ....
   test_mesh.py ..
   test_sim.py .......
   test_stt.py .

   ========================== 20 passed in 2.31 seconds ===========================
