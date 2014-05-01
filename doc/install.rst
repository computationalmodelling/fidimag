

How to install 
===============

Install external libraries.
---------------------------------------
Run the install.sh by using ::

   bash install.sh

in the pccp/libs folder.

Install pccp
---------------------------------------
In the pccp folder type ::

   make

to install the pccp.

Add the pccp to python path
---------------------------------------
Add the following to your .bashrc file::
   
   export PYTHONPATH="path_of_pccp_folder:$PYTHONPATH"
   export PYTHONPATH="path_of_pccp_parent_folder:$PYTHONPATH"

for example, suppose pccp is in the directory of ~/work, then::  

   export PYTHONPATH="~/work/pccp:$PYTHONPATH"
   export PYTHONPATH="~/work:$PYTHONPATH"



Quick test
===============
Go to the pccp/tests folder and run "py.test", a similar result is expected ::

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



