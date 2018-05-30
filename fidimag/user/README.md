# User Extensions

The user extensions directory is here to allow for 
the user to be able to straightforwardly run compiled
code within Fidimag for performance reasons. We
consider this an advanced feature and do not recommend
trying this unless you have experience writing and
building C/Cython programs.

Some of the energy classes perform callbacks to 
user-supplied functions. Performance for this 
is generally poor, as there is an overhead to 
calling Python functions repeatedly. Hence,
we place this folder here to allow you to expose
functions written in Cython/C conveniently.

An example has been supplied. We suggest copying
the folder and modifying each of the files in 
this. Please note that we have automated the
building of the extensions, but you can only
have a single Cython .pyx file per directory, 
because a single Cython module is created in
each folder. The module that is created will have 
the name of this file.

You do not explicitly need to write an __init__.py file;
your extension will be importable immediately from
fidimag.extensions.user.$FOLDERNAME

The __init__.py file lets you do the slightly shorter:
from fidimag.user.$FOLDERNAME import *

