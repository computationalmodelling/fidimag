=============
API reference
=============

This section is generated from the docstrings in the source, and lists the
classes and functions of the three Fidimag packages, with their arguments. Use
the search box to find something by name: searching for ``minimise``,
``NEBM_Geodesic`` or ``Demag`` will take you to the page for that object.

``fidimag.common``
    Everything that does not belong to a particular model: the meshes, the
    drivers and minimisers, the chain methods, the integrators, and the data
    and VTK saving.

``fidimag.atomistic``
    The discrete spin model, its ``Sim`` class and its energy classes.

``fidimag.micro``
    The micromagnetic model, its ``Sim`` class and its energy classes.

The physics behind the methods is described in the
:doc:`../physics_num_methods/index` section, which is the better place to
start; this one is for looking up the arguments of something you already know
you want.

.. autosummary::
   :toctree: generated
   :recursive:

   fidimag.common
   fidimag.atomistic
   fidimag.micro
