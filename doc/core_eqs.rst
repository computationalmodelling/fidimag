Core equations
===============

Fidimag can simulate systems using either a discrete spin formalism or a
continuum approximation of the material, i.e. micromagnetism. Spins are
described at a semi-classical level.

* Atomistic

  We describe the material by a lattice of magnetic moments
  :math:`\vec{\mu}_i=\mu_{s}\vec{S}_{i}`, with :math:`\vec{S}` as the spin
  direction (unit vector), whose ordering is given by the crystal structure of
  the magnetic solid. Currently, it is possible to specify a 2D/3D square
  lattice or a 2D hexagonal lattice. The magnetic moment is defined as
  :math:`\mu_{s}=g \mu_{B} S`, where `g` is the Landé g-factor, :math:`\mu_{B}`
  is the Bohr magneton and :math:`S` is the average spin (angular momentum)
  magnitude.

  The interactions of the magnetic moments are specified using the Heisenberg
  formalism.


* Micromagnetics

  In the continuum limit, we discretise the material as a mesh whose nodes are
  arranged in a cubic lattice and we use finite differences to evaluate the
  interactions. Instead of discrete spins, now we have a coordinate dependent
  magnetisation field whose magnitude is the magnetic moment per unit volume
  :math:`\vec{M}(\vec{r})=\mu_{s}\vec{m}/V` (:math:`\vec{m}` as a unit vector).
  Accordingly, every mesh node has assigned a magnetisation vector, which is
  the magnetisation field evaluated at that point. Because we are considering
  systems at zero temperature, it is more common to use the saturation
  magnetisation :math:`M_{s}` to describe the magnetisation field magnitude,
  i.e. :math:`\vec{M}=M_{s}\vec{m}`, where :math:`M_{s}` has units of `A/m`.

  The interactions under this approximation can be computed by taking
  the continuum limit of the interactions from the Heisenberg Hamiltonian.



Interactions
-----------------

At the level of atomic moments, :math:`\vec{\mu}_s`, originated from the
angular momentum of electrons in atoms, could be employed as the base unit in
simulations. There are several typical interactions between magnetic moments.
The total Hamiltonian is the summation of them.

.. math::
   \mathcal{H} = \mathcal{H}_{ex} + \mathcal{H}_{an} + \mathcal{H}_d + \mathcal{H}_{ext}


Exchange interaction
~~~~~~~~~~~~~~~~~~~~  
The classical Heisenberg Hamiltonian with the nearest-neighbor exchange interaction, 

.. math::
   \mathcal{H}_{ex} = -J \sum_{<i,j>}\vec{S}_i \cdot \vec{S}_j

where the summation is taken only once for each pair, so the effective field is 

.. math::
   \vec{H}_{i,ex} = \frac{J}{\mu_s} \sum_{<i,j>} \vec{S}_j


In the continuum limit the exchange energy could be written, 

.. math::
   E_{ex} = \int_\Omega A (\nabla \vec{m})^2 dx

so the corresponding effective field is

.. math::
   \vec{H} = \frac{2 A}{\mu_0 M_s} \nabla^2 \vec{m}

Once we implemented the Heisenberg exchange interaction, the effective field in the continuum case
can be computed by the same codes with 

.. math::
  J_x = 2A \frac{\Delta y \Delta z}{\Delta x}

Note that we needs the factor of :math:`\mu_0` to convert the units from T to A/m.

Anisotropy 
~~~~~~~~~~~
The Hamiltonian for uniaxial anisotropy with easy axis in x direction is expressed as,

.. math::
   \mathcal{H}_{an} = - D \sum_i (\vec{S}_{x,i})^2

and the corresponding field is

.. math::
   \vec{H}_{i,an} = \frac{2 D}{\mu_s} S_{x,i} \vec{e}_x


UniaxialAnisotropy 
~~~~~~~~~~~~~~~~~~~
The UniaxialAnisotropy energy of the magnetic system is defined as

.. math::
   E_{anis} = \int_\Omega K [ 1 - (\vec{m} \cdot \vec{u})^2 ] dx

the the effective field is,

.. math::
   \vec{H}=\frac{2 K}{\mu_0 M_s} (\vec{m} \cdot \vec{u}) \vec{u}

Dipolar interaction
~~~~~~~~~~~~~~~~~~~
The Hamiltonian for dipolar interaction is,

.. math::
   \mathcal{H}_{d}=-\frac{\mu_0}{4\pi}\sum_{i<j}\frac{3 (\vec{\mu}_i\cdot \vec{e}_{ij})(\vec{\mu}_j\cdot \vec{e}_{ij}) - \vec{\mu}_i \cdot \vec{\mu}_j}{r_{ij}^3} 

.. math::
   \vec{H}_{i,d} =\frac{\mu_0}{4\pi}\sum_{i \neq j}\frac{3 \vec{e}_{ij} (\vec{\mu}_j\cdot \vec{e}_{ij}) - \vec{\mu}_j}{r_{ij}^3}
   :label: eq_h_d


Dzyaloshinskii-Moriya interaction (DMI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DMI is an antisymmetric, anisotropic exchange coupling beteen spins (magnetic moments), 

.. math::
   \mathcal{H}_{dmi}= \sum_{<i,j>} \vec{D}_{ij}\cdot [\vec{S}_i \times \vec{S}_j]

Note that :math:`\vec{a}\cdot(\vec{b}\times\vec{c})=(\vec{a}\times\vec{b})\cdot\vec{c}`, the effective field can be computed by

.. math::
   \vec{H}_i = - \frac{1}{\mu_s} \frac{\partial \mathcal{H}}{\partial \vec{S}_i} = \frac{1}{\mu_s}  \sum_{<i,j>} \vec{D}_{ij}\times\vec{S}_j

For bulk materials :math:`\vec{D}_{ij} = D \vec{r}_{ij}` and for interfacial DMI one has :math:`\vec{D}_{ij} = D \vec{r}_{ij} \times \vec{e}_z`, in both cases the vector :math:`\vec{D}_{ij}` such that :math:`\vec{D}_{ij}=-\vec{D}_{ji}`.


In the continuum limit the bulk DMI energy could be written, 

.. math::
   E_{dmi} = \int_\Omega D_a \vec{m} \cdot (\nabla \times \vec{m}) dx

where :math:`D_a = -D/a^2` and the effective field is

.. math::
   \vec{H}=-\frac{2 D_a}{\mu_0 M_s} (\nabla \times \vec{m})



For the interfacial case, the effective field thus becomes,

.. math::
   \vec{H}=\frac{2 D}{M_s a^2} (\vec{e}_x \times \frac{\partial \vec{m}}{\partial y} - \vec{e}_y \times \frac{\partial \vec{m}}{\partial x} )

Compared with the effective field [PRB 88 184422]

.. math::
   \vec{H}=\frac{2 D_a}{\mu_0 M_s} ((\nabla \cdot \vec{m}) \vec{e}_z - \nabla m_z)

we have :math:`D_a = D/a^2`, note that there is no negative sign for the interfacial case.


.. Similar to the exchange case, the effective field in the continuum case
.. can be computed by the same codes with 

.. .. math::
..  D_x = D \Delta y \Delta z

.. Also, note that we needs the factor of :math:`\mu_0` to convert the units from T to A/m.

Zeeman energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The zeeman energy is,

.. math::
   \mathcal{H}_{dmi}= - \sum_{i} \mu_s \vec{H}_{ext}\cdot  \vec{S}_i


Basically, we will follow the above equations to write codes.


Landau-Lifshitz-Gilbert (LLG) equation
---------------------------------------

* Atomistic

For the discrete theory, the dynamics of the magnetic moments is governed by
the LLG equation,

.. math::
   \frac{\partial \vec{S}_i}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \vec{S}_i \times (\vec{H}_i + \alpha \vec{S}_i \times \vec{H}_i) ]

where :math:`\vec{\mu}_s = |\vec{\mu}_i|`, :math:`0\leq\alpha\leq 1` is the
Gilbert damping constant, :math:`\gamma` is the Gilbert gyromagnetic ratio
(which sets the time scale) and the effective field :math:`\vec{H}_i` is
defined using the Hamiltonian :math:`\mathcal{H}` as

.. math::
   \vec{H}_i = - \frac{1}{\mu_s} \frac{\partial \mathcal{H}}{\partial \vec{S}_i}.

The gyromagnetic ratio of a free electron is :math:`\gamma = 1.76e11\,\text{Hz T}^{-1}`.

* Micromagnetics

In the micromagnetic limit, the equation has a similar structure

.. math::
   \frac{\partial \vec{m}}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \vec{m} \times (\vec{H} + \alpha \vec{m} \times \vec{H}) ]

where :math:`0\leq\alpha\leq 1` is the Gilbert damping constant and
:math:`\gamma` is the Gilbert gyromagnetic ratio (which sets the time scale).
The effective field :math:`\vec{H}` for this case is defined as

.. math::
   \vec{H} = - \frac{1}{\mu_{0}M_{s}} \frac{\partial \mathcal{H}}{\partial \vec{m}}.

The gyromagnetic ratio of a free electron is :math:`\gamma = 2.21e5\,\text{Hz T}^{-1}`.
