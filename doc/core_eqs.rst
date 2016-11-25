Core equations
===============

Fidimag can simulate systems using either a discrete spin formalism or a
continuum approximation of the material, i.e. micromagnetism. Spins are
described at a semi-classical level.

* Atomistic

  We describe the material by a lattice of magnetic moments
  :math:`\vec{\mu}_i=\mu_{s}\vec{S}_{i}`, with :math:`\vec{S}` as the spin
  direction (unit vector). The ordering of the atoms or molecules is given by
  the crystal structure of the magnetic solid. Currently, it is possible to
  specify a 2D/3D square lattice or a 2D hexagonal lattice. The magnetic moment
  is defined as :math:`\mu_{s}=g \mu_{B} S`, where `g` is the Landé g-factor,
  :math:`\mu_{B}` is the Bohr magneton and :math:`S` is the average spin
  (angular momentum) magnitude.

  Interactions between magnetic moments are specified using the Heisenberg
  formalism.  
|

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
------------

At the atomic level, the magnetic moment originates from the total angular
momentum of electrons in the material atoms. In ferromagnets, most of the
angular momentum comes from the spin, thus we normally just refer to this
quantity. There are multiple interactions between electrons that we can
describe using a semi-classical approximation, where the spin is treated as a
pseudo vector per every lattice site of the material. In this approximation,
magnetic interactions can be described using Heisenberg's formalism for the
exchange interaction. The total Hamiltonian for a magnetic system is the
sum of all these magnetic interactions:

.. math::
   \mathcal{H} = \mathcal{H}_{\text{ex}} + \mathcal{H}_{\text{an}} 
   + \mathcal{H}_{\text{d}} + \mathcal{H}_{\text{DMI}} + \mathcal{H}_{\text{Zeeman}}

where we have *exchange*, *anisotropy*, *dipolar interactions*,
*Dzyaloshinskii-Moriya interactions* and *Zeeman interaction*.


Exchange interaction
~~~~~~~~~~~~~~~~~~~~  

The classical Heisenberg Hamiltonian with the nearest-neighbor exchange
interaction reads

.. math::
   \mathcal{H}_{ex} = -J \sum_{<i,j>}\vec{S}_i \cdot \vec{S}_j

where the summation is performed only once for every pair of spins. The
effective field is

.. math::
   \vec{H}_{i,ex} = \frac{J}{\mu_s} \sum_{<i,j>} \vec{S}_j

|

In the continuum limit the exchange energy can be written as

.. math::
   E_{ex} = \int_{V} A (\nabla \vec{m})^2 dx

with :math:`V` as the volume of the system and :math:`A` the anisotropy constant
in :math:`\text{J m}^{-1}`. Correspondingly, the effective
field is

.. math::
   \vec{H} = \frac{2 A}{\mu_0 M_s} \nabla^2 \vec{m}

The effective field in the continuum approximation can be computed as

.. math::
  J_x = 2A \frac{\Delta y \Delta z}{\Delta x}

Note that we need the :math:`\mu_0` factor to convert units from T to A/m.

Anisotropy 
~~~~~~~~~~~

The Hamiltonian for uniaxial anisotropy with easy axis along the unitary
:math:`\hat{u}` direction is expressed as,

.. math::
   \mathcal{H}_{\text{an}} = - \mathval{K}_{u} \sum_i \left(\vec{S}_{i}\cdot\hat{u}\right)^2

with :math:`\mathval{K}_{u}` as the anisotropy constant in eV. The
corresponding field is

.. math::
   \vec{H}_{i,\text{an}} = \frac{2 \mathval{K}_{u}}{\mu_s} \left(\vec{S}_{i}\cdot\hat{u}\right)\hat{u}

|

In micromagnetics, the uniaxial anisotropy energy of the system is defined as

.. math::
   E_{anis} = \int_{V} K_{u} [ 1 - (\vec{m} \cdot \hat{u})^2 ]\, dV

with :math:`K_{u}` as the anisotropy constant in :math:`\text{J m}^{-3}`. The
effective field reads

.. math::
   \vec{H}=\frac{2 K_{u}}{\mu_0 M_s} \left(\vec{m} \cdot \hat{u}\right) \hat{u}

Dipolar interaction
~~~~~~~~~~~~~~~~~~~

The Hamiltonian for dipolar interactions is defined as

.. math::
   \mathcal{H}_{d}=-\frac{\mu_0 \mu_{s}^{2}}{4\pi} \sum_{i<j}
   \frac{3 (\vec{S}_i\cdot \hat{r}_{ij})(\vec{S}_j\cdot \hat{r}_{ij}) - \vec{S}_i \cdot \vec{S}_j}{\vec{r}_{ij}^3} 

with :math:`\vec{r}_{ij}` the spatial vector pointing from the :math:`i`-th to
the :math:`j`-th lattice site.  The effective field is

.. math::
   \vec{H}_{i,d} =\frac{\mu_0 \mu_{s}}{4\pi}\sum_{i \neq j}\frac{3 \hat{r}_{ij} (\vec{S}_j\cdot \hat{r}_{ij}) 
   - \vec{S}_j}{\vec{r}_{ij}^3}


Dzyaloshinskii-Moriya interaction (DMI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DMI is an antisymmetric, anisotropic exchange coupling between spins (magnetic moments), 

.. math::
   \mathcal{H}_{\text{DMI}}= \sum_{<i,j>} \vec{D}_{ij}\cdot [\vec{S}_i \times \vec{S}_j]

Noting that
:math:`\vec{a}\cdot(\vec{b}\times\vec{c})=(\vec{a}\times\vec{b})\cdot\vec{c}`,
the effective field can be computed as

.. math::
   \vec{H}_i = - \frac{1}{\mu_s} \frac{\partial \mathcal{H}}{\partial \vec{S}_i} = \frac{1}{\mu_s}  \sum_{<i,j>} \vec{D}_{ij}\times\vec{S}_j

For bulk materials :math:`\vec{D}_{ij} = D \vec{r}_{ij}` and for interfacial DMI one has :math:`\vec{D}_{ij} = D \vec{r}_{ij} \times \vec{e}_z`, in both cases the vector :math:`\vec{D}_{ij}` such that :math:`\vec{D}_{ij}=-\vec{D}_{ji}`.

|

In the continuum limit the bulk DMI energy is written as 

.. math::
   E_{dmi} = \int_\Omega D_a \vec{m} \cdot (\nabla \times \vec{m}) dx

where :math:`D_a = -D/a^2` and the effective field is

.. math::
   \vec{H}=-\frac{2 D_a}{\mu_0 M_s} (\nabla \times \vec{m})



For the interfacial case, the effective field becomes,

.. math::
   \vec{H}=\frac{2 D}{M_s a^2} (\vec{e}_x \times \frac{\partial \vec{m}}{\partial y} - \vec{e}_y \times \frac{\partial \vec{m}}{\partial x} )

Compared with the effective field [PRB 88 184422]

.. math::
   \vec{H}=\frac{2 D_a}{\mu_0 M_s} ((\nabla \cdot \vec{m}) \vec{e}_z - \nabla m_z)

where :math:`D_a = D/a^2`. Notice that there is no negative sign for the interfacial case.


.. Similar to the exchange case, the effective field in the continuum case
.. can be computed by the same codes with 

.. .. math::
..  D_x = D \Delta y \Delta z

.. Also, note that we needs the factor of :math:`\mu_0` to convert the units from T to A/m.

Zeeman energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Zeeman energy is,

.. math::
   \mathcal{H}_{\text{Zeeman}}= - \sum_{i} \mu_s \vec{H}_{ext}\cdot  \vec{S}_i



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

The gyromagnetic ratio of a free electron is :math:`\gamma = 1.76\times10^{11}\,\text{rad Hz T}^{-1}`.

* Micromagnetics

In the micromagnetic limit, the equation has a similar structure

.. math::
   \frac{\partial \vec{m}}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \vec{m} \times (\vec{H} + \alpha \vec{m} \times \vec{H}) ]

where :math:`0\leq\alpha\leq 1` is the Gilbert damping constant and
:math:`\gamma` is the Gilbert gyromagnetic ratio (which sets the time scale).
The effective field :math:`\vec{H}` for this case is defined as

.. math::
   \vec{H} = - \frac{1}{\mu_{0}M_{s}} \frac{\partial \mathcal{H}}{\partial \vec{m}}.

The Gilbert gyromagnetic ratio of a free electron is :math:`\gamma = 2.21\times10^{5}\,\text{Hz T}^{-1}`.
