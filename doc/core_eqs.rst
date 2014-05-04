

Core equations
===============

Landau-Lifshitz-Gilbert (LLG) equation
---------------------------------------
The dynamic of magnetic moment :math:`\vec{\mu}_i` is governed by LLG equation,

.. math::
   \frac{\partial \vec{S}_i}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \vec{S}_i \times (\vec{H}_i + \alpha \vec{S}_i \times \vec{H}_i) ]

where :math:`\vec{S}_i` is the unit vector of magnetic moment, 

.. math::
   \vec{S}_i=\frac{\vec{\mu}_i}{\mu_s},

and :math:`\vec{\mu}_s = |\vec{\mu}_i|`,  the effective field :math:`\vec{H}_i` is defined as

.. math::
   \vec{H}_i = - \frac{1}{\mu_s} \frac{\partial \mathcal{H}}{\partial \vec{S}_i}.



Interactions
-----------------
In the level of atomic moments, :math:`\vec{\mu}_s`, originated from the angular momentum of electrons in atoms, could be employed as the base unit in simulations. There are several typical interactions between magnetic moments. The total Hamiltonian is the summation of them.

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

Anisotropy 
~~~~~~~~~~~
The Hamiltonian for uniaxial anisotropy with easy axis in x direction is expressed as,

.. math::
   \mathcal{H}_{an} = - D \sum_i (\vec{S}_{x,i})^2

and the corresponding field is

.. math::
   \vec{H}_{i,an} = \frac{2 D}{\mu_s} S_{x,i} \vec{e}_x

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

In the continuum limit the bulk dmi energy could be written, 

.. math::
   E_{dmi} = \int_\Omega D \vec{m} \cdot (\nabla \times \vec{m}) dx

and the effective field is

.. math::
   \vec{H}=-\frac{2D}{\mu_0 M_s} (\nabla \times \vec{m})


Zeeman energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The zeeman energy is,

.. math::
   \mathcal{H}_{dmi}= - \sum_{i} \mu_s \vec{H}_{ext}\cdot  \vec{S}_i


Basically, we will follow the above equations to write codes.
