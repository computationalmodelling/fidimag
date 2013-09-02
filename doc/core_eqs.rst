

Core equations
===============

Interactions
-----------------
In the level of atomic moments, :math:`\vec{\mu}_s`, originated from the angular momentum of electrons in atoms, could be employed as the base unit in simulations. There are several typical interactions between magnetic moments. The total Hamiltonian is the summation of them.

.. math::
   \mathcal{H} = \mathcal{H}_{ex} + \mathcal{H}_{an} + \mathcal{H}_d


Exchange interaction
~~~~~~~~~~~~~~~~~~~~  
The classical Heisenberg Hamiltonian with the nearest-neighbor exchange interaction, 

.. math::
   \mathcal{H}_{ex} = -J \sum_{<i,j>}\vec{S}_i \vec{S}_j

where the summation is taken only once for each pair, and :math:`\vec{S}` is the unit vector of magnetic moment,

.. math::
   \vec{S}_i=\frac{\vec{\mu}_i}{\mu_i}.

Anisotropy 
~~~~~~~~~~~
The Hamiltonian for uniaxial anisotropy with easy axis in x direction is expressed as,

.. math::
   \mathcal{H}_{an} = - D \sum_i (\vec{S}_{x,i})^2

Dipolar interaction
~~~~~~~~~~~~~~~~~~~
The Hamiltonian for dipolar interaction is,

.. math::
   \mathcal{H}_{d}=-\frac{\mu_0}{4\pi}\sum_{i<j}\frac{3 (\vec{\mu}_i\cdot \vec{e}_{ij})(\vec{\mu}_j\cdot \vec{e}_{ij}) - \vec{\mu}_i \cdot \vec{\mu}_j}{r_{ij}^3} 

Dzyaloshinskii-Moriya interaction (DMI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DMI is an antisymmetric, anisotropic exchange coupling beteen spins (magnetic moments), 

.. math::
   \mathcal{H}_{dmi}=\sum_{<i,j>} \vec{D}_{ij}\cdot [\vec{S}_i \times \vec{S}_j]

Note that :math:`\vec{a}\cdot(\vec{b}\times\vec{c})=(\vec{a}\times\vec{b})\cdot\vec{c}`, the effective field can be computed by

.. math::
   \vec{H}_i = - \frac{1}{\mu_i}\frac{\partial \mathcal{H}}{\partial \vec{S}_i} =\sum_{<i,j>} \frac{1}{\mu_i} \vec{D}_{ij}\times\vec{S}_j

since :math:`\vec{D}_{ij}` is a vector such that :math:`\vec{D}_{ij}=-\vec{D}_{ji}`, we need to take care of the sign of the field according to the relative position of spins in implementation. 

Landau-Lifshitz-Gilbert (LLG) equation
---------------------------------------
The dynamic of magnetic moments are governed by LLG equation,

.. math::
   \frac{\partial \vec{S}_i}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \vec{S}_i \times (\vec{H}_i + \alpha \vec{S}_i \times \vec{H}_i) ]

where :math:`\vec{H}_i` is the effective field,

.. math::
   \vec{H}_i = - \frac{1}{\mu_i}\frac{\partial \mathcal{H}}{\partial \vec{S}_i}

Therefore, the corresponding fields are,

.. math::
   \vec{H}_{i,ex} =\frac{J}{\mu_i} \sum_{<i,j>} \vec{S}_j

.. math::
   \vec{H}_{i,an} = \frac{2 D}{\mu_i} S_{x,i} \vec{e}_x

.. math::
   \vec{H}_{i,d} =\frac{\mu_0}{4\pi}\sum_{i \neq j}\frac{3 \vec{e}_{ij} (\vec{\mu}_j\cdot \vec{e}_{ij}) - \vec{\mu}_j}{r_{ij}^3}
   :label: eq_h_d

Basically, we will follow the above equations to write codes.
