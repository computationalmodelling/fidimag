Monte Carlo Simulation
=======================

In the atomistic part of Fidimag, Monte Carlo based on Metropolis algorithm is integrated. 
The supported interactions include the exchange interaction, the bulk DMI, the external field and
cubic anisotropy. The total Hamiltonian of the system is therefore given by

.. math::
   \mathcal{H} = \mathcal{H}_{ex} + \mathcal{H}_{dmi} + \mathcal{H}_{ext} + \mathcal{H}_{c}

where :math:`\mathcal{H}_{ex}` is the nearest-neighbor exchange interaction, 

.. math::
   \mathcal{H}_{ex} = -J \sum_{<i,j>}\vec{m}_i \cdot \vec{m}_j

Note that the summation is taken only once for each pair and :math:`\vec{m}_i` is the unit vector of spin :math:`\vec{S}` at site :math:`i`. 
  
The Hamiltonian of bulk DMI can be expressed as, 

.. math::
   \mathcal{H}_{dmi}= \sum_{<i,j>} \vec{D}_{ij}\cdot [\vec{m}_i \times \vec{m}_j]

where :math:`\vec{D}_{ij} = D \vec{e}_{ij}` with :math:`\vec{e}_{ij}` is the unit vector between :math:`\vec{S}_{i}` and :math:`\vec{S}_{j}`.

The Hamiltonian of external field is

.. math::
   \mathcal{H}_{dmi}= - \sum_{i} \mu_s \vec{H}\cdot  \vec{m}_i

where :math:`\mu_s = g \mu_B S`.

The cubic anisotropy implemented in Fidimag is, 

.. math::
   \mathcal{H}_{c}= - \sum_{i} (K_c/2)  (m_{i,x}^4 + m_{i,y}^4 + m_{i,z}^4)

which is equivalent to the form 

.. math::
   \mathcal{H}_{c}= \sum_{i} K_c (m_{i,x}^2 m_{i,y}^2 + m_{i,y}^2 m_{i,z}^2 + m_{i,z}^2 m_{i,x}^2)


