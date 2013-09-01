

Extended equations
===================

Including the susceptibility 
-------------------------------
The dynamic equation we used including the susceptibility is,

.. math::
   \frac{\partial \vec{S}_i}{\partial t} = - \gamma \vec{S}_i \times \vec{H}_i + \alpha \gamma \vec{H}_i + \frac{1}{\chi} (1-S_i^2)\vec{S}_i
   :label: eq_llg_s

where :math:`\vec{S}_i` is defined by the ratio of magnetic moment :math:`\vec{\mu}_i` and the equalibrium magnitude of the effective magnetic moment :math:`\mu_s=M_0 a^3`, 

.. math::
   \vec{S}_i=\frac{\vec{\mu}_i}{\mu_s}.

The effective field is defined by,

.. math::
   \vec{H}_i = - \frac{1}{\mu_s}\frac{\partial \mathcal{H}}{\partial \vec{S}_i}

Therefore, the corresponding fields are,

.. math::
   \vec{H}_{i,ex} =\frac{J}{\mu_s} \sum_{<i,j>} \vec{S}_j

.. math::
   \vec{H}_{i,an} = \frac{2 D}{\mu_s} S_{x,i} \vec{e}_x



