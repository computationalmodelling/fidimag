

Developing note
===============

Data structure
---------------------------------------
The first thing we may want to know is that how we label each spin in space,
which can be peered from the index method in Mesh class ::

    def index(self, i, j, k):
        idx = i*self.nyz + j*self.nz + k
        return idx

So we label the cells along z axis first then along y axis and the last is x, therefore we using the following code to loop all cells ::

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				index = nyz * i + nz * j + k;
				//do something ...
			}
		}
	}

Exchange and DMI energy
--------------------------------------------------
Although we can compute the total energy according to the direct definitions, however, 
we also can compute them using the corresponding effective fields, the energy density is

.. math::
   \mathcal{H}_i = - \vec{S}_i \cdot \vec{H}_i

where 

.. math::
   \vec{H}_i = - \frac{\partial \mathcal{H}}{\partial \vec{S}_i}


and thus the total energy is

 .. math::
   \mathcal{H} = \frac{1}{2}\sum_i \mathcal{H}_i


Anisotropy and external field energy
--------------------------------------------------
We compute them directly from the direct definitions.