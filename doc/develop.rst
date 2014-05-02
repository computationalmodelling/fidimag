

Develop note
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





