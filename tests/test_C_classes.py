from fidimag.extensions.c_clib import PyExchangeEnergy
import numpy as np

A = np.ones(9, dtype=np.double)
print(A)
Exch = PyExchangeEnergy(A)
