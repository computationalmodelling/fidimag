import numpy
cimport numpy as np
np.import_array()


cdef extern from "fidimag_random.h":
    ctypedef struct mt19937_state:
        pass
    void initial_rng_mt19973(mt19937_state *state)