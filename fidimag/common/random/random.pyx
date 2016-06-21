#import numpy
cimport numpy as np
np.import_array()


cdef extern from "fidimag_random.h":
    ctypedef struct mt19937_state:
        pass
    
    mt19937_state *create_mt19937_state()
    void finalize_mt19937_state(mt19937_state *state)

    void initial_rng_mt19973(mt19937_state *state, int seed)
    double random_double(mt19937_state *state)
    void gauss_random_vector(mt19937_state *state, double *x, int n)

cdef extern from "time.h":
    ctypedef int time_t
    time_t time(time_t *timer)

cdef class rng_mt19937(object):
    cdef mt19937_state *_c_state
    cdef public int seed
    def __init__(self, seed=None):
        if seed:
            self.seed = int(seed)
        else:
            self.seed = time(NULL)
                
        self._c_state = create_mt19937_state()
        if self._c_state is NULL:
            raise MemoryError()
        
        initial_rng_mt19973(self._c_state, self.seed)

    def set_seed(self, seed):
        self.seed = int(seed)
        initial_rng_mt19973(self._c_state, self.seed)

    def random(self):
        """
            return a random number in [0,1]
        """
        return random_double(self._c_state)
    
    def fill_vector_gaussian(self, np.ndarray[np.float64_t, ndim=1] vector):
        gauss_random_vector(self._c_state, &vector[0], vector.shape[0])
    


    def __dealloc__(self):
        if self._c_state is not NULL:
            finalize_mt19937_state(self._c_state)
            self._c_state = NULL