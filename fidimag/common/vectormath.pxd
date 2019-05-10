cdef extern from "c_vectormath.h":
     void normalise(double *m, int *pins, int n)
     void normalise(double * a, int n)