import numpy
cimport numpy as np
np.import_array()

     

cdef extern from "dipolar.h":
    
    # used for demag
    ctypedef struct fft_demag_plan:
        pass

    fft_demag_plan * create_plan()
    void finalize_plan(fft_demag_plan * plan)
    void init_plan(fft_demag_plan * plan, double dx, double dy, double dz, int nx,int ny, int nz, int oommf)
    void compute_fields(fft_demag_plan * plan, double *spin, double *mu_s, double *field)
    void exact_compute(fft_demag_plan * plan, double *spin, double *mu_s, double *field)
    double compute_demag_energy(fft_demag_plan *plan, double *spin, double *mu_s, double *field)


cdef class FFTDemag(object):
    cdef fft_demag_plan *_c_plan

    def __cinit__(self, dx, dy, dz, nx, ny, nz, oommf):
        self._c_plan = create_plan()
        if self._c_plan is NULL:
            raise MemoryError()
        init_plan(self._c_plan, dx, dy, dz, nx, ny, nz, oommf)

    def free(self):
        self.__dealloc__()

    def __dealloc__(self):
        if self._c_plan is not NULL:
            finalize_plan(self._c_plan)
            self._c_plan = NULL


    def compute_field(self,np.ndarray[double, ndim=1, mode="c"] spin,
                        np.ndarray[double, ndim=1, mode="c"] mu_s,
                        np.ndarray[double, ndim=1, mode="c"] field):
        compute_fields(self._c_plan, &spin[0], &mu_s[0], &field[0])

    def compute_exact(self,
                      np.ndarray[double, ndim=1, mode="c"] spin,
                      np.ndarray[double, ndim=1, mode="c"] mu_s,
                      np.ndarray[double, ndim=1, mode="c"] field):
        exact_compute(self._c_plan, &spin[0], &mu_s[0], &field[0])
        
    def compute_energy(self,
                      np.ndarray[double, ndim=1, mode="c"] spin,
                      np.ndarray[double, ndim=1, mode="c"] mu_s,
                      np.ndarray[double, ndim=1, mode="c"] field):
        
        return compute_demag_energy(self._c_plan, &spin[0], &mu_s[0], &field[0])


cdef extern from "demagcoef.h":
    double CalculateSDA00(double x, double y, double z, double dx,double dy,double dz)
    double DemagNxxAsymptotic(double x, double y, double z, double dx,double dy,double dz)
    
def compute_Nxx(x, y, z, dx, dy, dz):
    return CalculateSDA00(x,y,z,dx,dy,dz)

def compute_Nxx_asy(x, y, z, dx, dy, dz):
    return DemagNxxAsymptotic(x,y,z,dx,dy,dz)
