import numpy
cimport numpy as np
np.import_array()

cdef extern from "dipolar.h":
    
    # used for demag
    ctypedef struct fft_demag_plan:
        pass

    fft_demag_plan * create_plan()
    void finalize_plan(fft_demag_plan * plan)
    void init_plan(fft_demag_plan * plan, double dx, double dy, double dz, int nx,int ny, int nz)
    void compute_fields(fft_demag_plan * plan, double *spin, double *mu_s, double *field)
    void exact_compute(fft_demag_plan * plan, double *spin, double *mu_s, double *field)
    double compute_demag_energy(fft_demag_plan *plan, double *spin, double *mu_s, double *field)
    void compute_dipolar_tensors(fft_demag_plan *plan)
    void compute_demag_tensors(fft_demag_plan *plan)
    void create_fftw_plan(fft_demag_plan *plan)
    void compute_demag_tensors_2dpbc(fft_demag_plan *plan, double pbc_2d_error, int sample_repeat_nx, int sample_repeat_ny, int dipolar_radius_sq)

cdef class FFTDemag(object):
    cdef fft_demag_plan *_c_plan
    
    #tensor_type could be 'dipolar', 'demag' or '2d_pbc'
    def __cinit__(self, dx, dy, dz, nx, ny, nz, tensor_type='dipolar', options={}):
        self._c_plan = create_plan()
        if self._c_plan is NULL:
            raise MemoryError()
        init_plan(self._c_plan, dx, dy, dz, nx, ny, nz)
	
        if tensor_type == 'dipolar':
            compute_dipolar_tensors(self._c_plan)
        elif tensor_type == 'demag':
            compute_demag_tensors(self._c_plan)
        elif tensor_type == '2d_pbc':
            pbc_2d_error = 1e-10
            sample_repeat_nx = -1
            sample_repeat_ny = -1
            asymptotic_radius = 32.0
            dipolar_radius = 10000.0
            if 'sample_repeat_nx' in options:
                sample_repeat_nx = options['sample_repeat_nx']
                sample_repeat_ny = options['sample_repeat_ny']
            if 'relative_tensor_error' in options:
                pbc_2d_error = options['relative_tensor_error']
            if 'asymptotic_radius' in options:
                asymptotic_radius = options['asymptotic_radius']
            if 'dipolar_radius' in options:
                dipolar_radius = options['dipolar_radius']
            
            compute_demag_tensors_2dpbc(self._c_plan, pbc_2d_error, sample_repeat_nx, sample_repeat_ny, dipolar_radius)
        
        else:
            raise Exception("Only support options 'dipolar', 'demag' and '2d_pbc'.")

        create_fftw_plan(self._c_plan)
	   

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
