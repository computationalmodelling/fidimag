import numpy as np
import numpy
import ctypes
cimport numpy as np
np.import_array()

cdef extern from "dipolar.h":

    # used for demag
    ctypedef struct fft_demag_plan:
        int nx, ny, nz
        double dx, dy, dz
        int lenx, leny, lenz
        int total_length
        double *tensor_xx
        double *tensor_xy
        double *tensor_xz
        double *tensor_yy
        double *tensor_yz
        double *tensor_zz
        double *hx
        double *hy
        double *hz
        complex *Nxx
        complex *Nxy
        complex *Nxz
        complex *Nyy
        complex *Nyz
        complex *Nzz
        complex *Hx
        complex *Hy
        complex *Hz
        complex *Mx
        complex *My
        complex *Mz

    fft_demag_plan * create_plan()
    void finalize_plan(fft_demag_plan * plan)
    void init_plan(fft_demag_plan * plan, double dx, double dy, double dz, int nx,int ny, int nz)
    void compute_fields(fft_demag_plan * plan, double *spin, double *mu_s, double *field)
    void exact_compute(fft_demag_plan * plan, double *spin, double *mu_s, double *field)
    double compute_demag_energy(fft_demag_plan *plan, double *spin, double *mu_s, double *field)
    void compute_dipolar_tensors(fft_demag_plan *plan)
    void compute_demag_tensors(fft_demag_plan *plan)
    void create_fftw_plan(fft_demag_plan *plan)
    void compute_demag_tensors_2dpbc(fft_demag_plan *plan, double *tensors, double pbc_2d_error, int sample_repeat_nx, int sample_repeat_ny, double dipolar_radius)
    void fill_demag_tensors_c(fft_demag_plan *plan, double *tensors)

cdef class FFTDemag(object):
    cdef fft_demag_plan *_c_plan
    cdef int total_length
    cdef np.float64_t[:] tensor_xx_p, tensor_xy_p, tensor_xz_p, tensor_yy_p, \
                         tensor_yz_p, tensor_zz_p, hx_p, hy_p, hz_p
    cdef np.complex128_t[:] Nxx_p, Nxy_p, Nxz_p, Nyy_p, Nyz_p, Nzz_p, Hx_p, \
                            Hy_p, Hz_p, Mx_p, My_p, Mz_p
    cdef public np.ndarray tensor_xx, tensor_xy, tensor_xz, tensor_yy, \
                         tensor_yz, tensor_zz, Nxx, Nxy, Nxz, Nyy, Nyz, Nzz, \
                         Mx, My, Mz, Hx, Hy, Hz, hx, hy, hz
    #tensor_type could be 'dipolar', 'demag' or '2d_pbc'
    def __cinit__(self, dx, dy, dz, nx, ny, nz, tensor_type='dipolar'):
        self._c_plan = create_plan()
        if self._c_plan is NULL:
            raise MemoryError()
        init_plan(self._c_plan, dx, dy, dz, nx, ny, nz)

        if tensor_type == 'dipolar':
            compute_dipolar_tensors(self._c_plan)
            create_fftw_plan(self._c_plan)
        elif tensor_type == 'demag':
            compute_demag_tensors(self._c_plan)
            create_fftw_plan(self._c_plan)
        elif tensor_type == '2d_pbc':
            pass

        else:
            raise Exception("Only support options 'dipolar', 'demag' and '2d_pbc'.")

        self.total_length = int(self._c_plan.total_length)
        self.tensor_xx_p = <np.float64_t[:self.total_length]> self._c_plan.tensor_xx
        self.tensor_xx = np.asarray(self.tensor_xx_p)

        self.Nxx_p = <np.complex128_t[:self.total_length]> self._c_plan.Nxx
        self.Nxx = np.asarray(self.Nxx_p)
        self.Nxy_p = <np.complex128_t[:self.total_length]> self._c_plan.Nxy
        self.Nxy = np.asarray(self.Nxy_p)
        self.Nxz_p = <np.complex128_t[:self.total_length]> self._c_plan.Nxz
        self.Nxz = np.asarray(self.Nxz_p)
        self.Nyy_p = <np.complex128_t[:self.total_length]> self._c_plan.Nyy
        self.Nyy = np.asarray(self.Nyy_p)
        self.Nyz_p = <np.complex128_t[:self.total_length]> self._c_plan.Nyz
        self.Nyz = np.asarray(self.Nyz_p)
        self.Nzz_p = <np.complex128_t[:self.total_length]> self._c_plan.Nzz
        self.Nzz = np.asarray(self.Nzz_p)

        self.Mx_p = <np.complex128_t[:self.total_length]> self._c_plan.Mx
        self.Mx = np.asarray(self.Mx_p)
        self.My_p = <np.complex128_t[:self.total_length]> self._c_plan.My
        self.My = np.asarray(self.My_p)
        self.Mz_p = <np.complex128_t[:self.total_length]> self._c_plan.Mz
        self.Mz = np.asarray(self.Mz_p)

        self.Hx_p = <np.complex128_t[:self.total_length]> self._c_plan.Hx
        self.Hx = np.asarray(self.Hx_p)
        self.Hy_p = <np.complex128_t[:self.total_length]> self._c_plan.Hy
        self.Hy = np.asarray(self.Hy_p)
        self.Hz_p = <np.complex128_t[:self.total_length]> self._c_plan.Hz
        self.Hz = np.asarray(self.Hz_p)

        self.hx_p = <np.float64_t[:self.total_length]> self._c_plan.hx
        self.hx = np.asarray(self.hx_p)
        self.hy_p = <np.float64_t[:self.total_length]> self._c_plan.hy
        self.hy = np.asarray(self.hy_p)
        self.hz_p = <np.float64_t[:self.total_length]> self._c_plan.hz
        self.hz = np.asarray(self.hz_p)


        self.tensor_xy_p = <np.float64_t[:self.total_length]> self._c_plan.tensor_xy
        self.tensor_xy = np.asarray(self.tensor_xy_p)
        self.tensor_xz_p = <np.float64_t[:self.total_length]> self._c_plan.tensor_xz
        self.tensor_xz = np.asarray(self.tensor_xz_p)
        self.tensor_yy_p = <np.float64_t[:self.total_length]> self._c_plan.tensor_yy
        self.tensor_yy = np.asarray(self.tensor_yy_p)
        self.tensor_yz_p = <np.float64_t[:self.total_length]> self._c_plan.tensor_yz
        self.tensor_yz = np.asarray(self.tensor_yz_p)
        self.tensor_zz_p = <np.float64_t[:self.total_length]> self._c_plan.tensor_zz
        self.tensor_zz = np.asarray(self.tensor_zz_p)


    def print_tensor(self):
        for k in range(self._c_plan.lenz):
            for j in range(self._c_plan.leny):
                for i in range(self._c_plan.lenx):
                    id = self._c_plan.nx*self._c_plan.ny
                    # print(tensor[])


    def free(self):
        self.__dealloc__()

    def __dealloc__(self):
        if self._c_plan is not NULL:
            finalize_plan(self._c_plan)
            self._c_plan = NULL

    def compute_tensors_2dpbc(self, np.ndarray[double, ndim=1, mode="c"] tensors,
                            pbc_2d_error, sample_repeat_nx, sample_repeat_ny, dipolar_radius):
        compute_demag_tensors_2dpbc(self._c_plan, &tensors[0], pbc_2d_error,
                            sample_repeat_nx, sample_repeat_ny, dipolar_radius)

    def fill_demag_tensors(self, np.ndarray[double, ndim=1, mode="c"] tensors):
        fill_demag_tensors_c(self._c_plan, &tensors[0])
        create_fftw_plan(self._c_plan)

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
