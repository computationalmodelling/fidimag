# distutils: language = c++
import numpy
cimport numpy as np
np.import_array()

cdef extern from "m_baryakhtar_clib.h":
    void compute_laplace_m(double *m, double *field, double *Ms, double dx, double dy, double dz,
        int nx, int ny, int nz)

    void compute_relaxation_field_c(double *m, double *field, double *Ms, double chi_inv, int n)

    void compute_perp_field_c(double *m, double *field, double *field_p, int n)

    void llg_rhs_baryakhtar(double *dm_dt, double *m, double *h, double *delta_h,
        double *alpha, double beta, int *pins,
        double gamma, int nxyz, int do_precession)

    void llg_rhs_baryakhtar_reduced(double *dm_dt, double *m, double *hp, double *delta_hp,
                        double *alpha, double beta, int *pins,
                        double gamma, int nxyz, int do_precession, double default_c)

def compute_laplace_field(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            np.ndarray[double, ndim=1, mode="c"] Ms,
                            dx, dy, dz,
                            nx, ny, nz):

    compute_laplace_m(&spin[0], &field[0], &Ms[0], dx, dy, dz, nx, ny, nz)

def compute_relaxation_field(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            np.ndarray[double, ndim=1, mode="c"] Ms,
                            chi_inv,n):

    compute_relaxation_field_c(&spin[0], &field[0], &Ms[0], chi_inv, n)


def compute_perp_field(np.ndarray[double, ndim=1, mode="c"] m,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            np.ndarray[double, ndim=1, mode="c"] field_p,
                            n):

    compute_perp_field_c(&m[0], &field[0], &field_p[0], n)


def compute_llg_rhs_baryakhtar(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                               np.ndarray[double, ndim=1, mode="c"] m,
                            np.ndarray[double, ndim=1, mode="c"] h,
                            np.ndarray[double, ndim=1, mode="c"] delta_h,
                            np.ndarray[double, ndim=1, mode="c"] alpha,
                            beta,
                            np.ndarray[int, ndim=1, mode="c"] pins,
                            gamma, nxyz,
                            do_precession):

    llg_rhs_baryakhtar(&dm_dt[0], &m[0], &h[0], &delta_h[0], &alpha[0], beta, &pins[0], gamma, nxyz, do_precession)


def compute_llg_rhs_baryakhtar_reduced(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                               np.ndarray[double, ndim=1, mode="c"] m,
                            np.ndarray[double, ndim=1, mode="c"] h,
                            np.ndarray[double, ndim=1, mode="c"] delta_h,
                            np.ndarray[double, ndim=1, mode="c"] alpha,
                            beta,
                            np.ndarray[int, ndim=1, mode="c"] pins,
                            gamma, nxyz,
                            do_precession, default_c):

    llg_rhs_baryakhtar_reduced(&dm_dt[0], &m[0], &h[0], &delta_h[0], &alpha[0], beta, &pins[0], gamma, nxyz, do_precession, default_c)
