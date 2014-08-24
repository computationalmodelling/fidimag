import numpy
cimport numpy as np
np.import_array()

cdef extern from "baryakhtar_clib.h":
    void compute_laplace_m(double *m, double *field, double A, double dx, double dy, double dz,
        int nx, int ny, int nz)
    
    void compute_exch_field_baryakhtar(double *m, double *field, double Me,
                        double chi_inv, double A, double dx, double dy, double dz,
                        int nx, int ny, int nz)
    
    void llg_rhs_baryakhtar(double *dm_dt, double *m, double *h, double *delta_h,
        double *alpha, double beta, int *pins,
        double gamma, int nxyz, int do_procession)

def compute_laplace_field(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            dx, dy, dz,
                            nx, ny, nz):

    compute_laplace_m(&spin[0], &field[0], 1.0, dx, dy, dz, nx, ny, nz)

def compute_exchange_field_baryakhtar(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            Me, chi_inv, 
                            A, 
                            dx, dy, dz,
                            nx, ny, nz):

    compute_exch_field_baryakhtar(&spin[0], &field[0], Me, chi_inv, A, dx, dy, dz, nx, ny, nz)


def compute_llg_rhs_baryakhtar(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                               np.ndarray[double, ndim=1, mode="c"] m,
                            np.ndarray[double, ndim=1, mode="c"] h,
                            np.ndarray[double, ndim=1, mode="c"] delta_h,
                            np.ndarray[double, ndim=1, mode="c"] alpha,
                            beta, 
                            np.ndarray[int, ndim=1, mode="c"] pins,
                            gamma, nxyz, 
                            do_procession):

    llg_rhs_baryakhtar(&dm_dt[0], &m[0], &h[0], &delta_h[0], &alpha[0], beta, &pins[0], gamma, nxyz, do_procession)