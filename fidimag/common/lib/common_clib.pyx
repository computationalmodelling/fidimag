import numpy
cimport numpy as np
np.import_array()

# -----------------------------------------------------------------------------

cdef extern from "common_clib.h":

    void llg_rhs(double * dm_dt, double * spin,
                 double *h, double *alpha, int *pins,
                 double gamma, int n, int do_precession, double default_c)


    void llg_rhs_jtimes(double *jtn, double *m, double *h,
                        double *mp, double *hp, double *alpha, int *pins,
                        double gamma, int n,
                        int do_precession, double default_c)

# -----------------------------------------------------------------------------

def compute_llg_rhs(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] alpha,
                np.ndarray[int, ndim=1, mode="c"] pins,
                gamma, n, do_precession, default_c):
    llg_rhs(&dm_dt[0], &spin[0], &field[0], &alpha[0], &pins[0],
            gamma, n, do_precession, default_c)


def compute_llg_jtimes(np.ndarray[double, ndim=1, mode="c"] jtn,
                np.ndarray[double, ndim=1, mode="c"] m,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] mp,
                np.ndarray[double, ndim=1, mode="c"] field_p,
                np.ndarray[double, ndim=1, mode="c"] alpha,
                np.ndarray[int, ndim=1, mode="c"] pins,
                gamma, n, do_precession, default_c):
    llg_rhs_jtimes(&jtn[0], &m[0], &field[0], &mp[0], &field_p[0],
                   &alpha[0], &pins[0], gamma, n, do_precession, default_c)
