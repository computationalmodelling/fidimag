import numpy
cimport numpy as np
np.import_array()

# -----------------------------------------------------------------------------

cdef extern from "common_clib.h":

    # From: llg.c
    void llg_rhs(double * dm_dt, double * spin,
                 double *h, double *alpha, int *pins,
                 double gamma, int n, int do_precession, double default_c)


    void llg_rhs_jtimes(double *jtn, double *m, double *h,
                        double *mp, double *hp, double *alpha, int *pins,
                        double gamma, int n,
                        int do_precession, double default_c)

    void compute_stt_field_c(double *spin, double *field,
                             double *jx, double *jy, double *jz,
                             double dx, double dy, double dz, int *ngbs, int n)

    # From: stt.c
    void llg_stt_rhs(double *dm_dt, double *m, double *h, double *h_stt,
                     double *alpha,double beta, double u0, double gamma, int n)

    void llg_stt_cpp(double *dm_dt, double *m, double *h, double *p,
                     double *alpha, int *pins, double *a_J,
                     double beta, double gamma, int n)

# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------

def compute_stt_field(np.ndarray[double, ndim=1, mode="c"] spin,
                      np.ndarray[double, ndim=1, mode="c"] field,
                      np.ndarray[double, ndim=1, mode="c"] jx,
                      np.ndarray[double, ndim=1, mode="c"] jy,
                      np.ndarray[double, ndim=1, mode="c"] jz,
                      dx, dy, dz,
                      np.ndarray[int, ndim=2, mode="c"] ngbs,
                      n
                      ):
    compute_stt_field_c(&spin[0], &field[0], &jx[0], &jy[0],&jz[0],
                        dx, dy, dz, &ngbs[0, 0], n)

def compute_llg_stt_rhs(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                        np.ndarray[double, ndim=1, mode="c"] spin,
                        np.ndarray[double, ndim=1, mode="c"] field,
                        np.ndarray[double, ndim=1, mode="c"] field_stt,
                        np.ndarray[double, ndim=1, mode="c"] alpha,
                        beta, u0, gamma, n):
    llg_stt_rhs(&dm_dt[0], &spin[0], &field[0], &field_stt[0],
                &alpha[0], beta, u0, gamma, n)



def compute_llg_stt_cpp(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                        np.ndarray[double, ndim=1, mode="c"] spin,
                        np.ndarray[double, ndim=1, mode="c"] field,
                        np.ndarray[double, ndim=1, mode="c"] p,
                        np.ndarray[double, ndim=1, mode="c"] alpha,
       	                np.ndarray[int, ndim=1, mode="c"] pin,
                        np.ndarray[double, ndim=1, mode="c"] a_J,
                        beta, gamma, n):
    llg_stt_cpp(&dm_dt[0], &spin[0], &field[0], &p[0],
                &alpha[0], &pin[0], &a_J[0], beta, gamma, n)
