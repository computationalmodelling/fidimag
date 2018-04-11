import numpy
# cimport numpy as np
# np.import_array()

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

    # -------------------------------------------------------------------------
    # From steepest_descent.c

    void sd_update_spin (double *spin, double *spin_last,
                         double *mxH, double *mxmxH, double *mxmxH_last, double *tau,
                         int* pins, int n)

    void sd_compute_step (double *spin, double *spin_last, double *field, double *scale,
                          double *mxH, double *mxmxH, double *mxmxH_last, double *tau,
                          int *pins, int n, int counter, double tmin, double tmax)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def compute_llg_rhs(double [:] dm_dt,
                    double [:] spin,
                    double [:] field,
                    double [:] alpha,
                    int [:] pins,
                    gamma, n, do_precession, default_c):
    llg_rhs(&dm_dt[0], &spin[0], &field[0], &alpha[0], &pins[0],
            gamma, n, do_precession, default_c)


def compute_llg_jtimes(double [:] jtn,
                       double [:] m,
                       double [:] field,
                       double [:] mp,
                       double [:] field_p,
                       double [:] alpha,
                       int [:] pins,
                       gamma, n, do_precession, default_c):
    llg_rhs_jtimes(&jtn[0], &m[0], &field[0], &mp[0], &field_p[0],
                   &alpha[0], &pins[0], gamma, n, do_precession, default_c)

# -----------------------------------------------------------------------------

def compute_stt_field(double [:] spin,
                      double [:] field,
                      double [:] jx,
                      double [:] jy,
                      double [:] jz,
                      dx, dy, dz,
                      int [:, :] ngbs,
                      n
                      ):
    compute_stt_field_c(&spin[0], &field[0], &jx[0], &jy[0],&jz[0],
                        dx, dy, dz, &ngbs[0, 0], n)

def compute_llg_stt_rhs(double [:] dm_dt,
                        double [:] spin,
                        double [:] field,
                        double [:] field_stt,
                        double [:] alpha,
                        beta, u0, gamma, n):
    llg_stt_rhs(&dm_dt[0], &spin[0], &field[0], &field_stt[0],
                &alpha[0], beta, u0, gamma, n)



def compute_llg_stt_cpp(double [:] dm_dt,
                        double [:] spin,
                        double [:] field,
                        double [:] p,
                        double [:] alpha,
                        int [:] pin,
                        double [:] a_J,
                        beta, gamma, n):
    llg_stt_cpp(&dm_dt[0], &spin[0], &field[0], &p[0],
                &alpha[0], &pin[0], &a_J[0], beta, gamma, n)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def compute_sd_spin(double [:] spin,
                    double [:] spin_last,
                    double [:] mxH,
                    double [:] mxmxH,
                    double [:] mxmxH_last,
                    double [:] tau,
                    int [:] pins,
                    n):

    sd_update_spin(&spin[0], &spin_last[0], &mxH[0],
                   &mxmxH[0], &mxmxH_last[0], &tau[0], &pins[0], n
                   )

def compute_sd_step(double [:] spin,
                    double [:] spin_last,
                    double [:] field,
                    double [:] scale,
                    double [:] mxH,
                    double [:] mxmxH,
                    double [:] mxmxH_last,
                    double [:] tau,
                    int [:] pins,
                    n, counter, tmin, tmax):

    sd_compute_step(&spin[0], &spin_last[0], &field[0], &scale[0], &mxH[0],
                    &mxmxH[0], &mxmxH_last[0], &tau[0], &pins[0],
                    n, counter, tmin, tmax
                    )
